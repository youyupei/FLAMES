# take a bam file as input and output read coverage for each position in the genome

import sys
import pysam
import numpy as np
import argparse
from parse_gene_anno import parseGFF3
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed


import helper

def parse_arg():
    parser = argparse.ArgumentParser(description='BAM file and GTF/GFF file parser')
    parser.add_argument('--bam_file', type=str, help='Path to the BAM file')
    parser.add_argument('--gtf_file', type=str, help='Path to the GTF/GFF file')
    parser.add_argument('--genome_ref', type=str, help='Path to the FASTA file for genome reference')
    parser.add_argument('--processes', type=int, default=1, help='Number of processes to use')
    args = parser.parse_args()


    # check if the file exists
    helper.check_exist([args.bam_file, args.gtf_file, args.genome_ref])

    # check the processes is valid (it should be a positive integer)
    if args.processes <= 0:
        helper.error_msg(f"Number of processes should be a positive integer, but got {args.processes}")
        sys.exit(1)

    return args


def merge_gene_group(gene_group):
    # make sure the strand are identical for all rows, else raise an warning
    if gene_group.Strand.nunique() != 1:
        helper.warning_msg(f"Strand information is not consistent for gene {gene_group.Gene_id.iloc[0]}"
                           "This gene will be skipped.")
        return None
    else:
        strand = gene_group.Strand.iloc[0]
        if strand not in ["+", "-"]:
            helper.warning_msg(f"Strand information is not valid for gene {gene_group.Gene_id.iloc[0]}"
                               "This gene will be skipped.")
            return None
    # make sure only 1 row in the group has Feature == "gene"
    if sum(gene_group.Feature == "gene") != 1:
        helper.warning_msg(f"In GTF file, Expect 1 gene with gene id {gene_group.Gene_id.iloc[0]}, but got {sum(gene_group.Feature == 'gene')}"
                           "This gene will be skipped.")
        return None
    else:
        rst_df = gene_group[gene_group.Feature == "gene"]
        # deep copy the gene information
        rst_df = rst_df.copy(deep=True)
        # add a new column to store the number of transcripts (concatenate the numbers separated by ",")
        # 
        if strand == "+":
            rst_df["TSS"] = ','.join([str(x) for x in gene_group[gene_group.Feature == "transcript"].Start])
        else:
            rst_df["TSS"] = ','.join([str(x) for x in gene_group[gene_group.Feature == "transcript"].End])
        return rst_df




# find plyA site
def find_polyA_site(gene, fasta_handle, window_size=10, 
                    minAorT=8, merge_dist=3, flank_butter=50):
    """
    find the polyA site for the gene
    return a list of tuple, each tuple is the start and end of the polyA site
    ordered by the position in the gene (decending  for + strand, acending for - strand)
    """
    seq = fasta_handle.fetch(gene.Chromosome, gene.Start, gene.End)
    if gene.Strand == "+":
        seq = np.array([int(x=='A') for x in seq])
    else:
        seq = np.array([int(x=='T') for x in seq])
    
    # get sliding window (size = 10) sum 
    window_sum = np.convolve(seq, np.ones(window_size), 'valid')
    # find location of consecutive window sum >= minAT
    polyA_site = np.where(window_sum >= minAorT)[0]
    if len(polyA_site):
        # merge consecutive polyA when the difference is <=3
        polyA_site = np.split(polyA_site, np.where(np.diff(polyA_site) > merge_dist)[0]+1)
        polyA_site = \
            [(gene.Start+x[0]-flank_butter, gene.Start+x[-1]+window_size+flank_butter) for x in polyA_site]
        
        if gene.Strand == "+":
            polyA_site = polyA_site[::-1]

        return polyA_site
    else:
        return []
def process_gene(grp, bam_file_path, fasta_file_path):
    """
    process the gene group
    """
    # read the  bam file
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    # read the fasta file
    fasta_file = pysam.FastaFile(fasta_file_path)
    # get the gene information
    gene = merge_gene_group(grp).iloc[0]
    # get strand
    reverse_strand = gene.Strand=='-'
    # find the polyA site
    polyA_site = find_polyA_site(gene, fasta_file)
    # filter the polyA site that overlap with a list of location
    filter_loc = [int(x) for x in gene.TSS.split(",")]
    polyA_site = \
        [x for x in polyA_site if not any([x[0]<y<x[1] for y in filter_loc])]

    # fetch the read for each site
    int_prim_list = []
    normal_list = []

    for read in bam_file.fetch(gene.Chromosome, gene.Start, gene.End):
        int_prim_flag = False
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue    

        # get length of soft clip
        soft_clip_start = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
        soft_clip_end = read.cigartuples[-1][1] if read.cigartuples[0][0] == 4 else 0
        
        for start,end in polyA_site:

            if reverse_strand:
                if start <= read.reference_start-soft_clip_start <= end:
                    int_prim_flag = True
                    int_prim_list.append(read.to_string())
                    break
                elif read.reference_start-soft_clip_start < start:
                    break

            else:
                
                if start <= read.reference_end+soft_clip_end <= end:
                    int_prim_flag = True
                    int_prim_list.append(read.to_string())
                    break
                elif end < read.reference_end+soft_clip_end:
                    break
        if not int_prim_flag:
            normal_list.append(read.to_string())

    bam_file.close()
    fasta_file.close()

    return int_prim_list, normal_list


def main():
    
    args = parse_arg()
    # Read the GTF file into a pandas DataFrame
    gtf_data = pd.read_csv(args.gtf_file, sep="\t", comment="#", header=None)
    # remove unneeded information
    gtf_data = gtf_data[gtf_data[2].isin(["gene", "transcript"])]
    gtf_data = gtf_data.iloc[:, [0,2,3,4,6,8]]
    # Extract the gene name from the last column
    gtf_data[8] = gtf_data[8].str.extract('gene_id "(.+?)";', expand=False)
    # add column names
    gtf_data.columns = ["Chromosome", "Feature", "Start", "End", "Strand", "Gene_id"]


    # group by gene name
    gtf_gene_group = gtf_data.groupby('Gene_id')
    template = pysam.AlignmentFile(args.bam_file, "rb")

    # output_normal_bam = pysam.AlignmentFile("normal.bam", "wb", 
    #                                         template=template)

    # output_int_prim = pysam.AlignmentFile("int_prim.bam", "wb", 
    #                                         template=template)
    output_normal_bam = open("normal.sam", "w")
    output_normal_bam.write(template.text)
    output_int_prim = open("int_prim.sam", "w")
    output_int_prim.write(template.text)
    if args.processes == 1:
        for i,(_,grp) in tqdm(enumerate(gtf_gene_group)):
            # if i > 10: # for testing
            #     break
            int_prim_list, normal_list = \
                process_gene(grp, args.bam_file, args.genome_ref)
            [output_normal_bam.write(x+'\n') for x in normal_list]
            [output_int_prim.write(x+'\n') for x in int_prim_list]
        # close the file
        output_normal_bam.close()
        output_int_prim.close()
    else:
        for future in helper.multiprocessing_submit(
                            process_gene, 
                            (x for _,x in gtf_gene_group),
                            pbar=True,
                            pbar_unit = "gene_group",
                            preserve_order = False,
                            n_process=args.processes, 
                            bam_file_path=args.bam_file, 
                            fasta_file_path=args.genome_ref):
        
            int_prim_list, normal_list = future.result()
            [output_normal_bam.write(x+'\n') for x in normal_list]
            [output_int_prim.write(x+'\n') for x in int_prim_list]
        # close the file
        output_normal_bam.close()
        output_int_prim.close()
    
    # Finished
    print("Finished")
        


if __name__ == "__main__":
    main()