# take a bam file as input and output read coverage for each position in the genome

import sys
import pysam
import numpy as np
import argparse
import pandas as pd
from tqdm import tqdm
import helper


# TODO:
# 1. remove potential internal priming site on 3' UTR


def parse_arg():
    parser = argparse.ArgumentParser(description='BAM file and GTF/GFF file parser')
    parser.add_argument('--bam_file', type=str, help='Path to the BAM file')
    parser.add_argument('--gtf_file', type=str, help='Path to the GTF/GFF file')
    parser.add_argument('--genome_ref', type=str, help='Path to the FASTA file for genome reference')
    parser.add_argument('--processes', type=int, default=1, help='Number of processes to use')
    # add a output file prefix
    parser.add_argument('--output_prefix', type=str, default="", help='Prefix for the output file')

    args = parser.parse_args()


    # check if the file exists
    helper.check_exist([args.bam_file, args.gtf_file, args.genome_ref])

    # check the processes is valid (it should be a positive integer)
    if args.processes <= 0:
        helper.error_msg(f"Number of processes should be a positive integer, but got {args.processes}")
        sys.exit(1)

    return args

def get_utr_from_transcript_group(grp):
    if not "UTR" in grp.Feature.values:
        return []
    # if the strand columns a uniquely +
    rst_3_utr = []
    if len(grp.Strand.unique()) == 1 and grp.Strand.unique()[0] == "+":
        for row in grp.sort_values('Start', ascending=False).itertuples():
            if row.Feature == "UTR":
                rst_3_utr.append((row.Start, row.End))
                break
            if row.Feature == "CDS":
                break
        return sorted(rst_3_utr)
    elif len(grp.Strand.unique()) == 1 and grp.Strand.unique()[0] == "-":
        for row in grp.sort_values('Start').itertuples():
            if row.Feature == "UTR":
                rst_3_utr.append((row.Start, row.End))
                break
            if row.Feature == "CDS":
                break
        return rst_3_utr
    else:
        helper.warning_msg(f"Strand information is not consistent for transcript in GTF."
                           "3' and 5' can not be distinguished, used both transcript_id: {grp.transcript_id.unique()}")
        rst_3_utr = list(zip(grp[grp.Feature == "UTR"].Start, grp[grp.Feature == "UTR"].End))
        if len(rst_3_utr) == 1:
            return rst_3_utr
        else:
            rst_3_utr = sorted(rst_3_utr)
            return [rst_3_utr[0], rst_3_utr[-1]]

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
            rst_df["TTS"] = ','.join([str(x) for x in gene_group[gene_group.Feature == "transcript"].End])
        else:
            rst_df["TTS"] = ','.join([str(x) for x in gene_group[gene_group.Feature == "transcript"].Start])
        return rst_df




# find plyA site
def find_polyA_site(gene, fasta_handle, window_size=10, 
                    minAorT=8, merge_dist=3, flank_butter=25):
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
    process the gene group (single process)
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
    filter_loc = [int(x) for x in gene.TTS.split(",")]

    # # number of polyA site before filtering
    # num_polyA_site = len(polyA_site)

    # filter the polyA site that overlap with the Known TTS
    polyA_site = \
        [x for x in polyA_site if not any([x[0]<y<x[1] for y in filter_loc])]


    # filter the polyA site that overlap with the 3' UTR
    polyA_site = \
        [x for x in polyA_site if not any(x[0]<y[1] and y[0]<x[1] for y in gene.UTR3)]

    # fetch the read for each site
    int_prim_sam_text = ''
    normal_sam_text = ''
    int_prim_counter = 0
    normal_read_counter = 0
    for read in bam_file.fetch(gene.Chromosome, gene.Start, gene.End):
        int_prim_flag = False
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue    

        # get length of soft clip or hard clip
        # clip_start = read.cigartuples[0][1] if read.cigartuples[0][0] in [4,5] else 0
        # clip_end = read.cigartuples[-1][1] if read.cigartuples[0][0] in [4,5] else 0

        for start,end in polyA_site:

            if reverse_strand:
                #if start <= read.reference_start-clip_start <= end or \
                if start <= read.reference_start <= end:
                    int_prim_flag = True
                    int_prim_sam_text += read.to_string() + '\n'
                    int_prim_counter += 1
                    break
                #elif read.reference_start-clip_start < start:
                elif read.reference_start < start:
                    break

            else:
                
                #if start <= read.reference_end+clip_end <= end:
                if start <= read.reference_end <= end:
                    int_prim_flag = True
                    int_prim_sam_text+= read.to_string() + '\n'
                    int_prim_counter += 1
                    break
                #elif end < read.reference_end+clip_end:
                elif end < read.reference_end:
                    break
        if not int_prim_flag:
            normal_sam_text += read.to_string() + '\n'
            normal_read_counter += 1

    bam_file.close()
    fasta_file.close()

    return int_prim_sam_text, normal_sam_text, int_prim_counter, normal_read_counter


def main():
    
    args = parse_arg()
    # Read the GTF file into a pandas DataFrame
    gtf_data = pd.read_csv(args.gtf_file, sep="\t", comment="#", header=None)
    
    # remove unneeded information
    gtf_data = gtf_data[gtf_data[2].isin(["gene", "transcript", "CDS", "UTR"])]
    gtf_data = gtf_data.iloc[:, [0,2,3,4,6,8]]
    
    # Extract the gene name from the last column
    gtf_data[9] = gtf_data[8].str.extract('transcript_id "(.+?)";', expand=False)
    gtf_data[8] = gtf_data[8].str.extract('gene_id "(.+?)";', expand=False)

    # add column names
    gtf_data.columns = ["Chromosome", "Feature", "Start", "End", "Strand", "Gene_id", "transcript_id"]

    # get 3' UTR
    utr3_df = gtf_data[gtf_data.Feature.isin(["CDS", "UTR"])].groupby(["Gene_id","transcript_id"])
    utr3_df = utr3_df.apply(get_utr_from_transcript_group)
    utr3_df = utr3_df.droplevel("transcript_id") 
    utr3_df= utr3_df.groupby("Gene_id").apply(lambda y: sorted(list(set().union(*y))))
    utr3_df.name = 'UTR3'

    # merge the 3' UTR information to the gtf_data
    # drop the transcript_id column
    gtf_data = gtf_data[gtf_data.Feature.isin(["gene", "transcript"])]
    gtf_data = gtf_data.drop("transcript_id", axis=1)
    gtf_data = gtf_data.merge(utr3_df, on="Gene_id", how="left")
    gtf_data['UTR3'] = gtf_data['UTR3'].fillna("").apply(list)
    

    # group by gene name
    gtf_gene_group = gtf_data.groupby('Gene_id')
    template = pysam.AlignmentFile(args.bam_file, "rb")

    # output_normal_bam = pysam.AlignmentFile("normal.bam", "wb", 
    #                                         template=template)

    # output_int_prim = pysam.AlignmentFile("int_prim.bam", "wb", 
    #                                         template=template)
    header = template.text

    out_fn_prefix = args.output_prefix
    if out_fn_prefix and not out_fn_prefix.endswith("_"):
        out_fn_prefix += "_"

    output_normal_bam = open(out_fn_prefix+"normal.sam", "w")
    output_normal_bam.write(header)
    output_int_prim = open(out_fn_prefix+"int_prim.sam", "w")
    output_int_prim.write(header)
    rst_int_prim_tot = 0
    rst_normal_read_tot = 0
    
    if args.processes == 1:
        for i,(_,grp) in tqdm(enumerate(gtf_gene_group)):
            int_prim_text, normal_text, int_prim_c, normal_read_c = \
                process_gene(grp, args.bam_file, args.genome_ref)
            output_normal_bam.write(normal_text)
            output_int_prim.write(int_prim_text)
            rst_int_prim_tot += int_prim_c
            rst_normal_read_tot += normal_read_c

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
        
            int_prim_text, normal_text, \
            int_prim_c, normal_read_c = future.result()

            rst_int_prim_tot += int_prim_c
            rst_normal_read_tot += normal_read_c
            output_normal_bam.write(normal_text)
            output_int_prim.write(int_prim_text)
            # [output_normal_bam.write(x+'\n') for x in normal_list]
            # [output_int_prim.write(x+'\n') for x in int_prim_list]
        # close the file
        output_normal_bam.close()
        output_int_prim.close()
    
    # Finished
    print("Finished")
    # 
    with open(out_fn_prefix+"summary.txt", 'w') as f:
        f.write(
            f"Total number of reads: {rst_int_prim_tot + rst_normal_read_tot}\n"
            f"Internal priming reads: {rst_int_prim_tot}\n"
            f"Other reads: {rst_normal_read_tot}\n"
            f"Proportion of internal priming reads: {rst_int_prim_tot/(rst_int_prim_tot+rst_normal_read_tot)}\n"
        )
        

if __name__ == "__main__":
    main()