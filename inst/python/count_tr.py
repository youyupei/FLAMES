# quantify transcript

import os
import sys
import gzip
import numpy as np
import fast_edit_distance 
from itertools import groupby
from collections import Counter, defaultdict, namedtuple
import multiprocessing as mp
import concurrent.futures
import pysam as ps
from datetime import datetime

from parse_gene_anno import parse_gff_tree
from filter_gff import annotate_filter_gff, annotate_full_splice_match_all_sample

import helper


# log errors when calling via R/reticulate:
# import traceback
# sys.stderr = open('stderr.out', 'a')
# try:
#     print(2/0) # code with error
# except:
#     traceback.print_exception(*sys.exc_info())


def umi_dedup(l, has_UMI, max_ed=1):
    """
    Deduplicate UMI by allowing edit distance <= max_ed
    Input:
        l: list of UMI
        has_UMI: whether the input has UMI
        max_ed: maximum edit distance allowed for UMI deduplication

    Return:
        tuple of (read count, dedup UMI count) 
            Note: if no UMI, the dedup UMI count is the same as read count
    """
    read_cnt = len(l)
    
    if has_UMI:
        dup_cnt = Counter(l)
        l_cnt = sorted(dup_cnt.most_common(), key=lambda x: (x[1], x[0]), reverse=True)
        if len(l_cnt) == 1:
            return 1, 1
        
        # remove UMI with edit distance <= max_ed
        rm_umi = {}
        for ith in range(len(l_cnt)-1):
            for jth in range(len(l_cnt)-1, ith, -1):  # first assess the low abundant UMI
                if l_cnt[jth][0] not in rm_umi:
                    if fast_edit_distance.edit_distance(
                            l_cnt[ith][0],l_cnt[jth][0], max_ed) <= max_ed:
                        rm_umi[l_cnt[jth][0]] = 1
                        l_cnt[ith] = (l_cnt[ith][0], l_cnt[ith][1]+l_cnt[jth][1])
        # return read count, dedup UMI count
        return read_cnt, len(l_cnt)-len(rm_umi)
    else:
        return read_cnt, read_cnt


def wrt_tr_to_csv(bc_tr_count_dict, transcript_dict, csv_f, 
                    transcript_dict_ref=None, has_UMI=False,print_saturation=False):
    """
    Write transcript count to csv file
    """
    f = gzip.open(csv_f, "wt")
    all_tr = set()
    for bc in bc_tr_count_dict:
        all_tr.update(list(bc_tr_count_dict[bc].keys()))
    all_tr = list(all_tr)
    
    # write header
    f.write("transcript_id,gene_id," +
            ",".join([x for x in bc_tr_count_dict])+"\n")
    for k, v in bc_tr_count_dict.items():
        if v is None:
            print (k,v)
    tr_cnt = {}
    dup_count = ()
    for tr in all_tr:
        read_counts, dup_counts = zip(
            *[umi_dedup(bc_tr_count_dict[bc][tr], has_UMI)
                 if tr in bc_tr_count_dict[bc] else (0,0) for bc in bc_tr_count_dict]
                 )
        
        tr_cnt[tr] = sum(dup_counts)

        if tr in transcript_dict:
            f.write(
                "{},{},".format(tr, transcript_dict[tr].parent_id))
        elif (transcript_dict_ref is not None) and (tr in transcript_dict_ref):
            f.write(
                "{},{},".format(tr, transcript_dict_ref[tr].parent_id))
        else:
            print("cannot find transcript in transcript_dict:", tr)
            exit(1)
        f.write(",".join([str(x) for x in dup_counts])+"\n")
    f.close()
    
    # print saturation if requested
    if print_saturation and has_UMI and sum(dup_counts):
        helper.green_msg(f"The estimated saturation is {1-len(dup_count)/sum(dup_count)}")

    return tr_cnt


def make_bc_dict(bc_anno):
    """
    Creates a dictionary of barcodes and their corresponding samples from a barcode annotation file.

    Args:
    bc_anno (str): The path to the barcode annotation file.

    Returns:
    dict: A dictionary with barcodes as keys and their corresponding samples as values.
    """
    with open(bc_anno) as f:
        # skip header line
        f.readline()

        bc_dict = dict()
        for line in f:
            line_vals = line.rstrip().split(',')
            sample = line_vals[0]
            bc = line_vals[1]

            bc_dict[bc] = sample

    return(bc_dict)


def query_len(cigar_string, hard_clipping=False):
    """
    https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
    Given a CIGAR string, return the number of bases consumed from the
    query sequence.
    CIGAR is a sequence of the form <operations><operator> such that operations is an integer giving 
    the number of times the operator is used
    M = match
    I = Insertion
    S = Soft clipping
    = = sequence match
    X = sequence mismatch
    """
    if (not hard_clipping):
        read_consuming_ops = ("M", "I", "S", "=", "X")
    else:
        read_consuming_ops = {"M", "I", "S", "H", "=", "X"}
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    this_len = 0
    for is_digit, val in cig_iter:
        if is_digit:
            this_len = int(''.join(val))
        else:
            if ''.join(val) in read_consuming_ops:
                result += this_len
    return result


def parse_realigned_bam(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, bc_file = False):
    """
    """

    fa_idx = dict((it.strip().split()[0], int(
        it.strip().split()[1])) for it in open(fa_idx_f))
    tr_count_dict = defaultdict(int)
    bc_tr_count_dict = {}
    bc_tr_badcov_count_dict = {}
    #tr_cov_dict = {}
    read_dict = {}
    aligned_read = namedtuple('aligned_read', ['tr', 'AS', 'tr_cov', 'read_cov', 'mapq'])
    cnt_stat = Counter()
    bamfile = ps.AlignmentFile(bam_in, "rb")
    # gene_names = bamfile.references

    if bc_file:
        bc_dict = make_bc_dict(bc_file)
    for rec in bamfile.fetch(until_eof=True):
        # update: treat mapping quality 0 as unmapped
        if rec.is_unmapped or rec.mapping_quality==0:
            cnt_stat["unmapped"] += 1
            continue

        # create a aligned_read namedtuple
        rec_align_info = aligned_read(
            tr = rec.reference_name, 
            AS = rec.get_tag("AS"), 
            tr_cov = (rec.reference_end-rec.reference_start)/fa_idx[rec.reference_name],
            read_cov = float(rec.query_alignment_length)/rec.infer_read_length(), 
            mapq = rec.mapping_quality)

        # map_st = rec.reference_start
        # map_en = rec.reference_end
        # tr = rec.reference_name
        # tr_cov = float(map_en-map_st)/fa_idx[tr]
        #tr_cov_dict.setdefault(tr, []).append(tr_cov)
        tr_count_dict[rec_align_info.tr] += 1

        if rec.query_name not in read_dict:
            read_dict.setdefault(rec.query_name, []).append(rec_align_info)
    
        else:
            if rec_align_info.AS > read_dict[rec.query_name][0].AS:
                read_dict[rec.query_name].insert(0, rec_align_info)
    
            # if same alignment score choose the one with higher transcript coverage, might be internal TSS
            elif rec_align_info.AS == read_dict[rec.query_name][0].AS and \
                    rec_align_info.read_cov == read_dict[rec.query_name][0].read_cov:
                
                if rec_align_info.tr_cov > read_dict[rec.query_name][0].tr_cov:
                    read_dict[rec.query_name].insert(0, rec_align_info)
            else:
                read_dict[rec.query_name].append(rec_align_info)
        
        if rec_align_info.tr not in fa_idx:
            # This shouldn't happen but worth checking
            cnt_stat["not_in_annotation"] += 1
            print("\t" + str(rec_align_info.tr), "not in annotation ???")
    
    # TODO: The original code defined that a transcript is kept if it has > min_sup_reads reads with > 90% coverage
    # However, it may result in a situation where a transcript with alternative TSS/TTS is discarded. I don't think 
    # it is a good idea to discard such transcripts. The alternative TSS/TTS should be considered at the isoform
    # identification step but not the quantification step.
    # tr_kept = dict((tr, tr) for tr in tr_cov_dict if len(
    #     [it for it in tr_cov_dict[tr] if it > 0.9]) > min_sup_reads)

    # TEMP_FIX: keep all transcripts
    tr_kept = set(tr for tr, cnt in tr_count_dict.items() if cnt > min_sup_reads)
    #unique_tr_count = Counter(read_dict[r][0][0]
    #                          for r in read_dict if read_dict[r][0][2] > 0.9)

    for r in read_dict:
        # filter out aligned transcripts that are not in tr_kept
        filtered_alignment = [it for it in read_dict[r] if it.tr in tr_kept]
        if len(filtered_alignment) > 0:
            hit = filtered_alignment[0]  # namedtuple('aligned_read', ['tr', 'AS', 'tr_cov', 'read_cov', 'mapq'])
        else:
            cnt_stat["no_good_match"] += 1
            continue
        # below line creates issue when header line has more than one _.
        # in this case, umi is assumed to be delimited from the barcode by the last _
        # bc, umi = r.split("#")[0].split("_")  # assume cleaned barcode
        try:
            bc, umi = r.split("#")[0].split("_")  # assume cleaned barcode
        except ValueError as ve:
            print(ve, ": ", ve.args, ".")
            raise ValueError(
                "Please check if barcode and UMI are delimited by \"_\"")

        if bc_file:
            bc = bc_dict[bc]
        
        if hit.tr_cov > min_tr_coverage and hit.read_cov > min_read_coverage:
            if len(filtered_alignment) > 1 and \
                    hit.AS == filtered_alignment[1].AS and \
                    hit.tr_cov == filtered_alignment[1].tr_cov:
                cnt_stat["ambigious_reads"] += 1
                if bc not in bc_tr_badcov_count_dict:
                    bc_tr_badcov_count_dict[bc] = {}
                bc_tr_badcov_count_dict[bc].setdefault(hit.tr, []).append(umi)

            else:
                if bc not in bc_tr_count_dict:
                    bc_tr_count_dict[bc] = {}
                bc_tr_count_dict[bc].setdefault(hit.tr, []).append(umi)
                cnt_stat["counted_reads"] += 1
        else:
            cnt_stat["not_enough_coverage"] += 1
            if bc not in bc_tr_badcov_count_dict:
                bc_tr_badcov_count_dict[bc] = {}
            bc_tr_badcov_count_dict[bc].setdefault(hit.tr, []).append(umi)
    print(("\t" + str(cnt_stat)))
    return bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept
    
def realigment_min_sup_reads_filter(bam_list, fa_idx_f, min_sup_reads):
    fa_idx = dict((it.strip().split()[0], int(
        it.strip().split()[1])) for it in open(fa_idx_f))
    tr_cov_dict = {}
    bamfiles = [ps.AlignmentFile(i, "rb") for i in bam_list]

    for bamfile in bamfiles:
        for rec in bamfile.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            tr = rec.reference_name
            tr_cov = float(rec.reference_end-rec.reference_start)/fa_idx[tr]
            tr_cov_dict[tr] = tr_cov_dict.setdefault(tr, 0) + (1 if tr_cov > 0.9 else 0)

    return set(tr for tr in tr_cov_dict if tr_cov_dict[tr] >= min_sup_reads)


def parse_realigned_bam_sc_multi_sample(bam_in, fa_idx_f, tr_kept, min_tr_coverage, min_read_coverage):
    """
    """
    fa_idx = dict((it.strip().split()[0], int(
        it.strip().split()[1])) for it in open(fa_idx_f))
    bc_tr_count_dict = {}
    bc_tr_badcov_count_dict = {}
    read_dict = {}
    cnt_stat = Counter()
    bamfile = ps.AlignmentFile(bam_in, "rb")

    for rec in bamfile.fetch(until_eof=True):
        if rec.is_unmapped:
            cnt_stat["unmapped"] += 1
            continue
        elif rec.reference_name not in tr_kept:
            continue
        tr = rec.reference_name
        tr_cov = float(rec.reference_end-rec.reference_start)/fa_idx[tr]
        if rec.query_name not in read_dict:
            read_dict.setdefault(rec.query_name, []).append((tr, rec.get_tag("AS"), tr_cov, float(
                rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
        else:
            if rec.get_tag("AS") > read_dict[rec.query_name][0][1]:
                read_dict[rec.query_name].insert(0, (tr, rec.get_tag("AS"), tr_cov, float(
                    rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
            # same aligned sequence
            elif rec.get_tag("AS") == read_dict[rec.query_name][0][1] and float(rec.query_alignment_length)/rec.infer_read_length() == read_dict[rec.query_name][0][3]:
                # choose the one with higher transcript coverage, might be internal TSS
                if tr_cov > read_dict[rec.query_name][0][2]:
                    read_dict[rec.query_name].insert(0, (tr, rec.get_tag("AS"), tr_cov, float(
                        rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
            else:
                read_dict[rec.query_name].append((tr, rec.get_tag("AS"), tr_cov, float(
                    rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
        if tr not in fa_idx:
            cnt_stat["not_in_annotation"] += 1
            print("\t" + str(tr), "not in annotation ???")
    for r in read_dict:
        tmp = read_dict[r]
        tmp = [it for it in tmp if it[0] in tr_kept]
        if len(tmp) > 0:
            hit = tmp[0]  # transcript_id, pct_ref, pct_reads
        else:
            cnt_stat["no_good_match"] += 1
            continue
        # below line creates issue when header line has more than one _.
        # in this case, umi is assumed to be delimited from the barcode by the last _
        # bc, umi = r.split("#")[0].split("_")  # assume cleaned barcode
        try:
            bc, umi = r.split("#")[0].split("_")  # assume cleaned barcode
        except ValueError as ve:
            print(ve, ": ", ve.args, ".")
            raise ValueError(
                "Please check if barcode and UMI are delimited by \"_\"")

        if len(tmp) == 1 and tmp[0][4] > 0:
            if bc not in bc_tr_count_dict:
                bc_tr_count_dict[bc] = {}
            bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
            cnt_stat["counted_reads"] += 1
        elif len(tmp) > 1 and tmp[0][1] == tmp[1][1] and tmp[0][3] == tmp[1][3]:
            if hit[1] > 0.8:
                if bc not in bc_tr_count_dict:
                    bc_tr_count_dict[bc] = {}
                bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
                cnt_stat["counted_reads"] += 1
            else:
                cnt_stat["ambigious_reads"] += 1
                if bc not in bc_tr_badcov_count_dict:
                    bc_tr_badcov_count_dict[bc] = {}
                bc_tr_badcov_count_dict[bc].setdefault(hit[0], []).append(umi)
        elif hit[2] < min_tr_coverage or hit[3] < min_read_coverage:
            cnt_stat["not_enough_coverage"] += 1
            if bc not in bc_tr_badcov_count_dict:
                bc_tr_badcov_count_dict[bc] = {}
            bc_tr_badcov_count_dict[bc].setdefault(hit[0], []).append(umi)
        else:
            if bc not in bc_tr_count_dict:
                bc_tr_count_dict[bc] = {}
            bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
            cnt_stat["counted_reads"] += 1
    print(("\t" + str(cnt_stat)))
    return bc_tr_count_dict, bc_tr_badcov_count_dict

def parse_realigned_bam_bulk(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage):
    """
    """
    print("Inputs: ", bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage)
    fa_idx = dict((it.strip().split()[0], int(
        it.strip().split()[1])) for it in open(fa_idx_f))
    bc_tr_count_dict = {}
    bc_tr_badcov_count_dict = {}
    tr_cov_dict = {}
    read_dict = {}
    cnt_stat = Counter()
    bamfiles = [ps.AlignmentFile(i, "rb") for i in bam_in]

    for bamfile in bamfiles:
        sample_name = os.path.basename(bamfile.filename.decode()).replace('_realign2transcript.bam','')
        for rec in bamfile.fetch(until_eof=True):
            if rec.is_unmapped:
                cnt_stat["unmapped"] += 1
                continue
            map_st = rec.reference_start
            map_en = rec.reference_end
            tr = rec.reference_name
            tr_cov = float(map_en-map_st)/fa_idx[tr]
            tr_cov_dict.setdefault(tr, []).append(tr_cov)
            if sample_name + "_" + rec.query_name not in read_dict:
                read_dict.setdefault(sample_name + "_" + rec.query_name, []).append((tr, rec.get_tag("AS"), tr_cov, float(
                    rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
            else:
                if rec.get_tag("AS") > read_dict[sample_name + "_" + rec.query_name][0][1]:
                    read_dict[sample_name + "_" + rec.query_name].insert(0, (tr, rec.get_tag("AS"), tr_cov, float(
                        rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
                # same aligned sequence
                elif rec.get_tag("AS") == read_dict[sample_name + "_" + rec.query_name][0][1] and float(rec.query_alignment_length)/rec.infer_read_length() == read_dict[sample_name + "_" + rec.query_name][0][3]:
                    # choose the one with higher transcript coverage, might be internal TSS
                    if tr_cov > read_dict[sample_name + "_" + rec.query_name][0][2]:
                        read_dict[sample_name + "_" + rec.query_name].insert(0, (tr, rec.get_tag("AS"), tr_cov, float(
                            rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
                else:
                    read_dict[sample_name + "_" + rec.query_name].append((tr, rec.get_tag("AS"), tr_cov, float(
                        rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
            if tr not in fa_idx:
                cnt_stat["not_in_annotation"] += 1
                print("\t" + str(tr), "not in annotation ???")
    tr_kept = dict((tr, tr) for tr in tr_cov_dict if len(
        [it for it in tr_cov_dict[tr] if it > 0.9]) > min_sup_reads)
    #unique_tr_count = Counter(read_dict[r][0][0]
    #                          for r in read_dict if read_dict[r][0][2] > 0.9)
    for r, v in read_dict.items():
        tmp = [it for it in v if it[0] in tr_kept]
        if len(tmp) > 0:
            hit = tmp[0]  # transcript_id, pct_ref, pct_reads
        else:
            cnt_stat["no_good_match"] += 1
            continue
        # below line creates issue when header line has more than one _.
        # in this case, umi is assumed to be delimited from the barcode by the last _
        # bc, umi = r.split("#")[0].split("_")  # assume cleaned barcode
        r_split = r.split("#")[0].split("_")
        if len(r_split) == 2:
            bc, umi = r_split
        elif len(r_split) > 2: # when '_' in file names
            umi = r_split[-1]
            bc = "_".join(r_split[:-1])
        else:
            raise ValueError("Please check if barcode and UMI are delimited by \"_\":\n" + r_split)

        if len(tmp) == 1 and tmp[0][4] > 0:
            if bc not in bc_tr_count_dict:
                bc_tr_count_dict[bc] = {}
            bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
            cnt_stat["counted_reads"] += 1
        elif len(tmp) > 1 and tmp[0][1] == tmp[1][1] and tmp[0][3] == tmp[1][3]:
            if hit[1] > 0.8:
                if bc not in bc_tr_count_dict:
                    bc_tr_count_dict[bc] = {}
                bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
                cnt_stat["counted_reads"] += 1
            else:
                cnt_stat["ambigious_reads"] += 1
                if bc not in bc_tr_badcov_count_dict:
                    bc_tr_badcov_count_dict[bc] = {}
                bc_tr_badcov_count_dict[bc].setdefault(hit[0], []).append(umi)
        elif hit[2] < min_tr_coverage or hit[3] < min_read_coverage:
            cnt_stat["not_enough_coverage"] += 1
            if bc not in bc_tr_badcov_count_dict:
                bc_tr_badcov_count_dict[bc] = {}
            bc_tr_badcov_count_dict[bc].setdefault(hit[0], []).append(umi)
        else:
            if bc not in bc_tr_count_dict:
                bc_tr_count_dict[bc] = {}
            bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
            cnt_stat["counted_reads"] += 1
    print(("\t" + str(cnt_stat)))
    return bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept

def tr_len_range(l):
    """
    [0,500] -> 0
    [500,1000] -> 1
    [1000,1500] -> 2
    [1500,2000] -> 3
    [2000,Inf] -> 4
    """
    return(min(4, l/500))


def realigned_bam_coverage(bam_in, fa_idx_f, coverage_dir):
    fa_idx = dict((it.strip().split()[0], int(
        it.strip().split()[1])) for it in open(fa_idx_f))
    left_clip_count = Counter()
    right_clip_count = Counter()
    tr_strand = Counter()
    bc_pct = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}}
    bc_cov_pct = {0: [], 1: [], 2: [], 3: [], 4: []}
    gene_pct = {0: [], 1: [], 2: [], 3: [], 4: []}
    bamfile = ps.AlignmentFile(bam_in, "rb")
    for rec in bamfile.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_supplementary or rec.is_secondary:
            continue
        bc, umi = rec.query_name.split("#")[0].split(
            "_")  # assume cleaned barcode
        map_st = rec.reference_start
        map_en = rec.reference_end
        tr = rec.reference_name
        if float(map_en-map_st)/fa_idx[tr] < 0.3:
            continue
        if rec.cigar[0][0] == 4:  # BAM_CSOFT_CLIP
            left_clip_count[rec.cigar[0][1]] += 1
        if rec.cigar[-1][0] == 4:  # BAM_CSOFT_CLIP
            right_clip_count[rec.cigar[-1][1]] += 1
        tr_strand[rec.is_reverse] += 1
        if not rec.is_reverse:
            pass
        gene_pct[tr_len_range(fa_idx[tr])].append(
            float(map_en-map_st)/fa_idx[tr])
        bc_pct[tr_len_range(fa_idx[tr])].setdefault(
            bc, []).append(float(map_st-0)/fa_idx[tr])
        bc_pct[tr_len_range(fa_idx[tr])].setdefault(
            bc, []).append(float(map_en-0)/fa_idx[tr])
        bc_cov_pct[tr_len_range(fa_idx[tr])].append(
            float(map_en-map_st)/fa_idx[tr])
    print(left_clip_count.most_common(30))
    print(right_clip_count.most_common(30))
    print(tr_strand)
    print(np.histogram(
        bc_pct[0][list(bc_pct[0].keys())[0]], bins=200, range=(0, 1)))
    for i in bc_pct:
        coverage_f = open(os.path.join(
            coverage_dir, "transcript_cov_per_cell.{}.csv".format(i)), "w")
        for bc in bc_pct[i]:
            lhi, _ = np.histogram(bc_pct[i][bc], bins=200, range=(0, 1))
            coverage_f.write("{},".format(bc)+",".join(str(it)
                             for it in lhi)+"\n")
        coverage_f.close()
    tr_cov_f = open(os.path.join(coverage_dir, "transcript_cov.csv"), "w")
    for i in gene_pct:
        lhi, _ = np.histogram(gene_pct[i], bins=200, range=(0, 1))
        tr_cov_f.write("{},".format(i)+",".join(str(it) for it in lhi)+"\n")
    tr_cov_f.close()

# this is the main function
def quantification(config_dict, annotation, outdir, pipeline):
    transcript_fa_idx = os.path.join(outdir, "transcript_assembly.fa.fai")
    isoform_gff3 = os.path.join(outdir, "isoform_annotated.gtf" if os.path.isfile(os.path.join(outdir, "isoform_annotated.gtf")) else "isoform_annotated.gff3")
    isoform_gff3_f = os.path.join(outdir, "isoform_annotated.filtered.gff3")
    FSM_anno_out = os.path.join(outdir, "isoform_FSM_annotation.csv")

    # parsing gff file in background (concurrent future)
    executor = concurrent.futures.ProcessPoolExecutor(5)
    futures = {}
    futures['parse_gff_tree_anno'] = executor.submit(parse_gff_tree, annotation)
    futures['parse_gff_tree_iso'] = executor.submit(parse_gff_tree, isoform_gff3)

    if pipeline == "sc_single_sample":
        realign_bam = os.path.join(outdir, "realign2transcript.bam")
        tr_cnt_csv = os.path.join(outdir, "transcript_count.csv.gz")
        tr_badcov_cnt_csv = os.path.join(outdir, "transcript_count.bad_coverage.csv.gz")
        bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept = parse_realigned_bam(
            realign_bam,
            transcript_fa_idx,
            config_dict["isoform_parameters"]["min_sup_cnt"],
            config_dict["transcript_counting"]["min_tr_coverage"],
            config_dict["transcript_counting"]["min_read_coverage"],
            bc_file = False)
        
        chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = futures['parse_gff_tree_anno'].result()
        chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = futures['parse_gff_tree_iso'].result()

        tr_cnt = wrt_tr_to_csv(bc_tr_count_dict, transcript_dict_i, tr_cnt_csv,
                            transcript_dict, "UMI" in config_dict["barcode_parameters"]["pattern"].keys())
        wrt_tr_to_csv(bc_tr_badcov_count_dict, transcript_dict_i, tr_badcov_cnt_csv,
                    transcript_dict, "UMI" in config_dict["barcode_parameters"]["pattern"].keys(),
                    print_saturation = False)
        annotate_filter_gff(isoform_gff3, annotation, isoform_gff3_f, FSM_anno_out,
                        tr_cnt, config_dict["isoform_parameters"]["min_sup_cnt"], verbose=False)
        return

    elif pipeline == "bulk":
        realign_bam = [os.path.join(outdir, f) for f in os.listdir(outdir) if f[-22:] == "realign2transcript.bam"]
        tr_cnt_csv = os.path.join(outdir, "transcript_count.csv.gz")
        tr_badcov_cnt_csv = os.path.join(outdir, "transcript_count.bad_coverage.csv.gz")
        bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept = parse_realigned_bam_bulk(
            realign_bam,
            transcript_fa_idx,
            config_dict["isoform_parameters"]["min_sup_cnt"],
            config_dict["transcript_counting"]["min_tr_coverage"],
            config_dict["transcript_counting"]["min_read_coverage"])


        chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = futures['parse_gff_tree_anno'].result()
        chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = futures['parse_gff_tree_iso'].result()

        tr_cnt = wrt_tr_to_csv(
            bc_tr_count_dict, transcript_dict_i, tr_cnt_csv,
                transcript_dict, has_UMI=False, print_saturation = False)
        wrt_tr_to_csv(
            bc_tr_badcov_count_dict, transcript_dict_i, tr_badcov_cnt_csv,
                      transcript_dict, has_UMI=False, print_saturation = False)
        annotate_filter_gff(
            isoform_gff3, annotation, isoform_gff3_f, FSM_anno_out, tr_cnt, 
                config_dict["isoform_parameters"]["min_sup_cnt"], verbose=False)
        return

    elif pipeline == "sc_multi_sample":
        realign_bam = [os.path.join(outdir, f) for f in os.listdir(outdir) if f[-23:] == "_realign2transcript.bam"]
        tr_kept = realigment_min_sup_reads_filter(realign_bam, transcript_fa_idx, config_dict["isoform_parameters"]["min_sup_cnt"])
        chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = futures['parse_gff_tree_anno'].result()
        chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = futures['parse_gff_tree_iso'].result()

        for sample_bam in  realign_bam:
            sys.stderr.write("parsing " + sample_bam + "...\n")
            sample = os.path.basename(sample_bam).replace('_realign2transcript.bam','')
            tr_cnt_csv = os.path.join(outdir, sample+ "_"+"transcript_count.csv.gz")
            tr_badcov_cnt_csv = os.path.join(outdir, sample+ "_"+"transcript_count.bad_coverage.csv.gz")
            bc_tr_count_dict, bc_tr_badcov_count_dict  = parse_realigned_bam_sc_multi_sample(
                sample_bam,
                transcript_fa_idx,
                tr_kept,
                config_dict["transcript_counting"]["min_tr_coverage"],
                config_dict["transcript_counting"]["min_read_coverage"])
            sys.stderr.write("parsing " + sample_bam + "done\n")
            tr_cnt = wrt_tr_to_csv(bc_tr_count_dict, transcript_dict_i, tr_cnt_csv,
                                   transcript_dict, "UMI" in config_dict["barcode_parameters"]["pattern"].keys())
            sys.stderr.write("wrt_tr_to_csv for" + sample_bam + "done\n")
            wrt_tr_to_csv(bc_tr_badcov_count_dict, transcript_dict_i, tr_badcov_cnt_csv,
                          transcript_dict, "UMI" in config_dict["barcode_parameters"]["pattern"].keys(),
                          print_saturation = False)
            del bc_tr_count_dict, bc_tr_badcov_count_dict, tr_cnt
            ##gc.collect()

        sys.stderr.write("annotate_full_splice_match_all_sample...\n")
        annotate_full_splice_match_all_sample(FSM_anno_out, isoform_gff3, annotation)

        return

    else:
        raise ValueError(f"Unknown pipeline type {pipeline}")


