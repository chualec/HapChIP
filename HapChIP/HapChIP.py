import argparse
import os
import operator
import sys
from itertools import combinations
import sys, getopt
from datetime import datetime
import pysam
from collections import Counter


def welcome(options):
    """
    This function displays a welcomd message
    and prints out the detected input files for:
      bam file
      VCF file
      output directory
      optional name
    """
    print("")
    print("*************")
    print("thank you for using the HapChIP")
    print("")
    if options.reads is not None:
        print("Detected bam file path: %s" % (options.reads))
    if options.variants is not None:
        print("Detected vcf file path: %s" % (options.variants))
    if options.out_dir is not None:
        print("Output directory: %s" % (options.out_dir))
    if options.name is not None:
        print("Detected output file name: %s" % (options.name))


def progressBar(current, total, end=False, text="NA"):
    """
    This function helps display a progress bar on the bottom of the screen while running.
    This helps users know whether the progress has stopped or not.
    """
    barLength = 20
    percent = float(current) * 100 / total
    raw_num = str(current) + "/" + str(total) + "\r"
    arrow = "-" * int(percent / 100 * barLength - 1) + ">"
    spaces = " " * (barLength - len(arrow))
    if not end:
        print("Progress: [%s%s] %d %% %s" % (arrow, spaces, percent, raw_num), end="\r")
    else:
        print("Progress: [%s%s] %d %% %s" % (arrow, spaces, percent, raw_num))
        print(text)
        print("")


# the cigar seq function expands the sequence information to account for soft clipping, insertions, and deletions
def cigarseq(bam_read):
    """
    This function reads the bam read's information from the CIGAR seq and alignment sequence.
    And returns a full-length sequence of of the bam read containing clipping, insertion, and deletion information.

    input:
      bam read object
    output:
      string
    """
    read_seq = bam_read.query_alignment_sequence
    quality_seq = bam_read.query_alignment_qualities
    cig_seq = ""
    # qual_seq = []
    for cigs in bam_read.cigar:
        if cigs[0] == 0:  # these indicate match
            cig_seq += read_seq[0 : cigs[1]]
            read_seq = read_seq[cigs[1] :]
            # qual_seq.extend(quality_seq[0:cigs[1]])
            # quality_seq = quality_seq[cigs[1]:]
        elif cigs[0] == 1:  # these indicate insertion
            pass
            read_seq = read_seq[cigs[1] :]
            # quality_seq = quality_seq[cigs[1]:]
        elif cigs[0] == 2:  # these indicate deletion
            cig_seq += "D" * cigs[1]
            # qual_seq.extend(["D"]*cigs[1])
        elif cigs[0] == 3:  # these indicate skipped
            cig_seq += "N" * cigs[1]
            # qual_seq.extend(["N"]*cigs[1])
    return cig_seq


def get_seq_base(bam_read, var_pos):
    """
    This function returns the base of the bam read based on position of variant.
    Input:
      bam read
      variant position (int)
    output:
      string
    """
    cigar_seq = cigarseq(bam_read)
    ref_base = var_pos - bam_read.pos - 1
    # return(cigar_seq[ref_base], str(qual_seq[ref_base]))
    return cigar_seq[ref_base]


def grab_genotype(rec):
    """
    This function grabs genotype information from the VCF file
    input:
      variant information
    output:
      string (genotype)
    """
    format_entry = str(rec).split()[-1]
    genotype = format_entry.split(":")[0]
    return genotype


def check_var(rec):
    """
    This function performs a QC check for the variant
    The criterias for exclusion are:
      multiple genotypes
      indel
      non-standard chromosome
    input:
      variant information
    output:
      True or False
    """
    checkok = True
    alleles = rec.alleles
    # check if indel
    if len(alleles[0]) > 1 or len(alleles[1]) > 1:
        checkok = False
    # check if het or phased
    elif grab_genotype(rec) == "1/2":
        checkok = False
    elif len(rec.chrom) > 5:
        checkok = False
    return checkok


def grab_haplotype(var_ref, var_alt, var_base, GT):
    """
    This function takes haplotype information and variant information
    to assign the genotype information to the variant
    input:
      variant reference (ATCG)
      variant alt (ATCG)
      variant base (ATCG)
      GT (genotype information)
    output:
      haplotype (represented as int)
        0 = unknown
        1 = 0|1
        2 = 1|0
    """
    hap = 0
    if var_base == var_ref and var_base != var_alt:
        if GT == "0|1":
            hap = 1
        elif GT == "1|0":
            hap = 2
    elif var_base != var_ref and var_base == var_alt:
        if GT == "0|1":
            hap = 2
        elif GT == "1|0":
            hap = 1
    return hap


def check_read(sam_read):
    """
    performs a QC check on the bam read
    removes if mapping quality is 0
    removes if no aligned read in bam file
    input:
      bam read
    output:
      True or False
    """
    if int(sam_read.mapping_quality) == 0:
        return False
    for cig in sam_read.cigar:
        if cig[0] == 0 and cig[1] > 1:
            return True
    return False


def reduce_hap(hap_list):
    """
    This function reduces haplotype lists to one output via hierchy
    0 = unknown
    1 = 0|1
    2 = 1|0
    3 = conflict
    input:
      list of strings representing haplotye
    output:
      int (0~3)
    """
    if 3 in hap_list:
        return 3
    elif 1 in hap_list and 2 in hap_list:
        return 3
    elif 1 in hap_list:
        return 1
    elif 2 in hap_list:
        return 2
    else:
        return 0


def parse_hichip(bam_in, vcf):
    """
    This is the main looping function.
    This loops through each variant position
    grab all bam reads for each variant position
    determines bam haplotype information for each variant-read pair
    and adds them to a dictionary
    input:
      bam file
      vcf file
    output:
      dictionary of read info
    """
    print("starting parse hichip...")
    var_count = 0
    num_lines = len(list(vcf.fetch()))
    readstring_dict = {}
    # for each line in the vcf file
    RX_library = {}
    for rec in vcf.fetch():
        var_count += 1
        if var_count % 1000 == 0:
            progressBar(var_count, num_lines)
        var_chr = rec.chrom
        var_pos = rec.pos
        var_ref = rec.ref
        GT = grab_genotype(rec)
        # check and do only het variants
        if check_var(rec) and (GT == "1|0" or GT == "0|1"):
            # for each read that overlaps the variant position
            for bam_read in bam_in.fetch(var_chr, var_pos - 1, var_pos):
                # print(bam_read.query_name, var_chr, var_pos)
                if check_read(
                    bam_read
                ):  # this is to help avoid weird situations where there is no cigar string.
                    RX = str(bam_read.query_name)
                    var_base = get_seq_base(bam_read, var_pos)
                    if str(var_base) != "D":
                        var_alt = str(rec.alts[0])
                        # block_num = get_block(block,var_count)
                        hap = grab_haplotype(var_ref, var_alt, var_base, GT)
                        try:
                            RX_library[RX].append(hap)
                            # RX_library[RX] = reduce_hap(RX_library[RX])
                        except:
                            RX_library[RX] = [hap]

    progressBar(num_lines, num_lines, end=True, text="split contig done")
    return RX_library


def vote_hap_majority(RX_dict, log):
    """
    This takes the read haplotype library and votes on which haplotype each read should belong in
    A read is determined to be a certain haplotype if:
      no conflicting haplotype is found
      have a 80% representation for one haplotype when more than 5 phased variants are found

    input:
      dictionary of reads and haplotype
    output:
      simplified phased dictionary of reads and haplotype
    """
    for RX in RX_dict:
        haplist = RX_dict[RX]
        len0 = haplist.count(0)
        len3 = haplist.count(3)
        len1 = haplist.count(1)
        len2 = haplist.count(2)
        length = len1 + len2
        # log.write("\t".join[str(RX) ,str(len0),str(len1),str(len2),str(len3), "\r"])
        if length <= 5:
            if len1 != 0 and len2 == 0:
                RX_dict[RX] = 1
            elif len1 == 0 and len2 != 0:
                RX_dict[RX] = 2
            else:
                RX_dict[RX] = 0
        else:
            if (len1 != 0 and len2 == 0) or (len1 / len2 > 0.8):
                RX_dict[RX] = 1
            elif len1 == 0 and len2 != 0 or (len2 / len1 > 0.8):
                RX_dict[RX] = 2
            else:
                RX_dict[RX] = 3
    return RX_dict


def write_summary(RX_dict, log):
    c = Counter(RX_dict.values())
    log.write("Haplotype Summary:\n")
    for i in range(0, 4):
        log.write("{} : {} \n".format(i, c[i]))


# def write_library(RX_library):
# print(RX_library)
# test = open("/mctp/share/users/chualec/hichip/vif/files/test.txt", "w")
# for i in RX_library:
#   strings = "\t".join[str(i), RX_library[i], "\r"]
#   test.write(strings)


def parse_hichip2(bam_in, RX_library, bam_out):
    """
    This functions takes the read-bam library and writes the bam read to different outputs files based on the haplotype

    input:
      bam file
      bam libary
      output files
    output:
      none
    """
    for bam_read in bam_in:
        if check_read(bam_read):
            RX = str(bam_read.query_name)
            try:
                hap = RX_library[RX]
            except:
                hap = 0
            bam_out[hap].write(bam_read)


def main(options):
    """
    This is the main function to read in information from the options given
    input:
      parsed options object
    output:
      none
    """
    welcome(options)
    print("")
    print("starting vif converter")
    # load the files given from the options using pysam
    startTime = datetime.now()
    # cleaning the file names and making the directory if needed
    out_dir_path = options.out_dir
    if out_dir_path[-1] != "/":
        out_dir_path += "/"

    if options.name is not None:
        output_name = out_dir_path + options.name
    else:
        output_name = out_dir_path + os.path.basename(options.reads)
        # print(os.path.basename(options.reads))

    if not os.path.isdir(out_dir_path):
        os.makedirs(out_dir_path)

    bam_in = pysam.AlignmentFile(options.reads)
    vcf = pysam.VariantFile(options.variants)
    log = open(output_name + ".log", "w")

    bam_out_0 = pysam.AlignmentFile(
        output_name + ".hap0.bam", template=bam_in, mode="wb"
    )
    bam_out_1 = pysam.AlignmentFile(
        output_name + ".hap1.bam", template=bam_in, mode="wb"
    )
    bam_out_2 = pysam.AlignmentFile(
        output_name + ".hap2.bam", template=bam_in, mode="wb"
    )
    bam_out_s = pysam.AlignmentFile(
        output_name + ".haps.bam", template=bam_in, mode="wb"
    )

    bam_outs = [bam_out_0, bam_out_1, bam_out_2, bam_out_s]

    RX_dict = parse_hichip(bam_in, vcf)
    print("finished parsing, starting voting...")
    RX_dict = vote_hap_majority(RX_dict, log)

    # Add a function here to write a output summary
    write_summary(RX_dict, log)

    bam_in.close()
    bam_in = pysam.AlignmentFile(options.reads)
    print("finished voting, starting parsing again...")
    parse_hichip2(bam_in, RX_dict, bam_outs)

    print("finished parsing...")
    print("runtime = ", (datetime.now() - startTime))
    print("")


if __name__ == "__main__":
    """
    This function parses options from the command line
    input:
      none
    output:
      parsed options
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--reads", type=str, required=True, help="bam file")
    parser.add_argument(
        "-v",
        "--variants",
        type=str,
        required=True,
        help="Phased variant calls to split reads",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        type=str,
        required=True,
        help="Directory to output split BAMs",
    )
    parser.add_argument("--name", type=str, required=False, help="Optional output name")
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
