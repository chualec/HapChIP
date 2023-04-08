import argparse
import os
import operator
import sys
from itertools import combinations
import sys, getopt
from datetime import datetime
import pysam

# the cigar seq function expands the sequence information to account for soft clipping, insertions, and deletions


def progressBar(current, total, end=False, text="NA"):
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


def cigarseq(bam_read):
    read_seq = bam_read.query_alignment_sequence
    quality_seq = bam_read.query_alignment_qualities
    cig_seq = ""
    qual_seq = []
    for cigs in bam_read.cigar:
        if cigs[0] == 0:  # these indicate match
            cig_seq += read_seq[0 : cigs[1]]
            read_seq = read_seq[cigs[1] :]
            qual_seq.extend(quality_seq[0 : cigs[1]])
            quality_seq = quality_seq[cigs[1] :]
        elif cigs[0] == 1:  # these indicate insertion
            read_seq = read_seq[cigs[1] :]
            quality_seq = quality_seq[cigs[1] :]
        elif cigs[0] == 2:  # these indicate deletion
            cig_seq += "D" * cigs[1]
            qual_seq.extend(["D"] * cigs[1])
        elif cigs[0] == 3:  # these indicate skipped
            cig_seq += "N" * cigs[1]
            qual_seq.extend(["N"] * cigs[1])
        elif cigs[0] == 4:  # these indicate soft_clip
            cig_seq += "S" * cigs[1]
            qual_seq.extend(["S"] * cigs[1])
        elif cigs[0] == 5:  # these indicate hard_clip
            pass
    return (cig_seq, qual_seq)


def get_seq_base(bam_read, var_pos):
    cigar_seq, qual_seq = cigarseq(bam_read)
    ref_base = var_pos - bam_read.pos - 1
    return (cigar_seq[ref_base], str(qual_seq[ref_base]))


# def grab_genotype(rec):
#   try:
#     return(dict(rec.bamples[0])['PGT'])
#   except:
#     return("0/1")
def grab_genotype(rec):
    return str(rec).split()[-1].split(":")[0]


def check_var(rec):
    checkok = True
    alleles = rec.alleles
    # check if indel
    if len(alleles[0]) > 1 or len(alleles[1]) > 1:
        checkok = False
    # check if het or phased
    # elif grab_genotype(rec) not in ["0/1", "0|1", "1|0"]:
    # checkok = False
    elif len(rec.chrom) > 6:
        checkok = False
    return checkok


def grab_haplotype(var_ref, var_alt, var_base, GT):
    hap = 0
    if var_base == var_ref:
        if GT == "0|1":
            hap = 1
        elif GT == "1|0":
            hap = 2
    elif var_base == var_alt:
        if GT == "0|1":
            hap = 2
        elif GT == "1|0":
            hap = 1
    return hap


def get_block(block, var_num):
    try:
        blocknum = block[str(var_num)]
    except:
        blocknum = 0
    return blocknum


def reduce_hap(hap_list):
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


def get_bam_hap(bam_read, vcf):
    chrom = str(bam_read.reference_name)
    refstart = int(bam_read.reference_start + 1)
    refend = int(bam_read.reference_end)
    hap_list = []
    for vcf_i in vcf.fetch(chrom, refstart, refend):
        if check_var(vcf_i):
            var_pos = vcf_i.pos
            var_ref = vcf_i.ref
            var_alt = str(vcf_i.alts[0])
            var_base, baseq = get_seq_base(bam_read, var_pos)
            GT = grab_genotype(vcf_i)
            vcf_hap = grab_haplotype(var_ref, var_alt, var_base, GT)
            hap_list.append(vcf_hap)
    return reduce_hap(hap_list)


def bam_list_decomb(bam_list, vcf):
    if len(bam_list) <= 1:
        return 0
    hap_list = []
    for bam_read in bam_list:
        bam_hap = get_bam_hap(bam_read, vcf)
        hap_list.append(bam_hap)
    return reduce_hap(hap_list)


def bam_list_write(bam_list, vcf, bam_outs):
    if len(bam_list) > 1:
        haplist = []
        ids = bam_list[0].query_name
        for bam_read in bam_list:
            hap = get_bam_hap(bam_read, vcf)
            haplist.append(
                "\t".join(
                    [
                        bam_read.reference_name,
                        str(bam_read.reference_start),
                        str(bam_read.reference_end),
                        str(bam_read.mapping_quality),
                        str(hap),
                    ]
                )
            )
        comb = combinations(range(1, len(bam_list) + 1), 2)
        for i in comb:
            bam_outs[4].write(
                ids + "\t" + haplist[i[0] - 1] + "\t" + haplist[i[1] - 1] + "\n"
            )


def check_read(bam_read):
    if int(bam_read.mapping_quality) == 0:
        return False
    for cig in bam_read.cigar:
        if cig[0] == 0 and cig[1] > 1:
            return True
    return False


def parse_hichip(bamfile, bam_outs):
    # for each read that overlaps the variant position
    bam_library = {}
    for bam_read in bamfile:
        if check_read(bam_read) and bam_read.is_supplementary == False:
            if str(bam_read.query_name) in bam_library:
                read2 = bam_library[bam_read.query_name]
                chr1 = str(bam_read.reference_name)
                start1 = str(bam_read.reference_start + 1)
                end1 = str(bam_read.reference_end)

                chr2 = str(read2.reference_name)
                start2 = str(read2.reference_start + 1)
                end2 = str(read2.reference_end)
                if abs(int(start2) - int(start1)) > 1000 or chr1 != chr2:
                    bam_outs[3].write(
                        "\t".join(
                            [
                                str(bam_read.query_name),
                                chr1,
                                start1,
                                end1,
                                chr2,
                                start2,
                                end2,
                                "\n",
                            ]
                        )
                    )
                if chr1 != chr2:
                    bam_outs[2].write(bam_read)
                    bam_outs[2].write(read2)
                elif abs(int(start1) - int(start2)) > 1000:
                    bam_outs[1].write(bam_read)
                    bam_outs[1].write(read2)
                else:
                    bam_outs[0].write(bam_read)
                    bam_outs[0].write(read2)
            else:
                bam_library[str(bam_read.query_name)] = bam_read
        # else:
        #   bam_outs[3].write(bam_read)

    # bam_outs[0].write(bam_read)


def welcome(options):
    print("")
    print("*************")
    print("thank you for using the vif converter")
    print("")
    if options.bam is not None:
        print("Detected bam file path: %s" % (options.bam))
    if options.out_dir is not None:
        print("Output directory: %s" % (options.out_dir))
    else:
        print("No phased block file detected")
    if options.name is not None:
        print("Detected output file name: %s" % (options.name))


def main(options):
    welcome(options)
    print("")
    print("starting vif converter")
    startTime = datetime.now()

    # get output path
    out_dir_path = options.out_dir
    if out_dir_path[-1] != "/":
        out_dir_path += "/"

    # make output name/location
    if options.name is not None:
        output_name = out_dir_path + options.name
    else:
        output_name = out_dir_path + os.path.basename(options.bam)

    # make output directory
    if not os.path.isdir(out_dir_path):
        os.makedirs(out_dir_path)

    # load in input files
    bam_in = pysam.AlignmentFile(options.bam, "r")

    # make output files
    bam_out_0 = pysam.AlignmentFile(
        output_name + ".cis_short.bam", template=bam_in, mode="wb"
    )
    bam_out_1 = pysam.AlignmentFile(
        output_name + ".cis_dist.bam", template=bam_in, mode="wb"
    )
    bam_out_2 = pysam.AlignmentFile(
        output_name + ".trans.bam", template=bam_in, mode="wb"
    )
    bam_out_s = pysam.AlignmentFile(
        output_name + ".bad.bam", template=bam_in, mode="wb"
    )
    bam_out_t = open(output_name + ".pairs.txt", "w")
    # record = open(output_name+".decomb4", "w")
    bam_outs = [bam_out_0, bam_out_1, bam_out_2, bam_out_t]

    parse_hichip(bam_in, bam_outs)

    print("runtime = ", (datetime.now() - startTime))
    print("")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", type=str, required=True, help="bam file")
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

# python3 ../vif_converter_2.py -s test.namesort.bam -v ../HC-12v2/HC-12v2.phased.phased.VCF.gz -o . --name vif2
