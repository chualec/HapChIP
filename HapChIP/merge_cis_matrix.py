import argparse
import os
import operator
import sys
from itertools import combinations
import sys, getopt
from datetime import datetime
import pysam

# the cigar seq function expands the sequence information to account for soft clipping, insertions, and deletions

import pandas as pd
import math


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


def parse_hichip(infiles, outfiles):
    pdmatrix = pd.DataFrame(
        0, columns=["HC-9", "HC-10", "HC-12", "HC-13", "HC-21", "HC-25"], index=["Sum"]
    )
    standard = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
    ]
    print(pdmatrix)
    for files in infiles:
        for row in files:
            colID = row.split()[1]
            pdmatrix.loc["Sum", colID] += 1
            # if row.split()[0].split()[0] in standard and row.split()[2] in standard:

            rowID = row.split()[0]
            if rowID in pdmatrix.index:
                pdmatrix.loc[rowID, colID] += 1
            else:
                pdmatrix.loc[rowID, colID] = 1
    dropped_df = pdmatrix.loc[((pdmatrix > 50)).any(1)]
    dropped_df.to_csv(outfiles)


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
    # out_dir_path = options.out_dir
    # if out_dir_path[-1] != "/":
    #   out_dir_path += "/"

    # make output name/location
    # if options.name is not None:
    #   output_name = out_dir_path + options.name
    # else:
    #   output_name = out_dir_path + os.path.basename(options.bam)
    #
    # #make output directory
    # if not os.path.isdir(out_dir_path):
    #    os.makedirs(out_dir_path)

    # load in input files
    hc9 = open(
        "/mctp/share/users/chualec/hichip/vif/files/HC-9/cis_trans_split/HC-9.pairs.trans.round.txt",
        "r",
    )
    hc10 = open(
        "/mctp/share/users/chualec/hichip/vif/files/HC-10/cis_trans_split/HC-10.pairs.trans.round.txt",
        "r",
    )
    hc12 = open(
        "/mctp/share/users/chualec/hichip/vif/files/HC-12/cis_trans_split/HC-12.pairs.trans.round.txt",
        "r",
    )
    hc13 = open(
        "/mctp/share/users/chualec/hichip/vif/files/HC-13/cis_trans_split/HC-13.pairs.trans.round.txt",
        "r",
    )
    hc21 = open(
        "/mctp/share/users/chualec/hichip/vif/files/HC-21/cis_trans_split/HC-21.pairs.trans.round.txt",
        "r",
    )
    hc25 = open(
        "/mctp/share/users/chualec/hichip/vif/files/HC-25/cis_trans_split/HC-25.pairs.trans.round.txt",
        "r",
    )
    infiles = [hc9, hc10, hc12, hc13, hc21, hc25]

    outfile = (
        "/mctp/share/users/chualec/hichip/vif/files/trans.all.bins.filtered.matrix"
    )
    parse_hichip(infiles, outfile)

    print("runtime = ", (datetime.now() - startTime))
    print("")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", type=str, required=False, help="bam file")
    parser.add_argument(
        "-o",
        "--out_dir",
        type=str,
        required=False,
        help="Directory to output split BAMs",
    )
    parser.add_argument("--name", type=str, required=False, help="Optional output name")
    parsed, unparsed = parser.parse_known_args()
    main(parsed)

# python3 ../vif_converter_2.py -s test.namesort.bam -v ../HC-12v2/HC-12v2.phased.phased.VCF.gz -o . --name vif2
