import unittest


import argparse
import os
import operator
import sys
from itertools import combinations
import sys, getopt
from datetime import datetime
import pysam


def test_readvcf():
  assert pysam.VariantFile("../unittest/unittest.vcf.gz"), "No VCF file vcf"

def test_readbam():
  assert pysam.AlignmentFile("../unittest/unittest.bam"), "No Bam file detected"

      
if __name__ == "__main__":
  test_readvcf()
  test_readbam()
  print("everything passed")
     
  
