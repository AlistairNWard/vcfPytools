#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

if __name__ == "__main__":
  main()

def main():

# Parse the command line options

  usage = "Usage: vcfTools.py union [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="append", type="string",
                    dest="vcfFiles", help="input vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.

  if options.vcfFiles == None:
    parser.print_help()
    print >> sys.stderr, "\nTwo input vcf files (-i) are required for performing union."
    exit(1)
  elif len(options.vcfFiles) != 2:
    print >> sys.stderr, "Tow input vcf files are required for performing union."

# Set the output file to stdout if no output file was specified.

  if options.output == None:
    outputFile = sys.stdout
  else:
    outputFile = open(options.output, 'w')

  v1 = vcf() # Define vcf object.
  v2 = vcf() # Define vcf object.

# Open the vcf files.

  v1.openVcf(options.vcfFiles[0])
  v2.openVcf(options.vcfFiles[1])

# Read in the header information.

  v1.parseHeader(options.vcfFiles[0])
  v2.parseHeader(options.vcfFiles[1])

# Check that the header for the two files contain the same samples.

  if v1.samplesList != v2.samplesList:
    print "vcf files contain different samples (or sample order)."
    exit(1)
  else:
    print "vcf files contain the same samples."
    print "Continuing with union."

# Close the vcf files.

  v1.closeVcf(options.vcfFiles[0])
  v2.closeVcf(options.vcfFiles[1])
