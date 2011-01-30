#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

if __name__ == "__main__":
  main()

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py test [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that a single  vcf file is given.
  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nAt least one vcf file (--in, -i) is required for performing intersection."
    exit(1)

  outputFile, writeOut = setOutput(options.output)

  v = vcf() # Define vcf object.

# Open the vcf files.
  v.openVcf(options.vcfFile)

# Read in the header information.
  v.parseHeader(options.vcfFile, False)

# Perform testing.
  while v.getRecord():
    continue

# Close the vcf files.
  v.closeVcf(options.vcfFile)

# End the program.
  return 0
