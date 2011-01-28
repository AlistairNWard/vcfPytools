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
  usage = "Usage: vcfPytools.py test [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file")

  (options, args) = parser.parse_args()

# Check that a single  vcf file is given.
  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nAt least one vcf file (--in, -i) is required for performing intersection."
    exit(1)

  v = vcf() # Define vcf object.

# Open the vcf files.
  v.openVcf(options.vcfFile)

# Read in the header information.
  v.parseHeader(options.vcfFile, False, True)

# Perform testing.
  success = 0
  while success == 0:
    success = v.getRecord()

# Close the vcf files.
  v.closeVcf(options.vcfFile)

# End the program.
  return 0
