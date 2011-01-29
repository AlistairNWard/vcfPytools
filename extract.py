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
  usage = "Usage: vcfPytools.py extract [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file (stdin for piped vcf)")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output validation file")
  parser.add_option("-s", "--reference-sequence",
                    action="store", type="string",
                    dest="referenceSequence", help="extract records from this reference sequence")
  parser.add_option("-r", "--region",
                    action="store", type="string",
                    dest="region", help="extract records from this region")

  (options, args) = parser.parse_args()

# Check that a vcf file is given.
  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (--in, -i) is required."
    exit(1)

# Check that either a reference sequence or region is specified,
# but not both.
  if not options.referenceSequence and not options.region:
    parser.print_help()
    print >> sys.stderr, "\nA region (--region, -r) or reference sequence (--reference-sequence, -s) must be supplied."
    exit(1)
  if options.referenceSequence and options.region:
    parser.print_help()
    print >> sys.stderr, "\nEither a region (--region, -r) or reference sequence (--reference-sequence, -s) must be supplied, but not both."
    exit(1)

# If a region was supplied, check the format.
  if options.region:
    if options.region.find(":") == -1 or options.region.find("..") == -1:
      print >> sys.stderr, "\nIncorrect format for region string.  Required: ref:start..end."
      exit(1)
    regionList = options.region.split(":",1)
    referenceSequence = regionList[0]
    try: start = int(regionList[1].split("..")[0])
    except ValueError:
      print >> sys.stderr, "region start coordinate is not an integer"
      exit(1)
    try: end = int(regionList[1].split("..")[1])
    except ValueError:
      print >> sys.stderr, "region end coordinate is not an integer"
      exit(1)

# Set the output file to stdout if no output file was specified.
  outputFile, writeOut = setOutput(options.output)

  v = vcf() # Define vcf object.

# Open the file.
  v.openVcf(options.vcfFile)

# Read in the header information.
  v.parseHeader(options.vcfFile, writeOut, True)

# Read through all the entries and write out records in the correct
# reference sequence.
  if options.referenceSequence:
    while v.getRecord() == 0:
      if v.referenceSequence == options.referenceSequence: outputFile.write(v.record)

# Read through all the entries and write out records in the correct
# region.
  elif options.region:
    while v.getRecord() == 0:
      if v.referenceSequence == referenceSequence and v.position >= start and v.position <= end: outputFile.write(v.record)

# Close the file.
  v.closeVcf(options.vcfFile)

# Terminate the program cleanly.
  return 0
