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
  usage = "Usage: vcfPytools.py validate [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file (stdin for piped vcf)")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output validation file")

  (options, args) = parser.parse_args()

# Check that a vcf file is given.
  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (-i) is required."
    exit(1)

# Set the output file to stdout if no output file was specified.
  outputFile, writeOut = setOutput(options.output)

  v = vcf() # Define vcf object.

# Open the file.
  v.openVcf(options.vcfFile)

# Read in the header information.
  v.parseHeader(options.vcfFile, writeOut, True)
  v.processInfo = True
  v.processGenotypes = True
  parsedReferenceSequences = {}
  sortError = False
  duplicateError = False

# Read through all the entries.
  previousReference = ""
  previousPosition = 0
  success = 0
  while success == 0:
    success = v.getRecord()

# Check that the current position isn't before the previous
# position (if in the same reference sequence) or that the
# reference sequence is different from the previous record, but
# has been seen before.  In either of these cases, the vcf file
# is not correctly sorted.
    if previousReference != v.referenceSequence:
      if parsedReferenceSequences.has_key(v.referenceSequence): sortError = True
    elif v.position < previousPosition: sortError = True
    elif v.position == previousPosition: duplicateError = True
      
    previousReference = v.referenceSequence
    previousPosition = v.position
    parsedReferenceSequences[ v.referenceSequence ] = True

# For each field in the info string, check that there exists a
# line in the header explaining the field and that the field
# contains the correct number and type of entries.
    for tag in v.infoTags:
      number, tagType, value = v.getInfo(tag)

# Retrieve all genotype information.  Any errors in these fields
# will be found.
    for sample in v.samplesList:
      for tag in v.genotypeFormats:
        number, value = v.getGenotypeInfo(sample, tag)

# Close the file.
  v.closeVcf(options.vcfFile)

# Print that no error were found (if there had been errors, the
# program would have terminated).
  if sortError == True: print >> sys.stderr, "vcf file is not sorted."
  elif duplicateError == True: print >> sys.stderr, "vcf file contains duplicate entries."
  else: print >> sys.stderr, "No errors found with vcf file"

# Terminate the program cleanly.
  return 0
