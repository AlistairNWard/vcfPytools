#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

import intersect
from intersect import setVcfPriority
from intersect import writeVcfRecord

if __name__ == "__main__":
  main()

# Calculate the union of the two vcf files.
def unionVcf(v1, v2, priority, outputFile):
  success1 = v1.getRecord()
  success2 = v2.getRecord()
  currentReferenceSequence = v1.referenceSequence

# Finish when the end of either file has been reached.
  while success1 or success2:

# If the end of the first vcf file is reached, write out the
# remaining records from the second vcf file.
    if not success1:
      outputFile.write(v2.record)
      success2 = v2.getRecord()

# If the end of the second vcf file is reached, write out the
# remaining records from the first vcf file.
    if not success2:
      outputFile.write(v1.record)
      success1 = v1.getRecord()

    if v1.referenceSequence == v2.referenceSequence:
      if v1.position == v2.position:
        writeVcfRecord(priority, v1, v2, outputFile)
        success1 = v1.getRecord()
        success2 = v2.getRecord()
      elif v2.position > v1.position:
        success1 = v1.parseVcf(v2.referenceSequence, v2.position, True, outputFile)
      elif v1.position > v2.position: success2 = v2.parseVcf(v1.referenceSequence, v1.position, True, outputFile)
    else:
      if v1.referenceSequence == currentReferenceSequence: success1 = v1.parseVcf(v2.referenceSequence, v2.position, True, outputFile)
      elif v2.referenceSequence == currentReferenceSequence: success2 = v2.parseVcf(v1.referenceSequence, v1.position, True, outputFile)
      currentReferenceSequence == v1.referenceSequence

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py union [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="append", type="string",
                    dest="vcfFiles", help="input vcf files")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")
  parser.add_option("-p", "--priority-file",
                    action="store", type="string",
                    dest="priorityFile", help="output record from this file")

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.
  if options.vcfFiles == None:
    parser.print_help()
    print >> sys.stderr, "\nTwo input vcf files (-i) are required for performing union."
    exit(1)
  elif len(options.vcfFiles) != 2:
    print >> sys.stderr, "Two input vcf files are required for performing union."

# Set the output file to stdout if no output file was specified.
  outputFile, writeOut = setOutput(options.output)

# If no priority is given to either file (from the -p command line
# option), set priorityQuality to True.  In this case, the record
# written to the output file will be that with the higest quality.
# If a priority is given, check that the file is one of the input
# vcf files.
  priority = setVcfPriority(options.priorityFile, options.vcfFiles) # intersect.py

  v1 = vcf() # Define vcf object.
  v2 = vcf() # Define vcf object.

# Open the vcf files.
  v1.openVcf(options.vcfFiles[0])
  v2.openVcf(options.vcfFiles[1])

# Read in the header information.
  v1.parseHeader(options.vcfFiles[0], writeOut)
  v2.parseHeader(options.vcfFiles[1], writeOut)

# Check that the header for the two files contain the same samples.
  if v1.samplesList != v2.samplesList:
    print >> sys.stderr, "vcf files contain different samples (or sample order)."
    exit(1)
  else:
    writeHeader(outputFile, v1, False) # tools.py

# Calculate the union.
  unionVcf(v1, v2, priority, outputFile)

# Close the vcf files.
  v1.closeVcf(options.vcfFiles[0])
  v2.closeVcf(options.vcfFiles[1])

# Terminate the program cleanly.
  return 0
