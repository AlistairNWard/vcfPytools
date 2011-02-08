#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

#Calculate the unique fraction.
def uniqueVcf(v1, v2, outputFile):
  success1 = v1.getRecord()
  success2 = v2.getRecord()
  currentReferenceSequence = v1.referenceSequence

# If the end of the first vcf file is reached, it can
# have no more unique records, so terminate.
  while success1:

# If the end of the second file is reached, output all
# of the remaining records from the first vcf file as
# they must all be unique.
    if not success2:
      outputFile.write(v1.record)
      success1 = v1.getRecord()

    if v1.referenceSequence == v2.referenceSequence and v1.referenceSequence == currentReferenceSequence:
      if v1.position == v2.position:
        success1 = v1.getRecord()
        success2 = v2.getRecord()
      elif v2.position > v1.position:
        success1 = v1.parseVcf(v2.referenceSequence, v2.position, True, outputFile)
      elif v1.position > v2.position: success2 = v2.parseVcf(v1.referenceSequence, v1.position, False, None)
    else:
      if v1.referenceSequence == currentReferenceSequence: success1 = v1.parseVcf(v2.referenceSequence, v2.position, True, outputFile)
      elif v2.referenceSequence == currentReferenceSequence: success2 = v2.parseVcf(v1.referenceSequence, v1.position, False, None)

# If the last record for a reference sequence is the same for both vcf
# files, they will both have referenceSequences different from the
# current reference sequence.  Change the reference sequence to reflect
# this and proceed.
      else:
        if v1.referenceSequence != v2.referenceSequence:
          print >> sys.stderr, "ERROR: Reference sequences for both files are unexpectedly different."
          print >> sys.stderr, "Check that both files contain records for the following reference sequences:"
          print >> sys.stderr, "\t", v1.referenceSequence, " and ", v2.referenceSequence
          exit(1)
      currentReferenceSequence = v1.referenceSequence

if __name__ == "__main__":
  main()

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py unique [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="append", type="string",
                    dest="vcfFiles", help="input vcf files")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.
  if options.vcfFiles == None:
    parser.print_help()
    print >> sys.stderr, "\nTwo input vcf files (-i) are required for performing calculation of unique fraction."
    exit(1)
  elif len(options.vcfFiles) != 2:
    print >> sys.stderr, "Two input vcf files are required for performing calculation of unique fraction."

# Set the output file to stdout if no output file was specified.
  outputFile, writeOut = setOutput(options.output)

  v1 = vcf() # Define vcf object.
  v2 = vcf() # Define vcf object.

# Open the vcf files.
  v1.openVcf(options.vcfFiles[0])
  v2.openVcf(options.vcfFiles[1])

# Read in the header information.
  v1.parseHeader(options.vcfFiles[0], writeOut)
  v2.parseHeader(options.vcfFiles[1], writeOut)

# Make it clear to the user which unique fraction is being
# calculated.  It is always the first vcf file inputted.
  print >> sys.stderr, "\nGenerating records unique to:", options.vcfFiles[0]

# Check that the header for the two files contain the same samples.
  if v1.samplesList != v2.samplesList:
    print >> sys.stderr, "vcf files contain different samples (or sample order)."
    exit(1)
  else:
    taskDescriptor = "##vcfPytools=generate variants unique to " + options.vcfFiles[0] + " when compared to " + options.vcfFiles[1]
    if v1.hasHeader: writeHeader(outputFile, v1, False, taskDescriptor) # tools.py
    else: writeHeader(outputFile, v2, False, taskDescriptor) # tools.py

# Calculate the unique fraction.
  uniqueVcf(v1, v2, outputFile)

# Check that the input files had the same list of reference sequences.
# If not, it is possible that there were some problems.
  checkReferenceSequenceLists(v1.referenceSequenceList, v2.referenceSequenceList) # tools.py

# Close the vcf files.
  v1.closeVcf(options.vcfFiles[0])
  v2.closeVcf(options.vcfFiles[1])

# Terminate the program cleanly.
  return 0
