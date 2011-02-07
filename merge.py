#!/usr/bin/python

import sys
import os.path
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py merge [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="append", type="string",
                    dest="vcfFiles", help="input vcf files")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.
  if not options.vcfFiles or len(options.vcfFiles) <= 0:
    print >> sys.stderr, "Multiple vcf files must be included."
    exit(1)

# Open the output file for writing, or set the output file to
# stdout if an output file wasn't defined.
  outputFile, writeOut = setOutput(options.output) # tools.py

  for index, vcfFile in enumerate(options.vcfFiles):
    v = vcf() # Define vcf object
    v.openVcf(vcfFile) # Open the vcf file
    v.parseHeader(vcfFile, True)

# Store the header from the first vcf file.  The samplesList from 
# all other vcf files being merged will be checked against this.
# Also, print out the header.
    if index == 0:
      samples = list(v.samplesList)
      taskDescriptor = "##vcfPytools=merge "
      for file in options.vcfFiles: taskDescriptor += file + ", "
      taskDescriptor = taskDescriptor.rstrip(", ")
      writeHeader(outputFile, v, False, taskDescriptor) # tools.py
    else:
      if v.samplesList != samples:
        print >> sys.stderr, "WARNING: Different samples in file: ", vcfFile

# print out the records.
    while v.getRecord():
      outputFile.write(v.record)

    v.closeVcf(vcfFile) # Close the vcf file

# Merge is complete,  Close the output file.
  outputFile.close()

# Terminate the program cleanly.
  return 0
