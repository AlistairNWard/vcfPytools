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
  outputFile, writeOut = setOutput(options.output) # tools.py

# Open the output file for writing, or set the output file to
# stdout if an output file wasn't defined.
  if options.output == None:
    outputFile = sys.stdout
  else:
    outputFile = open(options.output,'w')

  for index, vcfFile in enumerate(options.vcfFiles):
    v = vcf() # Define vcf object
    v.openVcf(vcfFile) # Open the vcf file
    v.parseHeader(vcfFile, True)

# Store the header from the first vcf file.  The samplesList from 
# all other vcf files being merged will be checked against this.
# Also, print out the header.
    if index == 0:
      samples = v.samplesList
      outputFile.write( v.headerText ) if v.headerText != "" else None
      outputFile.write( v.headerInfoText ) if v.headerInfoText != "" else None
      outputFile.write( v.headerFormatText ) if v.headerFormatText != "" else None
      outputFile.write( v.headerTitles )
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
