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

  usage = "Usage: vcfTools.py validate [options]"
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

  if options.output == None:
    outputFile = sys.stdout
    writeOut = False
  else:
    outputFile = open(options.output, 'w')
    writeOut = True

  v = vcf() # Define vcf object.

# Open the file.

  v.openVcf(options.vcfFile)

# Read in the header information.

  v.parseHeader(options.vcfFile, writeOut, True)
  v.processInfo = True
  v.processGenotypes = True

# Read through all the entries.

  for line in v.filehandle:
    v.getRecord(line)

# For each field in the info string, check that there exists a
# line in the header explaining the field and that the field
# contains the correct number and type of entries.

    for tag in v.infoTags:
      number, value = v.getInfo(tag)

# Retrieve all genotype information.  Any errors in these fields
# will be found.

    for sample in v.samplesList:
      for tag in v.genotypeFormats:
        number, value = v.getGenotypeInfo(sample, tag)

# Close the file.

  v.closeVcf(options.vcfFile)

# Print that no error were found (if there had been errors, the
# program would have terminated).

  print "No errors found with vcf file"
