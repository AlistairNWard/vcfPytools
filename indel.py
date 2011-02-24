#!/usr/bin/python

import sys
import os.path
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

if __name__ == "__main__": main()

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py indel [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output statistics file")

  (options, args) = parser.parse_args()

# Check that a single  vcf file is given.
  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nAt least one vcf file (--in, -i) is required for performing intersection."
    exit(1)

# Set the output file to stdout if no output file was specified.
  outputFile, writeOut = setOutput(options.output) # tools.py
  v = vcf() # Define vcf object.

# Open the file.
  v.openVcf(options.vcfFile)
  v.processInfo = True

# Read in the header information.
  v.parseHeader(options.vcfFile, writeOut)

# Initialise some variables.
  insertionDistribution = {}
  deletionDistribution = {}
  noInsertions = 0
  noDeletions = 0
  totalInsertionLength = 0
  totalDeletionLength = 0
  noMNPs = 0
  noSNPs = 0

# Parse the file.
  while v.getRecord():
    if v.infoTags.has_key("INS"):
      tagNumber, tagType, tagValue = v.getInfo("INS")
      insertionDistribution[int(tagValue[0])] = insertionDistribution.get(int(tagValue[0]), 0) + 1
      noInsertions += 1
      totalInsertionLength += int(tagValue[0])
    elif v.infoTags.has_key("DEL"):
      tagNumber, tagType, tagValue = v.getInfo("DEL")
      deletionDistribution[int(tagValue[0])] = deletionDistribution.get(int(tagValue[0]), 0) + 1
      noDeletions += 1
      totalDeletionLength += int(tagValue[0])
    elif v.infoTags.has_key("MNP"):
      noMNPs += 1
    else:
      noSNPs += 1

# Print out the distributions.
  allKeys = {}
  for key in insertionDistribution.iterkeys(): allKeys[int(key)] = 1
  for key in deletionDistribution.iterkeys(): allKeys[int(key)] = 1
  allKeysList = allKeys.keys()
  allSortedKeys = sorted(list(allKeysList))
  print >> outputFile, "Indel statistics for vcf file:\t", options.vcfFile, "\n"
  print >> outputFile, "Number of SNPs:\t", noSNPs
  print >> outputFile, "Number of MNPs:\t", noMNPs
  print >> outputFile, "Number of insertions: \t", noInsertions, "\t(length: ", totalInsertionLength, "bp)"
  print >> outputFile, "Number of deletions: \t", noDeletions, "\t(length: ", totalDeletionLength, "bp)"
  print >> outputFile, "Total:\t\t\t", noInsertions + noDeletions, "\t(length: ", totalInsertionLength + totalDeletionLength, "bp)\n"
  print >> outputFile, "ins/del number ratio:\t", float(noInsertions)/float(noDeletions)
  print >> outputFile, "ins/del length ratio:\t", float(totalInsertionLength)/float(totalDeletionLength), "\n"
  print >> outputFile, "size\tins\tdel\tins/del ratio"
  for key in allSortedKeys:
    insertions = insertionDistribution.get(key, 0)
    deletions = deletionDistribution.get(key, 0)
    ratio = float(insertions)/float(deletions) if float(deletions) != 0 else 0
    print >> outputFile, key, "\t", insertions, "\t", deletions, "\t", ratio
