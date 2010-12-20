#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

class statistics:
  def __init__(self):
    self.referenceSequences = {}
    self.totalSnps = {}
    self.transitions = {}
    self.transversions = {}
    self.multiAllelic = {}
    self.distributions = {}

  def processStats(self, ref, alt, multiAllelic, referenceSequence):

# Increment the number of SNPs.

    self.totalSnps[referenceSequence] = self.totalSnps.get(referenceSequence, 0) + 1

# Determine if the SNP is a transition, transversion or multi-allelic.

    if multiAllelic == True:
      self.multiAllelic[referenceSequence] = self.multiAllelic.get(referenceSequence, 0) + 1
    else:
      if ref < alt:
        alleles = ref + alt
      else:
        alleles = alt + ref
    
      if alleles.lower() == "ag" or alleles.lower() == "ct":
        self.transitions[referenceSequence] = self.transitions.get(referenceSequence, 0) + 1
      elif alleles.lower() == "ac" or alleles.lower() == "at" or alleles.lower() == "cg" or alleles.lower() == "gt":
        self.transversions[referenceSequence] = self.transversions.get(referenceSequence, 0) + 1

  def printStats(self):
    for ref in sorted(self.referenceSequences):

      totalSnps = self.totalSnps[ref] if self.totalSnps.has_key(ref) else 0
      transitions = self.transitions[ref] if self.transitions.has_key(ref) else 0
      transversions = self.transversions[ref] if self.transversions.has_key(ref) else 0
      multiAllelic = self.multiAllelic[ref] if self.multiAllelic.has_key(ref) else 0

      print ref, totalSnps, transitions, transversions, multiAllelic

# Initialise data structures for the distributions.

  def initialiseDistributions(self):
    print "INITIALISE"

if __name__ == "__main__":
  main()

def main():

# Parse the command line options

  usage = "Usage: vcfTools.py stats [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file (stdin for piped vcf)")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")
  parser.add_option("-d", "--distributions",
                    action="append", type="string",
                    dest="distributions", help="plot distributions of variables in the info fields" + \
                    " (all includes all info fields in header)")

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

  stats = statistics() # Define statistics object

  v.parseHeader(options.vcfFile, writeOut, True)

# If distributions for all the info fields listed in the header are
# requested, populate options.distributions with these values.

  if options.distributions:
    v.processInfo = True
    if options.distributions[0].lower() == "all" and len(options.distributions) == 1:
      for tag in v.infoHeaderTags:
        options.distributions.append(tag)
      del(options.distributions[0])
    else:
      for tag in options.distributions:
        if tag.lower() == "all":
          print "If distributions for all info fields are required, include -d all only"
          exit(0)

# Check that the requested info fields exist in the vcf file and
# initialise statistics dictionaries.

  if options.distributions:
    for tag in options.distributions:
      v.checkInfoFields(tag)
      stats.initialiseDistributions()

# Read through all the entries.

  for line in v.filehandle:
    v.getRecord(line)
    if options.distributions:
      for tag in options.distributions:
        tagNumber, tagValue = v.getInfo(tag)

# Close the file.

  v.closeVcf(options.vcfFile)
