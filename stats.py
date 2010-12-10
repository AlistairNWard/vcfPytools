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

  def processStats(self, ref, alt, multiAllelic, chromosome):

# Increment the number of SNPs.

    self.totalSnps[chromosome] = self.totalSnps.get(chromosome, 0) + 1

# Determine if the SNP is a transition, transversion or multi-allelic.

    if multiAllelic == True:
      self.multiAllelic[chromosome] = self.multiAllelic.get(chromosome, 0) + 1
    else:
      if ref < alt:
        alleles = ref + alt
      else:
        alleles = alt + ref
    
      if alleles.lower() == "ag" or alleles.lower() == "ct":
        self.transitions[chromosome] = self.transitions.get(chromosome, 0) + 1
      elif alleles.lower() == "ac" or alleles.lower() == "at" or alleles.lower() == "cg" or alleles.lower() == "gt":
        self.transversions[chromosome] = self.transversions.get(chromosome, 0) + 1

  def printStats(self):
    for ref in sorted(self.referenceSequences):

      totalSnps = self.totalSnps[ref] if self.totalSnps.has_key(ref) else 0
      transitions = self.transitions[ref] if self.transitions.has_key(ref) else 0
      transversions = self.transversions[ref] if self.transversions.has_key(ref) else 0
      multiAllelic = self.multiAllelic[ref] if self.multiAllelic.has_key(ref) else 0

      print ref, totalSnps, transitions, transversions, multiAllelic

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

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.

  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (-i) is required."
    exit(1)

# Set the output file to stdout if no output file was specified.

  if options.output == None:
    outputFile = sys.stdout
  else:
    outputFile = open(options.output, 'w')

  v = vcf() # Define vcf object.
  stats = statistics() # Define statistics object

# Open the file.

  v.openVcf(options.vcfFile)

# Read in the header information.

  v.parseHeader(options.vcfFile)

# Read through all the entries.

  for line in v.filehandle:
    v.readRecord(line)
    v.resolveInfo(v.info)

    stats.referenceSequences[v.chromosome] = stats.referenceSequences.get(v.chromosome, 0) + 1
    stats.processStats(v.ref, v.alt, v.multiAllelic,v.chromosome)

# Close the file.

  v.closeVcf(options.vcfFile)

# Output the statistics

  stats.printStats()
