#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

class statistics:
  def __init__(self):
    self.referenceSequences = {}
    self.novelSnps = {}
    self.knownSnps = {}
    self.transitions = {}
    self.transversions = {}
    self.multiAllelic = {}
    self.distributions = {}

  def processGeneralStats(self, rsid, ref, alt, multiAllelic, referenceSequence):
    self.referenceSequences[referenceSequence] = True

# Increment the number of SNPs.

    if rsid == ".":
      self.novelSnps[referenceSequence] = self.novelSnps.get(referenceSequence, 0) + 1
    else:
      self.knownSnps[referenceSequence] = self.knownSnps.get(referenceSequence, 0) + 1

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

  def printGeneralStats(self):
    allNovelSnps = 0
    allKnownSnps = 0
    allTransitions = 0
    allTransversions = 0
    allMultiAllelic = 0

    print '%(text1)20s  %(text2)18s  %(text3)18s  %(text4)7s  %(text5)15s  %(text6)15s  %(text7)11s  %(text8)15s' % \
          {"text1": "reference sequence", "text2": "total # novel SNPS", "text3": "total # known SNPs", \
           "text4": "% dbSNP", "text5": "# transitions", "text6": "# transversions", "text7": "Ts/Tv ratio", "text8": "# multiAllelic"}
    for ref in sorted(self.referenceSequences):
      novelSnps = self.novelSnps[ref] if self.novelSnps.has_key(ref) else 0
      knownSnps = self.knownSnps[ref] if self.knownSnps.has_key(ref) else 0
      transitions = self.transitions[ref] if self.transitions.has_key(ref) else 0
      transversions = self.transversions[ref] if self.transversions.has_key(ref) else 0
      multiAllelic = self.multiAllelic[ref] if self.multiAllelic.has_key(ref) else 0

      allNovelSnps += novelSnps
      allKnownSnps += knownSnps
      allTransitions += transitions
      allTransversions += transversions
      allMultiAllelic += multiAllelic
      dbsnp = 100*knownSnps/(knownSnps + novelSnps)
      dbsnp = 100*float(knownSnps)/( float(knownSnps) + float(novelSnps) ) if (knownSnps + novelSnps) != 0 else 0
      tsTv = float(transitions)/float(transversions) if transversions != 0 else 0

      print '%(ref)20s  %(novelSnps)18d  %(knownSnps)18d  %(dbsnp)7.2f  %(transitions)15d  %(transversions)15d  %(tstv)11.2f  %(multiAllelic)15d' % \
            {"ref": ref, "novelSnps": novelSnps, "knownSnps": knownSnps, "dbsnp": dbsnp, "transitions": transitions, \
             "transversions": transversions, "tstv": tsTv, "multiAllelic": multiAllelic}

    dbsnp = 100*float(allKnownSnps)/( float(allKnownSnps) + float(allNovelSnps) ) if (allKnownSnps + allNovelSnps) != 0 else 0
    tsTv = float(allTransitions)/float(allTransversions) if allTransversions != 0 else 0

    print
    print '%(ref)20s  %(novelSnps)18d  %(knownSnps)18d  %(dbsnp)7.2f  %(transitions)15d  %(transversions)15d  %(tstv)11.2f  %(multiAllelic)15d' % \
          {"ref": "total", "novelSnps": allNovelSnps, "knownSnps": allKnownSnps, "dbsnp": dbsnp, "transitions": allTransitions, \
           "transversions": allTransversions, "tstv": tsTv, "multiAllelic": allMultiAllelic}

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
    stats.processGeneralStats(v.rsid, v.ref, v.alt, v.multiAllelic, v.referenceSequence)
    if options.distributions:
      for tag in options.distributions:
        tagNumber, tagValue = v.getInfo(tag)

# Close the file.

  v.closeVcf(options.vcfFile)

# Print out the stats.

  stats.printGeneralStats()
