#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

class statistics:
  def __init__(self):
    self.referenceSequences = {}
    self.novelTransitions = {}
    self.knownTransitions = {}
    self.novelTransversions = {}
    self.knownTransversions = {}
    self.multiAllelic = {}
    self.distributions = {}

  def processGeneralStats(self, referenceSequence, rsid, ref, alt, multiAllelic, filters):
    self.referenceSequences[referenceSequence] = True

# Determine if the SNP is a transition, transversion or multi-allelic.

    if multiAllelic == True:
      if referenceSequence not in self.multiAllelic:
        self.multiAllelic[referenceSequence] = {}

      self.multiAllelic[referenceSequence][filters] = self.multiAllelic[referenceSequence].get(filters, 0) + 1
    else:
      if ref < alt:
        alleles = ref + alt
      else:
        alleles = alt + ref
    
# Increment the number of transitions.  Keep track of whether they are novel or known.

      if alleles.lower() == "ag" or alleles.lower() == "ct":
        if rsid == ".":
          if referenceSequence not in self.novelTransitions:
            self.novelTransitions[referenceSequence] = {}

          self.novelTransitions[referenceSequence][filters] = self.novelTransitions[referenceSequence].get(filters, 0) + 1
        else:
          if referenceSequence not in self.knownTransitions:
            self.knownTransitions[referenceSequence] = {}

          self.knownTransitions[referenceSequence][filters] = self.knownTransitions[referenceSequence].get(filters, 0) + 1

# Increment the number of transitions.  Keep track of whether they are novel or known.

      elif alleles.lower() == "ac" or alleles.lower() == "at" or alleles.lower() == "cg" or alleles.lower() == "gt":
        if rsid == ".":
          if referenceSequence not in self.novelTransversions:
            self.novelTransversions[referenceSequence] = {}

          self.novelTransversions[referenceSequence][filters] = self.novelTransversions[referenceSequence].get(filters, 0) + 1
        else:
          if referenceSequence not in self.knownTransversions:
            self.knownTransversions[referenceSequence] = {}

          self.knownTransversions[referenceSequence][filters] = self.knownTransversions[referenceSequence].get(filters, 0) + 1

  def printGeneralStats(self):
    allNovelTransitions = {}
    allKnownTransitions = {}
    allNovelTransversions = {}
    allKnownTransversions = {}
    allMultiAllelic = {}
    allFilters = {}
    allFilters["total"] = True

    print '%(text1)20s  %(text2)58s  %(text3)7s  %(text4)22s' % \
          {"text1": "", \
           "text2": "--------------------------# SNPs--------------------------", \
           "text3": "", \
           "text4": "------ts/tv ratio-----"}
    print '%(text1)20s  %(text2)10s  %(text3)10s  %(text4)10s  %(text5)10s  %(text6)10s  %(text7)7s  %(text8)6s  %(text9)6s  %(text10)6s  %(text11)14s' % \
          {"text1" : "filter", \
           "text2" : "total ", \
           "text3" : "novel ts", \
           "text4" : "novel tv", \
           "text5" : "known ts", \
           "text6" : "known tv", \
           "text7" : "% dbSNP", \
           "text8" : "total", \
           "text9" : "novel", \
           "text10": "known", \
           "text11": "# multiAllelic"}

    for ref in sorted(self.referenceSequences):
      novelTransitions = {}
      knownTransitions = {}
      novelTransversions = {}
      knownTransversions = {}
      multiAllelic = {}

      print "\nreference sequence:", ref

# Count up the novel transitions for each filter.

      if ref in self.novelTransitions:
        for key, value in self.novelTransitions[ref].iteritems():
          novelTransitions["total"] = novelTransitions.get("total", 0) + value
          allNovelTransitions["total"] = allNovelTransitions.get("total", 0) + value
          filters = key.split(";")
          for filter in filters:
            novelTransitions[filter] = novelTransitions.get(filter, 0) + value
            allNovelTransitions[filter] = allNovelTransitions.get(filter, 0) + value
            allFilters[filter] = True

# Count up the known transitions for each filter.

      if ref in self.knownTransitions:
        for key, value in self.knownTransitions[ref].iteritems():
          knownTransitions["total"] = knownTransitions.get("total", 0) + value
          allKnownTransitions["total"] = allKnownTransitions.get("total", 0) + value
          filters = key.split(";")
          for filter in filters:
            knownTransitions[filter] = knownTransitions.get(filter, 0) + value
            allKnownTransitions[filter] = allKnownTransitions.get(filter, 0) + value
            allFilters[filter] = True

# Count up the novel transversions for each filter.

      if ref in self.novelTransversions:
        for key, value in self.novelTransversions[ref].iteritems():
          novelTransversions["total"] = novelTransversions.get("total", 0) + value
          allNovelTransversions["total"] = allNovelTransversions.get("total", 0) + value
          filters = key.split(";")
          for filter in filters:
            novelTransversions[filter] = novelTransversions.get(filter, 0) + value
            allNovelTransversions[filter] = allNovelTransversions.get(filter, 0) + value
            allFilters[filter] = True

# Count up the known transversions for each filter.

      if ref in self.knownTransversions:
        for key, value in self.knownTransversions[ref].iteritems():
          knownTransversions["total"] = knownTransversions.get("total", 0) + value
          allKnownTransversions["total"] = allKnownTransversions.get("total", 0) + value
          filters = key.split(";")
          for filter in filters:
            knownTransversions[filter] = knownTransversions.get(filter, 0) + value
            allKnownTransversions[filter] = allKnownTransversions.get(filter, 0) + value
            allFilters[filter] = True

# Count up the number of multi-allelic sites.

      if ref in self.multiAllelic:
        for key, value in self.multiAllelic[ref].iteritems():
          filters = key.split(";")
          for filter in filters:
            multiAllelic[filter] = multiAllelic.get(filter, 0) + value
            allMultiAllelic[filter] = allMultiAllelic.get(filter, 0) + value
            allFilters[filter] = True

# Create a list of the filters and put "total" and "PASS" at the end.

      filterList = []
      for filter, value in allFilters.iteritems():
        if filter != "total" and filter != "PASS":
          filterList.append(filter)

      filterList.append("total")
      if "PASS" in allFilters:
        filterList.append("PASS")

# Calculate the dbsnp fraction and Ts/Tv ratio for each filter.

      for index, filter in enumerate(filterList):
        novelTs = novelTransitions.get(filter, 0)
        novelTv = novelTransversions.get(filter, 0)
        knownTs = knownTransitions.get(filter, 0) 
        knownTv = knownTransversions.get(filter, 0)
        multi = multiAllelic.get(filter, 0)

        novel = novelTs + novelTv
        known = knownTs + knownTv
        totalSnp = novelTs + novelTv + knownTs + knownTv

        transitions = novelTs + knownTs
        transversions = novelTv + knownTv

        dbsnp = 100 * float(known) / ( float(known) + float(novel) ) if (known + novel) != 0 else 0
        noveltstv = float(novelTs) / float(novelTv) if novelTv != 0 else 0.
        knowntstv = float(knownTs) / float(knownTv) if knownTv != 0 else 0.
        tstv = float(transitions) / float(transversions) if transversions != 0 else 0.

        if filter == "total":
          print "\n               ---------------------------------------------------------------------------------------------------------------------------------"

        print '%(filter)20s  %(total)10d  %(novelTs)10d  %(novelTv)10d  %(knownTs)10d  %(knownTv)10d  %(dbsnp)7.2f  \
%(totaltstv)6.2f  %(noveltstv)6.2f  %(knowntstv)6.2f  %(multi)14d' % \
              {"filter" : filter, \
               "total" : totalSnp, \
               "novelTs" : novelTs, \
               "novelTv" : novelTv, \
               "knownTs" : knownTs, \
               "knownTv" : knownTv, \
               "dbsnp" : dbsnp, \
               "totaltstv" : tstv, \
               "noveltstv" : noveltstv, \
               "knowntstv": knowntstv, \
               "multi": multi}

        if len(filterList) - 1 == index:
          print "               ---------------------------------------------------------------------------------------------------------------------------------"

    print "\nTotal for all reference sequences"
    for index, filter in enumerate(filterList):
      novelTs = allNovelTransitions.get(filter, 0)
      novelTv = allNovelTransversions.get(filter, 0)
      knownTs = allKnownTransitions.get(filter, 0) 
      knownTv = allKnownTransversions.get(filter, 0)
      multi = multiAllelic.get(filter, 0)

      novel = novelTs + novelTv
      known = knownTs + knownTv
      totalSnp = novelTs + novelTv + knownTs + knownTv

      transitions = novelTs + knownTs
      transversions = novelTv + knownTv

      #novel = allNovelTransitions.get(filter, 0) + allNovelTransversions.get(filter, 0)
      #known = allKnownTransitions.get(filter, 0) + allKnownTransversions.get(filter, 0)
      #transitions = allNovelTransitions.get(filter, 0) + allKnownTransitions.get(filter, 0)
      #transversions = allNovelTransversions.get(filter, 0) + allKnownTransversions.get(filter, 0)
      #multi = multiAllelic.get(filter, 0)

      dbsnp = 100 * float(known) / ( float(known) + float(novel) ) if (known + novel) != 0 else 0
      noveltstv = float(novelTs) / float(novelTv) if novelTv != 0 else 0
      knowntstv = float(knownTs) / float(knownTv) if knownTv != 0 else 0
      tstv = float(transitions) / float(transversions) if transversions != 0 else 0

      if filter == "total":
        print "\n               ---------------------------------------------------------------------------------------------------------------------------------"

      print '%(filter)20s  %(total)10d  %(novelTs)10d  %(novelTv)10d  %(knownTs)10d  %(knownTv)10d  %(dbsnp)7.2f  \
%(totaltstv)6.2f  %(noveltstv)6.2f  %(knowntstv)6.2f  %(multi)14d' % \
            {"filter" : filter, \
             "total" : totalSnp, \
             "novelTs" : novelTs, \
             "novelTv" : novelTv, \
             "knownTs" : knownTs, \
             "knownTv" : knownTv, \
             "dbsnp" : dbsnp, \
             "totaltstv" : tstv, \
             "noveltstv" : noveltstv, \
             "knowntstv": knowntstv, \
             "multi": multi}

      if len(filterList) - 1 == index:
        print "               ---------------------------------------------------------------------------------------------------------------------------------"

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
    stats.processGeneralStats(v.referenceSequence, v.rsid, v.ref, v.alt, v.multiAllelic, v.filters)
    if options.distributions:
      for tag in options.distributions:
        tagNumber, tagValue = v.getInfo(tag)

# Close the file.

  v.closeVcf(options.vcfFile)

# Print out the stats.

  stats.printGeneralStats()
