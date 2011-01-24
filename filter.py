#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

if __name__ == "__main__":
  main()

def filterFail(text, file):
  print >> sys.stderr, text
  if file != None:
    os.remove(file)
  exit(1)

def main():

# Parse the command line options

  usage = "Usage: vcfPytools.py filter [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")
  parser.add_option("-q", "--quality",
                    action="store", type="int",
                    dest="quality", help="filter out SNPs with qualities lower than selected value")
  parser.add_option("-n", "--info",
                    action="append", type="string", nargs=2,
                    dest="infoFilters", help="filter based on entries in the info string")
  parser.add_option("-r", "--remove-genotypes",
                    action="store_true", default=False,
                    dest="removeGeno", help="remove the genotype strings from the vcf file")

  (options, args) = parser.parse_args()

# Check that a single vcf file is given.

  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (-i, --input) is required for vcf filtering."
    exit(1)

# Set the output file to stdout if no output file was specified.

  if options.output == None:
    outputFile = sys.stdout
    writeOut = False
  else:
    outputFile = open(options.output, 'w')
    writeOut = True

  v = vcf() # Define vcf object.

# Open the vcf file.

  v.openVcf(options.vcfFile)

# Read in the header information.

  v.parseHeader(options.vcfFile, writeOut, True)

# Check that specified filters from the info field are either integers or floats.

  if options.infoFilters:
    v.processInfo = True # Process the info string
    filters = {}
    for filter, value in options.infoFilters:
      if v.infoHeaderTags.has_key(filter):
        if v.infoHeaderTags[filter][1].lower() == "integer":
          try:
            filters[filter] = int(value)
          except ValueError:
            text = "Filter " + filter + " requires an integer entry, not " + str(type(value))
            filterFail(text, options.output)

        if v.infoHeaderTags[filter][1].lower() == "float":
          try:
            filters[filter] = float(value)
          except ValueError:
            text = "Filter " + filter + " requires an float entry, not " + str(type(value))
            filterFail(text, options.output)

      else:
        text = "Filter " + filter + " has no explanation in the header.  Unknown type for the entry."
        filterFail(text, options.output)

# Parse the vcf file and check if any of the filters are failed.  If
# so, build up a string of failed filters.

  writeHeader(outputFile, v, options.removeGeno)
  for line in v.filehandle:
    filterString = ""
    v.getRecord(line)

# Check for quality filtering.

    if options.quality != None:
      if int(v.quality) < options.quality:
        filterString = filterString + ";" + "Q" + str(options.quality) if filterString != "" else "Q" + str(options.quality)

# Check for filtering on info string filters.

    if options.infoFilters:
      for filter, value in filters.iteritems():
        if v.infoTags.has_key(filter):
          if type(value) == int:
            if int(v.infoTags[filter]) < value:
              filterString = filterString + ";" + filter + str(value) if filterString != "" else filter + str(value)
          elif type(value) == float:
            if float(v.infoTags[filter]) < value:
              filterString = filterString + ";" + filter + str(value) if filterString != "" else filter + str(value)

    filterString = "PASS" if filterString == "" else filterString
    v.filters = filterString
    newRecord = v.buildRecord(options.removeGeno)
    outputFile.write( newRecord )

# Close the vcf files.

  v.closeVcf(options.vcfFile)
  exit(0)
