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

# Check that the reference and alternate in the dbsnp vcf file match those
# from the input vcf file.
def checkRefAlt(vcfRef, vcfAlt, dbsnpRef, dbsnpAlt, ref, position):
  text = "WARNING: ref and alt alleles differ between vcf and dbSNP. " + ref + ":" + str(position) + " vcf: " + \
         vcfRef + "/" + vcfAlt + ", dbsnp: " + dbsnpRef + "/" + dbsnpAlt

  if vcfRef.lower() != dbsnpRef.lower():
    if vcfRef.lower() != dbsnpAlt.lower(): print >> sys.stderr, text
  else:
    if vcfAlt.lower() != dbsnpAlt.lower(): print >> sys.stderr, text

# Intersect two vcf files.  It is assumed that the two files are
# sorted by genomic coordinates and the reference sequences are
# in the same order.
def annotateVcf(v, d, outputFile):
  success1 = v.getRecord()
  success2 = d.getRecord()
  currentReferenceSequence = v.referenceSequence

# Finish when the end of the first file has been reached.
  while success1:

# If the end of the dbsnp vcf file is reached, write out the
# remaining records from the vcf file.
    if not success2:
      outputFile.write(v.record)
      success1 = v.getRecord()

    if v.referenceSequence == d.referenceSequence:
      if v.position == d.position:
        v.rsid = d.getDbsnpInfo()
        checkRefAlt(v.ref, v.alt, d.ref, d.alt, v.referenceSequence, v.position)
        record = v.buildRecord(False)
        outputFile.write(record)

        success1 = v.getRecord()
        success2 = d.getRecord()
      elif d.position > v.position: success1 = v.parseVcf(d.referenceSequence, d.position, True, outputFile)
      elif v.position > d.position: success2 = d.parseVcf(v.referenceSequence, v.position, False, None)
    else:
      if v.referenceSequence == currentReferenceSequence: success1 = v.parseVcf(d.referenceSequence, d.position, True, outputFile)
      elif d.referenceSequence == currentReferenceSequence: success2 = d.parseVcf(v.referenceSequence, v.position, False, None)
      currentReferenceSequence = v.referenceSequence

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py dbsnp [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf files")
  parser.add_option("-d", "--dbsnp",
                    action="store", type="string",
                    dest="dbsnpFile", help="input dbsnp vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that a single  vcf file is given.
  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (--in, -i) is required for dbsnp annotation."
    exit(1)

# Check that a dbsnp vcf file is included.
  if options.dbsnpFile == None:
    parser.print_help()
    print >> sys.stderr, "\ndbSNP vcf file is required (-d, --dbsnp)."
    exit(1)

# Set the output file to stdout if no output file was specified.
  outputFile, writeOut = setOutput(options.output) # tools.py

  v = vcf() # Define vcf object.
  d = vcf() # Define dbsnp vcf object.
  d.processInfo = True
  d.dbsnpVcf = True

# Open the vcf files.
  v.openVcf(options.vcfFile)
  d.openVcf(options.dbsnpFile)

# Read in the header information.
  v.parseHeader(options.vcfFile, writeOut)
  d.parseHeader(options.dbsnpFile, writeOut)

# Add an extra line to the vcf header to indicate the file used for
# performing dbsnp annotation.
  taskDescriptor = "##vcfPytools=annotated vcf file with dbSNP file " + options.dbsnpFile
  writeHeader(outputFile, v, False, taskDescriptor) # tools.py

# Annotate the vcf file.
  annotateVcf(v, d, outputFile)

# Check that the input files had the same list of reference sequences.
# If not, it is possible that there were some problems.
  checkReferenceSequenceLists(v.referenceSequenceList, d.referenceSequenceList) # tools.py

# Close the vcf files.
  v.closeVcf(options.vcfFile)
  d.closeVcf(options.dbsnpFile)

# End the program.
  return 0
