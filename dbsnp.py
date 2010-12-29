#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

if __name__ == "__main__":
  main()

# Check that the reference and alternate in the dbsnp vcf file match those
# from the input vcf file.

def checkRefAlt(vcfRef, vcfAlt, dbsnpRef, dbsnpAlt, ref, position):
  if vcfRef.lower() != dbsnpRef.lower() or vcfAlt.lower() != dbsnpAlt.lower():
    text = "WARNING: " + ref + ":" + str(position) + \
           " has different bases than the dbsnp entry\n\tref: " + dbsnpRef + \
           "(" + vcfRef + "), alt: " + dbsnpAlt + "(" + vcfAlt + ")\n"

    print >> sys.stderr, text

def main():

# Parse the command line options

  usage = "Usage: vcfTools.py intersect [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file")
  parser.add_option("-d", "--dbsnp",
                    action="store", type="string",
                    dest="dbsnpFile", help="input dbsnp vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that a single vcf file is given.

  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (-i, --input) is required for dbsnp annotation."
    exit(1)

# Check that a dbsnp vcf file is included.

  if options.dbsnpFile == None:
    parser.print_help()
    print >> sys.stderr, "\ndbSNP vcf file is required (-d, --dbsnp)."
    exit(1)

# Set the output file to stdout if no output file was specified.

  if options.output == None:
    outputFile = sys.stdout
    writeOut = False
  else:
    outputFile = open(options.output, 'w')
    writeOut = True

  v = vcf() # Define vcf object.
  dbsnp = vcf() # Define dbsnp vcf object.
  dbsnp.processInfo = True
  dbsnp.dbsnpVcf = True

# Read in the reference sequences present in the dbsnp vcf file.

  dbsnp.openVcf(options.dbsnpFile)
  dbsnp.parseHeader(options.dbsnpFile, False, False)
  for dbsnpLine in dbsnp.filehandle:
    dbsnp.getRecord(dbsnpLine)
    dbsnp.referenceSequences[ dbsnp.referenceSequence ] = False
  vcfReferenceSequences = dbsnp.referenceSequences
  dbsnp.closeVcf(options.dbsnpFile)

# Open the vcf files.

  v.openVcf(options.vcfFile)
  dbsnp.openVcf(options.dbsnpFile)

# Read in the header information.

  v.parseHeader(options.vcfFile, writeOut, True)
  dbsnp.parseHeader(options.dbsnpFile, writeOut, True)

# Write out the header information to the new output file.  Include an extra
# header line indicating the dbSNP file used for the annotation.

  dbsnpText = "##dbsnp=" + options.dbsnpFile + "\n"
  outputFile.write( v.headerText ) if v.headerText != "" else None
  outputFile.write( dbsnpText )
  outputFile.write( v.headerInfoText ) if v.headerInfoText != "" else None
  outputFile.write( v.headerFormatText ) if v.headerFormatText != "" else None
  outputFile.write( v.headerTitles)

# Get the first line of the second vcf file.

  for dbsnpLine in dbsnp.filehandle:
    dbsnp.getRecord(dbsnpLine)
    break

# Determine the intersection of the vcf and dbSNP vcf files.

  for vcfLine in v.filehandle:
    v.getRecord(vcfLine)

# Check if the records from the two vcf files correspond to the same
# reference sequence.  If so, search up to the same position and
# write out the record if it exists in the second vcf file.

    if v.referenceSequence == dbsnp.referenceSequence:
      vcfReferenceSequences[v.referenceSequence] = True
      if v.position == dbsnp.position:
        v.rsid = dbsnp.getDbsnpInfo(dbsnpLine)
        checkRefAlt(v.ref, v.alt, dbsnp.ref, dbsnp.alt, v.referenceSequence, v.position)
        newRecord = v.buildRecord()
        outputFile.write( newRecord )
      elif v.position > dbsnp.position:
        for dbsnpLine in dbsnp.filehandle:
          dbsnp.getRecord(dbsnpLine)
          if v.referenceSequence != dbsnp.referenceSequence:
            outputFile.write( vcfLine )
            break
          if dbsnp.position > v.position:
            outputFile.write( vcfLine )
            break
          elif v.position == dbsnp.position:
            v.rsid = dbsnp.getDbsnpInfo(dbsnpLine)
            checkRefAlt(v.ref, v.alt, dbsnp.ref, dbsnp.alt, v.referenceSequence, v.position)
            newRecord = v.buildRecord()
            outputFile.write( newRecord )
            break
      else:
        outputFile.write( vcfLine )

# Check if the dbsnp file includes the reference sequence.

    elif vcfReferenceSequences.has_key(v.referenceSequence):

# If the reference sequence in the record from the first vcf file exists
# in the second, but has not been read yet, parse through the second
# vcf file until this reference sequence is reached, then search for the
# same position.

      if vcfReferenceSequences[v.referenceSequence] == False:
        for dbsnpLine in dbsnp.filehandle:
          rsid = dbsnp.getRecord(dbsnpLine)
          vcfReferenceSequences[dbsnp.referenceSequence] = True
          if v.referenceSequence == dbsnp.referenceSequence:
            if v.position == dbsnp.position:
              v.rsid = dbsnp.getDbsnpInfo(dbsnpLine)
              checkRefAlt(v.ref, v.alt, dbsnp.ref, dbsnp.alt, v.referenceSequence, v.position)
              newRecord = v.buildRecord()
              outputFile.write( newRecord )
              break
            elif v.position < dbsnp.position:
              outputFile.write( vcfLine )
              break

# If the reference sequence in the record from the first vcf file exists
# in the second and has already been parsed, close and reopen the second
# vcf file, then allow the search to begin again from the beginning of the
# file.

      elif vcfReferenceSequences[v.referenceSequence] == True:
        dbsnp.closeVcf(options.dbsnpFile)
        dbsnp.openVcf(options.dbsnpFile)
        dbsnp.parseHeader(options.dbsnpFile, False, False)
        for ref in vcfReferenceSequences:
          vcfReferenceSequences[ref] = False
        for dbsnpLine in dbsnp.filehandle:
          dbsnp.getRecord(dbsnpLine)
          if v.referenceSequence == dbsnp.referenceSequence:
            if v.position == dbsnp.position:
              v.rsid = dbsnp.getDbsnpInfo(dbsnpLine)
              checkRefAlt(v.ref, v.alt, dbsnp.ref, dbsnp.alt, v.referenceSequence, v.position)
              newRecord = v.buildRecord()
              outputFile.write( newRecord )
              break
            elif v.position < dbsnp.position:
              outputFile.write( vcfLine )
              break

    else:
      outputFile.write( vcfLine )

# Close the vcf files.

  v.closeVcf(options.vcfFile)
  dbsnp.closeVcf(options.dbsnpFile)
  exit(0)
