#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

# Check that the reference and alternate in the dbsnp vcf file match those
# from the input vcf file.

def checkRefAlt(vcfRef, vcfAlt, dbsnpRef, dbsnpAlt, ref, position):
  if vcfRef.lower() != dbsnpRef.lower() or vcfAlt.lower() != dbsnpAlt.lower():
    text = "WARNING: " + ref + ":" + str(position) + \
           " has different bases than the dbsnp entry\n\tref: " + dbsnpRef + \
           "(" + vcfRef + "), alt: " + dbsnpAlt + "(" + vcfAlt + ")\n"

    print >> sys.stderr, text

def calculateIntersection(v, dbsnp, line, dbsnpLine, vcfReferenceSequences, outputFile):

# If the second vcf file is at a different reference sequence, parse
# through the file until records for this reference are found.

  currentReferenceSequence = v.referenceSequence
  if dbsnp.referenceSequence != v.referenceSequence:
    if vcfReferenceSequences[dbsnp.referenceSequence] == "unparsed": vcfReferenceSequences[dbsnp.referenceSequence] = "skipped"
    while dbsnp.referenceSequence != v.referenceSequence:
      dbsnpLine = dbsnp.filehandle.readline()
      if not dbsnpLine:
        print >> sys.stderr, "Error occurred in intersection calculation."
        print >> sys.stderr, "Couldn't locate reference sequence:", dbsnp.referenceSequence
        exit(1)
      dbsnp.getRecord(dbsnpLine)

  while v.referenceSequence == currentReferenceSequence:

# If all of the entries in the second vcf file have been processed,
# parse through the remaining records in v as no more intersections
# can be found for this reference sequence.

    if dbsnp.referenceSequence != currentReferenceSequence:
      while v.referenceSequence == currentReferenceSequence:
        outputFile.write(line)
        line = v.filehandle.readline()
        if not line: break
        v.getRecord(line)
      break

# If the position in the first vcf file is smaller than that in the
# second vcf file, move to the next record in the first vcf file as
# this is not a shared record.

    if v.position < dbsnp.position:
      outputFile.write(line)
      line = v.filehandle.readline()
      if not line: break
      v.getRecord(line)

# If the positions are equal, thie record is present in both vcf files
# and so is written to the output.

    elif v.position == dbsnp.position:
      v.rsid = dbsnp.getDbsnpInfo(dbsnpLine)
      checkRefAlt(v.ref, v.alt, dbsnp.ref, dbsnp.alt, v.referenceSequence, v.position)
      newRecord = v.buildRecord(False)
      outputFile.write( newRecord )

      line = v.filehandle.readline()
      if not line: break
      v.getRecord(line)

      dbsnpLine = dbsnp.filehandle.readline()
      if not dbsnpLine: break
      dbsnp.getRecord(dbsnpLine)

    else:
      if dbsnp.referenceSequence == currentReferenceSequence:
        dbsnpLine = dbsnp.filehandle.readline()

# If the second vcf file is exhausted, parse through the remaining
# records for this reference sequence in v as no more intersections
# can be found for it.

        if not dbsnpLine: 
          while v.referenceSequence == currentReferenceSequence:
            outputFile.write(line)
            line = v.filehandle.readline()
            if not line: break
            v.getRecord(line)
          break
        dbsnp.getRecord(dbsnpLine)

        while dbsnp.referenceSequence == currentReferenceSequence and dbsnp.position <= v.position:

# If dbsnp.position = v.position, also iterate the record in v.

          if v.position == dbsnp.position:
            v.rsid = dbsnp.getDbsnpInfo(dbsnpLine)
            checkRefAlt(v.ref, v.alt, dbsnp.ref, dbsnp.alt, v.referenceSequence, v.position)
            newRecord = v.buildRecord(False)
            outputFile.write( newRecord )

            line = v.filehandle.readline()
            if not line: break
            v.getRecord(line)

          dbsnpLine = dbsnp.filehandle.readline()
          if not dbsnpLine: 
            while v.referenceSequence == currentReferenceSequence:
              outputFile.write( line )
              line = v.filehandle.readline()
              if not line: break
              v.getRecord(line)
            break
          dbsnp.getRecord(dbsnpLine)

# If dbsnp has moved on to the next reference sequence, parse through
# the rest of the records in v until the end of this reference
# sequence as no more intersections can be found for this
# reference sequence.

      else:
        while v.referenceSequence == currentReferenceSequence:
          outputFile.write( line )
          line = v.filehandle.readline()
          if not line: break
          v.getRecord(line)

  return v, dbsnp, line, dbsnpLine, vcfReferenceSequences

if __name__ == "__main__":
  main()

def main():

# Parse the command line options

  usage = "Usage: vcfPytools.py dbsnp [options]"
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

# Check that multiple vcf files are given.

  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (--in, -i) is required for dbsnp calculation."
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
  dbsnp = vcf() # Define vcf object.
  dbsnp.processInfo = True
  dbsnp.dbsnpVcf = True

# Read in the reference sequences present in the second vcf file.

  dbsnp.openVcf(options.dbsnpFile)
  dbsnp.parseHeader(options.dbsnpFile, False, False)
  for dbsnpLine in dbsnp.filehandle:
    dbsnp.getRecord(dbsnpLine)
    dbsnp.referenceSequences[ dbsnp.referenceSequence ] = "unparsed"
  vcfReferenceSequences = dbsnp.referenceSequences.copy()
  dbsnp.closeVcf(options.dbsnpFile)

# Open the vcf files.

  v.openVcf(options.vcfFile)
  dbsnp.openVcf(options.dbsnpFile)

# Read in the header information.

  v.parseHeader(options.vcfFile, writeOut, True)
  dbsnp.parseHeader(options.dbsnpFile, writeOut, True)

# Check that the header for the two files contain the same samples.

  dbsnp.headerText = dbsnp.headerText + "##dbsnp=" + options.dbsnpFile + "\n"
  writeHeader(outputFile, dbsnp, False) # tools.py

# Get the first record from both vcf files.

  line = v.filehandle.readline()
  dbsnpLine = dbsnp.filehandle.readline()
  v.getRecord(line)
  dbsnp.getRecord(dbsnpLine)

# Calculate the intersection. Check if the records from the two vcf files 
# correspond to the same reference sequence.  If so, search up to 
# the same position and write out the record if it exists in the 
# second vcf file.

  while True:
    if vcfReferenceSequences.has_key(v.referenceSequence) and vcfReferenceSequences[v.referenceSequence] == "skipped":
      dbsnp.closeVcf(options.dbsnpFile)
      dbsnp.openVcf(options.dbsnpFile)
      dbsnp.parseHeader(options.dbsnpFile, False, False)
      dbsnpLine = dbsnp.filehandle.readline()
      dbsnp.getRecord(dbsnpLine)
      for key, value in vcfReferenceSequences.iteritems():
        if vcfReferenceSequences[key] == "skipped":
          vcfReferenceSequences[key] = "unparsed"

    if vcfReferenceSequences.has_key(v.referenceSequence) and vcfReferenceSequences[v.referenceSequence] != "completed":
      vcfReferenceSequences[v.referenceSequence] = "completed"
      v, dbsnp, line, dbsnpLine, vcfReferenceSequences = calculateIntersection(v, dbsnp, line, dbsnpLine, vcfReferenceSequences, outputFile)
    elif not v.referenceSequence in vcfReferenceSequences:
      currentReferenceSequence = v.referenceSequence
      while v.referenceSequence == currentReferenceSequence:
        outputFile.write( line )
        line = v.filehandle.readline()
        if not line: break
        v.getRecord(line)

# If the end of the first vcf file has been reached, there can be no
# more intersections, so the calculation is complete.

    if not line: break

# Close the vcf files.

  v.closeVcf(options.vcfFile)
  dbsnp.closeVcf(options.dbsnpFile)

  exit(0)
