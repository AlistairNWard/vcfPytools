#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

import tools
from tools import *

def calculateIntersection(v1, v2, line1, line2, vcfReferenceSequences, outputFile, priority):

# If the second vcf file is at a different reference sequence, parse
# through the file until records for this reference are found.

  currentReferenceSequence = v1.referenceSequence
  if v2.referenceSequence != v1.referenceSequence:
    if vcfReferenceSequences[v2.referenceSequence] == "unparsed": vcfReferenceSequences[v2.referenceSequence] = "skipped"
    while v2.referenceSequence != v1.referenceSequence:
      line2 = v2.filehandle.readline()
      if not line2:
        print >> sys.stderr, "Error occurred in intersection calculation."
        print >> sys.stderr, "Couldn't locate reference sequence:", v2.referenceSequence
        exit(1)
      v2.getRecord(line2)

  while v1.referenceSequence == currentReferenceSequence:

# If all of the entries in the second vcf file have been processed,
# parse through the remaining records in v1 as no more intersections
# can be found for this reference sequence.

    if v2.referenceSequence != currentReferenceSequence:
      while v1.referenceSequence == currentReferenceSequence:
        line1 = v1.filehandle.readline()
        if not line1: break
        v1.getRecord(line1)
      break

# If the position in the first vcf file is smaller than that in the
# second vcf file, move to the next record in the first vcf file as
# this is not a shared record.

    if v1.position < v2.position:
      line1 = v1.filehandle.readline()
      if not line1: break
      v1.getRecord(line1)

# If the positions are equal, thie record is present in both vcf files
# and so is written to the output.

    elif v1.position == v2.position:
      if priority == 0:
        if v1.quality >= v2.quality: outputFile.write(line1)
        else: outputFile.write(line2)
      elif priority == 1: outputFile.write(line1)
      elif priority == 2: outputFile.write(line2)

      line1 = v1.filehandle.readline()
      if not line1: break
      v1.getRecord(line1)

      line2 = v2.filehandle.readline()
      if not line2: break
      v2.getRecord(line2)

    else:
      if v2.referenceSequence == currentReferenceSequence:
        line2 = v2.filehandle.readline()

# If the second vcf file is exhausted, parse through the remaining
# records for this reference sequence in v1 as no more intersections
# can be found for it.

        if not line2: 
          while v1.referenceSequence == currentReferenceSequence:
            line1 = v1.filehandle.readline()
            if not line1: break
            v1.getRecord(line1)
          break
        v2.getRecord(line2)

        while v2.referenceSequence == currentReferenceSequence and v2.position <= v1.position:

# If v2.position = v1.position, also iterate the record in v1.

          if v1.position == v2.position:
            if priority == 0:
              if v1.quality >= v2.quality: outputFile.write(line1)
              else: outputFile.write(line2)
            elif priority == 1: outputFile.write(line1)
            elif priority == 2: outputFile.write(line2)

            line1 = v1.filehandle.readline()
            if not line1: break
            v1.getRecord(line1)

          line2 = v2.filehandle.readline()
          if not line2: 
            while v1.referenceSequence == currentReferenceSequence:
              line1 = v1.filehandle.readline()
              if not line1: break
              v1.getRecord(line1)
            break
          v2.getRecord(line2)

# If v2 has moved on to the next reference sequence, parse through
# the rest of the records in v1 until the end of this reference
# sequence as no more intersections can be found for this
# reference sequence.

      else:
        while v1.referenceSequence == currentReferenceSequence:
          line1 = v1.filehandle.readline()
          if not line1: break
          v1.getRecord(line1)

  return v1, v2, line1, line2, vcfReferenceSequences

if __name__ == "__main__":
  main()

def main():

# Parse the command line options

  usage = "Usage: vcfPytools.py intersection [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="append", type="string",
                    dest="vcfFiles", help="input vcf files")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")
  parser.add_option("-p", "--priority-file",
                    action="store", type="string",
                    dest="priorityFile", help="output record from this file")

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.

  if options.vcfFiles == None:
    parser.print_help()
    print >> sys.stderr, "\nTwo input vcf files (-i) are required for performing intersection."
    exit(1)
  elif len(options.vcfFiles) != 2:
    print >> sys.stderr, "Two input vcf files are required for performing intersection."

# Set the output file to stdout if no output file was specified.

  if options.output == None:
    outputFile = sys.stdout
    writeOut = False
  else:
    outputFile = open(options.output, 'w')
    writeOut = True

# If no priority is given to either file (from the -p command line
# option), set priorityQuality to True.  In this case, the record
# written to the output file will be that with the higest quality.
# If a priority is given, check that the file is one of the input
# vcf files.

  if options.priorityFile == None:
    priority = 0
  else:
    if options.priorityFile != options.vcfFiles[0] and options.priorityFile != options.vcfFiles[1]:
      print >> sys.stderr, "The file defined as having priority for writing out records must"
      print >> sys.stderr, "be one of the two input vcf files."
      exit(1)

    elif options.priorityFile == options.vcfFiles[0]:
      priority = 1
    
    elif options.priorityFile == options.vcfFiles[1]:
      priority = 2

  v1 = vcf() # Define vcf object.
  v2 = vcf() # Define vcf object.

# Read in the reference sequences present in the second vcf file.

  v2.openVcf(options.vcfFiles[1])
  v2.parseHeader(options.vcfFiles[1], False, False)
  for line2 in v2.filehandle:
    v2.getRecord(line2)
    v2.referenceSequences[ v2.referenceSequence ] = "unparsed"
  vcfReferenceSequences = v2.referenceSequences.copy()
  v2.closeVcf(options.vcfFiles[1])

# Open the vcf files.

  v1.openVcf(options.vcfFiles[0])
  v2.openVcf(options.vcfFiles[1])

# Read in the header information.

  v1.parseHeader(options.vcfFiles[0], writeOut, True)
  v2.parseHeader(options.vcfFiles[1], writeOut, True)

# Check that the header for the two files contain the same samples.

  if v1.samplesList != v2.samplesList:
    print >> sys.stderr, "vcf files contain different samples (or sample order)."
    exit(1)
  else:
    writeHeader(outputFile, v1) # tools.py

# Get the first record from both vcf files.

  line1 = v1.filehandle.readline()
  line2 = v2.filehandle.readline()
  v1.getRecord(line1)
  v2.getRecord(line2)

# Calculate the union. Check if the records from the two vcf files 
# correspond to the same reference sequence.  If so, search up to 
# the same position and write out the record if it exists in the 
# second vcf file.

  while True:
    if vcfReferenceSequences.has_key(v1.referenceSequence) and vcfReferenceSequences[v1.referenceSequence] == "skipped":
      v2.closeVcf(options.vcfFiles[1])
      v2.openVcf(options.vcfFiles[1])
      v2.parseHeader(options.vcfFiles[1], False, False)
      line2 = v2.filehandle.readline()
      v2.getRecord(line2)
      for key, value in vcfReferenceSequences.iteritems():
        if vcfReferenceSequences[key] == "skipped":
          vcfReferenceSequences[key] = "unparsed"

    if vcfReferenceSequences.has_key(v1.referenceSequence) and vcfReferenceSequences[v1.referenceSequence] != "completed":
      vcfReferenceSequences[v1.referenceSequence] = "completed"
      v1, v2, line1, line2, vcfReferenceSequences = calculateIntersection(v1, v2, line1, line2, vcfReferenceSequences, outputFile, priority)
    elif not v1.referenceSequence in vcfReferenceSequences:
      currentReferenceSequence = v1.referenceSequence
      while v1.referenceSequence == currentReferenceSequence:
        line1 = v1.filehandle.readline()
        if not line1: break
        v1.getRecord(line1)

# If the end of the first vcf file has been reached, there can be no
# more intersections, so the calculation is complete.

    if not line1: break

# Close the vcf files.

  v1.closeVcf(options.vcfFiles[0])
  v2.closeVcf(options.vcfFiles[1])
  exit(0)
