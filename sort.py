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

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py sort [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="store", type="string",
                    dest="vcfFile", help="input vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that a single vcf file is given.
  if options.vcfFile == None:
    parser.print_help()
    print >> sys.stderr, "\nInput vcf file (-i, --input) is required for vcf sorting."
    exit(1)

# Set the output file to stdout if no output file was specified.
  outputFile, writeOut = setOutput(options.output)

  v = vcf() # Define vcf object.

# Open the vcf file.
  v.openVcf(options.vcfFile)

# Read in the header information.
  v.parseHeader(options.vcfFile, writeOut)
  writeHeader(outputFile, v, False)

# Parse the vcf file and for each reference sequence, output
# the position to a temp file.  These files will allow each
# reference sequence to be dealt with individually, reducing
# the amount of information requiring storage in memory.
  tempFiles = {}
  tempPositionsFiles = {}
  while v.getRecord():
    tempPositionsFile = "positions." + v.referenceSequence + ".vcfPytools.tmp"
    tempFile = "records." + v.referenceSequence + ".vcfPytools.tmp"
    if tempFile not in tempFiles:
      tempFilehandle = open(tempFile,'w')
      tempPositionsFilehandle = open(tempPositionsFile,'w')
      tempPositionsFiles[tempPositionsFile] = tempPositionsFilehandle
      tempFiles[tempFile] = tempFilehandle

    tempPositionsFiles[tempPositionsFile].write( str(v.position) + "\n" )
    tempFiles[tempFile].write(v.record)

# Close the temp files.
  for tempFile in tempFiles:
    tempFiles[tempFile].close()
  for tempPositionsFile in tempPositionsFiles:
    tempPositionsFiles[tempPositionsFile].close()

# Read in the positions, sort them and replace the temp positions
# file with a sorted file.
  for tempPositionsFile in tempPositionsFiles:
    positionsList = []
    tempPositionsFilehandle = open(tempPositionsFile, 'r')
    for line in tempPositionsFilehandle:
      positionsList.append( int(line) )

    tempPositionsFilehandle.close()
    os.remove(tempPositionsFile)

# Sort the list.
    sortedPositionsList = sorted(positionsList)
    tempPositionsFilehandle = open(tempPositionsFile, 'w')
    for line in sortedPositionsList:
      tempPositionsFilehandle.write( str(line) + "\n" )
    tempPositionsFilehandle.close()

# Read through the sorted positions file and the file containing
# only the records for that reference sequence and write the
# records to the output file in order.
  for referenceSequence in v.referenceSequencesList:
    v1 = vcf()
    positionsFile = "positions." + referenceSequence + ".vcfPytools.tmp"
    recordsFile = "records." + referenceSequence + ".vcfPytools.tmp"
    positionsFilehandle = open(positionsFile, 'r')
    v1.openVcf(recordsFile)

# Get the first position from the positions file and the first
# record from the records file.
    success = v1.getRecord()
    position = int( positionsFilehandle.readline() )

# Create a dictionary containing the positions of records that
# have been written to their own temp file.
    storedRecords = {}

# Sort and output.
    while True:

# If this position has already been seen in the vcf file, then
# it was written to a temp file and a dictionary key set.  Now
# write this record to file and delete the dictionay key.
      if position in storedRecords:
        storedTemp = "stored." + str(position) + ".vcfPytools.tmp"
        storedTempFilehandle = open(storedTemp,'r')
        storedLine = storedTempFilehandle.readline()
        outputFile.write( storedLine )
        storedTempFilehandle.close()
        os.remove(storedTemp)
        del storedRecords[position]

        position = int( positionsFilehandle.readline() )
        if not position: break

# If the position in the positions file is less than that of the
# record from the vcf file, this record needs to be saved until
# the correct time.  Create a temp file to hold it and add the
# position to a dictionary.
      elif position < v1.position:
        storedTemp = "stored." + str(v1.position) + ".vcfPytools.tmp"
        storedTempFilehandle = open(storedTemp,'w')
        storedTempFilehandle.write(v1.record)
        storedTempFilehandle.close()
        storedRecords[v1.position] = True

        success = v1.getRecord()
        if not success: break

# If the position in the vcf file agrees with that in the positions
# file, then the record can be written to file.
      elif position == v1.position:
        outputFile.write(v1.record)

        success = v1.getRecord()
        if not success: break

        position = int( positionsFilehandle.readline() )
        if not position: break

      elif position > v1.position:
        print >> sys.stderr, "I don't know how this happened."

# Check if any records remain in stored files.  If so, write them out
# to the output file in order.
    if len( storedRecords ) > 0:
      for position in sorted(storedRecords):
        storedTemp = "stored." + str(position) + ".vcfPytools.tmp"
        storedTempFilehandle = open(storedTemp,'r')
        storedLine = storedTempFilehandle.readline()
        outputFile.write( storedLine )
        storedTempFilehandle.close()
        os.remove(storedTemp)
        del storedRecords[position]

# Close and delete the temp files for this reference sequence.
    positionsFilehandle.close()
    v1.closeVcf(recordsFile)
    os.remove(positionsFile)
    os.remove(recordsFile)

# Close the vcf file.
  v.closeVcf(options.vcfFile)

# Terminate the program cleanly.
  return 0
