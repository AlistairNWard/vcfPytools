#!/usr/bin/python

import os.path
import sys

# Determine whether to output to a file or stdout.
def setOutput(output):
  if output == None:
    outputFile = sys.stdout
    writeOut = False
  else:
    output = os.path.abspath(output)
    outputFile = open(output, 'w')
    writeOut = True

  return outputFile, writeOut

# Write the header to file.
def writeHeader (outputFile, v, removeGenotypes):
  outputFile.write( v.headerText ) if v.headerText != "" else None
  outputFile.write( v.headerInfoText ) if v.headerInfoText != "" else None
  outputFile.write( v.headerFormatText ) if v.headerFormatText != "" else None
  if removeGenotypes:
    line = v.headerTitles.split("\t")
    newHeaderTitles = line[0]
    for i in range(1,8):
      newHeaderTitles = newHeaderTitles + "\t" + line[i]
    newHeaderTitles = newHeaderTitles + "\n"
    outputFile.write( newHeaderTitles )
  else:
    outputFile.write( v.headerTitles )

# Check that the two reference sequence lists are identical.
# If there are a different number or order, the results may
# not be as expected.
def checkReferenceSequenceLists(list1, list2):
  if len(list1) != len(list2):
    print >> sys.stderr, "WARNING: Input files contain a different number of reference sequences."
    print >> sys.stderr, "Results may not be as expected."
    print >> sys.stderr, "Ensure that input files have the same reference sequences in the same order."
  elif list1 != list2:
    print >> sys.stderr, "WARNING: Input files contain different or differently ordered reference sequences."
    print >> sys.stderr, "Results may not be as expected."
    print >> sys.stderr, "Ensure that input files have the same reference sequences in the same order."
