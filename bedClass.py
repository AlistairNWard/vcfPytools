#!/usr/bin/python

import os.path
import sys

class bed:
  def __init__(self):
    self.numberTargets = 0

  def openBed(self, filename):
    if filename == "stdin":
      self.filehandle = sys.stdin
    else:
      exists = os.path.exists(filename)
      if exists == False:
        print >> sys.stderr, "Failed to find file: ",filename
        exit(1)

      self.filehandle = open(filename,"r")

# Get a bed record.

  def getRecord(self, line):
    self.numberTargets = self.numberTargets + 1
    self.ref = ""
    self.start = 0
    self.end = 0

# bed file should be 0-based, half-open, so the start coordinate
# must be that in the bed file plus one.

    entries = line.split("\t")
    self.ref = entries[0]
    try:
      self.start = int(entries[1]) + 1
    except:
      text = "start position need is not an integer"
      self.generalError(text, "start", entries[1])

    try:
      self.end = int(entries[2])
    except:
      text = "end position need is not an integer"
      self.generalError(text, "end", entries[2])

# Close the bed file.

  def closeBed(self, filename):
    self.filehandle.close()

# Define error messages for different handled errors.

  def generalError(self, text, field, fieldValue):
    print >> sys.stderr, "\nError encountered when attempting to read:"
    if field != "":
      print >> sys.stderr, "\t", field, ":             ", fieldValue
    print >> sys.stderr,  "\n", text
    exit(1)
