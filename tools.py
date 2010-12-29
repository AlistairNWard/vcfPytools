#!/usr/bin/python

def writeHeader (outputFile, v):
  outputFile.write( v.headerText ) if v.headerText != "" else None
  outputFile.write( v.headerInfoText ) if v.headerInfoText != "" else None
  outputFile.write( v.headerFormatText ) if v.headerFormatText != "" else None
  outputFile.write( v.headerTitles )
