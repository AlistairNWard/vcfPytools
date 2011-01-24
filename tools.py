#!/usr/bin/python

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
