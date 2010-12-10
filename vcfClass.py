#!/usr/bin/python

import os.path
import sys

class vcf:
  def __init__(self):
    self.header = ""
    self.genotypes = False
    self.infoField = {}
    self.referenceSequences = {}
    self.chromosome = ""
    self.position = -1

  def openVcf(self, filename):
    if filename == "stdin":
      self.filehandle = sys.stdin
    else:
      exists = os.path.exists(filename)
      if exists == False:
        print "Failed to find file: ",filename
        exit(1)

      self.filehandle = open(filename,"r")

# Parse the vcf header.

  def parseHeader(self, filename):
    for line in self.filehandle:
      if line.startswith("##"):
        self.header = self.header + line
        tagValue = line.split("=",2)
      elif line.startswith("#"):
        self.header = self.header + line
        infoFields = line.rstrip("\n").split("\t")

# Strip the end of line character from the last infoFields entry.

        numberInfoFields = len(infoFields)
        if numberInfoFields > 8:
          if numberInfoFields - 9 == 1:
            print numberInfoFields - 9, " sample present in vcf file in: ", filename
          else:
            print numberInfoFields - 9, " samples present in vcf file in: ", filename
          self.samplesList = infoFields[9:]
          self.samplesList
          self.genotypes = True
        else:
          print "No samples present in the header."
          print "No genotype information available."
        break
      else:
        print "No header lines present."
        print "Terminating program."
        exit(1)

# Get the next line of information from the vcf file.

  def readRecord(self,line):
    vcfEntries      = line.rstrip("\n").split("\t")
    self.chromosome = vcfEntries[0]
    self.position   = int(vcfEntries[1])
    self.rsid       = vcfEntries[2]
    self.ref        = vcfEntries[3]
    self.alt        = vcfEntries[4]
    self.quality    = (vcfEntries[5])
    self.filters    = vcfEntries[6]
    self.info       = vcfEntries[7]

# Check for multiple alternate alleles.

    checkAlt = self.alt.split(",")
    if len(checkAlt) > 1:
      self.multiAllelic = True
    else:
      self.multiAllelic = False

# Read in string of genotype information, if genotypes are
# present in the file.

    if self.genotypes:
      self.genotypeFormat = vcfEntries[8].split(":")
      self.genotypes = vcfEntries[9:]

# Resolve all of the information in the info field and add to a
# dictionary.

  def resolveInfo(self,info):
    self.infoField = {}

    fields = info.split(";")
    for element in fields:
      info = element.split("=") if element.find("=") else element
      if len(info) == 1:
        self.infoField[info[0]] = True
      elif len(info) == 2:
        self.infoField[info[0]] = info[1]
      else:
        print "Unknown info string: ",element
        exit(1)

# Resolve genotype information.

  def resolveGenotypes(self,genotypeFormat,genotypes):
    print genotypes
    exit(0)

# Close the vcf file.

  def closeVcf(self, filename):
    self.filehandle.close()
