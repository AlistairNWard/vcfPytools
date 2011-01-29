#!/usr/bin/python

import os.path
import sys

class vcf:
  def __init__(self):
    self.hasHeader = True
    self.headerText = ""
    self.headerInfoText = ""
    self.headerFormatText = ""
    self.headerTitles = ""
    self.infoHeaderTags = {}
    self.formatHeaderTags = {}
    self.genotypes = False
    self.infoField = {}
    self.referenceSequences = {}
    self.referenceSequencesList = []
    self.referenceSequence = ""
    self.position = -1
    self.samplesList = []

    self.processInfo = False
    self.processGenotypes = False
    self.dbsnpVcf = False

# Open a vcf file.
  def openVcf(self, filename):
    if filename == "stdin": self.filehandle = sys.stdin
    else:
      try: self.filehandle = open(filename,"r")
      except IOError:
        print >> sys.stderr, "Failed to find file: ",filename
        exit(1)

# Parse the vcf header.
  def parseHeader(self, filename, writeOut):
    while self.getHeaderLine(filename, writeOut):
      continue

# Determine the type of information in the header line.
  def getHeaderLine(self, filename, writeOut):
    self.headerLine = self.filehandle.readline().rstrip("\n")
    if self.headerLine.startswith("##INFO"): success = self.headerInfo(writeOut, "info")
    elif self.headerLine.startswith("##FORMAT"): success = self.headerInfo(writeOut, "format")
    elif self.headerLine.startswith("##"): success = self.headerAdditional()
    elif self.headerLine.startswith("#"): success = self.headerTitleString(filename, writeOut)
    else: success = self.noHeader(filename, writeOut)

    return success

# Read information on an info field from the header line.
  def headerInfo(self, writeOut, lineType):
    if lineType == "info": self.headerInfoText = self.headerInfoText + self.headerLine + "\n"
    elif lineType == "format": self.headerFormatText = self.headerFormatText + self.headerLine + "\n"
    tag = self.headerLine.split("=",1)
    tagID = (tag[1].split("ID=",1))[1].split(",",1)

# Check if this info field has already been defined.
    if (lineType == "info" and self.infoHeaderTags.has_key(tagID[0])) or (lineType == "format" and self.formatHeaderTags.has_key(tagID[0])):
      print >> sys.stderr, "Info tag \"", tagID[0], "\" is defined multiple times in the header."
      exit(1)

# Determine the number of entries, entry type and description.
    tagNumber = (tagID[1].split("Number=",1))[1].split(",",1)
    tagType = (tagNumber[1].split("Type=",1))[1].split(",",1)
    try: tagDescription = ( ( (tagType[1].split("Description=\"",1))[1] ).split("\">") )[0]
    except IndexError: tagDescription = ""
    tagID = tagID[0]; tagNumber = tagNumber[0]; tagType = tagType[0]

# Check that the number of fields associated with the tag is either
# an integer or a '.' to indicate variable number of entries.
    if tagNumber == ".": tagNumber = "variable"
    else:
      try: tagNumber = int(tagNumber)
      except ValueError:
        print >> sys.stderr, "\nError parsing header.  Problem with info tag:", tagID
        print >> sys.stderr, "Number of fields associated with this tag is not an integer or '.'"
        exit(1)

    if lineType == "info": self.infoHeaderTags[tagID] = tagNumber, tagType, tagDescription
    if lineType == "format": self.formatHeaderTags[tagID] = tagNumber, tagType, tagDescription

    return True

# Read additional information contained in the header.
  def headerAdditional(self):
    self.headerText = self.headerText + self.headerLine + "\n"

    return True

# Read in the column titles to check that all standard fields
# are present and read in all the samples.
  def headerTitleString(self, filename, writeOut):
    self.headerTitles = self.headerLine + "\n"

# Strip the end of line character from the last infoFields entry.
    infoFields = self.headerLine.split("\t")
    if len(infoFields) > 8:
      if len(infoFields) - 9 == 1 and writeOut: print >> sys.stderr, numberInfoFields - 9, " sample present in vcf file: ", filename
      elif writeOut: print >> sys.stderr, numberInfoFields - 9, " samples present in vcf file: ", filename
      self.samplesList = infoFields[9:]
      self.genotypes = True
    elif len(infoFields) == 8:
      if writeOut: print >> sys.stderr, "No samples present in the header.  No genotype information available."
    else:
      print >> sys.stderr, "Not all vcf standard fields are available."
      exit(1)

    return False

# If there is no header in the vcf file, close and reopen the
# file so that the first line is avaiable for parsing as a 
# vcf record.
  def noHeader(self, filename, writeOut):
    if writeOut: print >> sys.stderr, "No header lines present in", filename
    self.hasHeader = False
    self.closeVcf(filename)
    self.openVcf(filename)

    return False

# Check that info fields exist.
  def checkInfoFields(self, tag):
    if self.infoHeaderTags.has_key(tag) == False:
      print >> sys.stderr, "Info tag \"", tag, "\" does not exist in the header."
      exit(1)

# Get the next line of information from the vcf file.
  def getRecord(self):
    self.record = self.filehandle.readline()
    if not self.record: return False

    vcfEntries      = self.record.rstrip("\n").split("\t")
    self.referenceSequence = vcfEntries[0]
    self.position   = int(vcfEntries[1])
    self.rsid       = vcfEntries[2]
    self.ref        = vcfEntries[3]
    self.alt        = vcfEntries[4]
    self.quality    = vcfEntries[5]
    self.filters    = vcfEntries[6]
    self.info       = vcfEntries[7]
    self.hasInfo    = True
    self.infoTags   = {}

    if len(vcfEntries) > 8:
      self.genotypeFormatString = vcfEntries[8]
      self.genotypeFormats = {}
      self.genotypes = vcfEntries[9:]
      self.hasGenotypes = True
      self.genotypeFields = {}
      self.phased = False
    else:
      self.hasGenotypes = False
      self.processGenotypes = False

# Add the reference sequence to the dictionary.  If it didn't previously
# exist append the reference sequence to the end of the list as well. 
# This ensures that the order in which the reference sequences appeared
# in the header can be preserved.
    if self.referenceSequence not in self.referenceSequences:
      self.referenceSequences[self.referenceSequence] = True
      self.referenceSequencesList.append(self.referenceSequence)

# Check for multiple alternate alleles.
    checkAlt = self.alt.split(",")
    if len(checkAlt) > 1: self.multiAllelic = True
    else: self.multiAllelic = False

# Resolve all of the information in the info field and add to a
# dictionary, if the information is to be read.
    if self.processInfo:

# First break the info string into its constituent elements.
      infoEntries = self.info.split(";")

# As long as some info fields exist, place them into a dictionary.
      if len(infoEntries) == 1 and infoEntries[0] == "": self.hasInfo = False
      else:
        self.hasInfo = True
        for entry in infoEntries:
          infoEntry = entry.split("=")
          if len(infoEntry) == 1: self.infoTags[infoEntry[0]] = True
          elif len(infoEntry) > 1: self.infoTags[infoEntry[0]] = infoEntry[1]

# Read in string of genotype information, if genotypes are
# present in the file and are requested.
    if self.processGenotypes:
      self.genotypeFormats = self.genotypeFormatString.split(":")

# Check that the number of genotype fields is equal to the number of samples
      if len(self.samplesList) != len(self.genotypes):
        text = "The number of genotypes is different to the number of samples"
        self.generalError(text, "", "")

# Add the genotype information to a dictionary.
      for i in range( len(self.samplesList) ):
        genotypeInfo = self.genotypes[i].split(":")
        self.genotypeFields[ self.samplesList[i] ] = {}

# Check that there are as many fields as in the format field.  If not, this must
# be because the information is not known.  In this case, it is permitted that
# the genotype information is either . or ./.
        if genotypeInfo[0] == "./." or genotypeInfo[0] == "." and len(self.genotypeFormats) != len(genotypeInfo): 
          self.genotypeFields[ self.samplesList[i] ] = "."
        else:
          if len(self.genotypeFormats) != len(genotypeInfo):
            text = "The number of genotype fields is different to the number specified in the format string"
            self.generalError(text, "sample", self.samplesList[i])

          for j in range( len(self.genotypeFormats) ): self.genotypeFields[ self.samplesList[i] ][ self.genotypeFormats[j] ] = genotypeInfo[j]

    return True

# Parse through the vcf file until the correct reference sequence is
# encountered and the position is greater than or equal to that requested.
  def parseVcf(self, referenceSequence, position, writeOut, outputFile):
    success = True
    if self.referenceSequence != referenceSequence:
      while self.referenceSequence != referenceSequence and success:
        if writeOut: outputFile.write(self.record)
        success = self.getRecord()

    while self.referenceSequence == referenceSequence and self.position < position and success:
      if writeOut: outputFile.write(self.record)
      success = self.getRecord()

    return success

# Get the information for a specific info tag.  Also check that it contains
# the correct number and type of entries.
  def getInfo(self, tag):
    result = []

# Check if the tag exists in the header information.  If so,
# determine the number and type of entries asscoiated with this
# tag.
    if self.infoHeaderTags.has_key(tag):
      infoNumber = self.infoHeaderTags[tag][0]
      infoType = self.infoHeaderTags[tag][1]
      numberValues = infoNumber

# First check that the tag exists in the information string.  Then split
# the entry on commas.  For flag entries, do not perform the split.
      if self.infoTags.has_key(tag):
        if numberValues == 0 and type(self.infoTags[tag]) == bool: result = True
        elif numberValues != 0 and type(self.infoTags[tag]) == bool:
          print >> sys.stderr, "ERROR"
          exit(1)
        else:
          fields = self.infoTags[tag].split(",")
          if len(fields) != numberValues:
            text = "Unexpected number of entries"
            self.generalError(text, "information tag", tag)

          for i in range(infoNumber):
            try: result.append(fields[i])
            except IndexError:
              text = "Insufficient values. Expected: " + self.infoHeaderTags[tag][0]
              self.generalError(text, "tag:", tag)
      else: numberValues = 0

    else:
      text = "information field does not have a definition in the header"
      self.generalError(text, "tag", tag)

    return numberValues, infoType, result

# Get the genotype information.
  def getGenotypeInfo(self, sample, tag):
    result = []
    if self.formatHeaderTags.has_key(tag):
      infoNumber = self.formatHeaderTags[tag][0]
      infoType = self.formatHeaderTags[tag][1]
      numberValues = infoNumber

      if self.genotypeFields[sample] == "." and len(self.genotypeFields[sample]) == 1:
        numberValues = 0
        result = "."
      else:
        if self.genotypeFields[sample].has_key(tag):
          if tag == "GT":
            if len(self.genotypeFields[sample][tag]) != 3 and len(self.genotypeFields[sample][tag]) != 1:
              text = "Unexected number of characters in genotype (GT) field"
              self.generalError(text, "sample", sample)

# If a diploid call, check whether or not the genotype is phased.
            elif len(self.genotypeFields[sample][tag]) == 3:
              self.phased = True if self.genotypeFields[sample][tag][1] == "|" else False
              result.append( self.genotypeFields[sample][tag][0] )
              result.append( self.genotypeFields[sample][tag][2] )
            elif len(self.genotypeFields[sample][tag]) == 3:
              result.append( self.genotypeFields[sample][tag][0] )
          else:
            fields = self.genotypeFields[sample][tag].split(",")
            if len(fields) != numberValues:
              text = "Unexpected number of characters in " + tag + " field"
              self.generalError(text, "sample", sample)

            for i in range(infoNumber): result.append(fields[i])
    else:
      text = "genotype field does not have a definition in the header"
      self.generalError(text, "tag", tag)

    return numberValues, result

# Parse the dbsnp entry.  If the entry conforms to the required variant type,
# return the dbsnp rsid value, otherwise ".".
  def getDbsnpInfo(self):

# First check that the variant class (VC) is listed as SNP.
    vc = self.info.split("VC=",1)
    if vc[1].find(";") != -1: snp = vc[1].split(";",1) 
    else:
      snp = []
      snp.append(vc[1])

    if snp[0].lower() == "snp": rsid = self.rsid
    else: rsid = "."

    return rsid

# Build a new vcf record.
  def buildRecord(self, removeGenotypes):
    record = self.referenceSequence + "\t" + \
                str(self.position) + "\t" + \
                self.rsid + "\t" + \
                self.ref + "\t" + \
                self.alt + "\t" + \
                self.quality + "\t" + \
                self.filters + "\t" + \
                self.info

    if self.hasGenotypes == True and not removeGenotypes:
      record = record + "\t" + self.genotypeFormatString
      for genotype in self.genotypes: record = record + "\t" + genotype

    record = record + "\n"

    return record

# Close the vcf file.
  def closeVcf(self, filename):
    self.filehandle.close()

# Define error messages for different handled errors.
  def generalError(self, text, field, fieldValue):
    print >> sys.stderr, "\nError encountered when attempting to read:"
    print >> sys.stderr, "\treference sequence : ", self.referenceSequence
    print >> sys.stderr, "\tposition :           ", self.position
    if field != "": print >> sys.stderr, "\t", field, ":             ", fieldValue
    print >> sys.stderr,  "\n", text
    exit(1)
