#!/usr/bin/python

import os.path
import sys
import optparse

import vcfClass
from vcfClass import *

if __name__ == "__main__":
  main()

def main():

# Parse the command line options

  usage = "Usage: vcfTools.py union [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="append", type="string",
                    dest="vcfFiles", help="input vcf file")
  parser.add_option("-o", "--out",
                    action="store", type="string",
                    dest="output", help="output vcf file")

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.

  if options.vcfFiles == None:
    parser.print_help()
    print >> sys.stderr, "\nTwo input vcf files (-i) are required for performing union."
    exit(1)
  elif len(options.vcfFiles) != 2:
    print >> sys.stderr, "Two input vcf files are required for performing union."

# Set the output file to stdout if no output file was specified.

  if options.output == None:
    outputFile = sys.stdout
    writeOut = False
  else:
    outputFile = open(options.output, 'w')
    writeOut = True

  v1 = vcf() # Define vcf object.
  v2 = vcf() # Define vcf object.

# Read in the reference sequences present in the second vcf file.

  v2.openVcf(options.vcfFiles[1])
  v2.parseHeader(options.vcfFiles[1], False, False)
  for line2 in v2.filehandle:
    v2.getRecord(line2)
    v2.referenceSequences[ v2.referenceSequence ] = False
  vcfReferenceSequences = v2.referenceSequences
  v2.closeVcf(options.vcfFiles[1])

# Open the vcf files.

  v1.openVcf(options.vcfFiles[0])
  v2.openVcf(options.vcfFiles[1])

# Read in the header information.

  v1.parseHeader(options.vcfFiles[0], writeOut, True)
  v2.parseHeader(options.vcfFiles[1], writeOut, True)

# Check that the header for the two files contain the same samples.

  if v1.samplesList != v2.samplesList:
    print "vcf files contain different samples (or sample order)."
    exit(1)
  else:
    outputFile.write( v1.header )

# Get the first line of the second vcf file.

  for line2 in v2.filehandle:
    v2.getRecord(line2)
    break

# Calculate the intersection.

  for line1 in v1.filehandle:
    v1.getRecord(line1)

# Check if the records from the two vcf files correspond to the same
# reference sequence.  If so, search up to the same position and
# write out the record if it exists in the second vcf file.

    if v1.referenceSequence == v2.referenceSequence:
      vcfReferenceSequences[v1.referenceSequence] = True
      if v1.position == v2.position:
        outputFile.write( line1 )
      elif v1.position > v2.position:
        for line2 in v2.filehandle:
          v2.getRecord(line2)
          if v1.referenceSequence != v2.referenceSequence:
            outputFile.write( line2 )
            break
          if v2.position > v1.position:
            break
          elif v1.position == v2.position:
            outputFile.write( line2 )
            break
          else:
            outputFile.write( line2 )
            break

# If the reference sequence in the record from the first vcf file exists
# in the second, but has not been read yet, parse through the second
# vcf file until this reference sequence is reached, then search for the
# same position.

    elif vcfReferenceSequences[v1.referenceSequence] == False:
      for line2 in v2.filehandle:
        v2.getRecord(line2)
        vcfReferenceSequences[v2.referenceSequence] = True
        if v1.referenceSequence == v2.referenceSequence:
          if v1.position == v2.position:
            outputFile.write( line1 )
            break
          elif v1.position < v2.position:
            break
        else:
          outputFile.write( line2 )

# If the reference sequence in the record from the first vcf file exists
# in the second and has already been parsed, close and reopen the second
# vcf file, then allow the search to begin again from the beginning of the
# file.

    elif vcfReferenceSequences[v1.referenceSequence] == True:
      print >> sys.stderr, "WARNING: Hit reference sequence for second time."
      print >> sys.stderr, "Check that vcf file is sorted."
      exit(1)

# Close the vcf files.

  v1.closeVcf(options.vcfFiles[0])
  v2.closeVcf(options.vcfFiles[1])
