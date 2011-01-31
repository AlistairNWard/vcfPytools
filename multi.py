#!/usr/bin/python

import sys
import os.path
import optparse
import itertools
import subprocess

import vcfClass
from vcfClass import *

import tools
from tools import *

def main():

# Parse the command line options
  usage = "Usage: vcfPytools.py multi [options]"
  parser = optparse.OptionParser(usage = usage)
  parser.add_option("-i", "--in",
                    action="append", type="string", nargs=2,
                    dest="vcfFiles", help="input vcf files with identifier")

  (options, args) = parser.parse_args()

# Check that multiple vcf files are given.
  if not options.vcfFiles or len(options.vcfFiles) < 2:
    print >> sys.stderr, "At least two vcf files must be included."
    exit(1)

# Create a list of zero's with as many entries as there
# are vcf files.  This list will be used to determine the
# permutations of intersect and unique later in the program.
  rootPermutations = []
  outputFilenames = {}
  for i in range(0, len(options.vcfFiles)): rootPermutations.append(0)

# Calculate the different permutations and create scripts.
  for i in range(0, len(options.vcfFiles) ):
    uniquePermutations = {}
    permutations = itertools.permutations(rootPermutations, len(options.vcfFiles) )
    for permutation in permutations:
      permutationList = list(permutation)
      if not uniquePermutations.has_key(permutation):
        uniquePermutations[permutation] = 1

# Arrange the list so that the first element is a zero.  A one
# represents a NOT and thus the "unique" tool is required.  In
# using the "unique" tool, however, the first entry must the
# vcf file for whom the records are required, i.e. it cannot
# be the NOT file.  The actual output list is now the vcf files
# and not 1's and 0's.
        fileOrder = list(options.vcfFiles)
        if permutationList[0] == 1:
          original = fileOrder[0]
          for fileID in range(1, len(options.vcfFiles)):
            if permutationList[fileID] == 0:
              permutationList[fileID] = 1
              permutationList[0] = 0

              fileOrder[0] = fileOrder[fileID]
              fileOrder[fileID] = original
              break

# Now build the command.
        commandExec = "python " + sys.argv[0] + " "
        if permutationList[1] == 0: command = commandExec + "intersect --in " + fileOrder[0][0] + " --in " + fileOrder[1][0]
        else: command = commandExec + "unique --in " + fileOrder[0][0] + " --in " + fileOrder[1][0]
        if len(options.vcfFiles) > 2:
          for fileID in range(2, len(options.vcfFiles) ):
            if permutationList[fileID] == 0: command += " | " + commandExec + "intersect --in stdin --in " + fileOrder[fileID][0]
            else: command += " | " + commandExec + "unique --in stdin --in " + fileOrder[fileID][0]

# Generate the output filename.
        outputFile = "variantsFrom"
        fileID = 0
        for include in permutation:
          if include == 0: outputFile += "." + options.vcfFiles[fileID][1]
          fileID += 1
        command += " --out " + outputFile + ".vcf"

# Execute the command.
        success = subprocess.call(command, shell=True)
        if success != 0:
          print >> sys.stderr, "\nThe following command failed:\n"
          print >> sys.stderr, command
          exit(1)

        outputFilenames.setdefault(len(options.vcfFiles) - i, []).append(outputFile)

    rootPermutations[i] = 1

# Generate and execute commands to find the unions of each set of files.
# This means that there will exist a vcf file for each segment of the
# Wenn diagram as well as a vcf file including the union of all records
# present in two files etc.
  for key, value in outputFilenames.iteritems():
    unionCommand = "python " + sys.argv[0] + " union "
    fileID = 1
    for filename in value:
      if fileID < 3: unionCommand += "--in " + filename + ".vcf "
      else: unionCommand += " | python " + sys.argv[0] + " union --in stdin --in " + filename + ".vcf "
      fileID += 1

    outputFile = "variantsFrom" + str(key) + "only.vcf"
    unionCommand += "--out " + outputFile
    if key != len(options.vcfFiles):
      success = subprocess.call(unionCommand, shell=True)
      if success != 0:
        print >> sys.stderr, "\nThe following command failed:\n"
        print >> sys.stderr, command
        exit(1)

# Terminate the program cleanly.
  return 0
