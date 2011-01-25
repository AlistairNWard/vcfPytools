#!/usr/bin/python

import os.path
import sys

__author__ = "alistair ward"
__version__ = "version 0.11"
__date__ = "january 2011"

def main():

  usage = "Usage: vcfPytools.py [tool] [options]\n\n" + \
          "Available tools:\n" + \
          "  dbsnp:\n\tAnnotate the vcf file with dbsnp membership (requires dbsnp in vcf format).\n" + \
          "  filter:\n\tFilter the vcf file.\n" + \
          "  intersect:\n\tGenerate the intersection of two vcf files.\n" + \
          "  merge:\n\tMerge a list of vcf files.\n" + \
          "  sort:\n\tSort a vcf file.\n" + \
          "  stats:\n\tGenerate statistics from a vcf file.\n" + \
          "  union:\n\tGenerate the union of two vcf files.\n" + \
          "  unique:\n\tGenerate the unique fraction from two vcf files.\n" + \
          "  validate:\n\tValidate the input vcf file.\n\n" + \
          "vcfPytools.py [tool] --help for information on a specific tool."

# Determine the requested tool.

  if len(sys.argv) > 1:
    tool = sys.argv[1]
  else:
    print >> sys.stderr, usage
    exit(1)

  if tool == "dbsnp":
    import dbsnp
    dbsnp.main()
  if tool == "filter":
    import filter
    filter.main()
  if tool == "intersect":
    import intersect
    intersect.main()
  elif tool == "merge":
    import merge
    merge.main()
  elif tool == "sort":
    import sort
    sort.main()
  elif tool == "stats":
    import stats
    stats.main()
  elif tool == "union":
    import union
    union.main()
  elif tool == "unique":
    import unique
    unique.main()
  elif tool == "validate":
    import validate
    validate.main()
  elif tool == "--help" or tool == "-h" or tool == "?":
    print >> sys.stderr, usage
  else:
    print >> sys.stderr, "Unknown tool: ",tool
    print >> sys.stderr, "\n", usage
    exit(1)

if __name__ == "__main__":
  main()
