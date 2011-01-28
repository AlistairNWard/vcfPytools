#!/usr/bin/python

import os.path
import sys

# Create a script for generating plots in R.  log is a variable determining
# whether to plot on a log scale.  log=0 means linear scales, log=1 has a
# log x-scale, log=1 a log y-scale and log=2 has both x and y log scales.

def createRScript(output, xLabel, log):
  file = "vcfPytoolsRScript.R"
  R = open(file, 'w')
  text = "pdf(\"" + output + "\")\n"
  print >> R, text
  print >> R, "data <- scan( \"Rdata\", list(0,0) )"
  if log == 0:
    print >> R, "plot( data[[1]], data[[2]], xlab=\"", xLabel, "\",ylab=\"Number of SNPs\",type=\"l\",col=\"orange\",axes=T,lwd=3)"
  elif log == 1:
    print >> R, "plot( data[[1]], data[[2]], xlab=\"", xLabel, "\",ylab=\"Number of SNPs\",type=\"l\",col=\"orange\",axes=T,lwd=3,log=\"x\")"
  elif log == 2:
    print >> R, "plot( data[[1]], data[[2]], xlab=\"", xLabel, "\",ylab=\"Number of SNPs\",type=\"l\",col=\"orange\",axes=T,lwd=3,log=\"y\")"
  elif log == 3:
    print >> R, "plot( data[[1]], data[[2]], xlab=\"", xLabel, "\",ylab=\"Number of SNPs\",type=\"l\",col=\"orange\",axes=T,lwd=3,log=\"xy\")"
  print >> R, "box()"
  print >> R, "dev.off()"
  R.close()

  return file
