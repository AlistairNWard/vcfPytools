#!/usr/bin/python

import os.path
import sys

def createRHistScript(output):
  file = "vcfPytoolsRScript.R"
  R = open(file, 'w')
  text = "pdf(\"" + output + "\")\n"
  R.write( text )
  print >> R, "data <- scan( \"Rdata\", list(0,0) )"
  print >> R, "plot( data[[1]], data[[2]] )"
  print >> R, "dev.off()"

  return file
