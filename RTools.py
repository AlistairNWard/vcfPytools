#!/usr/bin/python

import os.path
import sys

# Create a script for generating plots in R.  log is a variable determining
# whether to plot on a log scale.  log=0 means linear scales, log=1 has a
# log x-scale, log=1 a log y-scale and log=2 has both x and y log scales.

def createRScript(inputFile, tag):
  script = inputFile.rsplit(".",1)[0] + ".R"
  R = open(script, 'w')

# Plot the number of SNPs as a function of the tag.
  pdf = inputFile.rsplit(".",1)[0] + ".pdf"
  text = "pdf(\"" + pdf + "\")"
  print >> R, text
  print >> R, "data <- scan( \"" + inputFile + "\", list(0,0,0,0,0,0,0,0,0,0,0,0,0) )"
  print >> R, "colours<-c(\"orange\",\"darkred\",\"darkgreen\",\"deepskyblue\",\"black\")"
  print >> R, "total=data[[2]] + data[[3]] + data[[4]] + data[[5]]"
  print >> R, "par(mar=c(5,4,4,4))"
  print >> R, "plot( data[[1]], data[[2]], xlab=\"", tag, "\",ylab=\"Number of SNPs\",ylim=c(0,max(total)),type=\"l\",col=\"orange\",axes=T,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], data[[3]], xlab=\"\",ylab=\"\",ylim=c(0,max(total)),type=\"l\",col=\"darkred\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], data[[4]], xlab=\"\",ylab=\"\",ylim=c(0,max(total)),type=\"l\",col=\"darkgreen\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], data[[5]], xlab=\"\",ylab=\"\",ylim=c(0,max(total)),type=\"l\",col=\"deepskyblue\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], total, xlab=\"\",ylab=\"\",ylim=c(0,max(total)),type=\"l\",col=\"black\",axes=F,lwd=3)"
  print >> R, "box()"
  print >> R, "legend(\"topright\",c(\"Novel transitions\",\"Novel transversions\",\"Known transitions\",\"Known transversions\",\"Total\"),fill=colours)"

# Plot the ts/tv ratio for known, novel and all SNPs.
  pdf = inputFile.rsplit(".",1)[0] + ".tstv.pdf"
  text = "pdf(\"" + pdf + "\")"
  print >> R, text
  print >> R, "par(mar=c(5,4,4,4))"
  print >> R, "noveltstv=data[[2]]/data[[3]]"
  print >> R, "knowntstv=data[[4]]/data[[5]]"
  print >> R, "totaltstv=(data[[2]]+data[[4]]()/(data[[3]]+data[[5]])"
  print >> R, "plot( data[[1]], noveltstv, xlab=\"", tag, "\",ylab=\"Transitions/Transversions\",ylim=c(1.0,3.0),type=\"l\",col=\"orange\",axes=T,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], knowntstv, xlab=\"\",ylab=\"\",ylim=c(1.0,3.0),type=\"l\",col=\"darkred\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], totaltstv, xlab=\"\",ylab=\"\",ylim=c(1.0,3.0),type=\"l\",col=\"darkgreen\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], total, xlab=\"\",ylab=\"\",ylim=c(1.0,max(total)),log=\"y\",,type=\"l\",col=\"deepskyblue\",axes=F,lwd=3)"
  print >> R, "axis(4)"
  print >> R, "mtext(\"Number of SNPs\",side=4,line=3)"
  print >> R, "dev.off()"
  print >> R, "box()"
  print >> R, "legend(\"bottomright\",c(\"Novel\",\"Known\",\"Total\",\"Number of SNPs\"),fill=colours)"
  print >> R, "dev.off()"
  R.close()

  return script
