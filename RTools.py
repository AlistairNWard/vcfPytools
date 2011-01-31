#!/usr/bin/python

import os.path
import sys

# Create a script for generating plots in R.  log is a variable determining
# whether to plot on a log scale.  log=0 means linear scales, log=1 has a
# log x-scale, log=1 a log y-scale and log=2 has both x and y log scales.

def createRScript(inputFile, tag):
  script = os.path.abspath(inputFile).rsplit(".",1)[0] + ".R"
  R = open(script, 'w')

# Plot the number of SNPs as a function of the tag.
  pdf = inputFile.rsplit(".",1)[0] + ".pdf"
  text = "pdf(\"" + pdf + "\")"
  print >> R, text
  print >> R, "data <- scan( \"" + inputFile + "\", list(0,0,0,0,0,0,0,0,0,0,0,0,0) )"
  print >> R, "colours<-c(\"orange\",\"darkred\",\"darkgreen\",\"deepskyblue\",\"black\")"
  print >> R, "total=data[[2]] + data[[3]] + data[[4]] + data[[5]]"
  print >> R, "maxNumber=max(total)"
  print >> R, "par(mar=c(5,4,4,4))"
  print >> R, "plot( data[[1]], data[[2]]/maxNumber, xlab=\"", tag, "\",ylab=\"Number of SNPs\",ylim=c(0,1),type=\"l\",col=\"orange\",axes=T,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], data[[3]]/maxNumber, xlab=\"\",ylab=\"\",ylim=c(0,1),type=\"l\",col=\"darkred\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], data[[4]]/maxNumber, xlab=\"\",ylab=\"\",ylim=c(0,1),type=\"l\",col=\"darkgreen\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], data[[5]]/maxNumber, xlab=\"\",ylab=\"\",ylim=c(0,1),type=\"l\",col=\"deepskyblue\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], total/maxNumber, xlab=\"\",ylab=\"\",ylim=c(0,1),type=\"l\",col=\"black\",axes=F,lwd=3)"
  print >> R, "box()"
  print >> R, "legend(\"topright\",c(\"Novel transitions\",\"Novel transversions\",\"Known transitions\",\"Known transversions\",\"Total\"),fill=colours)"
  print >> R, "dev.off()"

# Plot the ts/tv ratio for known, novel and all SNPs.
  pdf = inputFile.rsplit(".",1)[0] + ".less.tstv.pdf"
  text = "pdf(\"" + pdf + "\")"
  print >> R, text
  print >> R, "par(mar=c(5,4,4,4))"
  print >> R, "noveltstv=data[[6]]/data[[7]]"
  print >> R, "knowntstv=data[[8]]/data[[9]]"
  print >> R, "totaltstv=(data[[6]]+data[[8]])/(data[[7]]+data[[9]])"
  print >> R, "total=data[[6]]+data[[7]]+data[[8]]+data[[9]]"
  print >> R, "maxNumber=max(total)"
  print >> R, "plot( data[[1]], noveltstv, xlab=\"", tag, "\",ylab=\"Transitions/Transversions\",ylim=c(1.0,3.0),type=\"l\",col=\"orange\",axes=T,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], knowntstv, xlab=\"\",ylab=\"\",ylim=c(1.0,3.0),type=\"l\",col=\"darkred\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], totaltstv, xlab=\"\",ylab=\"\",ylim=c(1.0,3.0),type=\"l\",col=\"darkgreen\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], total/maxNumber, xlab=\"\",ylab=\"\",ylim=c(0.,1.),type=\"l\",col=\"deepskyblue\",axes=F,lwd=3)"
  print >> R, "axis(4)"
  print >> R, "mtext(\"Number of SNPs\",side=4,line=3)"
  print >> R, "box()"
  print >> R, "legend(\"bottomright\",c(\"Novel\",\"Known\",\"Total\",\"Percentage of SNPs\"),fill=colours)"
  print >> R, "dev.off()"

# Plot the ts/tv ratio for known, novel and all SNPs.
  pdf = inputFile.rsplit(".",1)[0] + ".greater.tstv.pdf"
  text = "pdf(\"" + pdf + "\")"
  print >> R, text
  print >> R, "par(mar=c(5,4,4,4))"
  print >> R, "noveltstv=data[[10]]/data[[11]]"
  print >> R, "knowntstv=data[[12]]/data[[13]]"
  print >> R, "totaltstv=(data[[10]]+data[[12]])/(data[[11]]+data[[13]])"
  print >> R, "total=data[[10]]+data[[11]]+data[[12]]+data[[13]]"
  print >> R, "maxNumber=max(total)"
  print >> R, "plot( data[[1]], noveltstv, xlab=\"", tag, "\",ylab=\"Transitions/Transversions\",ylim=c(1.0,3.0),type=\"l\",col=\"orange\",axes=T,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], knowntstv, xlab=\"\",ylab=\"\",ylim=c(1.0,3.0),type=\"l\",col=\"darkred\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], totaltstv, xlab=\"\",ylab=\"\",ylim=c(1.0,3.0),type=\"l\",col=\"darkgreen\",axes=F,lwd=3)"
  print >> R, "par(new=T)"
  print >> R, "plot( data[[1]], total/maxNumber, xlab=\"\",ylab=\"\",ylim=c(0.,1.),type=\"l\",col=\"deepskyblue\",axes=F,lwd=3)"
  print >> R, "axis(4)"
  print >> R, "mtext(\"Number of SNPs\",side=4,line=3)"
  print >> R, "box()"
  print >> R, "legend(\"bottomright\",c(\"Novel\",\"Known\",\"Total\",\"Percentage of SNPs\"),fill=colours)"
  print >> R, "dev.off()"
  R.close()

  return script
