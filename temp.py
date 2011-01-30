#!/usr/bin/python

import os
import sys
import itertools
from math import factorial

def binomial(n, k):
  coeff = factorial(n)/( (factorial(n - k)) * factorial(k) )

  return coeff

if __name__ == "__main__":
  n=4
  files = []
  files = ["A", "B", "C", "D"]

  notList = []
  for i in range(0, n):
    notList.append(0)

  changeElement = 0
  for comb in range(0, n + 1):
    uniquePerms = {}
    print "\nNew set"
    coeff = binomial(n, comb)

    perm = itertools.permutations(notList, n)
    for permValue in perm:
      if not uniquePerms.has_key(permValue): uniquePerms[permValue] = 1

    for perm in uniquePerms.iterkeys():
      text = ""
      for i in range(0, n):
        if perm[i] == 0: text = text + files[i]
        else: text = text + " NOT " + files[i]

      print text

    try: notList[ changeElement ] = 1
    except IndexError: break
    changeElement = changeElement + 1
