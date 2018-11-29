# Use family information to generate pseudocontrols which contain the untransmitted parental allele
# for every variant in the proband

# USAGE: ./generate_pseudocontrols [vcf file] [ped file] 

import hts
import os
import streams
import strutils

# pedigree info
var
  ped_filename = commandLineParams()[1] 
  ped_line:string = ""
  ped = newFileStream(ped_filename, fmRead)
  probands = newSeq[string]()

while ped.readline(ped_line):
  var toks = ped_line.split('\t')
  if parseInt(toks[5]) == 1:
    probands.add(toks[1])
  else:
    continue 

stderr.write_line $len(probands) & " probands identified by case status in ped file"

# define vcf input and check that it opens
var
  input_vcf:VCF

if not open(input_vcf, commandLineParams()[0], threads = 3):
  quit "Couldn't open vcf"

var
  samples = input_vcf.samples
  x = newSeq[int32]()

# loop through the variants
for variant in input_vcf:
  for sample in samples:
      if not (sample in probands):
        continue
      else:
        echo variant.genotypes
      break
  break



# open vcf for writing
#var wtr:VCF
#doAssert(open(wtr, "pseudocontrols_"&commandLineParams()[0], mode = "w"))
#wtr.header = vcf.header
#doAssert(wtr.write_header())

# var x = newSeq[int32]()
# alts = variant.format.genotypes(x).alts
