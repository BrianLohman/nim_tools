# Generate pseudocontrol F1s with the untransmitted parental allele
# 
# Input: VCF of trios, ped file
# Output: VCF of pseudocontrols

# Core functions provided by Brent Pedersen and modified for use here

import strutils
import tables
import hts
import os
import random 
randomize()

# define types
type
  Sample* = ref object
    family_id*: string
    id*: string
    paternal_id: string
    maternal_id: string
    mom*: Sample
    dad*: Sample
    sex: int
    affected: int
    kids: seq[Sample]
    i:int

proc `$`*(s:Sample): string =
  return format("Sample($1)", s.id)

proc parse_ped(path: string): seq[Sample] =
  result = new_seq_of_cap[Sample](10)

  var look = newTable[string,Sample]()

  for line in lines(path):
    if line[0] == '#':
      continue
    var toks = line.strip().split('\t')

    var s = Sample(family_id: toks[0], id: toks[1], kids:new_seq[Sample](), paternal_id: toks[2], maternal_id:toks[3], sex: parseInt(toks[4]), affected: parseInt(toks[5]))
    result.add(s)
    look[s.id] = s
  
  for s in result:
    if s.paternal_id in look:
      s.dad = look[s.paternal_id]
      s.dad.kids.add(s)
    if s.maternal_id in look:
      s.mom = look[s.maternal_id]
      s.mom.kids.add(s)

# define input variables
var
  vcf:VCF
  ped_path = commandLineParams()[1]
  probands = newSeq[Sample]()

if not open(vcf, commandLineParams()[0], threads = 3):
  quit "Couldn't open VCF"

# assign sample index
var si = newTable[string,int]()
for i, s in vcf.samples:
  si[s] = i

# get samples in same order as VCF
var samples = new_seq[Sample](len(si))
for sample in parse_ped(ped_path):
  sample.i = si[sample.id]
  samples[sample.i] = sample
  if sample.affected == 2:
    probands.add(sample)

stderr.write_line($len(probands) & " probands found")

# loop through the vcf
var
  g = newSeq[int32]()
  untransmitted = newSeq[Sample]()

for variant in vcf:
  var alts = variant.format.genotypes(g).alts
  
  for p in probands:
    var kid_a = alts[p.i]
    var mom_a = alts[p.mom.i]
    var dad_a = alts[p.dad.i]
    var un_a = mom_a + dad_a - kid_a
    
    if dad_a == -1 or mom_a == -1 or kid_a == -1:
      un_a = -1
    if dad_a == 1 and mom_a == 1 and kid_a == 1:
      var val = rand(3)
      if val == 0:
        un_a = 0
      elif val == 3:
        un_a = 2
      else:
        un_a = 1

    untransmitted.add(Sample(family_id: p.family_id, id: "un_"&p.family_id, kids:new_seq[Sample](), paternal_id: p.dad.id, maternal_id:p.mom.id, sex: 0, affected: 1))

#TODO: expand the hst-nim library to include writing new samples to an exisiting or new vcf file
# ??



