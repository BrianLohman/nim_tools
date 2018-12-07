# Generate pseudocontrol F1s with the untransmitted parental allele
# 
# Input: VCF of trios, ped file
# Output: VCF of pseudocontrols
# USAGE: pseudocontrols [ped] [vcf]
# 
# Core functions provided by Brent Pedersen and modified for use here
# Currently working. ~40 min for a single trio with 70M variants
# TODO: lines 130 and 141 produce warnings: cannot prove that type is initilized


import bpbiopkg/pedfile
import random
import hts
import os
import tables
randomize()

# define I/O
var
  args = commandLineParams()
  samples = parse_ped(args[0])
  vcf:VCF
  wtr:VCF

if not open(vcf, args[1]):
    quit "couldn't open VCF"

if not open(wtr, "added-control.vcf", mode="w"):
    quit "couldn't open VCF for writing"

# get samples in same order as vcf
samples = samples.match(vcf)
wtr.copy_header(vcf.header)

var sib_lookup = initTable[int,Sample]()
var kids = newSeq[Sample]()
for s in samples:
  if s.dad != nil and s.mom != nil:
    kids.add(s)
    wtr.add_sample("control_" & s.id)
    sib_lookup[wtr.n_samples - 1] = s

doAssert wtr.write_header
echo sib_lookup

# allocate memory
var
  x:seq[int32]
  ints = newSeq[int32]()
  floats = newSeq[float32]()
  strings = newSeq[string]()

# assign the number of alternate alleles for the pseudocontrols
proc sib_alts(v:Variant, i:int, sib:Sample, alts: seq[int8]): int =
  #@[0, 1, 0, 0] is of form: mom, dad, kid, pseudocontrol
  var kid_a = alts[sib.i]
  if kid_a == -1:
      return -1
  var mom_a = alts[sib.mom.i]
  var dad_a = alts[sib.dad.i]
  result = mom_a + dad_a - kid_a
  if result >= 2:
      return 2
  if result <= -2:
      return 0
  # special cases: missing genotypes or hets
  # randomly choose if all known genotypes are hets
  if dad_a == -1 or mom_a == -1 or kid_a == -1:
    result = -1
  if dad_a == 1 and mom_a == 1 and kid_a == 1:
    var val = rand(3)
    if val == 0:
      result = 0
    elif val == 3:
      result = 2
    else:
      result = 1

# table for genotype lookup
var gt_lookup = {0: @[2'i32, 2], 1: @[2'i32, 4], 2: @[4'i32, 4], -1: @[0'i32, 0'i32]}.toTable

# start with given kid and alter to pseudocontrol as necessary
proc make_from_sib(v:Variant, i:int, sib:Sample): bool =
  result = true
  if v.ALT.len > 1: return false
  var ad_type:BCF_TYPE
  var dp_type:BCF_TYPE
  for field in v.format.fields:
    if field.name == "AD":
        ad_type = field.vtype
    elif field.name == "DP":
        dp_type = field.vtype
    if field.vtype == BCF_TYPE.FLOAT:
      doAssert v.format.get(field.name, floats) == Status.OK
      for j in 0..<field.n_per_sample:
          floats[i*field.n_per_sample+j] = floats[sib.i*field.n_per_sample+j]
      doAssert v.format.set(field.name, floats) == Status.OK
    elif field.vtype == BCF_TYPE.CHAR:
      doAssert v.format.get(field.name, strings) == Status.OK
      strings[i] = ""
      doAssert v.format.set(field.name, strings) == Status.OK
    else:
      doAssert v.format.get(field.name, ints) == Status.OK
      for j in 0..<field.n_per_sample:
          ints[i*field.n_per_sample+j] = ints[sib.i*field.n_per_sample+j]
      doAssert v.format.set(field.name, ints) == Status.OK
      if field.name == "AD" and ints.len != 2*v.n_samples: return false

  var alts = v.format.genotypes(x).alts
  var kid_alts = v.sib_alts(i, sib, alts)
  alts[i] = kid_alts.int8
  var gt_ints = gt_lookup[kid_alts]
  doAssert v.format.get("GT", ints) == Status.OK
  ints[i*2] = gt_ints[0]
  ints[i*2+1] = gt_ints[1]
  doAssert v.format.set("GT", ints) == Status.OK

  doAssert v.format.get("AD", ints) == Status.OK
  if kid_alts == 0:
    ints[i*2] += ints[i*2+1]
    ints[i*2+1] = 0
  if kid_alts == 2:
    ints[i*2+1] += ints[i*2]
    ints[i*2] = 0
  if kid_alts == 1:
    if alts[sib.i] == 1:
        discard
    else:
      var dp = ints[i*2] + ints[i*2+1]
      ints[i*2] = int32(dp / 2)
      ints[i*2+1] = int32(dp / 2)
  if ad_type == INT8:
      for i, v in ints.mpairs:
          if v > 127: v = 127

  doAssert v.format.set("AD", ints) == Status.OK
  var dp = newSeq[int32](int(ints.len/2))
  #echo dp.len, " ", ints.len
  for i in 0..<dp.len:
      dp[i] = ints[2*i] + ints[2 * i + 1]
      if dp_type == INT8:
        if dp[i] > 127: dp[i] = 127

  var status = v.format.set("DP", dp)
  assert status == Status.OK

# loop through the ped and vcf applying the functions above
for v in vcf:
    #echo v
    v.vcf = wtr

    for i in samples.len..v.n_samples:
        if not make_from_sib(v, i, sib_lookup[i]): continue
    doAssert wtr.write_variant(v)
