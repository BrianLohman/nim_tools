# hash table based method of liftover
# use UCSC dictionary of snps (hg38) and rsID to update genomic coordinates
#  UCSC downlaod:
#  requires formatting prior to use here: cut -f 2,4,5

# USAGE: nim c -d:release -r bim_lift_hg37_hg38.nim
# RETURNS: plink bim format file with hg38 coordintes, marking snps that fail to liftover
# ASSUMES: 1) snp dictionary is gzipped and tab separated

import os
import tables
import strutils
import zip/gzipfiles
import streams

type position = tuple
  chrom: string
  position: int

var
  lineCnt = 0
  line:string = ""
  tbl = initTable[string, position]()
  snps_dict = newGZFileStream("/scratch/ucgd/lustre/work/u0806040/data/hg38_all_snps150_chr_pos_rsID.txt.gz")

if isNil(snps_dict):
  quit "could not open snp dictionary"

# build has table of rsID (key) and hg38 positions (value = (chrom, bp))
while snps_dict.readline(line):
  lineCnt += 1
  var toks = line.split('\t')
  tbl[toks[3]] = (chrom:toks[0], position:parseInt(toks[2]))
  #if lineCnt mod 1000000 == 0:
  #  echo lineCnt
close(snps_dict)

# query the hash table with the rsID and return the tuple of position, parsing to the chr and bp fields
var
  bim_filename = commandLineParams()[0]
  line2:string = ""
  o = open(bim_filename[0 .. ^9]&"hg38.bim", fmWrite)
  bim_file = newFileStream(bim_filename, fmRead)

if isNil(bim_file):
  quit "could not open " & bim_filename

while bim_file.readline(line2):
  var toks = line2.split('\t')
  var rsID = toks[1]

# write to file: mark key errors with "fail_" for donwstream removal with plink
  if rsID in tbl:
    o.write([tbl[rsID][0], toks[1], "0", $tbl[rsID][1], toks[4], toks[5]&'\n'].join("\t"))
  else:
    o.write([toks[0], "fail_"&toks[1], "0", $toks[3], toks[4], toks[5]&'\n'].join("\t"))
close(bim_file)
close(o)
