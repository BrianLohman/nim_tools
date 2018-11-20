# convert genomic position in a plink bim file from chr:BP to rsID
# use UCSC dictionary of snps (hg37): http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp150.txt.gz
#   requires reformatting to remove extra columns prior to use here: cut -f 2,4,5

# USAGE: ./bim_chrBP_rsID [bim file]
# RETURNS: plink bim format file with rsID in addtion to genomic coordinates
# ASSUMES: 1) snp dictionary is gzipped
import os
import tables
import strutils
import zip/gzipfiles
import streams
import times

let time = cpuTime()

type position = tuple
   chrom: string
   position: int

var
  snps_filename = "/scratch/ucgd/lustre/work/u0806040/data/hg37_all_snps150_chr_pos_rsID.txt.gz"
  lineCnt = 0
  line:string = ""
  tbl = initTable[position, string]()
  f = newGZFileStream(snps_filename)

if isNil(f):
  quit "couldn't open " & snps_filename

# build hash table of hg37 positions (key: chrom:bp) and rsID (value)
while f.readline(line):
  var toks = line.split('\t')
  tbl[(chrom:toks[0], position:parseInt(toks[2]))] = toks[3]
  lineCnt += 1
  if lineCnt mod 1000000 == 0:
    echo lineCnt
close(f)

# query the hash table with the tuple of chr:BP from the bim file
var
  bim_filename = commandLineParams()[0]
  lineCnt2 = 0
  line2:string = ""
  o = open(bim_filename[0 .. ^5]&"_rsID.bim", fmWrite)
  f2 = newFileStream(bim_filename, fmRead)

if isNil(f2):
  quit "couldn't open " & bim_filename

while f2.readline(line2):
  var toks = line2.split('\t')
  var bim_position:position = (toks[0], parseInt(toks[3]))
  
# write to file: mark key erros with "fail_" for downstream removal with plink
  if bim_position in tbl:
    o.write([toks[0], tbl[bim_position], "0", toks[3], toks[4], toks[5]&'\n'].join("\t"))
  else:
    o.write([toks[0], "fail_"&toks[0]&":"&toks[3], "0", toks[3], toks[4], toks[5]&'\n'].join("\t"))
  lineCnt2 += 1
close(f2)
close(o)

echo "Run time = ",cpuTime() - time
