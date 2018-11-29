# Build a regions file from a list of genes and a gene models dictionary
# Intended to use with bcftool view for extraction of subset of genes from large vcf

# USAGE: ./build_regions [gene_models] [list of genes]

# ASSUMPTIONS: 
#     1. the gene models are tab separated and contain only a single entry for each gene with columns: chrom, strand, start, end, gene (following UCSC convention)
#     2. the list of genes we are interested in is from SFARI base, comma separated, and gene is in field 2

import os
import streams
import strutils
import tables

# read in the gene models as a table:
type position = tuple
  chrom: string
  start: int
  stop: int

var
  gene_models = commandLineParams()[0]
  line:string = ""
  tbl = initTable[string, position]()
  g = newFileStream(gene_models, fmRead)

if isNil(g):
  quit "coudln't open " & gene_models

# build a hash table of the gene namnes and genomic coordinates
while g.readline(line):
  var toks = line.split('\t')
  tbl[toks[4]] = (chrom:toks[0], start:parseInt(toks[2]), stop:parseInt(toks[3]))
close(g)

# query the hash table
var
  genes_of_interest = commandLineParams()[1]
  line2:string = ""
  output = open(genes_of_interest[0 .. ^5]&"_regions.txt", fmWrite)
  gl = newFileStream(genes_of_interest, fmRead)

if isNil(gl):
  quit "couldn't open " & genes_of_interest

while gl.readLine(line2):
  var toks = line2.split(',')
  var gene:string = toks[1]

  # write to file
  if gene in tbl:
    output.write([tbl[gene][0], $tbl[gene][1], $tbl[gene][2]&'\n'].join("\t"))
  else:
    continue

close(gl)
close(output)

