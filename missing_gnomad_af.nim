# count the number of variants that are missing the gnomad allele frequency annotation
import hts
import os
import strutils

var 
  vcf:VCF
  no_csq = 0
  found = 0
  annotation:string = ""
  no_gnomad = 0
  valid_gnomad = 0

if not open(vcf, commandLineParams()[0], threads=3):
  quit "Couldn't open vcf"

for variant in vcf:
  if variant.info.get("CSQ", annotation) != Status.OK:
    no_csq += 1
    continue
  found += 1
  if found mod 10000 == 0:
    stderr.write_line $found & " variants sreened"

  var csq = annotation.split("|")
  if len(csq[45]) == 0:
    no_gnomad += 1
  if len(csq[45]) > 0:
    valid_gnomad += 1

echo "Variants with a CSQ field" & $found
echo "Variants with NO CSQ field" & $no_csq
echo "Variants with a valid gnomad allele frequency" & $valid_gnomad
echo "Variants with NO gnomad allele frequency" & $no_gnomad
