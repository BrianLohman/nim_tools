# collect the gnomAD AF of all variants in a VCF
# formatted for input to R for plotting a histogram
#
# USAGE: gnomAD_AF_hist [vep annotated vcf with gnomAD AF]

# ASSUMES THAT THE GNOMAD ALLELE FREQUENCY HAS BEEN ADDED BY VCFANNO AND IS IN THE INFO FIELD

import hts
import os
import strutils
import sets

var
  vcf:VCF
  annotation:string = ""
  no_csq = 0
  found = 0
  not_in_gnomad = 0
  in_gnomad = 0
  fh = commandLineParams()[0]
  high_out = open(fh[0..^5]&"_high_impact_gnomAD_AF.txt", fmWrite)
  med_out = open(fh[0..^5]&"_med_impact_gnomAD_AF.txt", fmWrite)
  low_out = open(fh[0..^5]&"_low_impact_gnomAD_AF.txt", fmWrite)

type Impact = enum  
  UNKNOWN = -1
  LOW = 0  
  MED = 1  
  HIGH = 2  

# define impact levels by VEP annotation as a set
var low_impacts = toSet(@["3_prime_UTR_truncation","non_canonical_start_codon","synonymous_variant","coding_sequence_variant","incomplete_terminal_codon_variant","stop_retained_variant","mature_miRNA_variant","5_prime_UTR_premature_start_codon_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","conserved_intron_variant","intron_variant","exon_variant","gene_variant","NMD_transcript_variant","non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_amplification","feature_elongation","miRNA","transcript_variant","start_retained","regulatory_region_variant","feature_truncation","non_coding_exon_variant","nc_transcript_variant","conserved_intergenic_variant","intergenic_variant","intergenic_region","intragenic_variant","non_coding_transcript_exon_variant","non_coding_transcript_variant","transcript","sequence_feature","non_coding"])

var med_impacts = toSet(@["disruptive_inframe_deletion","conservative_inframe_deletion","disruptive_inframe_insertion","conservative_inframe_insertion","duplication","inversion","exon_region","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","initiator_codon_variant","regulatory_region_ablation","5_prime_UTR_truncation","splice_region_variant"])

var high_impacts = toSet(@["chromosome_number_variation","transcript_ablation","exon_loss_variant","exon_loss","rare_amino_acid_variant","protein_protein_contact","structural_interaction_variant","feature_fusion","bidirectional_gene_fusion","gene_fusion","feature_ablation","splice_acceptor_variant","splice_donor_variant","start_retained_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification"])

var unknown_impacts = toSet(@["?","","UNKNOWN"])

if not open(vcf, commandLineParams()[0], threads = 3):
  quit "couldn't open vcf"

for variant in vcf:
  if variant.info.get("CSQ", annotation) != Status.OK:
    no_csq += 1
    continue
  found += 1
  if found mod 10000 == 0:
    stderr.write_line $found & " variants processed"
  
  # interuption for debugging
  #if found == 10:
  #  break

  # for use when gnomAD_AF is in the INFO field:
  # variants with missing values are assigned 0
  var
    gnomad_af = newSeq[float32]()
    status = variant.info.get("gnomAD_AF", gnomad_af)

  if status == Status.OK:
    in_gnomad += 1

  if status != Status.OK:
    not_in_gnomad += 1
    gnomad_af = @[0'f32]

  var max_impact = Impact.UNKNOWN
  for a in annotation.split(","):
    var asp = a.split("|")
    for imp in asp[1].split('&'):
      if imp in unknown_impacts:
        continue
      elif imp in low_impacts:
        max_impact = max(max_impact, Impact.LOW)
      elif imp in med_impacts:
        max_impact = max(max_impact, Impact.MED)
      elif imp in high_impacts:
        max_impact = Impact.HIGH
      else:
        quit "unknown impact: " & imp
  
  var impact = max_impact
  if impact == Impact.UNKNOWN:
    continue
    
  elif impact == Impact.LOW:
    continue
    
  elif impact == Impact.LOW:
    low_out.write_line(gnomad_af[0])
    continue
  
  elif impact == Impact.MED:
    med_out.write_line(gnomad_af[0])
    continue

  elif impact == Impact.HIGH:
    high_out.write_line(gnomad_af[0])
    continue
  
  else:
    quit "impact not medium or high " & $impact

echo "In gnomAD: " & $in_gnomad
echo "Not in gnomAD: " & $not_in_gnomad

high_out.close()
med_out.close()
