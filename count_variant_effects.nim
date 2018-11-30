# Use a VEP annotated VCF file to build a table where rows are individuals and columsn are ID,
# predicted impact, and score (sum(genotype * SFARI gene score (flipped to 1 as worst and 6 best)))

# Both gnomAD allele frequency and call rate filters are applied

# Requires "sfari_gene_score_dict.txt", a two column file tab delimited file of gene symbol and sfari score

# USAGE: count_variant_effects [annoated vcf] 

import strutils
import os
import hts
import sets
import tables
import strformat
import streams

# set the field number of the VEP annotation
const impact_field = 1

# set the threshold for gnomAD allele freq
var gnomAD_AF_threshold = 0.1

# set the threshold for call_rate
const call_threshold = 0.95

# enumerate of impact levels
type Impact = enum
  LOW = 0
  MED = 1
  HIGH = 2
  UNKNOWN = 3

# define output variable
type key = tuple[sample:string, impact:Impact, gt:int8]
type score_key = tuple[sample:string, impact:Impact]

# define impact levels by VEP annotation as a set
var low_impacts = toSet(@["3_prime_UTR_truncation","non_canonical_start_codon","synonymous_variant","coding_sequence_variant","incomplete_terminal_codon_variant","stop_retained_variant","mature_miRNA_variant","5_prime_UTR_premature_start_codon_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","conserved_intron_variant","intron_variant","exon_variant","gene_variant","NMD_transcript_variant","non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_amplification","feature_elongation","miRNA","transcript_variant","start_retained","regulatory_region_variant","feature_truncation","non_coding_exon_variant","nc_transcript_variant","conserved_intergenic_variant","intergenic_variant","intergenic_region","intragenic_variant","non_coding_transcript_exon_variant","non_coding_transcript_variant","transcript","sequence_feature","non_coding"])
var med_impacts = toSet(@["disruptive_inframe_deletion","conservative_inframe_deletion","disruptive_inframe_insertion","conservative_inframe_insertion","duplication","inversion","exon_region","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","initiator_codon_variant","regulatory_region_ablation","5_prime_UTR_truncation","splice_region_variant"])
var high_impacts = toSet(@["chromosome_number_variation","transcript_ablation","exon_loss_variant","exon_loss","rare_amino_acid_variant","protein_protein_contact","structural_interaction_variant","feature_fusion","bidirectional_gene_fusion","gene_fusion","feature_ablation","splice_acceptor_variant","splice_donor_variant","start_retained_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification"])
var unknown_impacts = toSet(@["?","","UNKNOWN"])

# build dictionary for gene score lookup
var
  gene:string
  score:int
  gene_score_dict = initTable[gene, score]()
  gene_score_file = newFileStream("./sfari_gene_score_dict.txt", fmRead)
  dict_line:string = ""

if gene_score_file == nil:
  quit "couldn't find gene score file"

while gene_score_file.readline(dict_line):
  var toks = dict_line.split('\t')
  try:
    gene_score_dict[toks[0]] = parseInt(toks[1])
  except ValueError:
      gene_score_dict[toks[0]] = 6

close(gene_score_file)

# define input as vcf and check to see if it opens
var vcf:VCF
if not open(vcf, commandLineParams()[0], threads=3):
    quit "Couldn't open vcf"

# define output files
var output = open(commandLineParams()[0]&"_variant_counts_result.txt", fmWrite)
var score_output = open(commandLineParams()[0]&"_variant_scores_result.txt", fmWrite)

# initizlie tables
var tbl = newCountTable[key](8192)
var score_tbl = initTable[score_key, int](8192)

# define list of samples
var samples = vcf.samples

# for progress purposes
var found = 0
var missed = 0

# allocate memory
var ann:string = ""
var x = newSeq[int32]()

# count the number of variants without a gnomAD allele frequency
var no_gnomad:int = 0

# loop over variants in vcf file
for variant in vcf:
    if variant.info.get("CSQ", ann) != Status.OK:
        missed += 1
        if missed mod 500 == 0:
            stderr.write_line "didn't find CSQ field ", " missed:", missed, " found:", found
        continue
    found += 1
    # print status
    if found mod 10000 == 0:
        stderr.write_line "at:", $variant.CHROM, ":", $(variant.start + 1)
    
    # interruption for testing purposes, remove when format check complete
    #if found == 20000:
    #  break

    # check to see if a variant passes gnomAD filter
    var csq = ann.split("|")
    if len(csq[45]) == 0:
      no_gnomad += 1
    if len(csq[45]) > 0 and parseFloat(ann.split("|")[45]) > gnomAD_AF_threshold:
      continue

    # loop through annotation in the annotation field
    var impact: Impact
    for a in ann.split(","):
        var asp = a.split("|")
        # deal with multiple annotations and querey sets above 
        for imp in asp[impact_field].split('&'):
            if imp in low_impacts:
                impact =Impact.LOW 
            elif imp in high_impacts:
                impact = Impact.HIGH
            elif imp in med_impacts:
                impact = Impact.MED
            elif imp in unknown_impacts:
                impact = Impact.UNKNOWN
            else:
                quit "unknown impact:" & imp
    
    # skip variants of unknown impact
    if impact == Impact.UNKNOWN:
      continue

    # get number of alternate alleles
    var alts = variant.format.genotypes(x).alts

    # filter based on call rate threshold
    var call_rate:float32
    var called_counter, missing_counter = 0

    for i in alts:
      called_counter += 1
      if i == -1:
        missing_counter += 1
    call_rate = 1-(missing_counter/called_counter)

    if call_rate < call_threshold:
       continue
    
    # get sfari gene score associated with a variant
    var
      symbol:string = csq[3]
      weighted_score:int
    
    if symbol in gene_score_dict:
      weighted_score = (6 - gene_score_dict[symbol])
    else:
      continue

    # loop through samples for each variant
    for i, nalts in alts:
      if nalts == 0 or nalts == -1: 
        continue
      # older appraoch without weighting
      var k:key = (samples[i], impact, nalts)
      
      if k in tbl:
        tbl[k] += 1
      else:
        tbl[k] = 1
      # newer approach with weighted scores 
      var sc:score_key = (samples[i], impact)
      
      if not (sc in score_tbl):
        score_tbl[sc] = 0
      score_tbl[sc] += weighted_score * nalts.int

# write to file, looping over keys in table
for k, cnt in tbl:
  output.write_line(fmt"{k.sample}	{k.impact}	{k.gt}	{cnt}")

for k, score in score_tbl:
  score_output.write_line(fmt"{k.sample}	{k.impact}	{score}")
 
output.close()
score_output.close()

echo "Variants without a gnomAD allele frequency added by VEP " & $no_gnomad
