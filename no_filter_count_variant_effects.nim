# Use a VEP annotated VCF file to build a table of the number of variants and thier impact for each individual
# count_variant_effects [annoated vcf] 
import strutils
import os
import hts
import sets
import tables
import strformat

# set the field number of the VEP annotation
const impact_field = 1

# set the threshold for gnomAD allele freq
#var gnomAD_AF_threshold = 0.01

# set the threshold for call_rate
const call_threshold = 0.95

# enumerate of impact levels
type Impact = enum
  UNKNOWN = -1
  LOW = 0
  MED = 1
  HIGH = 2

# define output variable
type key = tuple[sample:string, impact:Impact, gt:int8]

# define impact levels by VEP annotation as a set
var low_impacts = toSet(@["3_prime_UTR_truncation","non_canonical_start_codon","synonymous_variant","coding_sequence_variant","incomplete_terminal_codon_variant","stop_retained_variant","mature_miRNA_variant","5_prime_UTR_premature_start_codon_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","conserved_intron_variant","intron_variant","exon_variant","gene_variant","NMD_transcript_variant","non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_amplification","feature_elongation","miRNA","transcript_variant","start_retained","regulatory_region_variant","feature_truncation","non_coding_exon_variant","nc_transcript_variant","conserved_intergenic_variant","intergenic_variant","intergenic_region","intragenic_variant","non_coding_transcript_exon_variant","non_coding_transcript_variant","transcript","sequence_feature","non_coding"])
var med_impacts = toSet(@["disruptive_inframe_deletion","conservative_inframe_deletion","disruptive_inframe_insertion","conservative_inframe_insertion","duplication","inversion","exon_region","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","initiator_codon_variant","regulatory_region_ablation","5_prime_UTR_truncation","splice_region_variant"])
var high_impacts = toSet(@["chromosome_number_variation","transcript_ablation","exon_loss_variant","exon_loss","rare_amino_acid_variant","protein_protein_contact","structural_interaction_variant","feature_fusion","bidirectional_gene_fusion","gene_fusion","feature_ablation","splice_acceptor_variant","splice_donor_variant","start_retained_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification"])
var unknown_impacts = toSet(@["?","","UNKNOWN"])

# define input as vcf
var vcf:VCF

# check to see that vcf opens
if not open(vcf, commandLineParams()[0], threads=3):
    quit "Couldn't open vcf"

# define output file
var output = open(commandLineParams()[0]&"_variant_counts_result.txt", fmWrite)

# initizlie table
var tbl = newCountTable[key]()
var samples = vcf.samples

var ann : string = ""
var x = newSeq[int32]()

# for progress purposes
var found = 0
var missed = 0

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

    # check to see if passes gnomAD filter
    #var tmp = ann.split("|")
    #if len(tmp[45]) > 0 and  parseFloat(ann.split("|")[45]) > gnomAD_AF_threshold:
    #  continue

    # loop through annotation in the annotation field, 1
    var max_impact = Impact.UNKNOWN
    for a in ann.split(","):
        var asp = a.split("|")
        # deal with multiple annotations and quiery sets above 
        for imp in asp[impact_field].split('&'):
            if imp in low_impacts:
                max_impact =max(max_impact, Impact.LOW)
            elif imp in med_impacts:
                max_impact = max(max_impact, Impact.MED)
            elif imp in high_impacts:
                max_impact = Impact.HIGH
                break
            elif imp in unknown_impacts:
                continue
            else:
                quit "unknown impact:" & imp
    
    var impact = max_impact
    # skip variants of unknown impact
    if impact == Impact.UNKNOWN:
      continue
    var alts = variant.format.genotypes(x).alts

    # filter based on call rate threshold
    var call_rate: float32
    var called_counter, missing_counter = 0

    for i in alts:
      called_counter += 1
      if i == -1:
        missing_counter += 1
    call_rate = 1-(missing_counter/called_counter)

    if call_rate < call_threshold:
       continue

    for i, nalts in alts:
      if nalts == 0 or nalts == -1: continue

      var k:key = (samples[i], impact, nalts)
      if k in tbl:
        tbl[k].inc
      else:
        tbl[k] = 1

# write to file, looping over keys in table
for k, cnt in tbl:
  output.write_line(fmt"{k.sample}	{k.impact}	{k.gt}	{cnt}")
  
output.close()
