# build a table of variants for testing wihtout having to query vcf every time
# columns are: impact (from vep), sfari gene score, gnomAD_AF, number of alts for every sample in the vcf
#
# Assumptions:
#   variants must have > 95% call rate among all samples
#   variants which vep classifies as unknown are skipped
#   variants which vep calls as in a gene not present in SFARI are skipped
#   variants which have no sfari gene score are assigned a score of 6 (worst possible)
#   variants which have no gnomAD allele frequency are assigned a value of 0

import os
import strutils
import strformat
import hts
import tables
import sets
import streams

# enumerate of impact levels
type Impact = enum
  UNKNOWN = -1
  LOW = 0
  MED = 1
  HIGH = 2
  
# define I/O
var 
  fh = commandLineParams()[0]
  o = open(fh[0..^5]&"_variant_table.txt", fmWrite)
  metadata = open(fh[0..^5]&"_variant_metadata.txt", fmWrite)

# define impact levels by VEP annotation as a set
var
  low_impacts = toSet(@["3_prime_UTR_truncation","non_canonical_start_codon","synonymous_variant","coding_sequence_variant","incomplete_terminal_codon_variant","stop_retained_variant","mature_miRNA_variant","5_prime_UTR_premature_start_codon_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","conserved_intron_variant","intron_variant","exon_variant","gene_variant","NMD_transcript_variant","non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_amplification","feature_elongation","miRNA","transcript_variant","start_retained","regulatory_region_variant","feature_truncation","non_coding_exon_variant","nc_transcript_variant","conserved_intergenic_variant","intergenic_variant","intergenic_region","intragenic_variant","non_coding_transcript_exon_variant","non_coding_transcript_variant","transcript","sequence_feature","non_coding"])
  med_impacts = toSet(@["disruptive_inframe_deletion","conservative_inframe_deletion","disruptive_inframe_insertion","conservative_inframe_insertion","duplication","inversion","exon_region","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","initiator_codon_variant","regulatory_region_ablation","5_prime_UTR_truncation","splice_region_variant"])
  high_impacts = toSet(@["chromosome_number_variation","transcript_ablation","exon_loss_variant","exon_loss","rare_amino_acid_variant","protein_protein_contact","structural_interaction_variant","feature_fusion","bidirectional_gene_fusion","gene_fusion","feature_ablation","splice_acceptor_variant","splice_donor_variant","start_retained_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification"])
  unknown_impacts = toSet(@["?","","UNKNOWN"])

# build dictionary for gene score lookup
# missing values are assigned 6 (worst)
var
  gene:string
  score:int
  gene_score_dict = initTable[gene, score]()
  gene_score_file = newFileStream("/scratch/ucgd/lustre/work/u0806040/data/sfari_gene_score_dict.txt", fmRead)
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

# run stats
var
  no_csq = 0
  valid_csq = 0
  unknown_impact = 0
  call_rate_fail = 0
  total_count = 0
  #gene_score_key_error = 0
  in_gnomad = 0
  no_gnomad = 0

# allocate memory
var
  csq:string = ""
  x = newSeq[int32]()

# write header to file
var t_samples = join(vcf.samples, "\t")
o.write_line(&"chrom\tend\tref\talt\tgene\timpact\tsfari_score\tgnomAD_AF\t{t_samples}")

# loop over variants in vcf
for variant in vcf:
  total_count += 1
  if variant.info.get("CSQ", csq) != Status.OK:
    no_csq += 1
    if no_csq mod 1000 == 0:
      stderr.write_line "CSQ field missing ", $no_csq, ", and with CSQ: ", $valid_csq
    continue

  valid_csq += 1
  
  # record gnomad allele freq
  var
    gnomad_af = newSeq[float32]()
    status = variant.info.get("gnomAD_AF", gnomad_af)
  if status == Status.UndefinedTag:
    quit "unknown gnomAD_AF field"
  if status != Status.OK:
    gnomad_af = @[0'f32]
    no_gnomad += 1

  if status == Status.OK:
    in_gnomad += 1

  # get number of alternate genotypes for each sample
  var alts = variant.format.genotypes(x).alts

  # filter based on call rate threshold
  var
    call_rate:float32
    called_counter, missing_counter = 0
    
  for i in alts:
    called_counter += 1
    if i == -1:
      missing_counter += 1
  call_rate = 1-(missing_counter/called_counter)
  
  if call_rate < 0.95:
    call_rate_fail += 1
    continue

  # parse the csq field
  # get sfari gene score assocaited with variant
  # as written, the last entry in the CSQ field that is in SFARI will be used
  # deal with overlapping genes by looping over annotations looking for sfari gene
  # deal with multiple impact annotations by spllitting on & and checking for highest impact
  var
    sfari_score = -1
    max_impact = Impact.UNKNOWN
    symbol = "empty"
  for a in csq.split(','):
    var
      asp = a.split('|', maxsplit = 4)
      this_symbol = asp[3]
    if this_symbol in gene_score_dict:
      sfari_score = gene_score_dict[this_symbol]
    else:
      continue
    
    for imp in asp[1].split('&'):

      if imp in low_impacts:
        if Impact.LOW > max_impact:
          symbol = this_symbol
        max_impact = max(max_impact, Impact.LOW)

      elif imp in med_impacts:
        if Impact.MED > max_impact:
          symbol = this_symbol
        max_impact = max(max_impact, Impact.MED)

      elif imp in high_impacts:
        if Impact.HIGH > max_impact:
          symbol = this_symbol
        max_impact = Impact.HIGH
      
      elif imp in unknown_impacts:
        continue
      else:
        quit "unknown impact:" & imp

  if sfari_score == -1:
    stderr.write_line "line with unknown genes"
    for a in csq.split(','):
      stderr.write_line a
    quit "debug"
  
  var impact = max_impact
  
  # skip variants of unknown impact
  if impact == Impact.UNKNOWN:
    unknown_impact += 1
    continue
  
  # write to file, tab separated (convert genotypes to tab separated string)
  var s_alts = newSeq[string](alts.len)
  for i, a in alts:
    s_alts[i] = $a
  var t_alts = join(alts, "\t")
  o.write_line(&"{variant.CHROM}\t{variant.POS}\t{variant.REF[0]}\t{variant.ALT[0]}\t{symbol}\t{impact}\t{sfari_score}\t{gnomad_af[0]}\t{talts}")

o.close()

# record metatdata
metadata.write_line("Infile = " & commandLineParams()[0])
metadata.write_line("Total variant count = " & $total_count)
metadata.write_line("Variants with CSQ field = " & $valid_csq)
metadata.write_line("Variants without CSQ field = " & $no_csq)
metadata.write_line("Variants failing call rate filter = " & $call_rate_fail)
metadata.write_line("Variants of unknown impact = " & $unknown_impact)
#metadata.write_line("Variants not associated with a SFARI gene = " & $gene_score_key_error)
metadata.write_line("Variants in gnomAD = " & $in_gnomad)
metadata.write_line("Variants not in gnomAD = " & $no_gnomad)
metadata.close()
