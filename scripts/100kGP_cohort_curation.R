library(readr)
library(dplyr)
library(stringr)
library(tidyverse)
library(Rlabkey)

tiered <- read_tsv('/re_gecip/enhanced_interpretation/mrintari/data/tiered_spliceogenic_missense_vars.tsv')



high_splice <- tiered %>% filter(SpliceAI >= 0.8)
write_tsv(high_splice,
          '/home/mrintari/re_gecip/enhanced_interpretation/mrintari/tiered_spliceogenic_missense_vars_0.8.tsv')

## ------- Access LabKey Tables Using API --------- ##

labkey.setCurlOptions(ssl_verifyhost = 2, ssl_verifypeer = TRUE,
                      NETRC_FILE = "/home/mrintari/re_gecip/enhanced_interpretation/mrintari/.netrc")
options(labkey.wafEncoding = FALSE)

labkey_to_df <- function(sql_query, database, maxrows){
  
  labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")
  
  labkey.executeSql(folderPath = database,
                    schemaName = "lists",
                    colNameOpt = "rname",
                    sql = sql_query,
                    maxRows = maxrows) %>%
    mutate(across(everything(), as.character))
}


labkey.setDefaults(apiKey="00bbf2bc111883df5f7b820abbb1c7240a6f84f9368ab5b20271a86a8ef81df5")
labkey.setDefaults(baseUrl="https://labkey-embassy.gel.zone/labkey/")


database <- "/main-programme/main-programme_v19_2024-10-31"

sql <- "SELECT * 
FROM tiering_data
WHERE Chromosome = '10' AND
Position = '113582002' AND 
Reference = 'C' AND
Alternate = 'T' AND
Participant_Type = 'Proband';"

query <- labkey_to_df(sql, database, 100000) %>% 
  select(participant_id, rare_diseases_family_id, sample_id, phenotype, 
         chromosome, position, reference, alternate, genotype,
         mode_of_inheritance, segregation_pattern, penetrance, tier,
         genomic_feature_hgnc, father_affected, mother_affected,
         full_brothers_affected, full_sisters_affected,
         participant_phenotypic_sex)


## Too many variants to manually check, Need to create ranked-priority system
## First, access all LabKey rows matching my tiered variants

# 1. Create lookup keys
tiered$CHROM <- gsub('chr', '', tiered$CHROM)

tiered <- tiered %>%
  mutate(lookup_key = paste(CHROM, POS, REF, ALT, sep = ":"))

# 2. Define function to process chunks
get_labkey_batch <- function(variant_keys) {
  # Format keys for SQL: '1:123:A:T','1:456:G:C',...
  formatted_keys <- paste0("'", variant_keys, "'", collapse = ",")
  
  # Note: Using || for string concatenation in LabKey SQL
  sql <- paste0("
    SELECT *
    FROM tiering_data
    WHERE (Chromosome || ':' || Position || ':' || Reference || ':' || Alternate) 
    IN (", formatted_keys, ")
    AND Participant_Type = 'Proband'
  ")
  
  # Use existing helper logic
  labkey.executeSql(
    baseUrl = "https://labkey-embassy.gel.zone/labkey/",
    folderPath = database,
    schemaName = "lists",
    colNameOpt = 'rname',
    sql = sql,
    maxRows = 100000
  )
}

# 3. Split keys into chunks of 500 and run
all_keys <- tiered$lookup_key
key_batches <- split(all_keys, ceiling(seq_along(all_keys) / 500))

# This will return a list of dataframes
results_list <- lapply(key_batches, get_labkey_batch)

# 4. Combine everything into one final dataframe
final_variants_df <- bind_rows(results_list) %>%
  mutate(across(everything(), as.character))


write_tsv(final_variants_df,
          "/home/mrintari/re_gecip/enhanced_interpretation/mrintari/spliceogenic_labkey_rows.tsv")


final_variants_df <- read_tsv('/re_gecip/enhanced_interpretation/mrintari/data/spliceogenic_labkey_rows.tsv',
                           col_types = cols(chromosome = col_character()))
final_variants <- final_variants_df %>% 
  filter(consequence_type == 'missense_variant')

# Count occurrences of each unique variant to find duplicate variants - variant
# in multiple probands
variant_counts <- final_variants %>% 
  group_by(chromosome, position, reference, alternate) %>% 
  tally(name = "occurrence_count") %>% 
  arrange(desc(occurrence_count))

tier2 <- final_variants %>% 
  filter(tier == 'TIER2')

t2_variant_counts <- tier2 %>% 
  group_by(chromosome, position, reference, alternate) %>% 
  tally(name = "occurrence_count") %>% 
  arrange(desc(occurrence_count))

# Get Scores for Tier 2 variants
t2_variant_counts$position <- as.numeric(t2_variant_counts$position)
t2_variant_counts <- t2_variant_counts %>% 
  left_join(tiered %>%
              select(CHROM, POS, REF, ALT, SpliceAI, REVEL, AlphaMissense),
            by = c("chromosome" = "CHROM",
                   "position" = "POS", 
                   "reference" = "REF",
                   "alternate" = "ALT"))


# Assess how many of the solved cases have any 'Spliceogenic Missense' variants
exit <- read_tsv('/re_gecip/enhanced_interpretation/mrintari/data/gmc_exit_questionnaire_2026-01-27_15-37-53.tsv')

exit <- exit %>%
  mutate(lookup_key = paste(chromosome, position, reference,
                            alternate, sep = ":"))

solved_splice <- exit %>% 
  semi_join(tiered, by = 'lookup_key') %>% 
  inner_join(tiered %>% select(lookup_key, SpliceAI, REVEL, AlphaMissense),
             by = "lookup_key") %>% 
  select(participant_id, family_id, additional_comments, acmg_classification,
         chromosome, position, reference, alternate, gene_name, SpliceAI, 
         REVEL, AlphaMissense, lookup_key)

# Assess deNovo variants
final_variants <- final_variants %>%
  mutate(lookup_key = paste(chromosome, position, reference,
                            alternate, sep = ":"))
deNovos <- final_variants %>% 
  filter(segregation_pattern == 'deNovo') %>% 
  inner_join(tiered %>% select(lookup_key, SpliceAI, REVEL, AlphaMissense),
             by = 'lookup_key') %>% 
  select(participant_id, rare_diseases_family_id, sample_id, phenotype, 
         chromosome, position, reference, alternate, genomic_feature_hgnc, 
         SpliceAI, REVEL, AlphaMissense, lookup_key) %>% 
  arrange(desc(SpliceAI))
  
panelApp <- read_lines('/re_gecip/enhanced_interpretation/mrintari/data/green_PanelApp_genes.txt')

my_deNovos <- deNovos %>% 
  filter(genomic_feature_hgnc %in% panelApp)

# None of the deNovo's are in relevant monoallelic PanelApp Genes

# Assessing Variants Marked as Uncertain
vucs <- read_tsv('re_gecip/enhanced_interpretation/mrintari/data/gmc_exit_questionnaire_VUCS_2026-01-30.tsv')
vucs <- vucs %>% 
  mutate(lookup_key = paste(chromosome, position, reference,
                            alternate, sep = ":"))

my_vucs <- vucs %>% 
  filter(case_solved_family != 'yes') %>% 
  semi_join(tiered, by = 'lookup_key') %>% 
  inner_join(tiered %>% select(lookup_key, SpliceAI, REVEL, AlphaMissense),
             by = "lookup_key") %>% 
  select(participant_id, family_id, case_solved_family, additional_comments,
         acmg_classification, chromosome, position, reference, alternate,
         gene_name, SpliceAI, REVEL, AlphaMissense, lookup_key)

# Assess Alt Hom variants
exit$participant_id <- as.character(exit$participant_id)
final_variants$participant_id <- as.character(final_variants$participant_id)
my_homz <- final_variants %>% 
  filter(genotype == 'alternate_homozygous') %>% 
  anti_join(exit, by = "participant_id")

# Filter for PanelApp Genes
my_homz <- my_homz %>% 
  filter(genomic_feature_hgnc %in% panelApp)

# Get Tally of Variant Frequency
my_homz_variant_counts <- my_homz %>% 
  group_by(lookup_key) %>% 
  tally(name = "occurrence_count") %>% 
  arrange(desc(occurrence_count))

my_homz_variant_counts <- my_homz_variant_counts %>% 
  left_join(tiered %>%
              select(lookup_key, SpliceAI, REVEL, AlphaMissense))


# Check ClinVar Spliceogenic Vars
clinvar <- read_tsv('/re_gecip/enhanced_interpretation/mrintari/data/ClinVar_Spliceogenic_Variants.tsv',
                    col_types=cols(Chromosome=col_character())) %>% select(-'...1')

tiered_clinvar <- final_variants %>% 
  filter(lookup_key %in% clinvar$GEL_ID)

rna_seq <- read_tsv('/re_gecip/enhanced_interpretation/mrintari/data/rnaseq_qc_metrics_2026-02-26_15-16-43.tsv')
rna_seq$participant_id <- as.character(rna_seq$participant_id)

clinvar_rna <- tiered_clinvar %>% 
  filter(participant_id %in% rna_seq$participant_id) %>% 
  inner_join(tiered %>% select(lookup_key, SpliceAI, REVEL, AlphaMissense),
             by = "lookup_key") %>% 
  select(-lastindexed, -modifiedby, -modified, -createdby, -created,
         -diimporthash, -container, -assembly, -ensembl_id, -so_term,
         -entityid) %>% 
  distinct(participant_id, .keep_all=TRUE) %>% 
  arrange(genomic_feature_hgnc) %>%
  select(participant_id, sample_id, lookup_key, phenotype, genomic_feature_hgnc,
         chromosome, position, reference, alternate, tier, genotype,
         mode_of_inheritance, segregation_pattern, penetrance,
         participant_phenotypic_sex, SpliceAI, REVEL, AlphaMissense) %>% 
  inner_join(rna_seq %>% select(participant_id, rna_folder_path),
             by="participant_id")

write_tsv(clinvar_rna,
          "/home/mrintari/re_gecip/enhanced_interpretation/mrintari/clinvar_rna.tsv")
