
# DOAC DDI CODEBOOK - ICD Code Definitions and Mapping for Outcome Variables

# This script defines the ICD-9 and ICD-10 codes for bleeding outcomes, trauma exclusion codes, and HAS-BLED codes. It also includes functions for forward/backward mapping using CMS GEMs

# Author: Kent
# Date: 3.28.25

#######################################

#################

# #Object vector
# object_oac <- c("Dabigatran Etexilate Mesylate", "Apixaban", "Rivaroxaban", "Warfarin Sodium", "Edoxaban")

#cont_enroll_req day requirement
cont_enroll_req <- 365

##################

# Section 1: Bleed Codes (ICD-9)

##################

# Definite bleed codes (Must be principal diagnosis) - codes "indicating" bleeding
gib_icd9_ind <- c("^5310.*", "^5312.*", "^5314.*", "^5316.*", "^5320.*", "^5322.*", "^5324.*", "^5326.*", "^5330.*", "^5332.*", "^5334.*",
                  "^5336.*", "^5340.*", "^5342.*", "^5344.*", "^5346.*", "^53501", "^53511", "^53521", "^53531", "^53541",
                  "^53551", "^53561", "^53783", "^4560", "^45620", "^5307", "^53021", "^53082", "^5780", "^4552", "^4555",
                  "^4558", "^56202", "^56203", "^56212", "^56213", "^56881", "^5693", "^56985", "^5781", "^5789")

gu_icd9_ind <- c("^59381", "^5997.*", "^6238", "^6266")

cerebral_icd9_ind <- c("^430", "^431", "^4320", "^4321", "^4329")

other_icd9_ind <- c("^4230", "^4590", "^56881", "^7191.*", "^7847", "^7848", "^7863.*")

# Combine into a single regex string
all_icd9_bleeds_ind <- paste(c(gib_icd9_ind, gu_icd9_ind, cerebral_icd9_ind, other_icd9_ind), collapse = "|")


# Possible Bleed Codes (suggestive codes requiring confirmation)
gib_icd9_possible <- c("^5311.*", "^5313.*", "^5315.*", "^5317.*", "^5319.*", "^5321.*",
                       "^5329.*", "^5331.*", "^5333.*", "^5335.*", "^5337.*",
                       "^5339.*", "^5341.*", "^5343.*", "^5345.*", "^5347.*", "^5349.*", "^53500",
                       "^53510", "^53520", "^53530", "^53540", "^53550", "^53560", "^56200",
                       "^56201", "^56210", "^56211", "^4550", "^4551", "^4553", "^4554", "^4556",
                       "^4557", "^4559", "^5301.*", "^53020")

all_gib_icd9_possible <- paste(gib_icd9_possible, collapse = "|")

unspec_icd9_possible <- c("^2851", "^2800", "^2859", "^79092", "^28749", "^2875")
all_unspec_icd9_possible <- paste(unspec_icd9_possible, collapse = "|")

gu_icd9_possible <- "^6262"

all_icd9_bleeds_possible <- paste(c(gib_icd9_possible, unspec_icd9_possible, gu_icd9_possible), collapse = "|")


# Additional Confirmatory Codes for Secondary Diagnoses
anemia_icd9 <- c("^2800" ,"^2851", "^2859")
orthostasis_icd9 <- "^4580"
syncope_icd9 <- "^7802"
all_comb_sec_icd9 <- paste(c(anemia_icd9, orthostasis_icd9, syncope_icd9), collapse = "|")


####################

# Section 2: Bleed Codes (ICD-10)

###################

# Definite bleed codes (Must be principal diagnosis) - codes "indicating" bleeding 

gib_icd10_ind <- c("^K250", "^K252", "^K254", "^K256", "^K260", "^K262", "^K264", "^K266",
                   "^K270", "^K272", "^K274", "^K276","^K280", "^K282", "^K284", "^K286", 
                   "^K2901", "^K2921", "^K2931", "^K2941", "^K2951", "^K2961", "^K2971", 
                   "^K2981", "^K2991", "^K31811", "^I8501", "^I8511","^K2211", "^K226", "^K920",
                   "^K5701", "^K5711", "^K5713", "^K5721", "^K5731", "^K5733", "^K5741", 
                   "^K5751", "^K5753", "^K5781", "^K5791", "^K5793", "^K661", "^K625", "^K5521", "^K921", "^K922")

gu_icd10_ind <- c("^N280", "^R310", "^R311", "^R3121", "^R3129", "^R319", "^N898", "^N921")

cerebral_icd10_ind <- c("^I6000", "^I6001", "^I6002", "^I6010", "^I6011", "^I6012", "^I602", "^I6030", 
                        "^I6031", "^I6032", "^I604", "^I6050", "^I6051", "^I606", "^I607", "^I608", 
                        "^I609", "^I610", "^I611", "^I612", "^I613", "^I614", "^I615", "^I616", "^I618", 
                        "^I619", "^I621", "^I6200", "^I6201", "^I6202", "^I6203", "^I629")

other_icd10_ind <- c("^I312", "^R58", "^K661", "^M2500", "^M25011", "^M25012", "^M25019", "^M25021", "^M25022", "^M25029", 
                     "^M25031", "^M25032", "^M25039", "^M25041", "^M25042", "^M25049", "^M25051", 
                     "^M25052", "^M25059", "^M25061", "^M25062", "^M25069", "^M25071", "^M25072", 
                     "^M25073", "^M25074", "^M25075", "^M25076", "^M2508", "^R040", "^R041", "^R042", 
                     "^R0481", "^R0489", "^R049")

all_icd10_bleeds_ind <- paste(c(gib_icd10_ind, gu_icd10_ind, cerebral_icd10_ind, other_icd10_ind), collapse = "|")


# Possible Bleed Codes (suggestive codes requiring confirmation)

gib_icd10_possible <- c("^K251", "^K253", "^K255", "^K257", "^K259","^K261", "^K263",
                        "^K265", "^K267", "^K269", "^K271", "^K273", "^K275", "^K277",
                        "^K279", "^K281", "^K283", "^K285", "^K287", "^K289", "^K2900",
                        "^K2920", "^K2930", "^K2940", "^K2950", "^K2960", "^K2970",
                        "^K2980", "^K2990", "^K5700", "^K5710", "^K5712", "^K5720",
                        "^K5730", "^K5732", "^K5740", "^K5750", "^K5752", "^K5780",
                        "^K5792", "^K640", "^K641", "^K642", "^K643", "^K644", "^K645",
                        "^K648", "^K649", "^K200", "^K208", "^K209", "^K210", "^K2210")
all_gib_icd10_possible <- paste(gib_icd10_possible, collapse = "|")

unspec_icd10_possible <- c("^D62", "^D500", "^D649", "^R791", "^D6959", "^D696")
all_unspec_icd10_possible <- paste(unspec_icd10_possible, collapse = "|")

all_icd10_bleeds_possible <- paste(c(gib_icd10_possible, unspec_icd10_possible), collapse = "|")

##Additional dx codes for rule in diagnoses

anemia_icd10 <- c("^D500" ,"^D62", "^D649")
orthostasis_icd10 <- "^I951"
syncope_icd10 <- "^R55"
all_comb_sec_icd10 <- paste(c(anemia_icd10, orthostasis_icd10, syncope_icd10), collapse = "|")


#####################

# Section 3: Trauma Exclusion Codes (ICD-9)

#####################

# Load trauma lookup table from codebook (sheet "trauma_icd9")
# The trauma lookup should include columns: Code, pattern, GI, GU, CNS, Other

trauma_lookup <- read_excel("C:/Users/kahanso2/Documents/doac-ddi/codes/doac_ddi_codebook.xlsx", sheet = "trauma_icd9")

# Create unified regex strings for trauma codes by site (ICD-9)
trauma_check_icd9 <- paste(trauma_lookup$pattern, collapse = "|")
trauma_check_icd9_GI <- paste(trauma_lookup %>% filter(GI == TRUE) %>% pull(pattern), collapse = "|")
trauma_check_icd9_GU <- paste(trauma_lookup %>% filter(GU == TRUE) %>% pull(pattern), collapse = "|")
trauma_check_icd9_CNS <- paste(trauma_lookup %>% filter(CNS == TRUE) %>% pull(pattern), collapse = "|")
trauma_check_icd9_Other <- paste(trauma_lookup %>% filter(Other == TRUE) %>% pull(pattern), collapse = "|")


######################

# Section 4: Trauma HCPCS Codes

######################

trauma_hcpcs <- c("^62000", "^62005", "^62010")
trauma_hcpcs_all <- paste(trauma_hcpcs, collapse = "|")


######################

# Section 5: Trauma Forward/Backward Mapping for ICD-10

#####################

# 1. Load the GEM Files (no decimals in codes)

forward_gem <- read.table("C:/Users/kahanso2/Documents/doac-ddi/data/2018_I9gem.txt", 
                          header = FALSE, stringsAsFactors = FALSE)
colnames(forward_gem) <- c("ICD9", "ICD10", "FLAG")

backward_gem <- read.table("C:/Users/kahanso2/Documents/doac-ddi/data/2018_I10gem.txt", 
                           header = FALSE, stringsAsFactors = FALSE)
colnames(backward_gem) <- c("ICD10", "ICD9", "FLAG")


# 2. Define Mapping Functions (using regex matching with perl=TRUE)
map_icd9_to_icd10 <- function(icd9_pattern, gem) {
  gem$ICD10[grepl(icd9_pattern, gem$ICD9, perl = TRUE)]
}

map_icd10_to_icd9 <- function(icd10_code, gem) {
  gem$ICD9[grepl(icd10_code, gem$ICD10, perl = TRUE)]
}

# 3. Define iterative mapping function that retains category flags
iterative_mapping_with_categories <- function(lookup, forward_gem, backward_gem) {
  
  # Start with the lookup rows.
  # Retain both the human-readable Code (source_code) and the regex pattern.
  current_icd9 <- lookup %>%
    mutate(source_code = Code,  # human-readable trauma code from the codebook
           source_pattern = pattern)  # regex pattern to use for matching
  
  # Data frame to accumulate forward mapping results.
  mapping_results <- data.frame(
    source_code = character(),
    ICD10  = character(),
    GI     = logical(),
    GU     = logical(),
    CNS    = logical(),
    Other  = logical(),
    iteration = integer(),
    stringsAsFactors = FALSE
  )
  
  final_icd10 <- character(0)
  iteration <- 1
  
  repeat {
    # -- Forward Mapping: For each current ICD9 row, find matching ICD10 codes using source_pattern.
    new_forward <- do.call(rbind, lapply(seq_len(nrow(current_icd9)), function(i) {
      pat <- current_icd9$source_pattern[i]
      icd10_matches <- map_icd9_to_icd10(pat, forward_gem)
      if (length(icd10_matches) > 0) {
        data.frame(
          source_code = current_icd9$source_code[i],  # retain the code from lookup
          ICD10     = icd10_matches,
          GI        = current_icd9$GI[i],
          GU        = current_icd9$GU[i],
          CNS       = current_icd9$CNS[i],
          Other     = current_icd9$Other[i],
          iteration = iteration,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))
    
    if (is.null(new_forward)) {
      new_forward <- data.frame(
        source_code=character(), ICD10=character(),
        GI=logical(), GU=logical(), CNS=logical(), Other=logical(),
        iteration=integer(), stringsAsFactors=FALSE
      )
    }
    
    # Retain only new ICD10 codes (not seen before)
    newly_found_icd10 <- setdiff(unique(new_forward$ICD10), final_icd10)
    if (length(newly_found_icd10) == 0) break
    
    final_icd10 <- unique(c(final_icd10, newly_found_icd10))
    mapping_results <- rbind(mapping_results, new_forward[new_forward$ICD10 %in% newly_found_icd10, ])
    
    # -- Backward Mapping: For each new ICD10, get ICD9 codes from backward GEM.
    backward_icd9 <- unique(unlist(lapply(newly_found_icd10, function(code) {
      map_icd10_to_icd9(code, backward_gem)
    })))
    
    # Build a union of all regex patterns from the lookup.
    pattern_union <- paste0("(?:", lookup$pattern, ")", collapse = "|")
    backward_icd9 <- backward_icd9[grepl(pattern_union, backward_icd9, perl = TRUE)]
    if (length(backward_icd9) == 0) break
    
    # For each backward ICD9 code, determine trauma Code (source_code) and flags from lookup.
    new_icd9 <- do.call(rbind, lapply(backward_icd9, function(code) {
      matches <- lookup[grepl(code, lookup$pattern, perl = TRUE), ]
      if (nrow(matches) > 0) {
        data.frame(
          source_code = paste(sort(unique(matches$Code)), collapse = "; "),
          GI     = any(matches$GI),
          GU     = any(matches$GU),
          CNS    = any(matches$CNS),
          Other  = any(matches$Other),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))
    if (is.null(new_icd9) || nrow(new_icd9)==0) break
    
    # Use the newly determined ICD9 (and associated Code) for the next iteration.
    current_icd9 <- new_icd9 %>%
      mutate(source_pattern = NA)  # source_pattern is not needed now; we already have source_code.
    
    iteration <- iteration + 1
  }
  
  # Aggregate mapping_results by ICD10.
  # For each ICD10, we combine the source_code values (the original trauma codes) via concatenation.
  final_mapping <- mapping_results %>%
    group_by(ICD10) %>%
    summarise(
      GI     = any(GI),
      GU     = any(GU),
      CNS    = any(CNS),
      Other  = any(Other),
      ICD9_list = paste(sort(unique(source_code)), collapse = "; ")
    ) %>%
    ungroup()
  
  final_mapping
}

# 4. Run the Iterative Mapping
# Run the mapping, referencing your trauma_lookup from Excel
final_mapping <- iterative_mapping_with_categories(trauma_lookup, forward_gem, backward_gem)

# (Optional) Add a caret prefix if desired
final_mapping$ICD10 <- ifelse(grepl("^\\^", final_mapping$ICD10, perl=TRUE),
                              final_mapping$ICD10,
                              paste0("^", final_mapping$ICD10))

# Add ICD-10 to excel file
# write.xlsx(final_mapping, file = "C:/Users/kahanso2/Documents/doac-ddi/codes/doac_ddi_codebook.xlsx", sheetName = "trauma_icd10", append = TRUE)

# Create unified regex strings for trauma codes by site (ICD-10)
trauma_check_icd10 <- paste(final_mapping$ICD10, collapse = "|")
trauma_check_icd10_GI <- paste(final_mapping %>% filter(GI == TRUE) %>% pull(ICD10), collapse = "|")
trauma_check_icd10_GU <- paste(final_mapping %>% filter(GU == TRUE) %>% pull(ICD10), collapse = "|")
trauma_check_icd10_CNS <- paste(final_mapping %>% filter(CNS == TRUE) %>% pull(ICD10), collapse = "|")
trauma_check_icd10_Other <- paste(final_mapping %>% filter(Other == TRUE) %>% pull(ICD10), collapse = "|")


#######################

# Section 6. Forward/backward mapping for HAS-BLED

######################
# 1. Load the HAS-BLED codebook from the Excel file.
#    The sheet "HAS-BLED" should contain columns: ABBR, CODE, ICD
hasbled_lookup <- read_excel("C:/Users/kahanso2/Documents/doac-ddi/codes/doac_ddi_codebook.xlsx", 
                             sheet = "has-bled_icd9")

# Define a function to extract codes based on HAS-BLED component and code type
get_codes <- function(df, categorization, code_type = "DIAG") {
  df %>%
    filter(Categorization == categorization, `Code Type` == code_type) %>%
    select(`ICD-9`) %>%
    pull()
}

# Now, create vectors for each component using the function

# Hypertension (H)
hasbled_htn_icd9 <- get_codes(hasbled_lookup, "H-Hypertension")

# Liver dysfunction (A)
hasbled_liver_icd9 <- get_codes(hasbled_lookup, "A-Liver dysfunction")

# Kidney dysfunction (A)
hasbled_kidney_icd9 <- get_codes(hasbled_lookup, "A-kidney dysfunction")

# Stroke/Systemic embolism (S)
hasbled_stroke_icd9 <- get_codes(hasbled_lookup, "S-Systemic embolism/stroke/TIA")

# Bleeding predisposition (B)
hasbled_bleed_icd9 <- get_codes(hasbled_lookup, "B-Bleeding predisposition")

# Heavy alcohol use or other drugs predisposing to bleeding (D)
# ICD-9 codes
hasbled_alc_icd9 <- get_codes(hasbled_lookup, "D-Heavy alcohol use or other drugs predisposing to bleeding", "DIAG")
# HCPCS codes (PROC)
hasbled_alc_proc <- get_codes(hasbled_lookup, "D-Heavy alcohol use or other drugs predisposing to bleeding", "PROC")

# 2. Forward backward mapping for HAS-BLED

# Helper functions using exact matching:
map_icd9_to_icd10_exact <- function(icd9_code, gem) {
  gem$ICD10[gem$ICD9 == icd9_code]
}

map_icd10_to_icd9_exact <- function(icd10_code, gem) {
  gem$ICD9[gem$ICD10 == icd10_code]
}

# Iterative mapping function for HAS-BLED codes using exact matching.
iterative_mapping_HASBLED_exact <- function(lookup, forward_gem, backward_gem) {
  # 'lookup' should have columns: categorization, ICD-9, and Code Type.
  # We assume all codes in the lookup are ICD-9 and complete.
  # (Procedure codes are not mapped, so ensure lookup has only DIAG codes.)
  
  # Start with the lookup codes.
  current_icd9 <- lookup %>%
    mutate(source_code = `ICD-9`) %>%
    select(source_code)
  
  # Data frame to accumulate forward mapping results.
  mapping_results <- data.frame(
    source_code = character(),
    ICD10 = character(),
    iteration = integer(),
    stringsAsFactors = FALSE
  )
  
  final_icd10 <- character(0)
  iteration <- 1
  
  repeat {
    # Forward mapping: for each current ICD-9 code, find matching ICD-10 codes.
    new_forward_list <- lapply(current_icd9$source_code, function(code) {
      icd10_matches <- map_icd9_to_icd10_exact(code, forward_gem)
      if (length(icd10_matches) > 0) {
        data.frame(
          source_code = code,
          ICD10 = icd10_matches,
          iteration = iteration,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    })
    new_forward <- do.call(rbind, new_forward_list)
    
    if (is.null(new_forward) || nrow(new_forward) == 0) break
    
    # Retain only new ICD-10 codes not seen before.
    newly_found_icd10 <- setdiff(unique(new_forward$ICD10), final_icd10)
    if (length(newly_found_icd10) == 0) break
    
    final_icd10 <- unique(c(final_icd10, newly_found_icd10))
    mapping_results <- rbind(mapping_results, new_forward[new_forward$ICD10 %in% newly_found_icd10, ])
    
    # Backward mapping: for each new ICD-10 code, find ICD-9 codes.
    backward_list <- lapply(newly_found_icd10, function(code) {
      map_icd10_to_icd9_exact(code, backward_gem)
    })
    backward_icd9 <- unique(unlist(backward_list))
    
    # Retain only those backward ICD-9 codes that are in the original lookup.
    valid_backward_icd9 <- intersect(backward_icd9, lookup$`ICD-9`)
    if (length(valid_backward_icd9) == 0) break
    
    # Prepare new ICD-9 codes for the next iteration.
    new_icd9_df <- data.frame(source_code = valid_backward_icd9, stringsAsFactors = FALSE)
    current_icd9 <- new_icd9_df
    
    iteration <- iteration + 1
  }
  
  # Aggregate mapping results by ICD-10 code.
  final_mapping <- mapping_results %>%
    group_by(ICD10) %>%
    summarise(
      HASBLED_ICD9 = paste(sort(unique(source_code)), collapse = "; "),
      iterations = max(iteration)
    ) %>%
    ungroup()
  
  final_mapping
}

# Create a lookup subset and mapping for each HAS-BLED component.
# Note: Only DIAG codes are used for mapping.

# 1. Hypertension (H)
lookup_htn <- hasbled_lookup %>%
  filter(Categorization == "H-Hypertension", `Code Type` == "DIAG")
mapped_htn <- iterative_mapping_HASBLED_exact(lookup_htn, forward_gem, backward_gem)
hasbled_htn_icd10 <- unique(mapped_htn$ICD10)

# 2. Liver Dysfunction (A)
lookup_liver <- hasbled_lookup %>%
  filter(Categorization == "A-Liver dysfunction", `Code Type` == "DIAG")
mapped_liver <- iterative_mapping_HASBLED_exact(lookup_liver, forward_gem, backward_gem)
hasbled_liver_icd10 <- unique(mapped_liver$ICD10)

# 3. Kidney Dysfunction (A)
lookup_kidney <- hasbled_lookup %>%
  filter(Categorization == "A-kidney dysfunction", `Code Type` == "DIAG")
mapped_kidney <- iterative_mapping_HASBLED_exact(lookup_kidney, forward_gem, backward_gem)
hasbled_kidney_icd10 <- unique(mapped_kidney$ICD10)

# 4. Stroke/Systemic Embolism (S)
lookup_stroke <- hasbled_lookup %>%
  filter(Categorization == "S-Systemic embolism/stroke/TIA", `Code Type` == "DIAG")
mapped_stroke <- iterative_mapping_HASBLED_exact(lookup_stroke, forward_gem, backward_gem)
hasbled_stroke_icd10 <- unique(mapped_stroke$ICD10)

# 5. Bleeding Predisposition (B)
lookup_bleed <- hasbled_lookup %>%
  filter(Categorization == "B-Bleeding predisposition", `Code Type` == "DIAG")
mapped_bleed <- iterative_mapping_HASBLED_exact(lookup_bleed, forward_gem, backward_gem)
hasbled_bleed_icd10 <- unique(mapped_bleed$ICD10)

# 6. Heavy Alcohol Use/Drugs Predisposing to Bleeding (D) - DIAG codes only
lookup_alc <- hasbled_lookup %>%
  filter(Categorization == "D-Heavy alcohol use or other drugs predisposing to bleeding", `Code Type` == "DIAG")
mapped_alc <- iterative_mapping_HASBLED_exact(lookup_alc, forward_gem, backward_gem)
hasbled_alc_icd10 <- unique(mapped_alc$ICD10)







##############################

# Section 7. Exclusion Lists for drugs, devices, and vitamins taht shouldnt cause DDI

excluded_mastfrm <- c("DEV", "CRE", "OIN", "LOT", "GEL", "EMU", "WAX", "TIN", "EMO", "FOA", "PAD", "PAS")
excluded_thrdtds <- c("Bulk Compounding Ingredient", "Vitamins, Prenatal, Misc Preps", 
                      "Vitamins W/Iron, Misc Preps.", "Vitamins W/Minerals, Misc Prep", 
                      "Vitamins, Plain, Misc Preps.", "Vitamins, Prenatal", "Amino Acid Supplements & Comb.", 
                      "Sterile Water", "Vitamin C & Combinations", "Polysaccharide-Iron Cmplx&Comb", "Lactobacillus & Comb.", 
                      "B-Complex/Iron/Vit. C & Comb.", "Skin/Mucous Membrane, Misc.", "Dietary Supplements, Misc.", "Iron Carbonyl & Comb.", 
                      "Iron Heme Polypeptide & Comb.", "Nutritional Supplements", "Opium Preparations", "Calcium Phosphate Dibasic", "Acidophilus & Comb.", 
                      "Hormones, Misc.", "Vitamin B Complex & Comb.", "B-Complex/Vitamin C", "Fish Oil & Comb.", "Omega-3 Fatty Acids", "Psyllium & Comb.", "Artificial Saliva, EENT")
excluded_gennme <- c("Antibacterial/Analgesic Combination", "Cough/Cold Combination", "Cyclobenzaprine HCl;Cream, Multi Ingredient")

# List of drugs to control for
nsaids <- c("Celecoxib", "Diclofenac", "Diflunisal", "Etodolac", "Fenoprofen Calcium", "Flurbiprofen", "Ibuprofen", 
            "Indomethacin", "Ketoprofen", "Ketorolac", "Meclofenamate Sodium", "Mefenamic Acid", "Meloxicam", "Nabumetone", 
            "Naproxen", "Oxaprozin", "Piroxicam", "Salsalate", "Sulindac", "Tolmetin Sodium")

antiplatelet <- c("Abciximab", "Anagrelide Hydrochloride", "Aspirin", "Cilostazol", "Clopidogrel", "Dipyridamole", "Eptifibatide", 
                  "Prasugrel Hydrochloride", "Ticagrelor", "Ticlopidine Hydrochloride", "Tirofiban Hydrochloride")

other_anticoag <- c("Fondaparinux Sodium", "Enoxaparin Sodium", "Argatroban", "Dalteparin", "Tinzaparin", "Heparin", "Bivalirudin", "Desirudin", "Phenprocoumon")

ssri_snri <- c("Citalopram", "Escitalopram", "Fluoxetine", "Fluvoxamine", "Milnacipran", "Paroxetine", "Sertraline")

giprotect <- c("Sucralfate", "Cimetidine", "Dexlansoprazole", "Esomeprazole", "Famotidine", "Lansoprazole", "Misoprostol", "Nizatidine", "Omeprazole", "Pantoprazole", "Rabeprazole", "Ranitidine")

anticoagulant <- c("Dabigatran Etexilate Mesylate", "Apixaban", "Rivaroxaban", "Warfarin Sodium", "Edoxaban")

full_anticoag <- c("Fondaparinux Sodium", "Enoxaparin Sodium", "Argatroban", "Dalteparin", "Tinzaparin", "Heparin", "Bivalirudin", "Desirudin", "Phenprocoumon", "Dabigatran Etexilate Mesylate", "Apixaban", "Rivaroxaban", "Warfarin Sodium", "Edoxaban" )

