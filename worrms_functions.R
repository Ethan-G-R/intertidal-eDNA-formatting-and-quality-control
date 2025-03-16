#install.packages("taxize")

# 14.03.25
# Dina's script for pulling functional information for taxa from WORMS

list.of.packages <- c("dplyr", 
                      "tidyverse", 
                      "phyloseq", 
                      "seqinr", 
                      "dada2", 
                      "sjmisc", 
                      "worrms",
                      "taxize",
                      "tibble", 
                      "taxadb",
                      "reshape2",
                      "decontam",
                      "microViz",
                      "cowplot",
                      "devtools",
                      "qiime2R", 
                      "microbiome",
                      "vegan",
                      "ggforce",
                      "ANCOMBC")
                      

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

invisible(lapply(list.of.packages, library, character.only = TRUE))

#Step 0 lets make a mock list of species

mock_sp_data <- data.frame(
   sp_name  = c("Aeolidia papillosa", "Galathea squamifera", "Hippolyte varians", 
                "Platynereis dumerilii", "Nicolea zostericola", "Ascidiella aspersa", 
                "Ophiothrix", "Enoplidae", "Heteronemertea", "Polycladida"))

library(writexl)
library(readxl)

ref_lib_taxa <- read_xlsx("trad_vs_eDNA_both_markers_table.xlsx", sheet = 3)
ref_lib_taxa <- ref_lib_taxa %>% rename(sp_name = lca_taxa_name_num_rm)

#Step 1: The first function called get_wormsid_noerror extracts unique identifiers for all taxon used in WoRMS called AphiaIDs.(15-mins)

get_wormsid_noerror <- function(sp_name){
  wm_id <- try(taxize::get_wormsid(sp_name,
                                   accepted = TRUE, 
                                   searchtype= "scientific", 
                                   marine_only = FALSE,
                                   ask = FALSE,
                                   row = 1,
                                   message = FALSE),
               silent = TRUE)
  
  #remove end word to attempt matching genus when previous no match
  
  if(class(wm_id) == "try-error"){
    sp_name <- str_remove(sp_name, "(\\s+\\w+)") 
    wm_id <- try(taxize::get_wormsid(sp_name,
                                     accepted = TRUE, 
                                     searchtype= "scientific", 
                                     marine_only = FALSE,
                                     ask = FALSE,
                                     row = 1,
                                     message = FALSE),
                 silent = TRUE)
  } 
  
  #convert non-matched species to NA AphiaID
  
  if(class(wm_id) == "try-error"){
    wm_id <- NA
  } 
  tibble(sciname = sp_name, aphiaID = as.double(unlist(wm_id)))
}

# Input is sp_name
# call the function
aphiaID_worms_output <- ref_lib_taxa$sp_name %>%
  purrr::map(get_wormsid_noerror, .progress = TRUE) %>%
  enframe() %>%
  unnest(cols = everything()) %>%
  select(-name)

# Step 2: The second function called get_wormsmeta accesses and formats any taxonomic meta data based on Aphia IDs (couple of minutes)

get_wormsmeta <- function(aphia_input){
  
  #split into smaller chunks for wm_record() to work
  taxadf_split <- split(aphia_input, (seq(nrow(aphia_input))-1) %/% 50) #split into smaller groups for worms to work
  temp_df <- data.frame()
  
  #run wm_record() through split list
  for (i in taxadf_split) {
    taxa_temp <- wm_record(id = i$aphiaID)
    temp_df = rbind(temp_df, taxa_temp)
  }
  tibble(temp_df)
}

# Input is 
# Call the function aphiaID_worms_output
worms_meta_output <- get_wormsmeta(aphiaID_worms_output)


#Step 3: The third function called get_worms_fgrp (originally developed here) accesses and formats information on broad taxonomic groupings of taxa based on Aphia IDs. (around 30 mins)

get_worms_fgrp <- function(AphiaID) {
  
  #' First, test if the species has attribute data
  attr_dat <- try(wm_attr_data(
    AphiaID, include_inherited = TRUE), silent = TRUE)
  
  #' set up out as null for later use
  out <- NULL
  
  if(!identical(class(attr_dat), "try-error")){
    
    #' if attribute data exists, test if functional group is there
    if("Functional group" %in% attr_dat$measurementType){
      fg_dat <- attr_dat %>% filter(measurementType == "Functional group")
      
      #' Insert if statement here about $children empty
      #' assign the $children - so that it can be used 
      children <- fg_dat$children
      
      if(max(lengths(children)) > 0 ){
        #' Extract the life stage information from the $children field
        life_stage <- children %>% bind_rows() %>%
          dplyr::select(measurementValue) %>%
          dplyr::rename(stage = measurementValue)
        
        #' add in rows for instances missing children
        idx <- which(lengths(children) == 0)
        if(length(idx) > 0){
          life_stage_null <- tibble(stage = rep(as.character(NA), length(children)))
          idy <- (1:length(children))[-idx]
          life_stage_null[idy,] <- life_stage
          life_stage <- life_stage_null
        }
        life_stage <- life_stage %>% bind_cols(., fg_dat)
        
        #' deal with cases where multiple records are returned for the same life stage:
        #' add a suffix to subsequent identical stages (adult_2, etc.)        
        life_stage <- life_stage %>% group_by(stage) %>% dplyr::mutate(nth_stage_val = 1:n()) %>% 
          ungroup() %>% 
          mutate(stage = case_when(
            nth_stage_val == 1 ~ stage,
            TRUE ~ paste(stage, nth_stage_val, sep = "_")
          )
          ) %>% 
          select(-nth_stage_val)
        
        #' create the output to return
        out <- tibble(AphiaID = as.numeric(life_stage$AphiaID),
                      stage = life_stage$stage, fun_grp = life_stage$measurementValue)
        
      } else{
        #' If no life stage info, assume stage is adult
        out <- tibble(AphiaID = as.numeric(fg_dat$AphiaID),
                      stage = "adult", fun_grp = fg_dat$measurementValue)
        
        #' Deal with cases where multiple records are returned (i.e. 2 or more adult fun_grps)
        if(nrow(out) > 1){
          out <- out %>% group_by(stage) %>% dplyr::mutate(nth_stage_val = 1:n()) %>% 
            ungroup() %>% 
            mutate(stage = case_when(
              nth_stage_val == 1 ~ stage,
              TRUE ~ paste(stage, nth_stage_val, sep = "_")
            )
            ) %>% 
            select(-nth_stage_val)
        }
      } 
    }
    
    #' add Pisces, from paraphyletic group: this takes priority over the above
    #' (e.g. class something as Pisces if it is a fish even if it is also listed as benthos)
    if ("Paraphyletic group" %in% attr_dat$measurementType) {
      if(first(
        attr_dat$measurementValue[attr_dat$measurementType == "Paraphyletic group"] == "Pisces")){
        out <- tibble(AphiaID = AphiaID, stage = "adult",
                      fun_grp = first(attr_dat$measurementValue[
                        attr_dat$measurementType == "Paraphyletic group"]))
      }
    }
  }
  #' add other paraphyletic groups: Algae, Algae > Macroalgae, Mangroves this DOES NOT takes priority over the above; i.e. only use if a functional group (e.g. phytoplankton) is not available
  if(is.null(out)){
    if (!identical(class(attr_dat), "try-error") && "Paraphyletic group" %in% attr_dat$measurementType) {
      out <- tibble(AphiaID = AphiaID, stage = "adult",
                    fun_grp = first(attr_dat$measurementValue[
                      attr_dat$measurementType == "Paraphyletic group"]))
    }
  }
  
  #' check taxonomy for other groups
  if(is.null(out)){
    taxo_dat <- try(wm_classification(AphiaID), silent = TRUE)
    if(identical(class(taxo_dat), "try-error")){
      fg <- as.character(NA)
    } else {
      fg <- case_when(
        "Aves" %in% taxo_dat$scientificname ~ "birds",
        "Mammalia" %in% taxo_dat$scientificname ~ "mammals",
        "Reptilia" %in% taxo_dat$scientificname ~ "reptiles",
        TRUE ~ as.character(NA)
      )
    }
    out <- tibble(AphiaID = AphiaID, stage = "adult", fun_grp = fg)
  }
  
  #' what if there are duplicate rows?
  out <- out[!duplicated(out), ]
  out <- out[!(is.na(out$stage)),] #gets rid of NA stages
  
  #' Replace 'not applicable' with NA
  out <- out %>% mutate_all(~ifelse(. == 'not applicable', NA, .))
  
  #' Keep only the life stage with the highest number (if applicable)
  out <- out %>%
    separate(stage, into = c("stage", "number"), sep = "_", remove = FALSE) %>%
    group_by(AphiaID, stage) %>%
    arrange(desc(number)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-number)
  
  #' Convert specific life stages to 'larva_other'
  out <- out %>% mutate(stage = case_when(
    stage %in% c('cerinula', 'larva > planula', 'medusa', 'zoea', 'nauplius', 'polyp', 'medusoid', 'ephyra', 'colony') ~ 'larva_other',
    TRUE ~ stage
  ))
  
  #' do some tidying of the output: tidy up functional_group text
  #' and spread to give one column per life stage
  out <- out %>% 
    distinct(AphiaID, stage, .keep_all = TRUE) %>%  # Add this line to remove duplicates
    mutate(functional_group = case_when(
      str_detect(fun_grp, ">") ~ tolower(word(fun_grp, -1)),
      fun_grp == "Pisces" ~ "fish",
      fun_grp == "Algae > Macroalgae" ~ "macroalgae",
      TRUE ~ tolower(fun_grp)
    )) %>%
    dplyr::select(-fun_grp) %>%
    spread(stage, functional_group)
  
  #' Change column name NA to other
  colnames(out)[colnames(out) == "NA"] <- "other"
  
  #' output the fg data
  out
}

# Input is worms_meta_output
# call function
worms_fgrp_output <- worms_meta_output %>%
  group_by(AphiaID) %>%
  do(get_worms_fgrp(AphiaID = .$AphiaID))

worms_fgrp_output

worms_fgrp_taxa <- full_join(worms_meta_output %>%
                               select(AphiaID, valid_name), 
                             worms_fgrp_output, by = "AphiaID")

worms_fgrp_taxa


################################################################################

get_worms_attr <- function(AphiaID) {
  
  #' Test if the species has attribute data
  attr_dat <- try(wm_attr_data(AphiaID, include_inherited = TRUE), silent = TRUE)
  
  #' Set up default tibble output
  out <- tibble(AphiaID = AphiaID, fun_grp = NA_character_, body_size = NA_character_, 
                feeding_type = NA_character_, ambi_group = NA_character_, stage = NA_character_)  
  
  if (!inherits(attr_dat, "try-error") && !is.null(attr_dat)) {
    
    #' Extract relevant attributes
    fg_dat <- attr_dat %>% filter(measurementType == "Functional group")
    size_dat <- attr_dat %>% filter(measurementType == "Body size (qualitative)")
    feeding_dat <- attr_dat %>% filter(measurementType == "Feedingtype")
    dev_dat <- attr_dat %>% filter(measurementType == "Development")
    ambi_dat <- attr_dat %>% filter(measurementType == "AMBI ecological group")
    
    #' Assign attributes if present
    fun_grp <- ifelse(nrow(fg_dat) > 0, fg_dat$measurementValue, NA_character_)
    body_size <- ifelse(nrow(size_dat) > 0, size_dat$measurementValue, NA_character_)
    feeding_type <- ifelse(nrow(feeding_dat) > 0, feeding_dat$measurementValue, NA_character_)
    dev_group <- ifelse(nrow(dev_dat) > 0, dev_dat$measurementValue, NA_character_)
    ambi_group <- ifelse(nrow(ambi_dat) > 0, ambi_dat$measurementValue, NA_character_)
    
    #' Extract stage information from the children column in Functional group data
    if (nrow(fg_dat) > 0 && !is.null(fg_dat$children)) {
      children <- fg_dat$children
      if(max(lengths(children)) > 0) {
        stage <- children %>% bind_rows() %>%
          dplyr::select(measurementValue) %>%
          dplyr::rename(stage = measurementValue)
        
        # Handling missing stages and creating unique stage labels if there are multiple stages
        stage <- stage %>% group_by(stage) %>% mutate(nth_stage_val = 1:n()) %>%
          ungroup() %>%
          mutate(stage = case_when(
            nth_stage_val == 1 ~ stage,
            TRUE ~ paste(stage, nth_stage_val, sep = "_")
          )) %>%
          select(-nth_stage_val)
      } else {
        stage <- tibble(stage = NA_character_)
      }
    } else {
      stage <- tibble(stage = NA_character_)
    }
    
    # Add the stage information to the output tibble
    out <- tibble(AphiaID = as.numeric(AphiaID), 
                  fun_grp = fun_grp, 
                  body_size = body_size, 
                  feeding_type = feeding_type, 
                  dev_group = dev_group, 
                  ambi_group = ambi_group,
                  stage = ifelse(nrow(stage) > 0, stage$stage, NA_character_))
  }
  
  #' Taxonomic fallback classification if no attributes were found
  if (all(is.na(out[-1]))) {  # Check if all columns except AphiaID are NA
    taxo_dat <- try(wm_classification(AphiaID), silent = TRUE)
    
    if (!inherits(taxo_dat, "try-error") && !is.null(taxo_dat)) {
      taxo_group <- case_when(
        "Aves" %in% taxo_dat$scientificname ~ "birds",
        "Mammalia" %in% taxo_dat$scientificname ~ "mammals",
        "Reptilia" %in% taxo_dat$scientificname ~ "reptiles",
        TRUE ~ NA_character_
      )
      out <- tibble(AphiaID = AphiaID, fun_grp = taxo_group, body_size = NA, feeding_type = NA, dev_group = NA, ambi_group = NA, stage = NA)
    }
  }
  
  #' Clean missing values and return distinct results
  out <- out %>%
    distinct(AphiaID, .keep_all = TRUE) %>%
    mutate(across(everything(), ~ ifelse(. == 'not applicable', NA, .)))
  
  return(out)
}

# Input is worms_meta_output
# call function
worms_attr_output <- worms_meta_output %>%
  group_by(AphiaID) %>%
  do(get_worms_attr(AphiaID = .$AphiaID))

worms_attr_taxa <- full_join(worms_meta_output %>%
                               select(AphiaID, valid_name), 
                             worms_attr_output, by = "AphiaID")
