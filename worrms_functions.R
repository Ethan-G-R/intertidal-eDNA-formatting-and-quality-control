#install.packages("taxize")

# 14.03.25
# Dina's script for pulling functional information for taxa from WORMS

library(taxize)
library(worrms)

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

# call the function
aphiaID_worms_output <- YOUR_TAXA_LIST %>%
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

# Call the function
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

# call function
worms_fgrp_output <- YOUR_UNIQUE_APHIAID_LIST %>%
  group_by(AphiaID) %>%
  do(get_worms_fgrp(AphiaID = .$AphiaID))
