# DOCUMENTATION -----------------------------------------------------------

#   Author:     Jesse Bruhn
#   Contact:    jesse_bruhn@brown.edu

# PREAMBLE ----------------------------------------------------------------

#print date
Sys.Date()

#Clear workspace
rm(list=ls())

#Load Packages
pckgs <- c("tidyverse")
lapply(pckgs, library, character.only=TRUE)

#Set Options
options(scipen=100)

#Clean Workspace
rm(pckgs)

#Session Info for Reproducibility
sessionInfo()

# NOTES -------------------------------------------------------------------

#this script takes the item response data from DESE and cleans it. 

#INTERESTING: Some of htese files have survey responses at the end. See 
#             Q1 - Q25 in teh 2001 data. I think these are survey responses

# CREATE SOME USEFUL FUNCTIONS --------------------------------------------



dataLoader <- function(grade, year, demographics = "nd", readFunction = read_csv, extension="txt"){
  
  file_path <- str_glue("raw_data/mcas_anonymized/g{grade}{year}{demographics}.{extension}")
  
  return(readFunction(file_path,
                      col_types = cols(.default = col_character())
  )
  
  )
  
}



# LOAD DATA ---------------------------------------------------------------

#vectors to poulate loops
grades <- c("03", "04", "05", "06", "07", "08", "10")
years <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", 
           "12", "13", "14", "15", "16", "17", "18", "19")
demographics <- c("d", "nd")

#list to store data frames
for(d in demographics){
  
  
  for(y in years){
    
    y_d <- str_glue("{y}_{d}")
    print(y_d)
    
    for(g in grades){
      
      #make some changes to account for idiosyncratic file names
      if(y =="15"){
        ext <- "dat"
        r_func <- read_tsv
      }else if(y %in% c("01", "02", "03", "04", "05")){
        ext <- "txt"
        r_func <- read_csv
      }else{
        ext <- "txt"
        r_func <- read_tsv
      }
      
      
      # else if(y %in% c("12", "13", "14", "16") & d == "nd" & g == "10"){
      #   ext <- "txt"
      #   r_func <- read_tsv
      # }
      
      
      if(y %in% c("08", "09", "10", "11", 
                  "12", "13", "14", "15") & 
         g == "10"){
        g <- "hs"
        
      }else if(y == "07" & g == "10"){
        
        g <- "HS"
        
      }
      
      #for some reason I do not have grade 05 in 2001 or 2002. At some point
      #I should check that I didn't just forget to download these from the 
      #website
      if(y %in% c("01", "02") &  g == "05"){
        next
      }
      
      print(g)
      
      if(g == "03"){
        
        df <- dataLoader(grade = g, 
                         year = y, 
                         demographics = d,
                         readFunction = r_func,
                         extension = ext)
        
        # if(dim(dat)[2]==1){
        #   break
        # }
        
        next
        
      }else{
        
        df <- df %>%
          bind_rows(dataLoader(grade = g, 
                               year = y, 
                               demographics = d,
                               readFunction = r_func,
                               extension = ext) 
          )  
        
      }
      
      
      
      
      
      
    }
    
    file_path <- str_glue("temp_files/{d}_{y}.rda")
    save(df, file = file_path)
    
  }
  
}





# STANDARDIZE VARIABLE NAMES ----------------------------------------------

#function to help rename variables
cleanVarnames <- function(inpath, y){
  
  load(inpath)
  
  if(y %in% c("2001", "2002")){
    df <- df %>%
      rename_with(str_to_lower) %>%
      rename(year = adminyear, 
             grade = grade, 
             gender = gender, 
             spec_ed = sped_off, 
             lep = lep_off, 
             frpl = frpl_off, 
             race = race_off, 
             migrant = migr_off, 
             e_test_stat = eteststa, 
             e_raw_score = erawsc, 
             m_test_stat = mteststa, 
             m_raw_score = mrawsc,
      ) %>%
      select(year, grade, gender, spec_ed, lep, frpl, race, migrant, e_test_stat, e_raw_score, 
             m_test_stat, m_raw_score, contains("wp"), contains("item"))  
  }else if(y %in% c("2003", "2004")){
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off, 
        lep = lep_off, 
        frpl = frpl_off, 
        race = race_off, 
        migrant = migr_off, 
        e_test_stat = eteststa, 
        e_raw_score = erawsc, 
        m_test_stat = mteststa, 
        m_raw_score = mrawsc,
        s_test_stat = steststa, 
        s_raw_score = srawsc 
        
      ) %>%
      select(year, grade,  gender, spec_ed, lep, frpl, race, migrant, e_test_stat, e_raw_score, 
             m_test_stat, m_raw_score, s_test_stat, s_raw_score, contains("wp"), contains("item")) 
    
  }else if(y %in% c("2005")){
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off, 
        lep = lep_off, 
        frpl = freelunch, 
        race = race_off, 
        migrant = migrant_of, 
        e_test_stat = eteststat, 
        e_raw_score = erawsc, 
        m_test_stat = mteststat, 
        m_raw_score = mrawsc,
        s_test_stat = steststat, 
        s_raw_score = srawsc 
        
      ) %>%
      select(year, grade,  gender, spec_ed, lep, frpl, race, 
             migrant, e_test_stat, e_raw_score, 
             m_test_stat, m_raw_score, s_test_stat, 
             s_raw_score, contains("wp"), contains("item")) 
    
    
  }else if(y %in% c("2006")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off,
        lep = lep_off, 
        frpl = freelunch_off, 
        race = race_off, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        migrant = migrant_off,
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_raw_score = erawsc, 
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_raw_score = mrawsc,
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_raw_score = srawsc 
        
      ) %>%
      select(year, grade, gender,  gender, spec_ed, lep, frpl, race, 
             years_in_ma, plan_504, first_language, disability_type, 
             need_level, migrant, 
             e_test_stat, e_alt_assesment, e_raw_score, 
             m_test_stat, m_alt_assesment, m_raw_score,  
             s_test_stat, s_alt_assesment, s_raw_score, contains("wp"), contains("item")) 
    
  }else if(y %in% c("2007", "2008")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off,
        spec_ed_placement = spedplacement,
        lep = lep_off, 
        frpl = freelunch_off, 
        race = race_off, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        migrant = migrant_off,
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_raw_score = erawsc, 
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_raw_score = mrawsc,
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_raw_score = srawsc 
        
      ) %>%
      select(year, grade, gender,  gender, spec_ed, spec_ed_placement, lep, frpl, race, 
             years_in_ma, plan_504, first_language, disability_type, 
             need_level, migrant, 
             e_test_stat, e_alt_assesment, e_raw_score, 
             m_test_stat, m_alt_assesment, m_raw_score,  
             s_test_stat, s_alt_assesment, s_raw_score, contains("wp"), contains("item")) 
    
  }else if(y %in% c("2009")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off,
        spec_ed_placement = spedplacement,
        lep = lep_off, 
        frpl = freelunch_off, 
        race = race_off, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        migrant = migrant_off,
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_raw_score = erawsc, 
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_raw_score = mrawsc,
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_raw_score = srawsc, 
        
        e_scale_score_2006 = escaleds_2006,
        e_scale_score_2007 = escaleds_2007,
        e_scale_score_2008 = escaleds_2008,
        m_scale_score_2006 = mscaleds_2006,
        m_scale_score_2007 = mscaleds_2007,
        m_scale_score_2008 = mscaleds_2008,
        
      ) %>%
      select(year, grade, gender,  gender, spec_ed, spec_ed_placement, lep, frpl, race, 
             years_in_ma, plan_504, first_language, disability_type, 
             need_level, migrant, 
             
             e_test_stat, e_alt_assesment, e_raw_score, 
             m_test_stat, m_alt_assesment, m_raw_score,  
             s_test_stat, s_alt_assesment, s_raw_score,
             
             e_scale_score_2006, e_scale_score_2007, e_scale_score_2008,
             
             m_scale_score_2006, m_scale_score_2007, m_scale_score_2008,
             
             contains("wp"), contains("item")) 
    
  }else if(y %in% c("2010", "2011")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off,
        spec_ed_placement = spedplacement,
        lep = lep_off, 
        frpl = freelunch_off, 
        race = race_off, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        migrant = migrant_off,
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_raw_score = erawsc, 
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_raw_score = mrawsc,
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_raw_score = srawsc, 
        
      ) %>%
      rename_with(.cols = contains("escaleds20"), 
                  .fn = ~str_replace(., "escaleds", "e_scale_score_")
      ) %>%
      rename_with(.cols = contains("mscaleds20"), 
                  .fn = ~str_replace(., "mscaleds", "m_scale_score_")
      ) %>%
      select(year, grade, gender,  gender, spec_ed, spec_ed_placement, lep, frpl, race, 
             years_in_ma, plan_504, first_language, disability_type, 
             need_level, migrant, 
             
             e_test_stat, e_alt_assesment, e_raw_score, 
             m_test_stat, m_alt_assesment, m_raw_score,  
             s_test_stat, s_alt_assesment, s_raw_score,
             
             contains("e_scale_score"),
             
             contains("m_scale_score"),
             
             contains("wp"),
             
             contains("item")) 
    
    
    
  }else if(y %in% c("2012")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off,
        spec_ed_placement = spedplacement,
        lep = lep_off, 
        frpl = freelunch_off, 
        race = race_off, 
        years_in_ma = yrsinmass,
        plan_504 = plan504_off,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_raw_score = erawsc, 
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_raw_score = mrawsc,
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_raw_score = srawsc, 
        
      ) %>%
      rename_with(.cols = contains("escaleds20"), 
                  .fn = ~str_replace(., "escaleds", "e_scale_score_")
      ) %>%
      rename_with(.cols = contains("mscaleds20"), 
                  .fn = ~str_replace(., "mscaleds", "m_scale_score_")
      ) %>%
      select(year, grade, gender,  gender, spec_ed, spec_ed_placement, lep, frpl, race, 
             years_in_ma, plan_504, first_language, disability_type, 
             need_level,  
             
             e_test_stat, e_alt_assesment, e_raw_score, 
             m_test_stat, m_alt_assesment, m_raw_score,  
             s_test_stat, s_alt_assesment, s_raw_score,
             
             contains("e_scale_score"),
             
             contains("m_scale_score"),
             
             contains("wp"),
             
             contains("item")) 
    
    
    
  }else if(y %in% c("2013", "2014")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade, 
        gender = gender, 
        spec_ed = sped_off,
        spec_ed_placement = spedplacement,
        lep = lep_off, 
        frpl = freelunch_off, 
        race = race_off, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_raw_score = erawsc, 
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_raw_score = mrawsc,
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_raw_score = srawsc, 
        
      ) %>%
      rename_with(.cols = contains("escaleds20"), 
                  .fn = ~str_replace(., "escaleds", "e_scale_score_")
      ) %>%
      rename_with(.cols = contains("mscaleds20"), 
                  .fn = ~str_replace(., "mscaleds", "m_scale_score_")
      ) %>%
      select(year, grade, gender,  gender, spec_ed, spec_ed_placement, lep, frpl, race, 
             years_in_ma, plan_504, first_language, disability_type, 
             need_level,  
             
             e_test_stat, e_alt_assesment, e_raw_score, 
             m_test_stat, m_alt_assesment, m_raw_score,  
             s_test_stat, s_alt_assesment, s_raw_score,
             
             contains("e_scale_score"),
             
             contains("m_scale_score"),
             
             contains("wp"),
             
             contains("item")) 
    
  }else if(y %in% c("2013", "2014")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade,
        gender = gender, 
        spec_ed_placement = spedplacement,
        economically_disadvantaged = ecodis, 
        race = race, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_mode = emode, 
        e_raw_score = erawsc, 
        
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_mode = mmode,
        m_raw_score = mrawsc,
        
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_mode = smode,
        s_subject = ssubject,
        s_raw_score = srawsc, 
        
      ) %>%
      rename_with(.cols = contains("escaleds20"), 
                  .fn = ~str_replace(., "escaleds", "e_scale_score_")
      ) %>%
      rename_with(.cols = contains("mscaleds20"), 
                  .fn = ~str_replace(., "mscaleds", "m_scale_score_")
      ) %>%
      select(year, grade, urban, town, county, gender, spec_ed_placement, ell, economically_disadvantaged, race, 
             years_in_ma, iep, plan_504, first_language, disability_type, 
             need_level,  
             
             accom_e, accom_m, accom_s, 
             calculator, humanreadaloud, texttospeech,
             
             e_test_stat, e_alt_assesment, e_mode, e_raw_score, 
             m_test_stat, m_alt_assesment, m_mode, m_raw_score,  
             s_test_stat, s_alt_assesment, s_mode, s_subject, s_raw_score,
             
             contains("e_scale_score"),
             
             contains("m_scale_score"),
             
             contains("wp"),
             
             contains("item")) 
    
  }else if(y %in% c("2017")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade,
        gender = gender, 
        spec_ed_placement = spedplacement,
        economically_disadvantaged = ecodis, 
        race = race, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_mode = emode, 
        e_attempt = eattempt,
        e_raw_score = erawsc, 
        
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_mode = mmode,
        m_attempt = mattempt,
        m_raw_score = mrawsc,
        
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_mode = smode,
        s_attempt = sattempt,
        s_subject = ssubject,
        s_raw_score = srawsc, 
        
        essay_1_idea_development = idea1, 
        essay_1_conventions_score = conv1,
        
        essay_2_idea_development = idea2, 
        essay_2_conventions_score = conv2,
        
        essay_3_idea_development = idea3, 
        essay_3_conventions_score = conv3,
        
        
        
        
      ) %>%
      rename_with(.cols = contains("escaleds20"), 
                  .fn = ~str_replace(., "escaleds", "e_scale_score_")
      ) %>%
      rename_with(.cols = contains("mscaleds20"), 
                  .fn = ~str_replace(., "mscaleds", "m_scale_score_")
      ) %>%
      select(year, grade, urban, town, county, gender, spec_ed_placement, ell, economically_disadvantaged, race, 
             years_in_ma, iep, plan_504, first_language, disability_type, 
             need_level,  
             
             accom_e, accom_m, accom_s, 
             calculator, humanreadaloud, texttospeech,
             
             e_test_stat, e_alt_assesment, e_mode, e_attempt, e_raw_score, 
             m_test_stat, m_alt_assesment, m_mode, m_attempt, m_raw_score,  
             s_test_stat, s_alt_assesment, s_mode, s_attempt, s_subject, s_raw_score,
             
             contains("e_scale_score"),
             
             contains("m_scale_score"),
             
             contains("essay"), 
             
             contains("wp"),
             
             contains("item")) 
    
  }else if(y %in% c("2018")){
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade,
        gender = gender, 
        spec_ed_placement = spedplacement,
        ell = el,
        economically_disadvantaged = ecodis, 
        race = race, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        
        calculator = accom_calculator, 
        readaloud = accom_readaloud, 
        scribe = accom_scribe,
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_mode = emode, 
        e_attempt = eattempt,
        e_raw_score = erawsc, 
        
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_mode = mmode,
        m_attempt = mattempt,
        m_raw_score = mrawsc,
        
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_mode = smode,
        s_attempt = sattempt,
        s_subject = ssubject,
        s_raw_score = srawsc, 
        
        essay_1_idea_development = idea1, 
        essay_1_conventions_score = conv1,
        
        essay_2_idea_development = idea2, 
        essay_2_conventions_score = conv2,
        
        essay_3_idea_development = idea3, 
        essay_3_conventions_score = conv3,
        
        
        
        
      ) %>%
      rename_with(.cols = contains("escaleds20"), 
                  .fn = ~str_replace(., "escaleds", "e_scale_score_")
      ) %>%
      rename_with(.cols = contains("mscaleds20"), 
                  .fn = ~str_replace(., "mscaleds", "m_scale_score_")
      ) %>%
      select(year, grade, urban, town, county, gender, spec_ed_placement, ell, economically_disadvantaged, race, 
             years_in_ma, iep, plan_504, first_language, disability_type, 
             need_level,  
             
             accom_e, accom_m, accom_s, 
             calculator, readaloud, scribe,
             
             e_test_stat, e_alt_assesment, e_mode, e_attempt, e_raw_score, 
             m_test_stat, m_alt_assesment, m_mode, m_attempt, m_raw_score,  
             s_test_stat, s_alt_assesment, s_mode, s_attempt, s_subject, s_raw_score,
             
             contains("e_scale_score"),
             
             contains("m_scale_score"),
             
             contains("essay"), 
             
             contains("wp"),
             
             contains("item"), 
             
             scl18ms:envdis7 
             
      ) 
    
  }else{
    
    df <- df %>%
      rename_with(str_to_lower) %>%
      mutate(year = y) %>%
      rename(
        grade = grade,
        gender = gender, 
        spec_ed_placement = spedplacement,
        ell = el,
        economically_disadvantaged = ecodis, 
        race = race, 
        years_in_ma = yrsinmass,
        plan_504 = plan504,
        first_language = firstlanguage,
        disability_type = natureofdis,
        need_level = levelofneed, 
        
        calculator = accom_calculator, 
        readaloud = accom_readaloud, 
        scribe = accom_scribe,
        
        e_test_stat = eteststat, 
        e_alt_assesment = ealt,
        e_mode = emode, 
        e_attempt = eattempt,
        e_raw_score = erawsc, 
        
        m_test_stat = mteststat,
        m_alt_assesment = malt,
        m_mode = mmode,
        m_attempt = mattempt,
        m_raw_score = mrawsc,
        
        s_test_stat = steststat,
        s_alt_assesment = salt,
        s_mode = smode,
        s_attempt = sattempt,
        s_subject = ssubject,
        s_raw_score = srawsc, 
        
        essay_1_idea_development = idea1, 
        essay_1_conventions_score = conv1,
        
        essay_2_idea_development = idea2, 
        essay_2_conventions_score = conv2,
        
        essay_3_idea_development = idea3, 
        essay_3_conventions_score = conv3,
        
        
        
        
      ) %>%
      rename_with(.cols = contains("escaleds20"), 
                  .fn = ~str_replace(., "escaleds", "e_scale_score_")
      ) %>%
      rename_with(.cols = contains("mscaleds20"), 
                  .fn = ~str_replace(., "mscaleds", "m_scale_score_")
      ) %>%
      select(year, grade, urban, town, county, gender, spec_ed_placement, ell, economically_disadvantaged,
             military, homeless, foster, migrant, race, 
             years_in_ma, iep, plan_504, first_language, disability_type, 
             need_level,  
             
             accom_e, accom_m, accom_s, 
             calculator, readaloud, scribe,
             
             e_test_stat, e_alt_assesment, e_mode, e_attempt, e_raw_score, 
             m_test_stat, m_alt_assesment, m_mode, m_attempt, m_raw_score,  
             s_test_stat, s_alt_assesment, s_mode, s_attempt, s_subject, s_raw_score,
             
             contains("e_scale_score"),
             
             contains("m_scale_score"),
             
             contains("essay"), 
             
             contains("wp"),
             
             contains("item") 
             
      ) 
    
  }
  
  
  
  return(df)
  
}

#2001 
item_responses_01_d <- cleanVarnames(inpath = "temp_files/d_01.rda", "2001")
save(item_responses_01_d, file = "clean_data/item_responses_2001_with_demographics.rda")
rm(item_responses_01_d)

#2002
item_responses_02_d <- cleanVarnames(inpath = "temp_files/d_02.rda", "2002")
save(item_responses_02_d, file = "clean_data/item_responses_2002_with_demographics.rda")
rm(item_responses_02_d)

#2003
item_responses_03_d <- cleanVarnames(inpath = "temp_files/d_03.rda", "2003")
save(item_responses_03_d, file = "clean_data/item_responses_2003_with_demographics.rda")
rm(item_responses_03_d)

#2004
item_responses_04_d <- cleanVarnames(inpath = "temp_files/d_04.rda", "2004")
save(item_responses_04_d, file = "clean_data/item_responses_2004_with_demographics.rda")
rm(item_responses_04_d)

#2005
item_responses_05_d <- cleanVarnames(inpath = "temp_files/d_05.rda", "2005")
save(item_responses_05_d, file = "clean_data/item_responses_2005_with_demographics.rda")
rm(item_responses_05_d)

#2006
item_responses_06_d <- cleanVarnames(inpath = "temp_files/d_06.rda", "2006")
save(item_responses_06_d, file = "clean_data/item_responses_2006_with_demographics.rda")
rm(item_responses_06_d)

#2007
item_responses_07_d <- cleanVarnames(inpath = "temp_files/d_07.rda", "2007")
save(item_responses_07_d, file = "clean_data/item_responses_2007_with_demographics.rda")
rm(item_responses_07_d)

#2008
item_responses_08_d <- cleanVarnames(inpath = "temp_files/d_08.rda", "2008")
save(item_responses_08_d, file = "clean_data/item_responses_2008_with_demographics.rda")
rm(item_responses_08_d)

#2009
item_responses_09_d <- cleanVarnames(inpath = "temp_files/d_09.rda", "2009")
save(item_responses_09_d, file = "clean_data/item_responses_2009_with_demographics.rda")
rm(item_responses_09_d)

#2010
item_responses_10_d <- cleanVarnames(inpath = "temp_files/d_10.rda", "2010")
save(item_responses_10_d, file = "clean_data/item_responses_2010_with_demographics.rda")
rm(item_responses_10_d)

#2011
item_responses_11_d <- cleanVarnames(inpath = "temp_files/d_11.rda", "2011")
save(item_responses_11_d, file = "clean_data/item_responses_2011_with_demographics.rda")
rm(item_responses_11_d)

#2012
item_responses_12_d <- cleanVarnames(inpath = "temp_files/d_12.rda", "2012")
save(item_responses_12_d, file = "clean_data/item_responses_2012_with_demographics.rda")
rm(item_responses_12_d)

#2013
item_responses_13_d <- cleanVarnames(inpath = "temp_files/d_13.rda", "2013")
save(item_responses_13_d, file = "clean_data/item_responses_2013_with_demographics.rda")
rm(item_responses_13_d)

#2014
item_responses_14_d <- cleanVarnames(inpath = "temp_files/d_14.rda", "2014")
save(item_responses_14_d, file = "clean_data/item_responses_2014_with_demographics.rda")
rm(item_responses_14_d)

#2015
#NOTE: MCAS data this year does not include test items

#2016
#NOTE: MCAS data this year does not include test items

#2017
item_responses_17_d <- cleanVarnames(inpath = "temp_files/d_17.rda", "2017")
save(item_responses_17_d, file = "clean_data/item_responses_2017_with_demographics.rda")
rm(item_responses_17_d)

#2018
item_responses_18_d <- cleanVarnames(inpath = "temp_files/d_18.rda", "2018")
save(item_responses_18_d, file = "clean_data/item_responses_2018_with_demographics.rda")
rm(item_responses_18_d)

#2019
item_responses_19_d <- cleanVarnames(inpath = "temp_files/d_19.rda", "2019")
save(item_responses_19_d, file = "clean_data/item_responses_2019_with_demographics.rda")
rm(item_responses_19_d)


