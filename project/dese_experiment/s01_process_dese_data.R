library(pacman)
p_load("tidyverse", "argparser", "feather",  "magic" , "pander", "dplyr",  "plotmo", "caret",  "cvms", "magrittr", "arrow", "here")

# parse arguments
arguments <- arg_parser("Process and Sample DESE data in 2022") %>% 
  add_argument(
    "--data_input_dir", 
    help= "data_input_dir", 
    default= here("DESE_data", "raw")
  ) %>%
  add_argument(
    "--data_output_dir", 
    help= "data_output_dir", 
    default= here("DESE_data", "sample_processed", "2022")
  ) %>%
  add_argument(
    "--grades", 
    help= "grades_considered", 
    default= c("03", "04", "05", "06", "07", "08", "10")
  ) %>%
  add_argument(
    "--year", 
    help= "exam_year", 
    default= c("22")
  ) %>%
  add_argument(
    "--sample_size", 
    help= "number of students to sample", 
    default= 2000
  ) %>%
  parse_args()


# filter raw data such that
## 1: focus on students participating in both math and english exams
## 2: filter out nonbinary items 
## 3: replace NA response with 0
# (added student id in case we lose track)
filter_item_response <- function(df){
  df <- df %>%
    mutate(student_id = paste0("s_id", 1:nrow(df))) %>%
    rename_with(str_to_lower) %>%
    filter(eteststat=="T" | eteststat=="RT") %>%
    filter(mteststat=="T" | mteststat=="RT") %>%
    select_if(~!all(is.na(.)))  %>%
    select(-(contains("item") & where(~ any(. > 1, na.rm=TRUE)))) %>%
    mutate(across(contains("item"), ~ ifelse(is.na(.), 0, .))) %>% 
    mutate(across(contains("mitem"), as.numeric)) %>%
    mutate(across(contains("eitem"), as.numeric))
    return(df)
}


set.seed(1)

# process and sample raw DESE data
data <-  arguments$grades %>%
  set_names(., .) %>%
  map(function(.x) read_tsv(file.path(arguments$data_input_dir, paste0("20", arguments$year), paste0("g",.x, "22d.txt")),
                            col_types = cols(.default = col_character()))) %>%
  map(function(.x) filter_item_response(.x)) %>%
  map(function(.x) .x %>% sample_n(arguments$sample_size, replace = FALSE)) 


# save sampled data

paste0(paste0("grade_", arguments[["grades"]], "_sampled"), ".feather") %>%
     map2(data, function(.x, .y) write_feather(.y, file.path(arguments[["data_output_dir"]], .x)))


        




