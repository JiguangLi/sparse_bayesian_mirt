library(pacman)
p_load("tidyverse", "argparser", "feather",  "yaml")


# filter data by grade and subject, filter out students who haven't taken exam
get_filtered_df <- function(data_dir, filter_grade, subject){
  year_df <- read_feather(data_dir)
  year_df$student_id <- seq(1, dim(year_df)[1])
  test_stat <- paste(subject, "_test_stat", sep="")
  df <- filter(year_df, grade==filter_grade, .data[[test_stat]] =="T")
  res_cols<- colnames(df)[grepl(paste0("^", subject, "item*"), colnames(df))]
  df<- df[rowSums(is.na(df[,res_cols]))!=length(res_cols),]
  return(df)
}

# filter out item responses only
get_response_df <- function(df, subject){
  response_df <- df %>% select(starts_with(paste0(subject,"item")))
  for(i in 1: length(colnames(response_df))){
    if(length(intersect(names(table(response_df[, i])), c("A", "B", "C", "D")))==3){
      response_df[, i] <- as.numeric(((response_df[,i]=="+")& (!is.na(response_df[,i]))))
    }else{
      response_df[,i][is.na(response_df[,i])] <- "0"
      
      response_df[,i]<- as.numeric(as.character(response_df[[colnames(response_df[, i])]]))
      
    }
  }
  return(response_df)
}


config_filepath <- "/Users/jiguangli/IBP_IRT/config/data_config.yaml"
config <- read_yaml(config_filepath)
filename <- paste("item", "responses", config[["year"]], "d.feather", sep= "_")
data_dir <- file.path(config[["feather_data_dir"]], filename)
df <- get_filtered_df(data_dir, config[["grade"]], config[["subject"]])
response_df<- get_response_df(df, config[["subject"]])

# filter out non-binary responses for now
binary_response <-response_df %>% 
  mutate_all(as.numeric) %>% 
  select_if(function(col) max(col) == 1) %>%
  mutate(student_id = 1:dim(.)[1]) %>%
  write_feather(.,"file.feather")

# sample
set.seed(42)
data_filename_format <- paste(config[["year"]], config[["grade"]], config[["subject"]], "%s.feather", sep= "_")
data_filename_format <- paste(config[["output_data_dir"]], data_filename_format, sep= "/")

train_df <- binary_response %>% 
  sample_n(config[["train_sample"]]) %>%
  write_feather(.,sprintf(data_filename_format, "train")) 

test_df <- binary_response %>% 
  anti_join(train_df, by = 'student_id') %>%
  sample_n(config[["test_sample"]]) %>%
  write_feather(.,sprintf(data_filename_format, "test"))















