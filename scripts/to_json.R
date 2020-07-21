#This file is made to convert some of our data for creating some .json files needed for the setup of the Pheweb

install.packages("jsonify")
install.packages("icd.data")

library(tidyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(reshape)
library(jsonlite)
library(icd.data)
library(stringr)

chapters <- fread("I:/RDecid/Ana_Sanchez/GWAS/icd9_ch_subch.csv", colClasses = "character" ) 

data_table <- data.frame(icd9cm_billable[[1]]) %>% mutate(TIPO = substr(code,0,3)) %>% left_join(chapters, by = "TIPO")

disease_code <- fread("I:/RDecid/Ana_Sanchez/categories_grupsmalal.csv") %>% mutate(code = gsub("TIPO|GROUP|Tipo", "", ICD9))  %>% 
mutate(TIPO = substr(code,0,3)) %>% left_join(chapters, by ="TIPO") %>% 
  select(Chapter, ICD9, Descripci?) %>% `colnames<-`(c("category", "phenocode", "phenostring")) %>% data.frame() 
#to json
json_table <- (jsonlite::toJSON(disease_code %>% select(phenocode, category), pretty = TRUE,  simplifyVector = TRUE))
write_json(json_table, "I:/RDecid/Ana_Sanchez/GWAS/categories_provados.json", pretty = TRUE,  simplifyVector = TRUE)
#312 category. phenocode, phenostring


write.table(disease_code %>%select(phenocode, category),"I:/RDecid/Ana_Sanchez/categories.csv", row.names = FALSE, quote= FALSE, sep = ";")
read_json("I:/RDecid/Ana_Sanchez/GWAS/categories_provados.json")

#categories csv is a mandatory file
dataset <- read_delim("C:/Users/studentdnabank/Desktop/categories.csv", ";", escape_double = FALSE, trim_ws = TRUE)
json_table <- (jsonlite::toJSON(dataset %>% select(phenocode, category), pretty = TRUE,  simplifyVector = TRUE))
