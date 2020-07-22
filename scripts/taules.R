#Making tables for showing the diseases selected
#Taula de les 54 malalties
#We first load the description files
data_table <- data.frame(icd9cm_billable[[1]]) %>%  `colnames<-`(c("code", "desc", "long"))
icd_desc <- fread("I:/RDecid/Ana_Sanchez/GWAS/icd9_ch_subch.csv", sep = ";", h=T, na.strings = c("NA", "")) %>%  `colnames<-`(c("TIPO_3", "DESC", "Chapter", "Subchapter", "Subchapter2"))
icd_desc2 <- fread("I:/RDecid/Ana_Sanchez/GWAS/icd9_ch_subch.csv", sep = ";", h=T, na.strings = c("NA", "")) %>%  `colnames<-`(c("code_group", "DESC", "Chapter", "Subchapter", "Subchapter2"))

#We create the table for 54 diseases
tipos_54<- read_excel("I:/RDecid/Ana_Sanchez/tipo_freq.xlsx") %>% mutate(code = gsub("TIPO|Tipo", "", TIPO)) %>% 
  left_join(data_table, by="code")%>% mutate(TIPO_3 = substr(code,0,3)) %>% left_join(icd_desc, by ="TIPO_3")%>%
  mutate(code_group= gsub("GROUP", "" , TIPO))%>% mutate(code_group = substr(code_group,0,3))%>% left_join(icd_desc2, by = "code_group") %>%
  mutate(chapter = coalesce(Chapter.x, Chapter.y))  %>% select(TIPO,CASE, CONTROL, FREQ, long, chapter) %>% 
  `colnames<-`(c("TIPO", "Cases", "Controls", "Freq","Description","Chapter")) %>% filter(Cases >=15)

#Save the table
write.table(tipos_54, "I:/RDecid/Ana_Sanchez/TIPO_52.csv", row.names = FALSE, quote= FALSE, sep = ";")

#Tacle of the 10 selected:
selected_tipos <- c("TIPO2749","TIPO36504","TIPO3659","TIPO4139","TIPO43491","TIPO5569","TIPO6961","TIPOC44","TIPOC50","GROUP272")
tipos_10 <- tipos_52 %>% filter(TIPO %in% selected_tipos)

#Safe the table

write.table(tipos_10, "I:/RDecid/Ana_Sanchez/TIPO_10.csv", row.names = FALSE, quote= FALSE, sep = ";")

#Taula dels 350
data_table2 <- data.frame(icd9cm_billable[[1]])
icd_desc3 = fread("I:/RDecid/Ana_Sanchez/GWAS/icd9_ch_subch.csv", sep = ";", h=T, na.strings = c("NA", ""))
taula_350 <- fread("I:/RDecid/Ana_Sanchez/TIPO_350.csv") 

