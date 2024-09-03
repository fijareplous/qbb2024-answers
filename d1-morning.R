##Day 1 Exercises

library("tidyverse")

df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1 )

#Opened and viewed df, confirmed SUBJECT is first column 
#view(df)

df %>%
  group_by(SUBJECT) %>%
  summarize(SUBJECT_n = n()) %>%
  arrange(desc(SUBJECT_n))
view(df)
#Subjects K-562 and GTEX-NPJ8 have the most samples (217 and 72, respectively)

df %>%
  group_by(SUBJECT) %>%
  summarize(SUBJECT_n = n()) %>%
  arrange(SUBJECT_n)
view(df)
#Subjects GTEX-1JMI6 and GTEX-1PAR6 have the least samples (1 each)

df %>%
  group_by(SMTSD) %>%
  summarize(SUBJECT_n = n()) %>%
  arrange(desc(SUBJECT_n))
view(df)
#Whole Blood and Muscle - Skeletal are the most common tissue types

df %>%
  group_by(SMTSD) %>%
  summarize(SUBJECT_n = n()) %>%
  arrange(SUBJECT_n)
view(df)
#Kidney - Medulla and Cervix - Ectocervix are the least common tissue types 

df_npj8 <- subset(df, SUBJECT == "GTEX-NPJ8")
view(df_npj8)

df_npj8 %>%
  group_by(SMTSD) %>%
  summarize(SMTSD_n = n()) %>%
  arrange(desc(SMTSD_n))
#Whole Blood samples are the most common in individual GTEX-NPJ8

#Opened df_npj8, saw that the SMTSISCH and SMNABTCH numbers are different

nrow(df %>%
  filter( !is.na(SMATSSCR) ) %>%
  group_by(SUBJECT) %>%
  summarize(meanSMATSSCR = mean(SMATSSCR)) %>%
  filter(meanSMATSSCR==0))
#15 subjects have a mean SMATSSCR of 0

df %>%
  filter( !is.na(SMATSSCR) ) %>%
  group_by(SUBJECT) %>%
  summarize(meanSMATSSCR = mean(SMATSSCR)) %>%
  arrange(meanSMATSSCR)

hist(df$SMATSSCR)
#Most individuals have a SMATSSCR close to 1 (roughly 10,000 individuals). 0 is the lowest value, and 3 is the highest SMATSSCR seen in the data. 

`#I would probably make a histogram showing the number of subjects with various mean scores








