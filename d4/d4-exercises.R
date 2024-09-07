library(tidyverse)
library(broom)

'
wget https://raw.githubusercontent.com/bxlab/cmdb-quantbio/main/
assignments/bootcamp/linear_regression/assignment/aau1043_dnm.csv
wget https://raw.githubusercontent.com/bxlab/cmdb-quantbio/main/
assignments/bootcamp/linear_regression/assignment/aau1043_parental_age.csv'

#Part 1.1
dnm <- read_csv(file = "~/qbb2024-answers/d4/aau1043_dnm.csv")


#Part 1.2
dnm_summary <- dnm %>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm = sum(Phase_combined == "father", na.rm = TRUE),
            n_maternal_dnm = sum(Phase_combined == "mother", na.rm = TRUE))
#phase combined evaluates to true or false depending on whether or not it 
#comes from father/mother


#Part 1.3 
ages <- read_csv(file = "~/qbb2024-answers/d4/aau1043_parental_age.csv")

#Part 1.4
joint <- left_join(dnm_summary, ages, by = "Proband_id")
joint

#Part 2.1
ggplot(data = joint, mapping = aes(x= Mother_age, y= n_maternal_dnm)) +
         geom_point() + ggtitle("Count of de Novo Mutations Across Maternal
         Ages") + xlab("Maternal Age") + ylab("Number of Maternal de Novo
          Mutations") + geom_smooth(method = "lm")

ggplot(data = joint, mapping = aes(x= Father_age, y= n_paternal_dnm)) +
  geom_point() + ggtitle("Count of de Novo Mutations Across Paternal Ages") +
    xlab("Paternal Age") + ylab("Number of Paternal de Novo Mutations") +
  geom_smooth(method = "lm")


#Part 2.2
mother_model <- lm(data = joint, formula = n_maternal_dnm ~1 + Mother_age)
mother_model
'
The p value for the relationship between maternal age and the number of maternal 
de novo mutations is <2.2x10^-16, which is extremely small, demonstrating that 
the chances of the two variables being independent and producing similar results
is extremely small (thus being statistically significant). The chart of the data
shows a correlation between the two variables, although there is less than I 
would expect for such a low p value. 
'

#Part 2.3
father_model <- lm(data = joint, formula = n_paternal_dnm ~1 + Father_age)
'
As with the previous question, The p value for the relationship between maternal 
age and the number of maternal de novo mutations is <2.2x10^-16, which is extremely
small, demonstrating that the results are statistically significant. The trends 
seen in the chart match these results, as paternal age and number of DNM appear
highly correlated. 
'

#Part 2.4
#In 2.3, the intercept was reported to be 7.018, and the slope is 0.457
new_proband <- tibble(Father_age = 50.5)
predict(father_model, new_proband)

#The number of de novo mutations should be close to 78.70 (or 79, really).

#Part 2.5
ggplot(data = joint) + geom_histogram(mapping = aes(x= n_paternal_dnm, fill =
    "Paternal"), alpha = .7,) + geom_histogram(mapping = aes(x= n_maternal_dnm, 
    fill = "Maternal"), alpha = .7) + xlab("De Novo Mutations") + ylab("Number of 
    Probands")

' This was as far as I got on Friday; I tried working on it more after bootcamp 
ended, but I am pretty lost. As such, I am planning to attend office hours on 
Monday to ask questions about the setup of part 2.5 and how to approach 2.6. My
apologies that the last two exercises are incomplete; I will submit them as soon
as I can.'
