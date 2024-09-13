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

Additionally, the slope of the trendline is 0.378, which suggests that the number 
of maternal dnms increases slightly with maternal age (very roughly by 1 dnm every
3 years, or .378 per year). This increase with age, while slight, is still 
statistically significant given the p value. 
'

#Part 2.3
father_model <- lm(data = joint, formula = n_paternal_dnm ~1 + Father_age)
father_model
'
As with the previous question, The p value for the relationship between maternal 
age and the number of maternal de novo mutations is <2.2x10^-16, which is extremely
small, demonstrating that the results are statistically significant. The trends 
seen in the chart match these results, as paternal age and number of DNM appear
highly correlated. 

Furthermore, the slope of the trendline 1.35, which suggests that the number of 
paternal dnms increases by about 1.35 per year. This increase is more substantial
than seen in maternal dnms, although both are statistically significant.
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

#Part 2.6
t.test(joint$n_paternal_dnm, joint$n_maternal_dnm, 
  alternative = "two.sided")
'
A two sided t-test would be appropriate given that we want to examine whether 
there is a statistical difference between the means of the maternally vs paternally
dnms per proband. Because we are comparing two groups, (not 3+, in the case of 
ANOVA tests), a t test is appropriate. 

The test produces a p value of >2.2x10^-16, which is extremely small, and sample
means of 52.02 and 12.7 for paternal and maternal dnms, respectively. This suggests
that the two groups are statistically different from one another and that those
differences are unlikely due to chance.
'
