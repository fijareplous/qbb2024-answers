library(tidyverse)
library(palmerpenguins)
library(ggthemes)
?penguins

#Selects just the island column
penguins[ , "island"]

#Selects multiple columns
penguins[ , c("species", "island")]

#Selects row two, specific columns
penguins[2, c("species", "island")]

#Selects 2nd row, 1st column
penguins[2,1]

#ggplot, could use scale_color_manual(values = c("red", "blue", "green")) 
#instead of scale_color_colorblind
ggplot(data = penguins) + 
  geom_point(mapping = aes(x = flipper_length_mm, 
                           y = body_mass_g, 
                           color = species,
                           shape = species)) + 
  scale_color_colorblind() +
  geom_smooth(mapping = aes(x = flipper_length_mm, 
                            y = body_mass_g), method = "lm") +
  xlab("Flipper Length (mm)") +
  ylab("Body Mass (g)") +
  ggtitle("Relationship Between Flipper Length and Body Mass")

#Saves ggplot as file (knows file extensions)
ggsave(filename = "~/qbb2024-answers/d1-afternoon/penguin-plot.pdf")  

#Does bill length/depth correspond to sex of penguin?
#male_penguins = penguins %>% filter(sex == "male")
#Alpha is the opacity, binwidth is the width of the bars, fill is the way the bars
#are filled in, and facet_grid separates the sexes into different graphs
#Adding %>% filter*!is.na(sex) to the end of penguins in the beginning will get rid 
#of NA samples
ggplot(data = penguins, mapping = aes(x = bill_length_mm, fill = sex)) + 
  geom_histogram(binwidth = 2, position = "identity", alpha = .7) + facet_grid(. ~ sex)

#Has flipper length gone up or down over the years?
#Need categories with boxplot which is why x = factor(year)
ggplot(data = penguins, mapping = aes(x = factor(year), y = body_mass_g, fill = sex)) + 
  geom_boxplot() + facet_grid(island ~ species)
