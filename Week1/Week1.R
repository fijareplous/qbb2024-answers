library(tidyverse)
library(ggplot2)

#Question 1.3

#Import coverage file:
genome_coverage <- read.delim("~/qbb2024-answers/Week1/genome_coverage_3x.txt")

#Need a column name
colnames(genome_coverage) <- c("coverage")

#Calculate coverage frequencies
cov_freq = genome_coverage %>% group_by(coverage) %>% summarize(n = n())

#Poisson probability mass function (PMF); probability of getting certain coverage
poisson_pmf <- dpois(cov_freq$coverage, lambda =3)

#Normal probability density function (PDF); prob. of getting certain coverage
norm_pdf <- dnorm(cov_freq$coverage, mean = 3, sd = sqrt(3))

'
ggplot time! geom_histogram plots the distribution of coverage across the genome;
the first geom_line displays the Poisson distribution, and the second geom_line
function displays the normal distribution. 
'

cov_plot <- ggplot() + 
  geom_histogram(data = genome_coverage, mapping = aes(x = 
     coverage, fill = "Genome Coverage (3x)"), color = "black", binwidth = 1) + 
  geom_line(data = cov_freq, mapping = aes(x = coverage, y = 
     poisson_pmf*nrow(genome_coverage), color = "Poisson Distribution")) + 
  geom_line(data = cov_freq, mapping = aes(x = coverage, y = norm_pdf*nrow(genome_coverage), 
     color = "Normal Distribution")) + 
  ylab("Frequency") + xlab("Coverage") + 
  scale_color_manual(name = "Legend", values = c("#9d0fff", "#e90268")) + 
  scale_fill_manual(name = " ", values = c("#b6ebe4")) + 
  ggtitle("Poisson and Normal Distributions Against 3x Genome Coverage")

ggsave(filename = "~/qbb2024-answers/Week1/ex1_3x_cov.png", 
       plot = cov_plot, width = 10, height = 6, dpi = 300)


#Question 1.5

#Import coverage file:
genome_coverage <- read.delim("~/qbb2024-answers/Week1/genome_coverage_10x.txt")

#Need a column name
colnames(genome_coverage) <- c("coverage")

#Calculate coverage frequencies
cov_freq = genome_coverage %>% group_by(coverage) %>% summarize(n = n())

#Poisson probability mass function (PMF); probability of getting certain coverage
poisson_pmf <- dpois(cov_freq$coverage, lambda =10)

#Normal probability density function (PDF); prob. of getting certain coverage
norm_pdf <- dnorm(cov_freq$coverage, mean = 10, sd = 3.16)

'
ggplot time! geom_histogram plots the distribution of coverage across the genome;
the first geom_line displays the Poisson distribution, and the second geom_line
function displays the normal distribution. 
'

cov_plot <- ggplot() + 
  geom_histogram(data = genome_coverage, mapping = aes(x = 
     coverage, fill = "Genome Coverage (10x)"), color = "black", binwidth = 1) + 
  geom_line(data = cov_freq, mapping = aes(x = coverage, y = 
     poisson_pmf*nrow(genome_coverage), color = "Poisson Distribution")) + 
  geom_line(data = cov_freq, mapping = aes(x = coverage, y = norm_pdf*nrow(genome_coverage), 
      color = "Normal Distribution")) + 
  ylab("Frequency") + xlab("Coverage") + 
  scale_color_manual(name = "Legend", values = c("#9d0fff", "#e90268")) + 
  scale_fill_manual(name = " ", values = c("#b6ebe4")) + 
  ggtitle("Poisson and Normal Distributions Against 10x Genome Coverage")

ggsave(filename = "~/qbb2024-answers/Week1/ex1_10x_cov.png", 
       plot = cov_plot, width = 10, height = 6, dpi = 300)


#Question 1.6

#Import coverage file:
genome_coverage <- read.delim("~/qbb2024-answers/Week1/genome_coverage_30x.txt")

#Need a column name
colnames(genome_coverage) <- c("coverage")

#Calculate coverage frequencies
cov_freq = genome_coverage %>% group_by(coverage) %>% summarize(n = n())

#Poisson probability mass function (PMF); probability of getting certain coverage
poisson_pmf <- dpois(cov_freq$coverage, lambda =30)

#Normal probability density function (PDF); prob. of getting certain coverage
norm_pdf <- dnorm(cov_freq$coverage, mean = 30, sd = 5.47)

'
ggplot time! geom_histogram plots the distribution of coverage across the genome;
the first geom_line displays the Poisson distribution, and the second geom_line
function displays the normal distribution. 
'

cov_plot <- ggplot() + 
  geom_histogram(data = genome_coverage, mapping = aes(x = 
     coverage, fill = "Genome Coverage (30x)"), color = "black", binwidth = 1) + 
  geom_line(data = cov_freq, mapping = aes(x = coverage, y = 
     poisson_pmf*nrow(genome_coverage), color = "Poisson Distribution")) + 
  geom_line(data = cov_freq, mapping = aes(x = coverage, y = norm_pdf*nrow(genome_coverage), 
      color = "Normal Distribution")) + 
  ylab("Frequency") + xlab("Coverage") + 
  scale_color_manual(name = "Legend", values = c("#9d0fff", "#e90268")) + 
  scale_fill_manual(name = " ", values = c("#b6ebe4")) + 
  ggtitle("Poisson and Normal Distributions Against 30x Genome Coverage")

ggsave(filename = "~/qbb2024-answers/Week1/ex1_30x_cov.png", 
       plot = cov_plot, width = 10, height = 6, dpi = 300)

