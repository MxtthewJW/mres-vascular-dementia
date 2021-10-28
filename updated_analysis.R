# Vascular dementia diagnostic comparison analysis

# Packages ----

library(tidyverse)
library(MASS)
library(msm)
library(DescTools)

# Power analysis for time periods

power.chisq.test(n = 32, w = 0.3, df = 1, sig.level = 0.05) # 40% power

# Time period results tables and chisquare/fishers ----
## Fishers test used for expected counts < 5 

icd_pand <- matrix(c(13, 3, 13, 3),
                   nrow = 2,
                   dimnames = list(c('VaD','NoVaD'),
                                   c('Pre-covid','Mid-covid')))
chisq.test(icd_pand)
fisher.test(icd_pand, alternative = "two.sided")

dsm4_pand <- matrix(c(12, 4, 8, 8),
                    nrow = 2,
                    dimnames = list(c('VaD','NoVaD'),
                                    c('Pre-covid','Mid-covid')))
chisq.test(dsm4_pand)

dsm5_pand <- matrix(c(14, 2, 11, 5),
                    nrow = 2,
                    dimnames = list(c('VaD','NoVaD'),
                                    c('Pre-covid','Mid-covid')))
chisq.test(dsm5_pand)
fisher.test(dsm5_pand, alternative = "two.sided")

ninds_pand <- matrix(c(10, 6, 8, 8),
                     nrow = 2,
                     dimnames = list(c('VaD','NoVaD'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(ninds_pand)

# power analysis for criteria comparison over separate time periods 

power.chisq.test(n = 16, w = 0.3, df = 1, sig.level = 0.05) # 22% power

# Criteria comparison tables for separate time periods with McNemar test ----
## McNemar test used given the correlation between IVs

#pre
icd_dsm4_pre <- matrix(c(12, 0, 1, 3),
                   nrow = 2,
                   dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                   c('DSM-IVVaD','DSM-IVNoVaD')))
mcnemar.test(icd_dsm4_pre)

icd_dsm5_pre <- matrix(c(12,2,1,1),
                       nrow = 2,
                       dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                       c('DSM-VVaD','DSM-VNoVaD')))
mcnemar.test(icd_dsm5_pre)

icd_ninds_pre <- matrix(c(10,0,3,3),
                    nrow = 2,
                    dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                    c('NINDSVaD','NINDSNoVaD')))
mcnemar.test(icd_ninds_pre)
#post
icd_dsm4_post <- matrix(c(8,0,5,3),
                       nrow = 2,
                       dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                       c('DSM-IVVaD','DSM-IVNoVaD')))
mcnemar.test(icd_dsm4_post)

icd_dsm5_post <- matrix(c(10,1,3,2),
                       nrow = 2,
                       dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                       c('DSM-VVaD','DSM-VNoVaD')))
mcnemar.test(icd_dsm5_post)

icd_ninds_post <- matrix(c(8,0,5,3),
                        nrow = 2,
                        dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                        c('NINDSVaD','NINDSNoVaD')))
mcnemar.test(icd_ninds_post)

# post-hoc power analysis for criteria comparison

power.chisq.test(n = 32, w = 0.3, df = 1, sig.level = 0.05) # 40% power

# Criteria comparison tables for aggregated time periods with McNemar test ----

icd_dsm4 <- matrix(c(20, 0, 6, 6),
                    nrow = 2,
                    dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                    c('DSM-IVVaD','DSM-IVNoVaD')))
mcnemar.test(icd_dsm4)

icd_dsm5 <- matrix(c(22, 3, 4, 3),
                   nrow = 2,
                   dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                   c('DSM-VVaD','DSM-VNoVaD')))
mcnemar.test(icd_dsm5)

icd_ninds <- matrix(c(18, 0, 8, 6),
                   nrow = 2,
                   dimnames = list(c('ICDVaD', 'ICDNoVaD'),
                                   c('NINDSVaD','NINDSNoVaD')))
mcnemar.test(icd_ninds)

dsm4_dsm5 <- matrix(c(19, 1, 6, 6),
                   nrow = 2,
                   dimnames = list(c('DSM-IVVaD','DSM-IVNoVaD'),
                                   c('DSM-VVaD','DSM-VNoVaD')))
mcnemar.test(dsm4_dsm5)

dsm4_ninds <- matrix(c(17, 1, 3, 11),
                   nrow = 2,
                   dimnames = list(c('DSM-IVVaD','DSM-IVNoVaD'),
                                   c('NINDSVaD','NINDSNoVaD')))
mcnemar.test(dsm4_ninds)

dsm5_ninds <- matrix(c(16, 2, 9, 5),
                   nrow = 2,
                   dimnames = list(c('DSM-VVaD','DSM-VNoVaD'),
                                   c('NINDSVaD','NINDSNoVaD')))
mcnemar.test(dsm5_ninds)

# power analysis for brain imaging element due to smaller sample 

power.chisq.test(n = 25, w = 0.3, df = 1, sig.level = 0.05) # 32% power

# Element comparison criteria tables and chisq test ----

icd_cogdec <- matrix(c(14, 2, 14, 2),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(icd_cogdec)
fisher.test(icd_cogdec, alternative = "two.sided")

icd_imaging <- matrix(c(11, 4, 9, 1),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(icd_imaging)
fisher.test(icd_imaging, alternative = "two.sided")

icd_temp <- matrix(c(8, 8, 3, 13),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(icd_temp)

icd_park <- matrix(c(5, 11, 2, 14),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(icd_park)
fisher.test(icd_park, alternative = "two.sided")

dsm4_cogdec <- matrix(c(13,3,12,4),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm4_cogdec)
fisher.test(dsm4_cogdec, alternative = "two.sided")

dsm4_imaging <- matrix(c(11,4,9,1),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm4_imaging)
fisher.test(dsm4_imaging, alternative = "two.sided")

dsm4_temp <- matrix(c(8,8,3,13),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm4_temp)

dsm4_park <- matrix(c(5,11,2,14),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm4_park)

dsm5_cogdec <- matrix(c(15,1,13,3),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm5_cogdec)

dsm5_imaging <- matrix(c(10,5,9,1),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm5_imaging)

dsm5_temp <- matrix(c(8,8,3,13),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm5_temp)

dsm5_park <- matrix(c(5,11,2,14),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm5_park)

ninds_cogdec <- matrix(c(12,4,9,7),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(ninds_cogdec)

ninds_imaging <- matrix(c(5,10,4,6),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(ninds_imaging)

ninds_temp <- matrix(c(0,16,2,14),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(ninds_temp)

ninds_park <- matrix(c(0,16,2,14),
                     nrow = 2,
                     dimnames = list(c('Met','NoMet'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(ninds_park)

# power for brain imaging

power.chisq.test(n = 32, w = 0.5, df = 1, sig.level = 0.05)

# Brain imaging availability

imaging_avail <- matrix(c(15, 10, 1, 6),
                   nrow = 2,
                   dimnames = list(c('Available', 'Unavailable'),
                                   c('Pre-covid','Mid-covid')))
fisher.test(imaging_avail, alternative = 'greater')

# power for severity (for both as samples different)

power.chisq.test(n = 25, w = 0.3, df = 1, sig.level = 0.05) # 32%
power.chisq.test(n = 18, w = 0.3, df = 1, sig.level = 0.05) # 25%

# Severity results tables and chi square/fishers ----

dsm5_sev <- matrix(c(8, 6, 4, 7),
                     nrow = 2,
                     dimnames = list(c('Major', 'Mild'),
                                     c('Pre-covid','Mid-covid')))
chisq.test(dsm5_sev)

ninds_sev <- matrix(c(2, 8, 3, 5),
                   nrow = 2,
                   dimnames = list(c('Major', 'Mild'),
                                   c('Pre-covid','Mid-covid')))
chisq.test(ninds_sev)
fisher.test(ninds_sev, alternative = "two.sided")

# participant information ----

# simulate 100 samples with same group statistical properties

total_samples <- 100
sample_size <- 32
participant <- rep(1:sample_size)
condition <- c(rep("Pre-covid", times = sample_size/2),
               rep("Mid-covid", times = sample_size/2))
all_data <- NULL

for (i in 1:total_samples) {
  sample <- i
  set.seed(1233 + i)
  dv <- c(rnorm(sample_size/2, 82.8125, 7.833422), rnorm(sample_size/2, 78.5625, 6.470124))
  my_data <- as_tibble(cbind(participant, condition, dv, sample))
  all_data <- rbind(my_data, all_data)
}

all_tidied_data <- all_data %>%
  mutate(condition = factor(condition), dv = as.integer(dv))

# visualising the age group mean samples

png("simulated_age_groups.png", units="in", width=16, height=10, res=500)
all_tidied_data %>%
  group_by(condition, sample) %>%
  summarise(average = mean(dv)) %>%
  ggplot(aes(x = condition, y = average, group = condition,
             label = sample, fill = condition)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .5, colour = 'black') +
  stat_summary(fun.data = "mean_cl_boot", colour = "black", shape = "diamond",
               size = 1.5, fill = 'white') +
  theme(text=element_text(family="TT Arial")) +
  labs(x = "Pandemic group", y = "Age (sample mean)") +
  theme_minimal(base_size = 28) +
  coord_flip()
dev.off()

# t-tests with results

result <- NULL
for (i in 1:total_samples) {
  result <- rbind(tidy(t.test(filter(all_tidied_data, 
                                     condition == "Pre-covid" & sample == i)$dv,
                              filter(all_tidied_data, 
                                     condition == "Mid-covid" & sample == i)$dv,
                              paired = FALSE)), result)
}

result

result %>% 
  filter(p.value < .05) %>%
  count()



  