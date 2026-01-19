#Updated to work with GitHub on 1/6/26
#Load libraries
library(tidyverse)
library(car)
library(emmeans)
library(ggeffects)
library(reshape2)
library(Rcpp)
library(DHARMa)
library(colorspace)
library(ggpubr)
library(lme4)
library(TMB)
library(glmmTMB)
library(performance)
library(dplyr)
library(parsedate)
library(here)
library(ggtext)
library(scales)
library(patchwork)
library(dplyr)
library(forcats)
library(ggplot2)
library(FSA)
library(nlstools)
library(plotrix)
library(Matrix)
library(tidyr)
library(lmerTest)
library(survival)
library(survminer)
library(ggsurvfit)
library(survivalAnalysis)
library(broom)
library(forestmodel)
library(knitr)
library(kableExtra)
library(flextable)
library(modelsummary)

#set working directory to work within GitHub
setwd("~/SequentialCoinfection/data")
total.spores <- read.csv("Coinfection Total Spore Yield Data.csv")
View(total.spores)

#Replace the #DIV/0! with NA for the cells that contain that value
total.spores[total.spores == '#DIV/0!'] <- NA

#Split apart treatment number and rep number into separate columns.
total.spores[c('TreatmentGroup', 'Rep')] <- str_split_fixed(total.spores$sample.ID, '-', 2)
View(total.spores)

#Use if else statements to retrieve the infection and exposure status of animals in treatments 1-10.
#Create new columns that denote Past and Metsch infection status and tell R to assign either Yes or No depending on treatment group

#Past exposed column
total.spores <- mutate(total.spores, PastExposed = ifelse(TreatmentGroup == "T2" | TreatmentGroup == "T3" | TreatmentGroup == "T4" | TreatmentGroup == "T5" | TreatmentGroup == "T6", "Yes", "No"))

#Metsch Exposed column
total.spores <- mutate(total.spores, MetschExposed = ifelse(TreatmentGroup == "T1" | TreatmentGroup == "T2", "No", "Yes"))

#Create new columns denoting whether animals were infected with Past, Metsch, or Both
total.spores <- mutate(total.spores, PastInfected = ifelse(Past.spore.in.animal > 0, "Yes", "No"))
total.spores <- mutate(total.spores, MetschInfected = ifelse(Metsch.spore.in.animal > 0, "Yes", "No"))
total.spores <- mutate(total.spores, Coinfected = ifelse(PastInfected == "Yes" & MetschInfected == "Yes", "Yes", "No"))

#Create column that indicates whether an animal was unexposed, exposed and infected, or exposed and uninfected for Metsch
total.spores <- mutate(total.spores, MetExpStatus = ifelse(MetschExposed == "No", "Unexposed", ifelse(MetschInfected == "Yes", "Infected", "ExposedUninfected")))

#Create column that indicates whether an animal was unexposed, exposed and infected, or exposed and uninfected for Past
total.spores <- mutate(total.spores, PastExpStatus = ifelse(PastExposed == "No", "Unexposed", ifelse(PastInfected == "Yes", "Infected", "ExposedUninfected")))

#Add DayMet as a variable to the total.spores dataframe
#The day for treatments 1 and 2 will be input as "0" for now
total.spores <- mutate(total.spores, DayMetFactor = ifelse(TreatmentGroup == "T1" | TreatmentGroup == "T2", 0,
                                                           ifelse(TreatmentGroup == "T3" | TreatmentGroup == "T7", 5, 
                                                                  ifelse(TreatmentGroup == "T4" | TreatmentGroup == "T8", 10, 
                                                                         ifelse(TreatmentGroup == "T5" | TreatmentGroup == "T9", 15, 30)))))
str(total.spores)

#recode DayMet as a factor for this version of the data. 
total.spores$DayMetFactor <- as.factor(total.spores$DayMetFactor)

#Add a numeric DayMet column as well
total.spores <- mutate(total.spores, DayMetNum = ifelse(TreatmentGroup == "T1" | TreatmentGroup == "T2", 0,
                                                        ifelse(TreatmentGroup == "T3" | TreatmentGroup == "T7", 5, 
                                                               ifelse(TreatmentGroup == "T4" | TreatmentGroup == "T8", 10, 
                                                                      ifelse(TreatmentGroup == "T5" | TreatmentGroup == "T9", 15, 30)))))

total.spores$DayMetNum <- as.integer(total.spores$DayMetNum)

#Check to see if Treatment Group and Infection Status variables are factors and change them if not.
total.spores$TreatmentGroup <- as.factor(total.spores$TreatmentGroup)
total.spores$PastExposed <- as.factor(total.spores$PastExposed)
total.spores$PastInfected <- as.factor(total.spores$PastInfected)
total.spores$PastExpStatus <- as.factor(total.spores$PastExpStatus)
total.spores$MetschExposed <- as.factor(total.spores$MetschExposed)
total.spores$MetschInfected <- as.factor(total.spores$MetschInfected)
total.spores$MetExpStatus <- as.factor(total.spores$MetExpStatus)
total.spores$Coinfected <- as.factor(total.spores$Coinfected)

#Reorder treatment groups
total.spores$TreatmentGroup <- factor(total.spores$TreatmentGroup, c("T1","T2","T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"))
str(total.spores)

#Convert metsch/past average to numeric variables, and total spore counts to integers, rather than them being characters
total.spores$Metsch.avg <- as.numeric(total.spores$Metsch.avg)
total.spores$Past.avg <- as.numeric(total.spores$Past.avg)
total.spores$Metsch.spore.in.animal <- as.integer(total.spores$Metsch.spore.in.animal)
total.spores$Past.spore.in.animal <- as.integer(total.spores$Past.spore.in.animal)

str(total.spores)

#Subset the data to only include relevant columns from here on out
#Spore dataset with DayMet as a factor
spores.factor <- subset(total.spores, select = c("Metsch.spore.in.animal", "Past.spore.in.animal", "PastInfected", "PastExposed", "PastExpStatus", "MetschInfected", "MetschExposed", "MetExpStatus", "Coinfected", "TreatmentGroup", "Rep", "DayMetFactor"))

#Spore dataset with DayMet as numeric
spores.numeric <- subset(total.spores, select = c("Metsch.spore.in.animal", "Past.spore.in.animal", "PastInfected", "PastExposed", "PastExpStatus", "MetschInfected", "MetschExposed", "MetExpStatus", "Coinfected", "TreatmentGroup", "Rep", "DayMetNum"))

####Creating figure 2#####

#facet plot of infection prevalence and spores 
past <- subset(spores.factor, TreatmentGroup %in% c("T2", "T3", "T4", "T5", "T6"))
past$DayMetFactor <- as.character(past$DayMetFactor)
past$DayMetFactor[past$DayMetFactor == "0"] <- "never"
past$DayMetFactor <- as.factor(past$DayMetFactor)
past$DayMetFactor <- factor(past$DayMetFactor,levels = c("never", "5", "10", "15", "30"))

#first part of the facet plot
#finding prevalence of past infection 
past$PastInfected_numeric <- ifelse(past$PastInfected == "Yes", 1, 0)
past <- na.omit(past)

#I embrace our AI future because I could not have written this without chatGPT 
# Step 1: Create the base prevalence data
prevalence_by_coinfection <- past %>%
  filter(PastInfected_numeric == 1) %>%
  group_by(DayMetFactor, Coinfected) %>%
  summarise(
    infected_n = n(),
    .groups = "drop"
  ) %>%
  left_join(
    past %>%
      group_by(DayMetFactor) %>%
      summarise(total_n = n(), .groups = "drop"),
    by = "DayMetFactor"
  ) %>%
  mutate(
    prevalence = infected_n / total_n,
    se = sqrt(prevalence * (1 - prevalence) / total_n)
  )

# Step 2: Add a "sum" row per DayMetFactor
prevalence_sum <- prevalence_by_coinfection %>%
  group_by(DayMetFactor) %>%
  summarise(
    infected_n = sum(infected_n),
    total_n = first(total_n),  # same total for both groups
    prevalence = infected_n / total_n,
    se = sqrt(prevalence * (1 - prevalence) / total_n),
    Coinfected = "Sum",  # New level
    .groups = "drop"
  )

# Step 3: Combine original and sum data
prevalence_combined <- bind_rows(prevalence_by_coinfection, prevalence_sum)

# Optional: make Coinfected labels more readable
prevalence_combined$Coinfected <- factor(
  prevalence_combined$Coinfected,
  levels = c("No", "Yes", "Sum"),
  labels = c("single infection", "coinfection", "total"))


#Making this facet plot by stacked with coinfected vs singly infected 
f1 <-ggplot(prevalence_combined, aes(x = DayMetFactor, y = prevalence, color = Coinfected)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "single infection" = "cornflowerblue",
      "coinfection" = "#7f39d4",
      "total" = "black"
    ), labels = c(
      "single infection" = "single infection",
      "coinfection" = "coinfection",
      "total" = "overall")) +
  ylab("Infection prevalence of *P. ramosa*") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) + labs(color = NULL) +
  coord_cartesian(ylim = c(0, 1))

# Spore count summary
past_spores <- subset(past, PastInfected == "Yes")

#Creating a dataset that averages spores by treatment 
past_spores_summary <- past_spores %>%
  group_by(DayMetFactor, Coinfected) %>%
  summarise(
    n = n(),
    mean_spores = mean(Past.spore.in.animal, na.rm = TRUE),
    se_spores = sd(Past.spore.in.animal, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# Optional: make Coinfected labels more readable
past_spores_summary$Coinfected <- factor(
  past_spores_summary$Coinfected,
  levels = c("No", "Yes"),
  labels = c("single infection", "coinfection"))

#Past spores plotted 
f2 <-ggplot(past_spores_summary, aes(x = DayMetFactor, y = mean_spores, color = Coinfected)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = mean_spores - se_spores, ymax = mean_spores + se_spores), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "single infection" = "cornflowerblue",
      "coinfection" = "#7f39d4"), labels = c(
      "single infection" = "single infection",
      "coinfection" = "coinfection")) +
  ylab("Mean spores per host *P. ramosa*") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) + labs(color = NULL) + scale_x_discrete(drop = FALSE) + scale_y_continuous(labels = scales::label_comma())


#The same data as above,but in a violin plot 
#ggplot(past_spores, aes(x = Coinfected, y = Past.spore.in.animal, color = DayMetFactor, shape = MetschExposed)) +
#  geom_violin(color = "black", fill = NA, trim = FALSE) +
#  geom_jitter(width = 0.3, size = 3, alpha = 1, show.legend = TRUE) +
#  theme_classic() +
#  labs(
#    x = NULL,
#    y = "Number of *P. ramosa* spores per host",
#    color = NULL
#  ) +
#  theme(
#    axis.title.x = ggtext::element_markdown(), 
#    axis.title.y = ggtext::element_markdown(),
#    legend.text = ggtext::element_markdown(), 
#    legend.title = ggtext::element_markdown(),
#    legend.position = c(.01, 1),
#    legend.justification = c("left", "top"),
#    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
#  ) +
#  scale_color_manual(
#    values = c("never" = "grey50", "5" = "#c1b0d6", "10" = "#ac8cd4", "15" = "#7f39d4", "30" = "#640ecc")
#  ) +
#  scale_x_discrete(labels = c("No" = "single infection", "Yes" = "coinfection")) +
#  scale_y_continuous(limits = c(0, NA), labels = scales::comma)

#Same thing as the two plots above, but for metsch 
#finding prevalence of metsch infection 
metsch_with_singlyexposed <- subset(spores.factor, TreatmentGroup %in% c("T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"))
metsch_with_singlyexposed$DayMetFactor <- as.character(metsch_with_singlyexposed$DayMetFactor)
metsch_with_singlyexposed$DayMetFactor[metsch_with_singlyexposed$DayMetFactor == "0"] <- "never"
metsch_with_singlyexposed$DayMetFactor <- as.factor(metsch_with_singlyexposed$DayMetFactor)
metsch_with_singlyexposed$DayMetFactor <- factor(metsch_with_singlyexposed$DayMetFactor,levels = c("never", "5", "10", "15", "30"))

#third part of the facet plot
#finding prevalence of metsch infection 
metsch_with_singlyexposed$MetschInfected_numeric <- ifelse(metsch_with_singlyexposed$MetschInfected == "Yes", 1, 0)
metsch_with_singlyexposed <- na.omit(metsch_with_singlyexposed)

# Step 1: Create the base prevalence data
prevalence_by_coinfection_metsch_with_singlyexposed <- metsch_with_singlyexposed %>%
  filter(MetschInfected_numeric == 1) %>%
  group_by(DayMetFactor, Coinfected, PastExpStatus) %>%
  summarise(
    infected_n = n(),
    .groups = "drop"
  ) %>%
  left_join(
    past %>%
      group_by(DayMetFactor) %>%
      summarise(total_n = n(), .groups = "drop"),
    by = "DayMetFactor"
  ) %>%
  mutate(
    prevalence = infected_n / total_n,
    se = sqrt(prevalence * (1 - prevalence) / total_n)
  )

# Step 2: Add a "sum" row per DayMetFactor
prevalence_metsch_sum_with_singly_exposed <-
  prevalence_by_coinfection_metsch_with_singlyexposed %>%
  filter(PastExpStatus %in% c("ExposedUninfected", "Infected")) %>%
  group_by(DayMetFactor) %>%
  summarise(
    infected_n = sum(infected_n),
    total_n = first(total_n),
    prevalence = infected_n / total_n,
    se = sqrt(prevalence * (1 - prevalence) / total_n),
    Coinfected = "Sum",
    PastExpStatus = "Summed",
    .groups = "drop"
  )

# Step 3: Combine original and sum data
prevalence_combined_metsch_with_singly_exposed <- bind_rows(prevalence_by_coinfection_metsch_with_singlyexposed, prevalence_metsch_sum_with_singly_exposed)

# Optional: make Coinfected labels more readable
prevalence_combined_metsch_with_singly_exposed$Coinfected <- factor(
  prevalence_combined_metsch_with_singly_exposed$Coinfected,
  levels = c("No", "Yes", "Sum"),
  labels = c("single infection", "coinfection", "total"))

prevalence_combined_metsch_with_singly_exposed$PastExpStatus <- factor(
  prevalence_combined_metsch_with_singly_exposed$PastExpStatus,
  levels = c("ExposedUninfected", "Unexposed", "Infected", "Summed"),
  labels = c("exposed but uninfected", "unexposed", "infected", "overall (only coexposed treatments)"))

#Third plot for the facet plot
f3<-ggplot(prevalence_combined_metsch_with_singly_exposed, aes(x = DayMetFactor, y = prevalence, color = Coinfected, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "single infection" = "tomato",
      "coinfection" = "#7f39d4",
      "total" = "black"
    ), labels = c(
      "single infection" = "single infection",
      "coinfection" = "coinfection",
      "total" = "overall")) +
  ylab("Infection prevalence of *A. monospora*") +
  xlab("Day of exposure to *A. monospora*") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) + labs(color = NULL) +
  coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(drop = FALSE) + labs(shape = "Exposed to *P. ramosa*") + scale_shape_manual(
    values = c("exposed but uninfected" = 16, "unexposed" = 15, "infected" = 17, "overall (only coexposed treatments)" = 18))

#fourth part of the facet plot
# Spore count summary
metsch_spores <- subset(metsch_with_singlyexposed, MetschInfected == "Yes")

#Creating a dataset that averages spores by treatment 
metsch_spores_summary <- metsch_spores %>%
  group_by(DayMetFactor, Coinfected, PastExpStatus) %>%
  summarise(
    n = n(),
    mean_spores = mean(Metsch.spore.in.animal, na.rm = TRUE),
    se_spores = sd(Metsch.spore.in.animal, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# Optional: make Coinfected labels more readable
metsch_spores_summary$Coinfected <- factor(
  metsch_spores_summary$Coinfected,
  levels = c("No", "Yes"),
  labels = c("single infection", "coinfection"))

#Past spores plotted 
f4<-ggplot(metsch_spores_summary, aes(x = DayMetFactor, y = mean_spores, color = Coinfected, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = mean_spores - se_spores, ymax = mean_spores + se_spores), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "single infection" = "tomato",
      "coinfection" = "#7f39d4"), labels = c(
        "single infection" = "single infection",
        "coinfection" = "coinfection")) +
  ylab("Mean spores per host *A. monospora*") +
  xlab("Day of exposure to *A. monospora*") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, .4),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) + labs(color = NULL) + scale_x_discrete(drop = FALSE) + scale_y_continuous(labels = scales::label_comma()) + labs(shape = "Exposed to *P. ramosa*") +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16, "Unexposed" = 15, "Infected" = 17),
    labels = c(
      "ExposedUninfected" = "exposed but uninfected",
      "Unexposed" = "unexposed",
      "Infected" = "infected")) + guides(
        shape = guide_legend(order = 1),
        color = guide_legend(order = 2)) + ylim(0, 175000)

#same data as above, but in a violin plot 
#ggplot(metsch_spores, aes(x = Coinfected, y = Metsch.spore.in.animal, color = DayMetFactor)) +
#  geom_violin(color = "black", fill = NA, trim = FALSE) +
#  geom_jitter(width = 0.3, size = 3, alpha = 1, show.legend = TRUE, shape = 17) +
#  theme_classic() +
#  labs(
#    x = "Infection status",
#    y = "Number of *A. monospora* spores per host",
#    color = NULL
#  ) +
#  theme(
#    axis.title.x = ggtext::element_markdown(), 
#    axis.title.y = ggtext::element_markdown(),
#    legend.text = ggtext::element_markdown(), 
#    legend.title = ggtext::element_markdown(),
#    legend.position = c(.01, 1),
#    legend.justification = c("left", "top"),
#    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
#  ) +
#  scale_color_manual(
#    values = c("never" = "grey50", "5" = "#c1b0d6", "10" = "#ac8cd4", "15" = "#7f39d4", "30" = "#640ecc")
#  ) +
#  scale_x_discrete(labels = c("No" = "single infection", "Yes" = "coinfection")) +
#  scale_y_continuous(limits = c(0, NA), labels = scales::comma)

infectionplot <-
  (f1 | f2) / (f3 | f4) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = ggtext::element_markdown(size = 14),
    legend.title = ggtext::element_markdown(size = 14),
    strip.text = element_text(size = 18)
  )

infectionplot
ggsave("~/SequentialCoinfection/figures/infectionplot_withexposurelabels.png", infectionplot, dpi = 600, width = 15, height = 12, units = "in")


###Creating figure 3####
#FECUNDITY ANALYSIS
egg <- read.csv("Coinfection Egg Data.csv")
View(egg)

#Add column that totals lifetime reproductive output for each individual
egg$LifeRep <- rowSums(egg[ ,c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49)], na.rm = TRUE)

#Split the treatment group and rep, code lifetime reproductive output as an integer
egg[c('TreatmentGroup', 'Rep')] <- str_split_fixed(egg$sample.ID, '-', 2)

#Recode variables
egg$TreatmentGroup <- as.factor(egg$TreatmentGroup)
egg$LifeRep <- as.integer(egg$LifeRep)

#Add in exposure and infection status columns and daymet columns from total.spores dataset
total.egg <- cbind(egg, total.spores[c("PastExposed", "MetschExposed", "PastInfected", "MetschInfected", "PastExpStatus", "MetExpStatus", "DayMetNum", "DayMetFactor")])

#Data Set Up: manipulate the egg data to show survival times

#Here, I added in the dates that animals were exposed to Metsch and Past, and calculated the length of
#survival post exposure to each parasite. If an animal did not survive 10 days after being exposed to 
#Past, it was removed from the model analysis. If an animal did not survive 5 days after being exposed
#to Metsch, it was removed from the model analysis. This is because premature death makes 
#it hard/impossible to know if an animal was infected with anything or survived past the window
#of infection and was not infected. Only 5 animals were removed and all were uninfected according 
#to spore count data. Start date of the experiment was chosen as 10/14/2022, which was the day 
#before Past exposures when all animals were still living (this can be changed if necessary, 
#but it wasn't used the in the mortality analysis).

#Subset the data
egg.life <- subset(total.egg, select=c("TreatmentGroup", "Rep", "PastExposed", "MetschExposed", "PastInfected", "MetschInfected", "PastExpStatus", "MetExpStatus", "DayMetNum", "DayMetFactor", "LifeRep", "death.date"))

#Add in a censor column for the survival analysis, where 1 means censored and 2 means dead
egg.life <- egg.life %>% 
  mutate(censor = ifelse(death.date == "", 1, 2))

#Add in the censor date for animals that did not die before the last day of the experiment
egg.life$death.date[egg.life$death.date == ""] <- "12/19/2022"

#Convert death date to a date in r
egg.life <- egg.life %>%
  mutate(death.date = as.Date(death.date, format = "%m/%d/%Y"))
str(egg.life)

#Add a column with experiment start date, which I set as the day before past exposure for now
egg.life <- egg.life %>%
  mutate(egg.life, start.date = "10/14/2022")
str(egg.life)

#Convert experiment start date to a date in r
egg.life <- egg.life %>%
  mutate(start.date = as.Date(start.date, format = "%m/%d/%Y"))

#Add a column with length of survival time in days from start of experiment
egg.life <- egg.life %>%
  mutate(lifesurvival = as.duration(start.date %--% death.date) / ddays(1))

#Subset the data to just treatments 3-6 for now
#egg.life.3to6 <- subset(egg.life, PastExposed == "Yes" & DayMetNum != "0" & DayMetFactor != "0" & PastInfected != "NA")

#Add in Past exposure date column
egg.life <- egg.life %>%
  mutate(PastExposureDate = "10/15/2022")
str(egg.life)

#Make past exposure date a date in r
egg.life <- egg.life %>%
  mutate(PastExposureDate = as.Date(PastExposureDate, format = "%m/%d/%Y"))

#Add in Metsch exposure date columns
egg.life <- egg.life %>%
  mutate(MetschExposureDate = ifelse(DayMetNum == "5", "10/20/2022", 
                                     ifelse(DayMetNum == "10", "10/25/2022",
                                            ifelse(DayMetNum == "15", "10/30/2022", "11/14/2022"))))

#Make metsch exposure dates into dates in r
egg.life <- egg.life %>%
  mutate(MetschExposureDate = as.Date(MetschExposureDate, format = "%m/%d/%Y"))
str(egg.life)

#Add in days of survival after exposure to past column
egg.life <- egg.life %>%
  mutate(PostPastSurvival = as.duration(PastExposureDate %--% death.date) / ddays(1))

#Add in days of survival after exposure to metsch column
egg.life <- egg.life %>%
  mutate(PostMetSurvival = as.duration(MetschExposureDate %--% death.date) / ddays(1))

#Subset the data to only include animals that made it a certain number of days past exposure to Past and Metsch (e stands for egg, l stands for life in the dataframe name)
el.final <- subset(egg.life, lifesurvival > 10)

#Add a factor column that has infection status as either None, Past, Metsch, or Coinf. This will be useful for graphing and model analysis.
el.final <- el.final %>%
  mutate(InfStatusFactor = ifelse(PastInfected == "Yes" & MetschInfected == "Yes", "Coinf",
                                  ifelse(PastInfected == "Yes", "Past", 
                                         ifelse(MetschInfected == "Yes", "Metsch", "Uninf"))))
el.final$InfStatusFactor <- as.factor(el.final$InfStatusFactor)
el.final$InfStatusFactor <- factor(el.final$InfStatusFactor, c("Uninf", "Past", "Metsch", "Coinf"))

el.final$DayMetFactor <- factor(el.final$DayMetFactor)
el.final$DayMetFactor <- relevel(el.final$DayMetFactor, "0")

#Never exposed get called uninfected 
el.final <- el.final %>%
  mutate(
    InfStatusFactor = if_else(
      TreatmentGroup == "T1",
      fct_na_value_to_level(InfStatusFactor, "Uninf"),  # convert NAs to "Uninf" only for T1
      InfStatusFactor  # keep the original value for other treatments
    )
  )

el.final <- el.final %>%
  mutate(
    PastExpStatus = if_else(
      TreatmentGroup == "T1",
      fct_na_value_to_level(PastExpStatus, "Unexposed"),  # convert NAs to "Uninf" only for T1
      PastExpStatus  # keep the original value for other treatments
    )
  )


el.final <- el.final %>%
  mutate(
    MetExpStatus = if_else(
      TreatmentGroup == "T1",
      fct_na_value_to_level(MetExpStatus, "Unexposed"),  # convert NAs to "Uninf" only for T1
      MetExpStatus  # keep the original value for other treatments
    )
  )

#Finding the mean lifetime fecundity for each treatment 
el.final_summary <- el.final %>%
  group_by(DayMetFactor, InfStatusFactor, MetExpStatus, PastExpStatus) %>%
  summarise(
    n = n(),
    total_fecundity = mean(LifeRep, na.rm = TRUE),
    se_fecundity = sd(LifeRep, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

clean_el.final_summary <- na.omit(el.final_summary)

levels(clean_el.final_summary$DayMetFactor) <- 
  sub("^0$", "never", levels(clean_el.final_summary$DayMetFactor))

#Creating a facet plot based on terminal infection 
uninf_egg_summary <- clean_el.final_summary %>%
  filter(InfStatusFactor == "Uninf")

e1<-ggplot(uninf_egg_summary, aes(x = DayMetFactor, y = total_fecundity, color = MetExpStatus, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = total_fecundity - se_fecundity, ymax = total_fecundity + se_fecundity), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "Unexposed" = "black",
      "ExposedUninfected" = "green3"), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed")) +
  ylab("Mean total fecundity per host") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, .45),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    color = "Exposed to *A. monospora*",
    shape = "Exposed to *P. ramosa*",
    title = "Uninfected"   # <-- add this line
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed")) + theme(plot.title = element_text(hjust = 0.5)) + ylim(0, 150)

#Creating the second graph for the egg facet plot 
past_egg_summary <- clean_el.final_summary %>%
  filter(InfStatusFactor == "Past")

e2<-ggplot(past_egg_summary, aes(x = DayMetFactor, y = total_fecundity, shape = MetExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  geom_errorbar(aes(ymin = total_fecundity - se_fecundity, ymax = total_fecundity + se_fecundity), 
                width = 0.1, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    shape = "Exposed to *A. monospora*",
    title = "*P. ramosa*"   # <-- add this line
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed")) + ylim(0,150)

#Creating the third graph for the egg facet plot 
metsch_egg_summary <- clean_el.final_summary %>%
  filter(InfStatusFactor == "Metsch")

e3 <-ggplot(metsch_egg_summary, aes(x = DayMetFactor, y = total_fecundity, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "tomato") +
  geom_errorbar(aes(ymin = total_fecundity - se_fecundity, ymax = total_fecundity + se_fecundity), 
                width = 0.1, position = position_dodge(width = 0.4), color = "tomato") +
  ylab("Mean total fecundity per host") +
  xlab("Day of exposure to *A. monospora*") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
     plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    shape = "Exposed to *P. ramosa*",
    title = "*A. monospora*"   # <-- add this line
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma()) +
scale_shape_manual(
  values = c(
    "ExposedUninfected" = 16,
    "Unexposed" = 15), labels = c(
      "ExposedUninfected" = "exposed but uninfected",
      "Unexposed" = "unexposed")) + ylim(0,150)

#Creating the fourth graph for the egg facet plot 
coinf_egg_summary <- clean_el.final_summary %>%
  filter(InfStatusFactor == "Coinf")

e4<-ggplot(coinf_egg_summary, aes(x = DayMetFactor, y = total_fecundity)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "#7f39d4") +
  geom_errorbar(aes(ymin = total_fecundity - se_fecundity, ymax = total_fecundity + se_fecundity), 
                width = 0.1, position = position_dodge(width = 0.4), color = "#7f39d4") +
  ylab(NULL) +
  xlab("Day of exposure to *A. monospora*") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  )  + scale_x_discrete(drop = FALSE) + scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title = "Coinfected"   # <-- add this line
  ) + ylim(0,150)

eggplot<- (e1 | e2) / (e3 | e4) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = ggtext::element_markdown(size = 14),
    legend.title = ggtext::element_markdown(size = 14),
    strip.text = element_text(size = 18)
  )
eggplot
ggsave("~/SequentialCoinfection/figures/eggplot.png", eggplot, dpi = 600, width = 15, height = 12, units = "in")

####CREATING FIGURE 4#####
##Survival 
#Finding the mean lifetime fecundity for each treatment 
survival_summary <- el.final %>%
  group_by(DayMetFactor, InfStatusFactor, MetExpStatus, PastExpStatus) %>%
  summarise(
    n = n(),
    mean_survival = mean(lifesurvival, na.rm = TRUE),
    se_survival = sd(lifesurvival, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

clean_survival_summary <- na.omit(survival_summary)

levels(clean_survival_summary$DayMetFactor) <- 
  sub("^0$", "never", levels(clean_survival_summary$DayMetFactor))

#Creating a facet plot based on terminal infection 
uninf_survival_summary <- clean_survival_summary %>%
  filter(InfStatusFactor == "Uninf")

s1<-ggplot(uninf_survival_summary, aes(x = DayMetFactor, y = mean_survival, color = MetExpStatus, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "Unexposed" = "black",
      "ExposedUninfected" = "green3"), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed")) +
  ylab("Mean lifespan of host in days") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = c(.01, .4),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    color = "Exposed to *A. monospora*",
    shape = "Exposed to *P. ramosa*",
    title = "Uninfected"   # <-- add this line
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0,75)) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed"))

past_survival_summary <- clean_survival_summary %>%
  filter(InfStatusFactor == "Past")

s2<-ggplot(past_survival_summary, aes(x = DayMetFactor, y = mean_survival, shape = MetExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), 
                width = 0.1, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  ylab(NULL) +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    shape = "Exposed to *A. monospora*",
    title = "*P. ramosa*"   # <-- add this line
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0,75)) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed")) 

metsch_survival_summary <- clean_survival_summary %>%
  filter(InfStatusFactor == "Metsch")

s3<-ggplot(metsch_survival_summary, aes(x = DayMetFactor, y = mean_survival, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "tomato") +
  geom_errorbar(aes(ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), 
                width = 0.1, position = position_dodge(width = 0.4), color = "tomato") +
  ylab("Mean lifespan of host in days") +
  xlab("Day of exposure to *A. monospora*") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    shape = "Exposed to *P. ramosa*",
    title = "*A. monospora*"   # <-- add this line
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0,75)) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed"))

coinf_survival_summary <- clean_survival_summary %>%
  filter(InfStatusFactor == "Coinf")

s4<-ggplot(coinf_survival_summary, aes(x = DayMetFactor, y = mean_survival)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "#7f39d4") +
  geom_errorbar(aes(ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), 
                width = 0.1, position = position_dodge(width = 0.4), color = "#7f39d4") +
  ylab(NULL) +
  xlab("Day of exposure to *A. monospora*") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, .4),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    title = "Coinfected"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0,75)) +
 theme(plot.title = element_text(hjust = 0.5))

survplot <- (s1 | s2) / (s3 | s4) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = ggtext::element_markdown(size = 14),
    legend.title = ggtext::element_markdown(size = 14),
    strip.text = element_text(size = 18)
  )
survplot
ggsave("~/SequentialCoinfection/figures/survplot_test.png", survplot, dpi = 600, width = 15, height = 12, units = "in")

#Survivor analysis - Hazard Ratio

#subsetting for only infected
#infected_el.final <- subset(el.final, !PastExpStatus == "ExposedUninfected")
#infected_el.final <- subset(infected_el.final, !MetExpStatus == "ExposedUninfected")

#censor animals that did not die
#infected_el.final$survobject <- with(infected_el.final, Surv(lifesurvival, censor))
#seq_lifespan <- survfit(survobject ~ InfStatusFactor, data = infected_el.final, conf.type = "log-log")
#summary(seq_lifespan)

#Censor animals who lived to the end 
el.final$survobject <- with(el.final, Surv(lifesurvival, censor))
#seq_lifespan <- survfit(survobject ~ InfStatusFactor, data = el.final, conf.type = "log-log")
#summary(seq_lifespan)

#Run Cox proportional hazards with interaction variables 
cox.seqlifespan <- coxph(survobject ~ InfStatusFactor*DayMetFactor, data = el.final)
summary(cox.seqlifespan)
cox.test <- cox.zph(cox.seqlifespan)
cox.test

#I think this is showing that the hazards aren't exactly proportional across all times, but it should be fine? 
#Cox is supposed to be fairly robust
plot(cox.test)

#Plot it out, not final plot 
ggforest(cox.seqlifespan, data = el.final)

# Get tidy model results (HRs, CIs, p-values, etc.)
tidy_fit <- tidy(cox.seqlifespan, exponentiate = TRUE, conf.int = TRUE)

# View all terms
print(tidy_fit$term)

#This can't be the fastest way but it is MY way
subset_terms <- tidy_fit %>%
  dplyr::filter(term %in% c(
    "InfStatusFactorPast",
    "InfStatusFactorMetsch", "InfStatusFactorCoinf", "DayMetFactor5", "DayMetFactor10", "DayMetFactor15",
    "DayMetFactor30", "InfStatusFactorCoinf:DayMetFactor5", "InfStatusFactorCoinf:DayMetFactor10", "InfStatusFactorCoinf:DayMetFactor15",
    "InfStatusFactorCoinf:DayMetFactor30"
  ))

#Reordering the terms and adding a p-value atericks column 
subset_terms$term <- factor(subset_terms$term, levels = c(
  "InfStatusFactorPast",
  "InfStatusFactorMetsch",
  "InfStatusFactorCoinf",
  "DayMetFactor5",
  "DayMetFactor10",
  "DayMetFactor15",
  "DayMetFactor30",
  "InfStatusFactorCoinf:DayMetFactor5",
  "InfStatusFactorCoinf:DayMetFactor10",
  "InfStatusFactorCoinf:DayMetFactor15",
  "InfStatusFactorCoinf:DayMetFactor30"
))

subset_terms<- subset_terms %>%
  filter(term != "InfStatusFactorCoinf:DayMetFactor30")

subset_terms <- subset_terms |>
  dplyr::mutate(sig_label = dplyr::case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

#Making my forest plot, this gets a depreciation warning because I am not specifying orientation of bars, but the code is running fine currently  
ggplot(subset_terms, aes(y = fct_rev(term), x = estimate)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  
  # Adjusted text for log scale
  geom_text(aes(
    x = conf.high * 1.05,  # use multiplication for log scale
    label = ifelse(p.value < 0.05,
                   paste0(sig_label, "\n(p=", signif(p.value, 2), ")"),
                   "")
  ),
  size = 3.5, hjust = 0, vjust = 0.5
  ) +
  
  # Add annotations for interpretation
  annotate("text", x = 0.7, y = Inf, label = "increased lifespan",
           hjust = 1, vjust = 1.2, size = 4, fontface = "italic") +
  annotate("text", x = 1.3, y = Inf, label = "decreased lifespan",
           hjust = 0, vjust = 1.2, size = 4, fontface = "italic") +
  
  # Italic and custom y-axis labels
  scale_y_discrete(labels = c(
    "InfStatusFactorPast" = expression(italic("P. ramosa")),
    "InfStatusFactorMetsch" = expression(italic("A. monospora")),
    "InfStatusFactorCoinf" = "Coinfected",
    "DayMetFactor5" = expression(italic("A. monospora") * " added on day 5"),
    "DayMetFactor10" = expression(italic("A. monospora") * " added on day 10"),
    "DayMetFactor15" = expression(italic("A. monospora") * " added on day 15"),
    "DayMetFactor30" = expression(italic("A. monospora") * " added on day 30"),
    "InfStatusFactorCoinf:DayMetFactor5" = expression("Coinfected: " * italic("A. monospora") * " added on day 5"),
    "InfStatusFactorCoinf:DayMetFactor10" = expression("Coinfected: " * italic("A. monospora") * " added on day 10"),
    "InfStatusFactorCoinf:DayMetFactor15" = expression("Coinfected: " * italic("A. monospora") * " added on day 15"),
    "InfStatusFactorCoinf:DayMetFactor30" = expression("Coinfected: " * italic("A. monospora") * " added on day 30")
  )) +
  
  # Log scale for x-axis with slight left expansion
  scale_x_log10(
    expand = expansion(mult = c(0.05, 0))
  ) +
  
  labs(x = "Hazard Ratio (log scale)", y = NULL) +
  theme_classic(base_size = 14)

ggsave("~/SequentialCoinfection/figures/el.final.survival.tiff", dpi = 600, width = 10.5, height = 6, units = "in", compression="lzw")
ggsave("~/SequentialCoinfection/figures/el.final.survival.png", dpi = 600, width = 10.5, height = 6, units = "in")


#subset for just coinfected? 
#el.final_justcoinf$survobject <- with(el.final_justcoinf, Surv(lifesurvival, censor))
#seq_lifespan <- survfit(survobject ~ DayMetFactor, data = el.final_justcoinf, conf.type = "log-log")
#summary(seq_lifespan)

#cox.justcoinf <- coxph(survobject ~ DayMetFactor, data =el.final_justcoinf)
#summary(cox.justcoinf)

#Plot it out 
#ggforest(cox.justcoinf, data = el.final_justcoinf)

#ggsurvplot(seq_lifespan, data = el.final_justcoinf, xlab = "lifespan", facet.by = "DayMetFactor",
#                                   conf.type = "log-log", conf.int = T, pval = TRUE,
#                                   pval.method = TRUE)

#Subset for just exposed but uninfected 
#uninf <- subset(el.final, InfStatusFactor == "Uninf")

#uninf$survobject <- with(uninf, Surv(lifesurvival, censor))
#uninf_lifespan <- survfit(survobject ~ PastExposed + MetschExposed, data = uninf, conf.type = "log-log")
#summary(uninf_lifespan)

#cox.uninflifespan <- coxph(survobject ~ PastExposed * MetschExposed, data = uninf)
#summary(cox.uninflifespan)
#cox.test <- cox.zph(cox.uninflifespan)
#cox.test
#plot(cox.test)

##### Figure 5 ######
#Body size analysis 
#Read in the body size data.
length <- read.csv("Coinfection Body Size_without_continous.csv")

length[c('TreatmentGroup', 'Rep')] <- str_split_fixed(length$sample.ID, '-', 2)
length$TreatmentGroup < - as.factor(length$TreatmentGroup)
length <- length %>%
  mutate(date1 = as.Date(date1, format = "%m/%d/%Y"), 
         date2 = as.Date(date2, format = "%m/%d/%Y"),
         date3 = as.Date(date3, format = "%m/%d/%Y"),
         date4 = as.Date(date4, format = "%m/%d/%Y"),
         date5 = as.Date(date5, format = "%m/%d/%Y"),
         date6 = as.Date(date6, format = "%m/%d/%Y"),
         date7 = as.Date(date7, format = "%m/%d/%Y"),
         date8 = as.Date(date8, format = "%m/%d/%Y"),
         date9 = as.Date(date9, format = "%m/%d/%Y"),
         past_exposed = as.Date(past_exposed, format = "%m/%d/%Y"),
         metsch_exposed = as.Date(metsch_exposed, format = "%m/%d/%Y"),
         expstart = as.Date(expstart, format = "%m/%d/%Y"))


#Try sub-setting the data so that only the size measurements (and not the dates on which they were measured) 
#are shown and view it again. nd stands for "no dates".
length.nd <- length[, c(1, 2, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 24, 25, 26)]
#Try reformatting the body size data, so that the size measurements are all in one column labeled "measurement". 
#Long.data is the new name because the type of data we are converting to is called "long data".
long.data <- melt(data = length.nd, 
                  id.vars = c("sample.ID", "TreatmentGroup", "Rep", "Terminal.infection", "past_exposed", "metsch_exposed", "expstart"), 
                  measured.vars = c("size1", "size2", "size3", "size4", "size5", "size6", "size7", "size8","size9"),
                  variable.name = "date", value.name = "bodysize")

#Create a real date column so dates aren't listed as size1, size2, etc. 
long.data$realdate <- with(long.data, ifelse(date == "size1", "10/18/2022", ifelse(date == "size2", "10/25/2022", ifelse(date == "size3", "11/01/2022", ifelse(date == "size4", "11/08/2022", ifelse(date == "size5", "11/15/2022", ifelse(date == "size6", "11/22/2022", ifelse(date == "size7", "11/29/2022", ifelse(date == "size8", "12/06/2022", "12/13/2022")))))))))

#Change date to be recognized as a date and not a character 
long.data <- long.data %>%
  mutate(realdate = as.Date(realdate, format = "%m/%d/%Y"))

#Add DayMet and Exposure Status as columns to dataframe
long.data$Coexposed <- with(long.data, ifelse(TreatmentGroup == "T3" | TreatmentGroup == "T4" | TreatmentGroup == "T5" | TreatmentGroup == "T6", "Yes", "No"))
long.data$DayMet <- with(long.data, ifelse(TreatmentGroup == "T1" | TreatmentGroup == "T2", "None",
                                           ifelse(TreatmentGroup == "T3" | TreatmentGroup == "T7", 5, 
                                                  ifelse(TreatmentGroup == "T4" | TreatmentGroup == "T8", 10, 
                                                         ifelse(TreatmentGroup == "T5" | TreatmentGroup == "T9", 15, 30)))))

long.data$Coexposed <- as.factor(long.data$Coexposed)
long.data$DayMet <- as.factor(long.data$DayMet)
long.data$DayMet <- factor(long.data$DayMet, c("None", "5", "10", "15", "30"))
long.data$TreatmentGroup <- factor(long.data$TreatmentGroup, c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"))
str(long.data)

#Adding in CSCAR suggested things i.e. making treatments continuous variables and log transforming body size 
long.data$days_since_metsch_exposure <- (long.data$realdate - long.data$metsch_exposed)
long.data$days_since_metsch_exposure <- (as.numeric(long.data$days_since_metsch_exposure))
long.data$days_since_metsch_exposure[is.na(long.data$days_since_metsch_exposure)] <- 0
long.data$days_since_past_exposure <- (long.data$realdate - long.data$past_exposed)
long.data$days_since_past_exposure <- (as.numeric(long.data$days_since_past_exposure))
long.data$days_since_past_exposure[is.na(long.data$days_since_past_exposure)] <- 0
long.data$days_diff <- (long.data$days_since_past_exposure - long.data$days_since_metsch_exposure)
long.data$time <-(long.data$realdate - long.data$expstart)
long.data$time <- (as.numeric(long.data$time))
long.data[long.data<0] <- 0

#transforming into a linear relationship 
long.data$transformedbody <- (log10(1 -long.data$bodysize/2185.39))
long.data$transformedbody <- (long.data$transformedbody*-1)

#adding a unique identifier to each daphnia for random intercept 
long.data <- transform(long.data, uniqueident = as.numeric(factor(sample.ID)))

#drop any body size that is NA 
long.data <- long.data %>% drop_na(transformedbody)

#looking around for outliers 
ggplot(long.data, aes(x=time, y=bodysize)) +  geom_point()

#This transformation looks pretty good 
ggplot(long.data, aes(x=time, y=transformedbody)) +  geom_point()

ggplot(long.data, aes(x=days_since_metsch_exposure, y=transformedbody)) +  geom_point()
ggplot(long.data, aes(x=days_since_past_exposure, y=transformedbody)) +  geom_point()

ggplot(long.data, aes(x=days_since_past_exposure, y=transformedbody, shape=TreatmentGroup, color=TreatmentGroup)) +
  geom_point() + geom_smooth(se = FALSE)

long.data$Terminal.infection <- as.factor(long.data$Terminal.infection)
long.data$Terminal.infection <- relevel(long.data$Terminal.infection, "none")

#Adding breaks so the plots match 
common_breaks <- seq(0, 65, by = 15)

scale_x_continuous(
  limits = c(0, 30),
  breaks = common_breaks
)

#Body size vs exposure 
long.data_noinfections <- subset(long.data,Terminal.infection=="none")
long.data_noinfections <- long.data_noinfections %>%
  mutate(
    exposure = case_when(
      TreatmentGroup == "T1" ~ "never_exposed",
      TreatmentGroup == "T2" ~ "only_past",
      TreatmentGroup %in% c("T3", "T4", "T5", "T6") ~ "coexposed",
      TreatmentGroup %in% c("T7", "T8", "T9", "T10") ~ "only_metsch",
      TRUE ~ NA_character_
    )
  )

l1<-ggplot(long.data_noinfections, aes(x=time, y=transformedbody, color=exposure)) + geom_point(alpha = .4) + 
  geom_smooth(
    se = FALSE,
    method = "gam",
    formula = y ~ s(x, k = 4)
  ) +
  theme_minimal() + 
  xlab("") + ylab("Transformed body size") +
   scale_color_manual(
    name = "Exposure status",
    values = c(
      "never_exposed" = "green3",
      "only_past" = "cornflowerblue",
      "only_metsch" = "tomato",
      "coexposed" = "#7f39d4"),
    labels = c(
      "never_exposed" = "never exposed",
      "only_past" = "*P. ramosa* exposed",
      "only_metsch" = "*A. monospora* exposed",
      "coexposed" = "coexposed"),  breaks = c("never_exposed", "coexposed", "only_metsch", "only_past")) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.3)) +
scale_x_continuous(
  limits = c(0, 65),
  breaks = common_breaks
)

#Graph showing body size plotted against terminal infection 
l2<-ggplot(long.data,
       aes(x = time,
           y = transformedbody,
           color = Terminal.infection)) +
  geom_point(alpha = 0.4) +
  geom_smooth(
    se = FALSE,
    method = "gam",
    formula = y ~ s(x, k = 4)
  ) +
  scale_color_manual(
    name = "Infection status",
    values = c(
      "none" = "green3",
      "past" = "cornflowerblue",
      "metsch" = "tomato",
      "coinfected" = "#7f39d4"
    ),
    labels = c(
      "none" = "uninfected",
      "past" = "*P. ramosa*",
      "metsch" = "*A. monospora*",
      "coinfected" = "coinfected"
    )
  ) +
  xlab("Time in days") +
  ylab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.3
    )
  ) +
  scale_x_continuous(
    limits = c(0, 65),
    breaks = common_breaks
  )

#graph about does metsch through time change size in coinfected individuals 
coinfected_bodysize <- subset(long.data, Terminal.infection=="coinfected")

l3<-ggplot(coinfected_bodysize, aes(x=time, y=transformedbody, color=DayMet)) + geom_point(alpha = .4) +
  geom_smooth(
    se = FALSE,
    method = "gam",
    formula = y ~ s(x, k = 4)) + 
  xlab("Time in days") + ylab("Transformed body size")  +
  scale_color_manual(
    name = "Day of exposure to *A. monospora*", values = c("None" = "grey50", "5" = "mediumpurple1", "10" = "mediumpurple2", "15" = "mediumpurple3", "30" = "mediumpurple4")) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.3
    )) +
  scale_x_continuous(
    limits = c(0, 65),
    breaks = common_breaks
  )

l_blank <- patchwork::plot_spacer()

bodysize <- (l1 | l2) / (l3 | l_blank)+
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.text = ggtext::element_markdown(size = 14),
    legend.title = ggtext::element_markdown(size = 14),
    strip.text = element_text(size = 18)
  )

bodysize
ggsave("~/SequentialCoinfection/figures/bodysize.png", bodysize, dpi = 600, width = 15, height = 12, units = "in")

####Figure 6######
####A giant table of models

#Overdispersion check taken from Michelle's resource quality code 
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#Does infection prevalence of past in coinfected treatments differ in comparison to control? NOTE: coinfected are considered as 
#having past in this model 
#past_only.glm <- glm(PastInfected_numeric ~ DayMetFactor, family = binomial, data = past)
#summary(past_only.glm)

#Does the result change, as an ANOVA? - No 
#past_only_anova <- anova(past_only.glm, test = "Chisq")
#past_only_anova

#All checks taken from Michelle's data. Looks good. 
#overdisp_fun(past_only.glm)

#testDispersion(past_only.glm)
#testZeroInflation(past_only.glm)

#past_only.glm_simResid <- simulateResiduals(fittedModel = past_only.glm)
#plot(past_only.glm_simResid)

#Doing the same analysis, but with coinfected not considered past infected
#past_coinfectednotpast <- past %>%
#  mutate(
#    PastInfected_numeric = if_else(
#      Coinfected == "Yes",
#      0,
#      PastInfected_numeric
#    )
#  )

#past_coinfectednotpast.glm <- glm(PastInfected_numeric ~ DayMetFactor, family = binomial, data = past_coinfectednotpast)
#summary(past_coinfectednotpast.glm)

#Does the result change, as an ANOVA? - No 
#past_coinfectednotpast.glm_anova <- anova(past_coinfectednotpast.glm, test = "Chisq")
#past_coinfectednotpast.glm_anova

#All checks taken from Michelle's data. Looks good. 
#overdisp_fun(past_coinfectednotpast.glm)

#testDispersion(past_coinfectednotpast.glm)
#testZeroInflation(past_coinfectednotpast.glm)

#past_coinfectednotpast.glm_simResid <- simulateResiduals(fittedModel = past_coinfectednotpast.glm)
#plot(past_coinfectednotpast.glm)

past_prevalence.glm <- glm(PastInfected_numeric ~ DayMetFactor+Coinfected, family = binomial, data = past)
summary(past_prevalence.glm)

overdisp_fun(past_prevalence.glm)

testDispersion(past_prevalence.glm)
testZeroInflation(past_prevalence.glm)

past_prevalence.glm_simResid <- simulateResiduals(fittedModel = past_prevalence.glm)
plot(past_prevalence.glm_simResid)

past_prevalence_anova <- anova(past_prevalence.glm, test = "Chisq")
past_prevalence_anova

#Now making a model for past and spore yield 
#There is high variance leading to overdispersion
mean(past_spores$Past.spore.in.animal)
var(past_spores$Past.spore.in.animal)

past_spore.glm <- glm(Past.spore.in.animal ~ DayMetFactor*Coinfected, data = past_spores)
summary(past_spore.glm)

#Michelle's overdispersion function shows this as being overdispersed but the second overdispersion test doesn't. Is there a case to 
#be made to switch to quassipoisson? 
overdisp_fun(past_spore.glm)

testDispersion(past_spore.glm)
testZeroInflation(past_spore.glm)

past_spore.glm_simResid <- simulateResiduals(fittedModel = past_spore.glm)
plot(past_spore.glm_simResid)

#Adding post-hoc emmeans 
past_spore_emm <- emmeans(
  past_spore.glm,
  ~ Coinfected | DayMetFactor,   # Compare Coinfected vs not, *within each treatment*
  type = "response"              # Converts log-odds to probabilities (prevalence)
)

summary(past_spore_emm)
#Does infection prevalence of metsch in coinfected treatments differ in comparison to control (also adding in past exposure)?
metsch_prevalence.glm <- glm(MetschInfected_numeric ~ DayMetFactor+Coinfected+PastExposed, family = binomial, data = metsch_with_singlyexposed)
summary(metsch_prevalence.glm)

overdisp_fun(metsch_prevalence.glm)

testDispersion(metsch_prevalence.glm)
testZeroInflation(metsch_prevalence.glm)

metsch_prevalence.glm_simResid <- simulateResiduals(fittedModel = metsch_prevalence.glm)
plot(metsch_prevalence.glm_simResid)

#Adding post-hoc for past exposure 
emm_metsch_prevalence <- emmeans(
  metsch_prevalence.glm,
  ~ PastExposed | DayMetFactor,
  type = "response"
)

summary(emm_metsch_prevalence)

#Metsch spores analysis 
metsch_spore.glm <- glm(Metsch.spore.in.animal ~ DayMetFactor*Coinfected*PastExposed, family = gaussian, data = metsch_spores)
summary(metsch_spore.glm)

overdisp_fun(metsch_spore.glm)

testDispersion(metsch_spore.glm)
testZeroInflation(metsch_spore.glm)

metsch_spore.glm_simResid <- simulateResiduals(fittedModel = metsch_spore.glm)
plot(metsch_spore.glm_simResid)

#Adding post-hoc for past exposure 
emm_metsch_spores <- emmeans(
  metsch_spore.glm,
  ~ PastExposed | DayMetFactor | Coinfected,
  type = "response"
)

summary(emm_metsch_spores)

#Does the prevalence of coinfection change over time? 
onlycoinfectiontreatments <- subset(spores.factor, TreatmentGroup %in% c("T3", "T4", "T5", "T6"))
levels(onlycoinfectiontreatments$Coinfected)[levels(onlycoinfectiontreatments$Coinfected) == "Yes"] <- "1"
levels(onlycoinfectiontreatments$Coinfected)[levels(onlycoinfectiontreatments$Coinfected) == "No"] <- "0"

hist(onlycoinfectiontreatments$Metsch.spore.in.animal)
hist(onlycoinfectiontreatments$Past.spore.in.animal)

onlycoinfectiontreatments.glm <- glm(Coinfected ~ DayMetFactor, family = binomial, data = onlycoinfectiontreatments)
summary(onlycoinfectiontreatments.glm)

overdisp_fun(onlycoinfectiontreatments.glm)

testDispersion(onlycoinfectiontreatments.glm)
testZeroInflation(onlycoinfectiontreatments.glm)

onlycoinfectiontreatments.glm_simResid <- simulateResiduals(fittedModel = onlycoinfectiontreatments.glm)
plot(onlycoinfectiontreatments.glm_simResid)

coinf_anova <- anova(onlycoinfectiontreatments.glm, test = "Chisq")
coinf_anova

#Moving on to fecundity
totalegg.glm <- glmmTMB(LifeRep ~ InfStatusFactor + lifesurvival + DayMetFactor, family = nbinom1(), data = el.final)
summary(totalegg.glm)

testDispersion(totalegg.glm)
testZeroInflation(totalegg.glm)

totalegg.glm_simResid <- simulateResiduals(fittedModel = totalegg.glm)
plot(totalegg.glm_simResid) # there are some patterns in the residual, but this is the best model 

#I don't really understand emmeans but here goes nothing 
emm_egg <- emmeans(totalegg.glm, ~ DayMetFactor * InfStatusFactor * lifesurvival) 
summary(emm_egg)

#totalegg.glm.df <- emmeans(totalegg.glm, specs = pairwise ~ InfStatusFactor | DayMetFactor | lifesurvival, type = "response")
#totalegg.glm.df
#coinfmetspores.emm.df <- as.data.frame(coinfmetspores.glm.df$emmeans)

plot(emm_egg)

#What's is up with the 30 day metsch? 
#justmetsch <- subset(spores.factor, TreatmentGroup %in% c("T7" ,"T8", "T9", "T10"))
#levels(justmetsch$MetschInfected)[levels(justmetsch$MetschInfected) == "Yes"] <- "1"
#levels(justmetsch$MetschInfected)[levels(justmetsch$MetschInfected) == "No"] <- "0"

#justmet_prop.glm <- glm(MetschInfected ~ DayMetFactor, family = binomial, data = justmetsch)
#summary(justmet_prop.glm)

#plot(justmet_prop.glm)

#testDispersion(justmet_prop.glm)
#testZeroInflation(justmet_prop.glm)

#justmet_prop.glm_simResid <- simulateResiduals(fittedModel = justmet_prop.glm)
#plot(justmet_prop.glm_simResid)

#justmet_spores <- subset(justmetsch, Metsch.spore.in.animal > 0)

#justmet_spores.glm <- glm(Metsch.spore.in.animal ~ DayMetFactor, data = justmet_spores)
#summary(justmet_spores.glm)

#plot(justmet_spores.glm)

#testDispersion(justmet_spores.glm)
#testZeroInflation(justmet_spores.glm)

#justmet_spores.glm_simResid <- simulateResiduals(fittedModel = justmet_spores.glm)
#plot(justmet_spores.glm_simResid)

#modelmix <- lmer(transformedbody ~ time * days_since_past_exposure * days_since_metsch_exposure * Terminal.infection * Coexposed + (1 | uniqueident), data = long.data) 
#modelmix2_drop1 <- lmer(transformedbody ~ time * days_since_past_exposure * days_since_metsch_exposure + Terminal.infection + (1 | uniqueident), data = long.data)
#anova(modelmix, modelmix2_drop1)

#This one is the best one, but I'm worried about how to deal with controls, since time has been taken out? 
#modelmix3_drop2 <- lmer(transformedbody ~ days_since_past_exposure * days_since_metsch_exposure + Terminal.infection + time + (1 | uniqueident), data = long.data) 
#summary(modelmix3_drop2)

#modelmix3_notime <- lmer(transformedbody ~ days_since_past_exposure * days_since_metsch_exposure + Terminal.infection + (1 | uniqueident), data = long.data) 
#anova(modelmix3_drop2, modelmix3_notime) #time definitely makes things worse

#modelmix3_nointeraction <- lmer(transformedbody ~ days_since_past_exposure + days_since_metsch_exposure + Terminal.infection + (1 | uniqueident), data = long.data) 
#anova(modelmix3_notime, modelmix3_nointeraction)

#summary(modelmix3_nointeraction)
#plot(modelmix3_nointeraction)

#modelmix3_nointeraction_simResid <- simulateResiduals(fittedModel = modelmix3_nointeraction)
#plot(modelmix3_nointeraction_simResid)

#modelmix4_interaction <- lmer(transformedbody ~ days_diff * days_since_metsch_exposure + Terminal.infection + (1 | uniqueident), data = long.data) 
#anova(modelmix4_interaction, modelmix3_nointeraction)

#summary(modelmix4_interaction)
#plot(modelmix4_interaction)
#modelmix4_interaction_simResid <- simulateResiduals(fittedModel = modelmix4_interaction)
#plot(modelmix4_interaction_simResid)

#modelmix6 <- lmer(transformedbody ~ days_diff + days_since_metsch_exposure + Terminal.infection + (1 | uniqueident), data = long.data) 

modelmix7 <- lmer(transformedbody ~ time * DayMet * Terminal.infection + (1 | uniqueident), data = long.data) 

summary(modelmix7)
plot(modelmix7)

modelmix7_simResid <- simulateResiduals(fittedModel = modelmix7)
plot(modelmix7_simResid)

testDispersion(modelmix7)
testZeroInflation(modelmix7)

#Need to keep terminal infection even though nothing is significant 
#modelmix5 <- lmer(transformedbody ~ days_diff * days_since_metsch_exposure + (1 | uniqueident), data = long.data) 
#anova(modelmix4_interaction, modelmix5)

#Plotting models to make a table 
tab1 <- as.data.frame(past_spore_emm)
tab3 <- as.data.frame(emm_metsch_spores)

tab1$model <- "P. ramosa spore model"
tab3$model <- "A. monospora spore model"

all_tabs <- rbind(tab1, tab3)

pub_tab <- pub_tab %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  )

spore_yield_table<-gt(pub_tab) %>%
  tab_header(
    title = "Estimated Marginal Means Across Spore Yield Models"
  )
gtsave(spore_yield_table,"emmeans_spore_yield_table.docx")


