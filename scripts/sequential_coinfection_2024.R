#Updated to work with GitHub on 1/5/26
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

#set whatever working directory you want 
setwd("C:/Users/kiram/OneDrive/Desktop/Duffy Lab/Sequential conifection past vs metsch/Data")
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
f1 <- ggplot(prevalence_combined, aes(x = DayMetFactor, y = prevalence, color = Coinfected)) +
  geom_point(size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "single infection" = "cornflowerblue",
      "coinfection" = "#7f39d4",
      "total" = "black"
    )
  ) +
  ylab("Infection Prevalence of *P. ramosa*") +
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
  coord_cartesian(ylim = c(0, 0.8))

# Spore count summary
past_spores <- subset(past, PastInfected == "Yes")

f2 <- ggplot(past_spores, aes(x = Coinfected, y = Past.spore.in.animal, color = DayMetFactor)) +
  geom_violin(color = "black", fill = NA, trim = FALSE) +
  geom_jitter(width = 0.3, size = 3, alpha = 1, show.legend = TRUE) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of *P. ramosa* spores in animal",
    color = NULL
  ) +
  theme(
    axis.title.x = ggtext::element_markdown(), 
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(), 
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  scale_color_manual(
    values = c("never" = "grey50", "5" = "#c1b0d6", "10" = "#ac8cd4", "15" = "#7f39d4", "30" = "#640ecc")
  ) +
  scale_x_discrete(labels = c("No" = "single infection", "Yes" = "coinfection")) +
  scale_y_continuous(limits = c(0, NA), labels = scales::comma)

#Same thing as the two plots above, but for metsch 
#finding prevalence of metsch infection 
metsch <- subset(spores.factor, TreatmentGroup %in% c("T2", "T3", "T4", "T5", "T6"))
metsch$DayMetFactor <- as.character(metsch$DayMetFactor)
metsch$DayMetFactor[metsch$DayMetFactor == "0"] <- "never"
metsch$DayMetFactor <- as.factor(metsch$DayMetFactor)
metsch$DayMetFactor <- factor(metsch$DayMetFactor,levels = c("never", "5", "10", "15", "30"))

#third part of the facet plot
#finding prevalence of past infection 
metsch$MetschInfected_numeric <- ifelse(metsch$MetschInfected == "Yes", 1, 0)
metsch <- na.omit(metsch)

# Step 1: Create the base prevalence data
prevalence_by_coinfection_metsch <- metsch %>%
  filter(MetschInfected_numeric == 1) %>%
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
prevalence_metsch_sum <- prevalence_by_coinfection_metsch %>%
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
prevalence_combined_metsch <- bind_rows(prevalence_by_coinfection_metsch, prevalence_metsch_sum)

# Optional: make Coinfected labels more readable
prevalence_combined_metsch$Coinfected <- factor(
  prevalence_combined_metsch$Coinfected,
  levels = c("No", "Yes", "Sum"),
  labels = c("single infection", "coinfection", "total"))

#Making this facet plot by stacked with coinfected vs singly infected 
f3 <- ggplot(prevalence_combined_metsch, aes(x = DayMetFactor, y = prevalence, color = Coinfected)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), shape = 17) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), 
                width = 0.1, position = position_dodge(width = 0.4)) +
  scale_color_manual(
    values = c(
      "single infection" = "tomato",
      "coinfection" = "#7f39d4",
      "total" = "black"
    )
  ) +
  ylab("Infection Prevalence of *A. monospora*") +
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
  coord_cartesian(ylim = c(0, 0.8))

#fourth part of the facet plot
metsch_spores <- subset(metsch, MetschInfected == "Yes")

f4 <- ggplot(metsch_spores, aes(x = Coinfected, y = Metsch.spore.in.animal, color = DayMetFactor)) +
  geom_violin(color = "black", fill = NA, trim = FALSE) +
  geom_jitter(width = 0.3, size = 3, alpha = 1, show.legend = TRUE, shape = 17) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of *A. monospora* spores in animal",
    color = NULL
  ) +
  theme(
    axis.title.x = ggtext::element_markdown(), 
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(), 
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  scale_color_manual(
    values = c("never" = "grey50", "5" = "#c1b0d6", "10" = "#ac8cd4", "15" = "#7f39d4", "30" = "#640ecc")
  ) +
  scale_x_discrete(labels = c("No" = "single infection", "Yes" = "coinfection")) +
  scale_y_continuous(limits = c(0, NA), labels = scales::comma)

infectionplot <- (f1 | f2) / (f3 | f4)
infectionplot
ggsave("infectionplot.png", infectionplot, dpi = 600, width = 10, height = 9, units = "in")

####Veering off into a different direction 
#Testing if the liklihood of an individual to be coinfected changed depending on when metsch was added in 
onlycoinfectiontreatments <- subset(spores.factor, TreatmentGroup %in% c("T3", "T4", "T5", "T6"))
levels(onlycoinfectiontreatments$Coinfected)[levels(onlycoinfectiontreatments$Coinfected) == "Yes"] <- "1"
levels(onlycoinfectiontreatments$Coinfected)[levels(onlycoinfectiontreatments$Coinfected) == "No"] <- "0"

hist(onlycoinfectiontreatments$Metsch.spore.in.animal)
hist(onlycoinfectiontreatments$Past.spore.in.animal)

onlycoinfectiontreatments.glm <- glm(Coinfected ~ DayMetFactor, family = binomial, data = onlycoinfectiontreatments)
summary(onlycoinfectiontreatments.glm)

#All checks taken from Michelle's data. Looks good. 
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(onlycoinfectiontreatments.glm)
 
testDispersion(onlycoinfectiontreatments.glm)
testZeroInflation(onlycoinfectiontreatments.glm)

onlycoinfectiontreatments.glm_simResid <- simulateResiduals(fittedModel = onlycoinfectiontreatments.glm)
plot(onlycoinfectiontreatments.glm_simResid)

#plot it out
onlycoinfectiontreatments$Coinfected <- as.numeric(onlycoinfectiontreatments$Coinfected)
dotchart(onlycoinfectiontreatments$Coinfected,
         groups = factor(onlycoinfectiontreatments$DayMetFactor),
         main = "Cleveland dotplot")

chai <- anova(onlycoinfectiontreatments.glm, test = "Chisq")
chai

#Do it again, but for metsch 
coinfmetspores <- subset(metsch, Metsch.spore.in.animal > 0)
hist(coinfmetspores$Metsch.spore.in.animal)
hist(coinfmetspores$Past.spore.in.animal)

coinfmetspores.glm <- glm(Metsch.spore.in.animal ~  DayMetFactor*Coinfected, data = coinfmetspores)
summary(coinfmetspores.glm)

plot(coinfmetspores.glm)

testDispersion(coinfmetspores.glm)
testZeroInflation(coinfmetspores.glm)

coinfmetspores.glm_simResid <- simulateResiduals(fittedModel = coinfmetspores.glm)
plot(coinfmetspores.glm_simResid)

library(knitr)
library(kableExtra)

coinfmetspores.glm.df <- emmeans(coinfmetspores.glm, specs = pairwise ~ Coinfected | DayMetFactor, type = "response")
coinfmetspores.glm.df
coinfmetspores.emm.df <- as.data.frame(coinfmetspores.glm.df$emmeans)

#Extract contrasts (which are the Tukey comparisons between groups for p-values)
coinf_emmeans <- coinfmetspores.glm.df$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

coinf_emmeans <- coinf_emmeans |> 
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

coinf_emmeans <- coinf_emmeans %>%
  rename(
    Contrast = contrast,
    `Day A. monospora was added` = DayMetFactor,
    `Estimate` = estimate,
    `Standard Error` = SE,
    `Lower CL` = lower.CL,
    `Upper CL` = upper.CL,
    `t-value` = t.ratio,
    `p-value` = p.value)

kable(coinf_emmeans, caption = "Pairwise Comparisons of A. monospora spore yields") %>%
  kable_styling()

library(flextable)

# Format ANOVA table
ft_anova <- flextable(chai) |>
  autofit() |>
  set_caption("Analysis of Deviance (Likelihood Ratio Test)")

# Format emmeans table
ft_emm <- flextable(coinf_emmeans) |>
  autofit() |>
  set_caption("Estimated Marginal Means (on response scale)")

# View in RStudio
ft_anova
ft_emm

read_docx() |>
  body_add_flextable(ft_model) |>
  body_add_par("Estimated Marginal Means") |>
  body_add_flextable(ft_emm) |>
  print(target = "model_and_emmeans.docx")

  
# Get the formula as text
model_formula <- as.character(formula(onlycoinfectiontreatments.glm))
model_formula_text <- paste(model_formula, collapse = " ")

# Add to table caption
library(flextable)
ft_anova <- flextable(chai) |>
  autofit() |>
  set_caption(paste0("Analysis of Deviance for Model: ", model_formula_text))

#Make a plot showing how co-infection makeup changes through time
#coinfmelt <- subset(onlycoinfectionspores, select = c(Past.spore.in.animal, Metsch.spore.in.animal, DayMetFactor))
#coinfmelt$uniqident <- sample(1:31, replace=FALSE)
#coinfmelt_melt <- melt(data = coinfmelt, 
#                  id.vars = c("uniqident", "DayMetFactor"), 
#    variable.name = "sporetype", value.name = "sporeload")

#ggplot(coinfmelt_melt, aes(x=DayMetFactor, y=sporeload, group=sporetype, color=sporetype, fill=sporetype)) + geom_bar(position = "dodge", stat = "summary") +
#  theme_classic() +  stat_summary(fun = "mean", geom = "col") +
#  stat_summary(fun.data = "mean_se")
  
#Looking at total spore number in  past
coinfandcontrol <- subset(spores.factor, TreatmentGroup %in% c("T3", "T4", "T5", "T6"))
onlycoinfectionspores <- subset(coinfandcontrol, Past.spore.in.animal > 0)
hist(onlycoinfectionspores$Metsch.spore.in.animal)
hist(onlycoinfectionspores$Past.spore.in.animal)

onlycoinfectionspores.glm <- glm(Past.spore.in.animal ~ DayMetFactor*Coinfected, data = onlycoinfectionspores)
summary(onlycoinfectionspores.glm)

plot(onlycoinfectionspores.glm)

testDispersion(onlycoinfectionspores.glm)
testZeroInflation(onlycoinfectionspores.glm)

onlycoinfectionspores.glm_simResid <- simulateResiduals(fittedModel = onlycoinfectionspores.glm)
plot(onlycoinfectionspores.glm_simResid)

onlycoinfectionspores.glm.df <- emmeans(onlycoinfectionspores.glm, specs = pairwise ~ Coinfected | DayMetFactor, type = "response")
onlycoinfectionspores.glm.df
past.prod.m1.emm.df <- as.data.frame(onlycoinfectionspores.glm.df$emmeans)

#Extract contrasts (which are the Tukey comparisons between groups for p-values)
past_emmeans <- onlycoinfectionspores.glm.df$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()

past_emmeans <- past_emmeans |> 
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

past_emmeans <- past_emmeans %>%
  rename(
    Contrast = contrast,
    `Day A. monospora was added` = DayMetFactor,
    `Estimate` = estimate,
    `Standard Error` = SE,
    `Lower CL` = lower.CL,
    `Upper CL` = upper.CL,
    `t-value` = t.ratio,
    `p-value` = p.value)

kable(past_emmeans, caption = "Pairwise Comparisons of P. ramosa spore yields") %>%
  kable_styling()


#Make a plot showing how co-infection makeup changes through time
#coinfmelt <- subset(onlycoinfectionspores, select = c(Past.spore.in.animal, Metsch.spore.in.animal, DayMetFactor))
#coinfmelt$uniqident <- sample(1:16, replace=FALSE)
#coinfmelt_melt <- melt(data = coinfmelt, 
#                       id.vars = c("uniqident", "DayMetFactor"), 
#                       variable.name = "sporetype", value.name = "sporeload")

#ggplot(coinfmelt_melt, aes(x=DayMetFactor, y=sporeload, group=sporetype, color=sporetype, fill=sporetype)) + geom_bar(position = "dodge", stat = "summary") +
#  theme_classic() +  stat_summary(fun = "mean", geom = "col") +
#  stat_summary(fun.data = "mean_se")

#Messing with proportion of metsch infection 
metschprop <- subset(spores.factor, TreatmentGroup %in% c("T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"))

#focusing on past exposed but not infected 
metschprop_pastexposed <- subset(metschprop, PastInfected == "No")
levels(metschprop_pastexposed$MetschInfected)[levels(metschprop_pastexposed$MetschInfected) == "Yes"] <- "1"
levels(metschprop_pastexposed$MetschInfected)[levels(metschprop_pastexposed$MetschInfected) == "No"] <- "0"

#Metch infection prevalence, there are some influential outliers but they are all legit so they are kept in 
metschprop_pastexposed.glm2 <- glm(MetschInfected ~ PastExposed*DayMetFactor, family = binomial, data = metschprop_pastexposed)
summary(metschprop_pastexposed.glm2)
plot(metschprop_pastexposed.glm2)

testDispersion(metschprop_pastexposed.glm2)
testZeroInflation(metschprop_pastexposed.glm2)

metschprop_pastexposed.glm2_simResid <- simulateResiduals(fittedModel = metschprop_pastexposed.glm2)
plot(metschprop_pastexposed.glm2_simResid)

metschprop_pastexposed.glm2.df <- emmeans(metschprop_pastexposed.glm2, specs = pairwise ~ PastExposed | DayMetFactor, type = "response")
metschprop_pastexposed.glm2.df

#look at spores 
hist(metschprop_pastexposed$Metsch.spore.in.animal)

#looks good 
metschspores_pastexposed <- subset(metschprop_pastexposed, Metsch.spore.in.animal > 0)
metschspores_pastexposed.glm <- glm(Metsch.spore.in.animal ~ PastExposed+DayMetFactor, data = metschspores_pastexposed)
summary(metschspores_pastexposed.glm)
plot(metschspores_pastexposed.glm)

testDispersion(metschspores_pastexposed.glm)
testZeroInflation(metschspores_pastexposed.glm)

metschspores_pastexposed.glm_simResid <- simulateResiduals(fittedModel = metschspores_pastexposed.glm)
plot(metschspores_pastexposed.glm_simResid)

metschspores_pastexposed.glm.df <- emmeans(metschspores_pastexposed.glm, specs = pairwise ~ PastExposed | DayMetFactor, type = "response")
metschspores_pastexposed.glm.df

#Moving into proportion infected with metsch when animals are either past infected or not but exposed are taken out 
metschprop_pastinfected <- subset(metschprop, PastInfected == "Yes" | PastExposed == "No")
levels(metschprop_pastinfected$MetschInfected)[levels(metschprop_pastinfected$MetschInfected) == "Yes"] <- "1"
levels(metschprop_pastinfected$MetschInfected)[levels(metschprop_pastinfected$MetschInfected) == "No"] <- "0"

metschprop_pastinfected.glm <- glm(MetschInfected ~ PastInfected+DayMetFactor, family = binomial, data = metschprop_pastinfected)
summary(metschprop_pastinfected.glm)

testDispersion(metschprop_pastinfected.glm)
testZeroInflation(metschprop_pastinfected.glm)

metschprop_pastinfected.glm_simResid <- simulateResiduals(fittedModel = metschprop_pastinfected.glm)
plot(metschprop_pastinfected.glm_simResid)

#Moving on to spores
metschspores_pastinfected <- subset(metschprop_pastinfected, Metsch.spore.in.animal > 0)
metschspores_pastinfected.glm <- glm(Metsch.spore.in.animal ~ DayMetFactor + Past.spore.in.animal + PastInfected, data = metschspores_pastinfected)
summary(metschspores_pastinfected.glm)
plot(metschspores_pastinfected.glm)

testDispersion(metschspores_pastinfected.glm)
testZeroInflation(metschspores_pastinfected.glm)

metschspores_pastinfected.glm_simResid <- simulateResiduals(fittedModel = metschspores_pastinfected.glm)
plot(metschspores_pastinfected.glm_simResid)

#Now on to past infection 
pastprop <- subset(spores.factor, TreatmentGroup %in% c("T2" ,"T3", "T4", "T5", "T6"))

#focusing on metsch exposed but not infected 
pastprop_metexposed <- subset(pastprop, MetschInfected == "No")
levels(pastprop_metexposed$PastInfected)[levels(pastprop_metexposed$PastInfected) == "Yes"] <- "1"
levels(pastprop_metexposed$PastInfected)[levels(pastprop_metexposed$PastInfected) == "No"] <- "0"

pastprop_metexposed.glm <- glm(PastInfected ~ MetschExposed+DayMetFactor, family = binomial, data = pastprop_metexposed)
summary(pastprop_metexposed.glm)

plot(pastprop_metexposed.glm)

testDispersion(pastprop_metexposed.glm)
testZeroInflation(pastprop_metexposed.glm)

pastprop_metexposed.glm_simResid <- simulateResiduals(fittedModel = pastprop_metexposed.glm)
plot(pastprop_metexposed.glm_simResid)

### spores
pastspores_metexposed <- subset(pastprop_metexposed, Past.spore.in.animal > 0)
pastspores_metexposed.glm <- glm(Past.spore.in.animal ~ MetschExposed + DayMetFactor, data = pastspores_metexposed)
summary(pastspores_metexposed.glm)

plot(pastspores_metexposed.glm)

testDispersion(pastprop_metexposed.glm)
testZeroInflation(pastprop_metexposed.glm)

pastspores_metexposed.glm_simResid <- simulateResiduals(fittedModel = pastspores_metexposed.glm)
plot(pastspores_metexposed.glm_simResid)

#What's is up with the 30 day metsch? 
justmetsch <- subset(spores.factor, TreatmentGroup %in% c("T7" ,"T8", "T9", "T10"))
levels(justmetsch$MetschInfected)[levels(justmetsch$MetschInfected) == "Yes"] <- "1"
levels(justmetsch$MetschInfected)[levels(justmetsch$MetschInfected) == "No"] <- "0"

justmet_prop.glm <- glm(MetschInfected ~ DayMetFactor, family = binomial, data = justmetsch)
summary(justmet_prop.glm)

plot(justmet_prop.glm)

testDispersion(justmet_prop.glm)
testZeroInflation(justmet_prop.glm)

justmet_prop.glm_simResid <- simulateResiduals(fittedModel = justmet_prop.glm)
plot(justmet_prop.glm_simResid)

justmet_spores <- subset(justmetsch, Metsch.spore.in.animal > 0)

justmet_spores.glm <- glm(Metsch.spore.in.animal ~ DayMetFactor, data = justmet_spores)
summary(justmet_spores.glm)

plot(justmet_spores.glm)

testDispersion(justmet_spores.glm)
testZeroInflation(justmet_spores.glm)

justmet_spores.glm_simResid <- simulateResiduals(fittedModel = justmet_spores.glm)
plot(justmet_spores.glm_simResid)

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

el.final.nouninf <- el.final %>%
  filter(InfStatusFactor != "Uninf")


totalegg.glm <- glmmTMB(LifeRep ~ DayMetFactor * InfStatusFactor * lifesurvival, family = nbinom1(), data = el.final.nouninf)
summary(totalegg.glm)

testDispersion(totalegg.glm)
testZeroInflation(totalegg.glm)

totalegg.glm_simResid <- simulateResiduals(fittedModel = totalegg.glm)
plot(totalegg.glm_simResid)

#I don't really understand emmeans but here goes nothing 
emm_egg <- emmeans(totalegg.glm, ~ DayMetFactor * InfStatusFactor * lifesurvival) 
summary(emm_egg)
totalegg.glm.df <- emmeans(totalegg.glm, specs = pairwise ~ InfStatusFactor | DayMetFactor | lifesurvival, type = "response")
totalegg.glm.df
coinfmetspores.emm.df <- as.data.frame(coinfmetspores.glm.df$emmeans)

plot(emm_egg)

#Extract the means for Egg Production Based on Infection Status and Day of Metsch Exposure Numeric
egg.means.3to6.DMnum <- el.3to6.final %>%
  group_by(DayMetNum, InfStatusFactor) %>%
  summarize(N  = length(LifeRep),
            LifeRep.mean = mean(LifeRep, na.rm = T),
            LifeRep.sd = sd(LifeRep, na.rm = T),
            LifeRep.se = LifeRep.sd / sqrt(N))

##Survival 
#Load necessary libraries
library(survival)
library(survminer)
library(ggsurvfit)
library(survivalAnalysis)

#Survivor analysis
infected_el.final <- subset(el.final, !PastExpStatus == "ExposedUninfected")
infected_el.final <- subset(infected_el.final, !MetExpStatus == "ExposedUninfected")
infected_el.final$survobject <- with(infected_el.final, Surv(lifesurvival, censor))
seq_lifespan <- survfit(survobject ~ InfStatusFactor, data = infected_el.final, conf.type = "log-log")
summary(seq_lifespan)

el.final$survobject <- with(el.final, Surv(lifesurvival, censor))
seq_lifespan <- survfit(survobject ~ InfStatusFactor, data = el.final, conf.type = "log-log")
summary(seq_lifespan)

cox.seqlifespan <- coxph(survobject ~ InfStatusFactor*DayMetFactor, data = el.final)
summary(cox.seqlifespan)
cox.test <- cox.zph(cox.seqlifespan)
cox.test
plot(cox.test)

cox.seqlifespan <- coxph(survobject ~ InfStatusFactor*DayMetFactor, data = el.final)
summary(cox.seqlifespan)
cox.test <- cox.zph(cox.seqlifespan)
cox.test
plot(cox.test)

#Plot it out 
ggforest(cox.seqlifespan, data = el.final)
forest_model

library(forestmodel)

forest_model(cox.seqlifespan)

library(broom)

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

#Making my forest plot 
library(forcats)

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

ggsave("el.final.survival.tiff", dpi = 600, width = 10.5, height = 6, units = "in", compression="lzw")
ggsave("el.final.survival.png", dpi = 600, width = 10.5, height = 6, units = "in")


#subset for just coinfected? 
#el.final_justcoinf$survobject <- with(el.final_justcoinf, Surv(lifesurvival, censor))
#seq_lifespan <- survfit(survobject ~ DayMetFactor, data = el.final_justcoinf, conf.type = "log-log")
#summary(seq_lifespan)

#cox.justcoinf <- coxph(survobject ~ DayMetFactor, data =el.final_justcoinf)
#summary(cox.justcoinf)

#Plot it out 
#ggforest(cox.justcoinf, data = el.final_justcoinf)

#ggsurvplot(seq_lifespan, data = el.final_justcoinf, xlab = "lifespan", facet.by = "DayMetFactor",
                                   conf.type = "log-log", conf.int = T, pval = TRUE,
                                   pval.method = TRUE)

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
##### BODY SIZE STARTS HERE ######

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
library(glmmTMB)
library(performance)
library(dplyr)
library(parsedate)
library(ggplot2)
library(FSA)
library(nlstools)
library(plotrix)
library(Matrix)
library(tidyr)
library(lmerTest)

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

#Model 6 is the model going forward
long.data_nocontrols <- subset(long.data, TreatmentGroup!= "T1")
modelmix6 <- lmer(transformedbody ~ days_diff * days_since_metsch_exposure + Terminal.infection + (1 | uniqueident), data = long.data_nocontrols) 

summary(modelmix6)
plot(modelmix6)

modelmix6_simResid <- simulateResiduals(fittedModel = modelmix6)
plot(modelmix6_simResid)

testDispersion(modelmix6)
testZeroInflation(modelmix6)

#Need to keep terminal infection even though nothing is significant 
#modelmix5 <- lmer(transformedbody ~ days_diff * days_since_metsch_exposure + (1 | uniqueident), data = long.data) 
#anova(modelmix4_interaction, modelmix5)

#model validation - plotting normality of residuals 
e <- resid(modelmix6)
hist(e, xlab = "Residuals", main = "")

plot(long.data_nocontrols$days_diff, e, xlab = "days difference",
     ylab = "Residuals")
plot(long.data_nocontrols$days_since_metsch_exposure, e, xlab = "days since metsch exposure",
     ylab = "Residuals")

#Ask someone who knows whats up if this is fine? 

#Make a graph, I believe in you
p1 <- ggplot(long.data_nocontrols, aes(x=time, y=transformedbody, color=Terminal.infection)) +
  geom_point(alpha = .4) + geom_smooth(se = FALSE) + theme_minimal() + xlab("time in days") + ylab("") +
  scale_color_manual(name = "infection status", labels = c("uninfected", "coinfected", expression(italic("A. monospora")), expression(italic("P. ramosa"))), values = c("none" = "grey50", "past" = "cornflowerblue", "metsch" = "red1", "coinfected" = "mediumpurple")) + theme(legend.position = c(0.15, 0.8))

##PLACE THIS IS SUPPLEMENTAL 
#long.data_nocontrols_noinfections <- subset(long.data_nocontrols,Terminal.infection=="none")

#p2 <- ggplot(long.data_nocontrols, aes(x=time, y=transformedbody, color=Coexposed)) + geom_point(alpha = .4) + geom_smooth(se = FALSE) + theme_minimal() + 
 # xlab("") + ylab("transformed body size") + theme(legend.position = c(0.1, 0.9)) +
 # scale_color_manual(name = "coexposed?", labels = c("no", "yes"), values = c("No" = "grey50", "Yes" = "red1"))

#graph about does metsch through time change size in coinfected individuals 
coinfected_bodysize <- subset(long.data, Terminal.infection=="coinfected" | TreatmentGroup=="T1")

p3 <- ggplot(coinfected_bodysize, aes(x=time, y=transformedbody, color=DayMet)) + geom_point(alpha = .4) + geom_smooth(se = FALSE) + theme_minimal() + 
  xlab("time in days") + ylab("transformed body size") + theme(legend.position = c(0.35, 0.8)) +
  scale_color_manual(
    name = expression(
      "number of days after " * italic("P. ramosa") * " exposure\nthat " * italic("M. bicuspidata") * " was added in coinfected individuals"), values = c("None" = "grey50", "5" = "mediumpurple1", "10" = "mediumpurple2", "15" = "mediumpurple3", "30" = "mediumpurple4"))

p2_blank <- patchwork::plot_spacer()
spacer <- ggplot() + theme_void() + theme(plot.margin = margin(t = 10))

layout <- (p2 | p1) / spacer / (p3 | p2_blank)

three_panel <- layout + plot_layout(heights = c(1, 0.1, 0.9))

three_panel

ggsave("three_panel_figure.png", three_panel, width = 10, height = 6, dpi = 300)


