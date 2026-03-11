#Load libraries
library(tidyverse)
library(patchwork)
library(ggtext)
library(here)
library(reshape2)
library(survival)
library(performance)
library(DHARMa)
library(glmmTMB)
library(emmeans)
library(car)
# library(survminer)

# Overdispersion check taken from earlier code
# Will get called later, loading it now
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Tell R where files are stored
here::i_am("scripts/sequential_coinfection.R")

# Load spore data
total.spores <- read.csv("data/Coinfection Total Spore Yield Data.csv")
# View(total.spores)

#Replace the #DIV/0! with NA for the cells that contain that value
total.spores[total.spores == '#DIV/0!'] <- NA

#Split apart treatment number and rep number into separate columns.
total.spores[c('TreatmentGroup', 'Rep')] <- str_split_fixed(total.spores$sample.ID, '-', 2)
# View(total.spores)

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

#Remove animals that had only 1 or 2 total spores across four fields of counts; don't feel confident about these diagnoses
#This removes four animals based on Metsch and four based on Past
spores.numeric <- spores.numeric %>%
  subset(!(Metsch.spore.in.animal %in% c(250, 500)) &
           !(Past.spore.in.animal %in% c(250, 500)))

# Remove the animals that died on or before day 10 (and therefore couldn't be diagnosed accurately)
spores.numeric <- spores.numeric %>%
  filter(!(TreatmentGroup == "T2" & Rep == 3),
         !(TreatmentGroup == "T2" & Rep == 25),
         !(TreatmentGroup == "T4" & Rep == 5),
         !(TreatmentGroup == "T8" & Rep == 10),
         !(TreatmentGroup == "T9" & Rep == 25))

#### Note: figure 1 is an overview of the experimental design, was not created in R

#### Creating figure 2 #####
#### Infection prevalence & spore yield #####

#facet plot of infection prevalence and spores 
past <- subset(spores.numeric, TreatmentGroup %in% c("T3", "T4", "T5", "T6"))

#first part of the facet plot
#finding prevalence of past infection 
past$PastInfected_numeric <- ifelse(past$PastInfected == "Yes", 1, 0)
past <- subset(past, !DayMetNum == "0")

# Step 1: Create the base prevalence data
prevalence_by_coinfection <- past %>%
  filter(PastInfected_numeric == 1) %>%
  group_by(DayMetNum, Coinfected) %>%
  summarise(
    infected_n = n(),
    .groups = "drop"
  ) %>%
  left_join(
    past %>%
      group_by(DayMetNum) %>%
      summarise(total_n = n(), .groups = "drop"),
    by = "DayMetNum"
  ) %>%
  mutate(
    prevalence = infected_n / total_n,
    se = sqrt(prevalence * (1 - prevalence) / total_n)
  )

# Step 2: Add a "sum" row per DayMetFactor
prevalence_sum <- prevalence_by_coinfection %>%
  group_by(DayMetNum) %>%
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

# make Coinfected labels more readable
prevalence_combined$Coinfected <- factor(
  prevalence_combined$Coinfected,
  levels = c("No", "Yes", "Sum"),
  labels = c("single infection", "coinfection", "total"))

#Making this facet plot by stacked with coinfected vs singly infected 
f1<-ggplot(prevalence_combined, aes(x = DayMetNum, y = prevalence, color = Coinfected, shape = Coinfected)) +
  geom_point(size = 4, position = position_dodge(width = 1.75)) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), 
                width = 0.1, position = position_dodge(width = 1.75)) +
  scale_color_manual(
    values = c(
      "single infection" = "cornflowerblue",
      "coinfection" = "#7f39d4",
      "total" = "black"
    ), labels = c(
      "single infection" = "single infection",
      "coinfection" = "coinfection",
      "total" = "total *P. ramosa* infection"
      )
    ) +
  scale_shape_manual(                                       
    values = c(
      "single infection" = 16,                             
      "coinfection" = 17,                                   
      "total" = 18                                          
    ),
    labels = c(                                             
      "single infection" = "single infection",
      "coinfection" = "coinfection",
      "total" = "total *P. ramosa* infection"
    )
  ) +
  ylab("*P. ramosa* infection prevalence") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) + 
  labs(color = NULL, shape = NULL) +                        # Added shape = NULL to remove both legend titles
  coord_cartesian(ylim = c(0, 1)) +
  guides(
    color = guide_legend(                                   # Added guides() to merge legends
      title = NULL,
      title.theme = ggtext::element_markdown()
    ),
    shape = guide_legend(
      title = NULL,
      title.theme = ggtext::element_markdown()
    )
  )

# Spore count summary
past_spores <- subset(past, PastInfected == "Yes")

#Creating a dataset that averages spores by treatment 
past_spores_summary <- past_spores %>%
  group_by(DayMetNum, Coinfected) %>%
  summarise(
    n = n(),
    mean_spores = mean(Past.spore.in.animal, na.rm = TRUE),
    se_spores = sd(Past.spore.in.animal, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# make Coinfected labels more readable
past_spores_summary$Coinfected <- factor(
  past_spores_summary$Coinfected,
  levels = c("No", "Yes"),
  labels = c("single infection", "coinfection"))

#Past spores plotted 
f2<-ggplot(past_spores_summary, aes(x = DayMetNum, y = mean_spores, color = Coinfected, shape = Coinfected)) +
  geom_point(size = 4, position = position_dodge(width = 1.75)) +
  geom_errorbar(aes(ymin = mean_spores - se_spores, ymax = mean_spores + se_spores), 
                width = 0.1, position = position_dodge(width = 1.75)) +
  scale_color_manual(
    values = c(
      "single infection" = "cornflowerblue",
      "coinfection" = "#7f39d4"), labels = c(
      "single infection" = "single infection",
      "coinfection" = "coinfection")) +
  scale_shape_manual(                                       
    values = c(
      "single infection" = 16,                             # Circle
      "coinfection" = 17                                   # Triangle
    ),
    labels = c(                                          
      "single infection" = "single infection",
      "coinfection" = "coinfection"
    )
  ) +
  ylab("*P. ramosa* spores per host") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) + labs(color = NULL) + scale_y_continuous(labels = scales::label_comma()) + ylim(0, 900000) +
  guides(
    color = guide_legend(                               
      title = NULL,
      title.theme = ggtext::element_markdown()
    ),
    shape = guide_legend(
      title = NULL,
      title.theme = ggtext::element_markdown()
    )
  )


#Same thing as the two plots above, but for metsch 
#finding prevalence of metsch infection 
metsch_with_singlyexposed <- subset(spores.numeric, TreatmentGroup %in% c("T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"))

#third part of the facet plot
#finding prevalence of metsch infection 
metsch_with_singlyexposed$MetschInfected_numeric <- ifelse(metsch_with_singlyexposed$MetschInfected == "Yes", 1, 0)
#metsch_with_singlyexposed <- na.omit(metsch_with_singlyexposed)

# Step 1: Create the base prevalence data
prevalence_by_coinfection_metsch_with_singlyexposed <- metsch_with_singlyexposed %>%
  filter(MetschInfected_numeric == 1) %>%
  group_by(DayMetNum, Coinfected, PastExpStatus) %>%
  summarise(
    infected_n = n(),
    .groups = "drop"
  ) %>%
  left_join(
    past %>%
      group_by(DayMetNum) %>%
      summarise(total_n = n(), .groups = "drop"),
    by = "DayMetNum"
  ) %>%
  mutate(
    prevalence = infected_n / total_n,
    se = sqrt(prevalence * (1 - prevalence) / total_n)
  )

# Step 2: Add a "sum" row per DayMetFactor
prevalence_metsch_sum_with_singly_exposed <-
  prevalence_by_coinfection_metsch_with_singlyexposed %>%
  filter(PastExpStatus %in% c("ExposedUninfected", "Infected")) %>%
  group_by(DayMetNum) %>%
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

# make Coinfected labels more readable
prevalence_combined_metsch_with_singly_exposed$Coinfected <- factor(
  prevalence_combined_metsch_with_singly_exposed$Coinfected,
  levels = c("No", "Yes", "Sum"),
  labels = c("single infection", "coinfection", "total"))

prevalence_combined_metsch_with_singly_exposed$PastExpStatus <- factor(
  prevalence_combined_metsch_with_singly_exposed$PastExpStatus,
  levels = c("ExposedUninfected", "Unexposed", "Infected", "Summed"),
  labels = c("exposed but uninfected", "unexposed", "infected", "overall (only coexposed treatments)"))

fig3_data <- prevalence_combined_metsch_with_singly_exposed 

# %>%
#  filter(!(Coinfected == "single infection" & PastExpStatus == "exposed but uninfected"))

fig3_data <- fig3_data %>%
  mutate(group = paste(Coinfected, PastExpStatus, sep = " | "))

write.csv(fig3_data, file = "fig3_data.csv", row.names = FALSE)

#Third plot for the facet plot
f3<-ggplot(fig3_data, aes(x = DayMetNum, y = prevalence, color = group, shape = group)) +
  geom_point(size = 4, position = position_dodge(width = 1.75)) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), 
                width = 0.1, position = position_dodge(width = 1.75)) +
  scale_color_manual(
    values = c(
      "single infection | exposed but uninfected" = "tomato4",
      "single infection | unexposed" = "tomato",
      "coinfection | infected" = "#7f39d4",
      "total | overall (only coexposed treatments)" = "black"
    ),
    labels = c(
      "single infection | exposed but uninfected" = "single infection, exposed to *P. ramosa*",
      "single infection | unexposed" = "single infection, never exposed to *P. ramosa*",
      "coinfection | infected" = "coinfection",
      "total | overall (only coexposed treatments)" = "total *A. monospora* infection in coexposed treatments"
    )
  ) +
  scale_shape_manual(
    values = c(
      "single infection | exposed but uninfected" = 15,
      "single infection | unexposed" = 16,
      "coinfection | infected" = 17,
      "total | overall (only coexposed treatments)" = 18
    ),
    labels = c(
      "single infection | exposed but uninfected" = "single infection, exposed to *P. ramosa*",
      "single infection | unexposed" = "single infection, never exposed to *P. ramosa*",
      "coinfection | infected" = "coinfection",
      "total | overall (only coexposed treatments)" = "total *A. monospora* infection in coexposed treatments"
    )
  ) +
  ylab("*A. monospora* infection prevalence") +
  xlab("Day of *A. monospora* exposure") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  coord_cartesian(ylim = c(0, 0.99)) +
  guides(
    color = guide_legend(
      title = NULL,
      title.theme = ggtext::element_markdown()
    ),
    shape = guide_legend(
      title = NULL,
      title.theme = ggtext::element_markdown()
    )
  )

# Metsch spore yield
# Spore count summary
metsch_spores <- subset(metsch_with_singlyexposed, MetschInfected == "Yes")

#Creating a dataset that averages spores by treatment 
metsch_spores_summary <- metsch_spores %>%
  group_by(DayMetNum, Coinfected, PastExpStatus) %>%
  summarise(
    n = n(),
    mean_spores = mean(Metsch.spore.in.animal, na.rm = TRUE),
    se_spores = sd(Metsch.spore.in.animal, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# make Coinfected labels more readable
metsch_spores_summary$Coinfected <- factor(
  metsch_spores_summary$Coinfected,
  levels = c("No", "Yes"),
  labels = c("single infection", "coinfection"))

fig4_data <- metsch_spores_summary %>%
  filter(!(Coinfected == "single infection" & PastExpStatus == "exposed but uninfected"))

fig4_data <- fig4_data %>%
  mutate(group = paste(Coinfected, PastExpStatus, sep = " | "))

#Metsch spores plotted 
f4<-ggplot(fig4_data, aes(x = DayMetNum, y = mean_spores, color = group, shape = group)) +
  geom_point(size = 4, position = position_dodge(width = 1.75)) +
  geom_errorbar(aes(ymin = mean_spores - se_spores, ymax = mean_spores + se_spores), 
                width = 0.1, position = position_dodge(width = 1.75)) +
  scale_color_manual(
    values = c(
      "single infection | ExposedUninfected" = "tomato4",
      "coinfection | Infected" = "#7f39d4",
      "single infection | Unexposed" = "tomato"
    ),
    labels = c(
      "single infection | ExposedUninfected" = "single infection, exposed to *P. ramosa*",
      "coinfection | Infected" = "coinfection",
      "single infection | Unexposed" = "single infection, never exposed to *P. ramosa*"
    )
  ) +
  scale_shape_manual(
    values = c(
      "single infection | ExposedUninfected" = 15,
      "coinfection | Infected" = 17,
      "single infection | Unexposed" = 16
    ),
    labels = c(
      "single infection | ExposedUninfected" = "single infected, exposed to *P. ramosa*",
      "coinfection | Infected" = "coinfection",
      "single infection | Unexposed" = "single infected, never exposed to *P. ramosa*"
    )
  ) +
  ylab("*A. monospora* spores per host") +
  xlab("Day of *A. monospora* exposure") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, 1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      title.theme = ggtext::element_markdown(),
      override.aes = list(shape = c(17, 15, 16))
    ),
    shape = "none"
  ) + ylim(0, 250000)


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
ggsave("figures/infectionplot_withexposurelabels.png", infectionplot, dpi = 600, width = 15, height = 12, units = "in")


#### Analyses related to Figure 2 ####

# Analyses related to Past infection and spore yield
#Does the prevalence of any Past infection change with day of Metsch addition?
past2 <- subset(spores.numeric, TreatmentGroup %in% c("T3", "T4", "T5", "T6"))
past2$PastInfected_numeric <- ifelse(past2$PastInfected == "Yes", 1, 0)

past_prevalence.glm <- glm(PastInfected_numeric ~ DayMetNum, family = binomial, data = past2)
summary(past_prevalence.glm)

overdisp_fun(past_prevalence.glm)

past_prevalence.glm_simResid <- simulateResiduals(fittedModel = past_prevalence.glm)
plot(past_prevalence.glm_simResid)

past_prevalence_anova <- anova(past_prevalence.glm, test = "Chisq")
past_prevalence_anova # This is the result for whether day of Metsch exposure influences overall Past prev in T3-T6

#Does the prevalence of coinfection change with day of Metsch addition? 
onlycoinfectiontreatments <- subset(spores.numeric, TreatmentGroup %in% c("T3", "T4", "T5", "T6"))
levels(onlycoinfectiontreatments$Coinfected)[levels(onlycoinfectiontreatments$Coinfected) == "Yes"] <- "1"
levels(onlycoinfectiontreatments$Coinfected)[levels(onlycoinfectiontreatments$Coinfected) == "No"] <- "0"

hist(onlycoinfectiontreatments$Metsch.spore.in.animal)
hist(onlycoinfectiontreatments$Past.spore.in.animal)

onlycoinfectiontreatments.glm <- glm(Coinfected ~ DayMetNum, family = binomial, data = onlycoinfectiontreatments)
summary(onlycoinfectiontreatments.glm)

overdisp_fun(onlycoinfectiontreatments.glm)

onlycoinfectiontreatments.glm_simResid <- simulateResiduals(fittedModel = onlycoinfectiontreatments.glm)
plot(onlycoinfectiontreatments.glm_simResid)

coinf_anova <- anova(onlycoinfectiontreatments.glm, test = "Chisq")
coinf_anova # This is the result for whether day of Metsch exposure influences coinfection prev in T3-T6

#Now looking at Past spore yield
mean(past_spores$Past.spore.in.animal)
var(past_spores$Past.spore.in.animal) 

past_spore.glm <- glmmTMB(Past.spore.in.animal ~ DayMetNum*Coinfected, family = nbinom1, data = past_spores)
summary(past_spore.glm) # This shows a marginally significant interaction between day of Metsch exposure & coinfection on Past spore yield

overdisp_fun(past_spore.glm)

past_spore.glm_simResid <- simulateResiduals(fittedModel = past_spore.glm)
plot(past_spore.glm_simResid)


emm_past_spore <- emmeans(
  past_spore.glm,
  ~ Coinfected | DayMetNum,
  type = "response",  at = list(DayMetNum = c(5, 10, 15, 30)))

summary(emm_past_spore)
pairs(emm_past_spore) # Pairwise comparisons for each day
# Note that day 15 is driven by single coinfected animal -- not a robust test

#Does infection prevalence of metsch in coexposed treatments differ in comparison to control?
metsch_prevalence.glm <- glm(MetschInfected_numeric ~ DayMetNum*PastExposed, family = binomial, data = metsch_with_singlyexposed)
summary(metsch_prevalence.glm)

overdisp_fun(metsch_prevalence.glm)

metsch_prevalence.glm_simResid <- simulateResiduals(fittedModel = metsch_prevalence.glm)
plot(metsch_prevalence.glm_simResid)

metsch_prevalence_anova <- anova(metsch_prevalence.glm, test = "Chisq")
metsch_prevalence_anova # This is the main result for Metsch infection prevalence

#Adding test for past exposure 
emm_metsch_prevalence <- emmeans(
  metsch_prevalence.glm,
  ~ PastExposed * DayMetNum,  
  type = "response", at = list(DayMetNum = c(5, 10, 15, 30)))

summary(emm_metsch_prevalence)
pairs(emm_metsch_prevalence, by = "DayMetNum") 

#Metsch spores analysis

# First, just looking at individuals that were coexposed to see if coinfection affects spore yield
metsch_spores_T3toT6 <- metsch_spores %>%
  subset(TreatmentGroup %in% c("T3", "T4", "T5", "T6"))

metsch_spore.glm.1 <- glmmTMB(Metsch.spore.in.animal ~ DayMetNum*Coinfected, family = nbinom1, data = metsch_spores_T3toT6) 
summary(metsch_spore.glm.1)

overdisp_fun(metsch_spore.glm.1) # This analysis says things are fine, but keep reading.....

metsch_spore.glm_simResid <- simulateResiduals(fittedModel = metsch_spore.glm.1)
plot(metsch_spore.glm_simResid) # Residuals definitely look wonky

emm_metsch_spores.1 <- emmeans(
  metsch_spore.glm.1,
  ~ DayMetNum * Coinfected, 
  type = "response", at = list(DayMetNum = c(5, 10, 15, 30)))

summary(emm_metsch_spores.1)
pairs(emm_metsch_spores.1, by = "DayMetNum")

# write.csv(metsch_spores_T3toT6, file = "metsch_spores_T3toT6.csv", row.names = FALSE)

# Checking whether the problem is a singleton in T5 
metsch_spores_T3toT6 %>%
  group_by(TreatmentGroup, Coinfected) %>%
  summarize(count = n())
# This confirms that there was only one coinfected individual in T5, which is probably causing the problems. 

# First, just looking at individuals that were coexposed to see if coinfection affects spore yield
metsch_spores_T3toT6noT5 <- metsch_spores %>%
  subset((TreatmentGroup %in% c("T3", "T4", "T6")) & Metsch.spore.in.animal > 500)

metsch_spore.glm.2 <- glmmTMB(Metsch.spore.in.animal ~ DayMetNum*Coinfected, family = nbinom1, data = metsch_spores_T3toT6noT5) 
summary(metsch_spore.glm.2)

overdisp_fun(metsch_spore.glm.2)

metsch_spore.glm_simResid2 <- simulateResiduals(fittedModel = metsch_spore.glm.2)
plot(metsch_spore.glm_simResid2) # Definitely looks better

emm_metsch_spores.2 <- emmeans(
  metsch_spore.glm.2,
  ~ DayMetNum * Coinfected, 
  type = "response", at = list(DayMetNum = c(5, 10, 30)))

summary(emm_metsch_spores.2)
pairs(emm_metsch_spores.2, by = "DayMetNum")


# Also look at spore yield in singly infected individuals, comparing singly exposed with coexposed-but-singly-infected animals
metsch_spores_singleinf <- metsch_spores %>%
  subset(Coinfected %in% c("No"))

metsch_spore.glm.singleinf <- glmmTMB(Metsch.spore.in.animal ~ DayMetNum*PastExposed, family = nbinom1, data = metsch_spores_singleinf) 
summary(metsch_spore.glm.singleinf)

overdisp_fun(metsch_spore.glm.singleinf) 

metsch_spore.glm_simResid_singleinf <- simulateResiduals(fittedModel = metsch_spore.glm.singleinf)
plot(metsch_spore.glm_simResid_singleinf) # Also overdispersed

# Skipping comparison of single infection spores due to dispersion issues

#### FECUNDITY #####
egg <- read.csv("data/Coinfection Egg Data.csv")
# View(egg)

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

str(egg.life)

#Subset the data to only include animals that made it a certain number of days past exposure to Past and Metsch (e stands for egg, l stands for life in the dataframe name)
el.final <- subset(egg.life, lifesurvival > 10)
str(el.final)

#Summarizing who was removed -- that's a whole lot of T1 animals, suggesting something weird happened with that treatment
el.removed <- subset(egg.life, lifesurvival < 11)
# write.csv(el.removed, file = "el.removed.csv", row.names = FALSE)

# Remove the animals that had ambiguous infections (only 1 or 2 Metsch or Past spores across four fields)
el.final <- el.final %>%
  filter(!(TreatmentGroup == "T6" & Rep == 13),
         !(TreatmentGroup == "T10" & Rep == 21),
         !(TreatmentGroup == "T4" & Rep == 15),
         !(TreatmentGroup == "T6" & Rep == 12),
         !(TreatmentGroup == "T3" & Rep == 3),
         !(TreatmentGroup == "T4" & Rep == 10),
         !(TreatmentGroup == "T4" & Rep == 4),
         !(TreatmentGroup == "T2" & Rep == 7))

#Add a factor column that has infection status as either None, Past, Metsch, or Coinf. This will be useful for graphing and model analysis.
el.final <- el.final %>%
  mutate(InfStatusFactor = ifelse(PastInfected == "Yes" & MetschInfected == "Yes", "Coinf",
                                  ifelse(PastInfected == "Yes", "Past", 
                                         ifelse(MetschInfected == "Yes", "Metsch", "Uninf"))))
el.final$InfStatusFactor <- as.factor(el.final$InfStatusFactor)
el.final$InfStatusFactor <- factor(el.final$InfStatusFactor, c("Uninf", "Past", "Metsch", "Coinf"))

el.final$DayMetFactor <- factor(el.final$DayMetFactor)
el.final$DayMetFactor <- relevel(el.final$DayMetFactor, "0")

# Cutting T1 from the analyses because so many animals died
el.final <- el.final %>%
  filter(!(TreatmentGroup == "T1"))
# That removed the 13 T1 animals that survived past day 10

# Also cutting T2 from analysis since it doesn't relate to central question (time since Metsch exposure)
el.final <- el.final %>%
  filter(!(TreatmentGroup == "T2"))

#adding a column that is reproduction divided by lifespan 
el.final$egg_per_day <- el.final$LifeRep/el.final$lifesurvival

#Finding the mean lifetime fecundity for each treatment 
el.final_summary <- el.final %>%
  group_by(DayMetNum, InfStatusFactor, MetExpStatus, PastExpStatus) %>%
  summarise(
    n = n(),
    total_fecundity = mean(LifeRep, na.rm = TRUE),
    se_totalfecundity = sd(LifeRep, na.rm = TRUE) / sqrt(n),
    perday_fecundity = mean(egg_per_day, na.rm = TRUE),
    se_perdayfecundity = sd(egg_per_day, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

clean_el.final_summary <- na.omit(el.final_summary)

#Drop the 0 time point 
clean_el.final_summary <- subset(clean_el.final_summary, !DayMetNum == "0")

# Mean for reference line for mean total uninfected eggs 
el.final_uninf <- subset(el.final, InfStatusFactor == "Uninf")
mean(el.final_uninf$LifeRep)  # 114.1034

#Creating the Past panel for the host fitness combo plot 
past_egg_summary <- clean_el.final_summary %>%
  filter(InfStatusFactor == "Past")

e1<-ggplot(past_egg_summary, aes(x = DayMetNum, y = total_fecundity, shape = MetExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  geom_errorbar(aes(ymin = total_fecundity - se_totalfecundity, ymax = total_fecundity + se_totalfecundity), 
                width = 0.1, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  geom_hline(yintercept = 114.1034, linetype = "dashed", color = "gray40") +   
  ylab("Lifetime fecundity") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.position = "none"     # This hides the legend
  ) +
  labs(title = "*P. ramosa*") +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15
    ),
    labels = c(
      "ExposedUninfected" = "exposed but uninfected",
      "Unexposed" = "unexposed"
    )
  ) +
  ylim(0, 150)


#Creating the Metsch panel for the host fitness combo plot 
metsch_egg_summary <- clean_el.final_summary %>%
  filter(InfStatusFactor == "Metsch")

e2<-ggplot(metsch_egg_summary, aes(x = DayMetNum, y = total_fecundity, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 1), color = "tomato") +
  geom_errorbar(aes(ymin = total_fecundity - se_totalfecundity, ymax = total_fecundity + se_totalfecundity), 
                width = 0.1, position = position_dodge(width = 1), color = "tomato") +
  ylab("Lifetime fecundity") +
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
    shape = NULL,
    title = "*A. monospora*"   
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed to *P. ramosa* but uninfected",
        "Unexposed" = "unexposed to *P. ramosa*")) + ylim(0,150) +
  geom_hline(yintercept = 114.1034, linetype = "dashed", color = "gray40")


#Creating the coinfected panel for the host fitness combo plot 
coinf_egg_summary <- clean_el.final_summary %>%
  filter(InfStatusFactor == "Coinf")

e3<-ggplot(coinf_egg_summary, aes(x = DayMetNum, y = total_fecundity)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "#7f39d4") +
  geom_errorbar(aes(ymin = total_fecundity - se_totalfecundity, ymax = total_fecundity + se_totalfecundity), 
                width = 0.1, position = position_dodge(width = 0.4), color = "#7f39d4") +
  ylab("Lifetime fecundity") +
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
  )  + scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title = "Coinfected"  
  ) + ylim(0,150) + geom_hline(yintercept = 114.1034, linetype = "dashed", color = "gray40")


### Analysis of fecundity data

# Look at just infecteds 

# all infected hosts
el.final.inf <- el.final %>%
  subset(!(MetschInfected == "No" & PastInfected == "No"))


# totalegg.glm.inf_gaussian <- glmmTMB(LifeRep ~ InfStatusFactor * DayMetNum, family = gaussian, data = el.final.inf) #DHARMa check won't run on gaussian
# totalegg.glm.inf_nb2 <- glmmTMB(LifeRep ~ InfStatusFactor * DayMetNum, family = nbinom2, data = el.final.inf) # DHARMa indicates significant deviation on QQ plots
totalegg.glm.inf_nb1 <- glmmTMB(LifeRep ~ InfStatusFactor * DayMetNum, family = nbinom1, data = el.final.inf)

totalegg.glm.inf_simResid_nb1 <- simulateResiduals(fittedModel = totalegg.glm.inf_nb1)
plot(totalegg.glm.inf_simResid_nb1) 

summary(totalegg.glm.inf_nb1)
Anova(totalegg.glm.inf_nb1, type = "III") 

emm_trends <- emtrends(totalegg.glm.inf_nb1,
                       ~ InfStatusFactor,
                       var = "DayMetNum")
summary(emm_trends) # testing whether the slope of DayMetNum is significant within each infection class

# Compare slopes across infection classes
pairs(emm_trends) # not using this, though confirms what is obvious from looking at the plot and from emtrends

### Survival ###

#Finding the mean lifetime fecundity for each treatment 
survival_summary <- el.final %>%
  group_by(DayMetNum, InfStatusFactor, MetExpStatus, PastExpStatus) %>%
  summarise(
    n = n(),
    mean_survival = mean(lifesurvival, na.rm = TRUE),
    se_survival = sd(lifesurvival, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

clean_survival_summary <- na.omit(survival_summary)

#Drop never exposed 
clean_survival_summary<-subset(clean_survival_summary, !DayMetNum == "0")

# write.csv(clean_survival_summary, file = "clean_survival_summary.csv", row.names = FALSE)

#Creating a panels for combo plot 

#creating a reference line 
mean(el.final_uninf$lifesurvival)  # 60.36782

past_survival_summary <- clean_survival_summary %>%
  filter(InfStatusFactor == "Past")

s1<-ggplot(past_survival_summary, aes(x = DayMetNum, y = mean_survival, shape = MetExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), 
                width = 0.1, position = position_dodge(width = 0.4), color = "cornflowerblue") +
  ylab("Host lifespan (days)") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = "none",
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    shape = NULL,
    title = "*P. ramosa*" 
  ) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0,80)) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed but uninfected",
        "Unexposed" = "unexposed")) + geom_hline(yintercept = 60.36782, linetype = "dashed", color = "gray40")

metsch_survival_summary <- clean_survival_summary %>%
  filter(InfStatusFactor == "Metsch")

s2<-ggplot(metsch_survival_summary, aes(x = DayMetNum, y = mean_survival, shape = PastExpStatus)) +
  geom_point(size = 4, position = position_dodge(width = 1.25), color = "tomato") +
  geom_errorbar(aes(ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), 
                width = 0.1, position = position_dodge(width = 1.25), color = "tomato") +
  ylab("Host lifespan (days)") +
  xlab(NULL) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    plot.title = ggtext::element_markdown(hjust = 0.5),
    legend.position = c(.01,1),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    shape = NULL,
    title = "*A. monospora*"   
  ) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0,80)) +
  scale_shape_manual(
    values = c(
      "ExposedUninfected" = 16,
      "Unexposed" = 15), labels = c(
        "ExposedUninfected" = "exposed to *P. ramosa* but uninfected",
        "Unexposed" = "unexposed to *P. ramosa*")) + geom_hline(yintercept = 60.36782, linetype = "dashed", color = "gray40")


coinf_survival_summary <- clean_survival_summary %>%
  filter(InfStatusFactor == "Coinf")

s3<-ggplot(coinf_survival_summary, aes(x = DayMetNum, y = mean_survival)) +
  geom_point(size = 4, position = position_dodge(width = 0.4), color = "#7f39d4") +
  geom_errorbar(aes(ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), 
                width = 0.1, position = position_dodge(width = 0.4), color = "#7f39d4") +
  ylab("Host lifespan (days)") +
  xlab("Day of exposure to *A. monospora*") +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(.01, .4),
    legend.justification = "none",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  ) +
  labs(
    title = "Coinfected"
  ) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0,80)) +
 theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = 60.36782, linetype = "dashed", color = "gray40")



comboplot <- (e1 | s1) / (e2 | s2 ) / (e3 | s3) +
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
comboplot
ggsave("figures/comboplot.png", comboplot, dpi = 300, width = 15, height = 18, units = "in")


#Lifespan analysis 
#Survivor analysis - Hazard Ratio

#Censor animals who lived to the end 
el.final$survobject <- with(el.final, Surv(lifesurvival, censor))
seq_lifespan <- survfit(survobject ~ InfStatusFactor, data = el.final, conf.type = "log-log")
summary(seq_lifespan)

#Run Cox proportional hazards with interaction variables 
cox_model <- coxph(survobject ~ InfStatusFactor*DayMetNum, robust = TRUE, data = el.final)
summary(cox_model)
cox.test <- cox.zph(cox_model) 
cox.test # Hazards are not proportional 


#Not surprising that hazards aren't proportional. Will analyze separately for each panel of the figure

# just singly Past infected:
el.final.singlepast <- el.final %>%
  subset((MetschInfected == "No" & PastInfected == "Yes"))

cox_model.singlepast <- coxph(survobject ~ DayMetNum, robust = TRUE, data = el.final.singlepast)
summary(cox_model.singlepast)

cox.test <- cox.zph(cox_model.singlepast) 
cox.test # Hazards are proportional 

# just singly Metsch infected:
el.final.singlemet <- el.final %>%
  subset((MetschInfected == "Yes" & PastInfected == "No"))

cox_model.singlemet <- coxph(survobject ~ DayMetNum, robust = TRUE, data = el.final.singlemet)
summary(cox_model.singlemet)
cox.test <- cox.zph(cox_model.singlemet) 
cox.test # Hazards are not proportional 

# just coinfected
el.final.coinf <- el.final %>%
  subset((MetschInfected == "Yes" & PastInfected == "Yes"))

cox_model.coinf <- coxph(survobject ~ DayMetNum, robust = TRUE, data = el.final.coinf)
summary(cox_model.coinf)

cox.test <- cox.zph(cox_model.coinf)
cox.test  # Hazards are proportional 



# Old version analyzing lifespan by day of death -- these distributions are super wonky

# write.csv(el.final, file = "el.final.csv", row.names = FALSE)

# el.final_died <- subset(el.final, !censor == "1") # MAD doesn't think we want this censoring
lifespan.glm <- glmmTMB(lifesurvival ~ InfStatusFactor*DayMetNum, family = gaussian, data = el.final)
summary(lifespan.glm)


lifespan.glm_simResid <- simulateResiduals(fittedModel = lifespan.glm)
plot(lifespan.glm_simResid)  # looks wonky

# Will look at just each panel of the figure for a separate analysis

# just singly Metsch infected:
el.final.singlemet <- el.final %>%
  subset((MetschInfected == "Yes" & PastInfected == "No"))

lifespan.glm.singlemet <- glmmTMB(lifesurvival ~ DayMetNum, family = gaussian, data = el.final.singlemet)
summary(lifespan.glm.singlemet)
# significant effect of day of Metsch addition

# just coinfected
el.final.coinf <- el.final %>%
  subset((MetschInfected == "Yes" & PastInfected == "Yes"))

lifespan.glm.coinf <- glmmTMB(lifesurvival ~ DayMetNum, family = gaussian, data = el.final.coinf)
summary(lifespan.glm.coinf)
# significant effect of day of Metsch addition on coinfecteds


# just singly Past infected:
el.final.singlepast <- el.final %>%
  subset((MetschInfected == "No" & PastInfected == "Yes"))

lifespan.glm.singlepast <- glmmTMB(lifesurvival ~ DayMetNum, family = gaussian, data = el.final.singlepast)
summary(lifespan.glm.singlepast)
# no significant effect of day of Metsch addition on singly Past infected lifespan



##### Body size figure for supplement ######
#Body size analysis 
#Read in the body size data.
length <- read.csv("data/Coinfection Body Size_without_continous.csv")

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

#Removing T1 individuals because of high mortality
long.data <- long.data %>%
  subset(TreatmentGroup != "T1")

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

max(long.data$bodysize)

#transforming into a linear relationship 
long.data$transformedbody <- (log10(1 -long.data$bodysize/2185.39))
long.data$transformedbody <- (long.data$transformedbody*-1)

#adding a unique identifier to each daphnia for random intercept 
long.data <- transform(long.data, uniqueident = as.numeric(factor(sample.ID)))

#drop any body size that is NA 
long.data <- long.data %>% drop_na(transformedbody)

#looking for outliers 
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
ggsave("figures/bodysize.png", bodysize, dpi = 600, width = 15, height = 12, units = "in")

# Redoing with untransformed body sizes
l1v2<-ggplot(long.data_noinfections, aes(x=time, y=bodysize, color=exposure)) + geom_point(alpha = .4) + 
  geom_smooth(
    se = FALSE,
    method = "gam",
    formula = y ~ s(x, k = 4)
  ) +
  theme_minimal() + 
  xlab("") + ylab("Body size (microns)") +
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
    legend.position = c(0.5, 0.3),
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
l2v2<-ggplot(long.data,
           aes(x = time,
               y = bodysize,
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

l3v2<-ggplot(coinfected_bodysize, aes(x=time, y=bodysize, color=DayMet)) + geom_point(alpha = .4) +
  geom_smooth(
    se = FALSE,
    method = "gam",
    formula = y ~ s(x, k = 4)) + 
  xlab("Time in days") + ylab("Body size (microns)")  +
  scale_color_manual(
    name = "Day of exposure to *A. monospora*", values = c("None" = "grey50", "5" = "mediumpurple1", "10" = "mediumpurple2", "15" = "mediumpurple3", "30" = "mediumpurple4")) +
  theme_classic() +
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(),
    legend.position = c(0.5, 0.3),
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

bodysizev2 <- (l1v2 | l2v2) / (l3v2 | l_blank) +
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

bodysizev2
ggsave("figures/bodysizev2.png", bodysizev2, dpi = 300, width = 15, height = 12, units = "in")

