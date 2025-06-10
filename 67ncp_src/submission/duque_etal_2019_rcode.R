####
### duque_etal_2019_rcode.R
### Created by Jeffrey R. Stevens on 9 Jan 2019 (jeffrey.r.stevens@gmail.com),
###	 finalized on 21 Mar 2019
### Summary: This script analyzes and produces plots for mesotocin effects on social bonding analysis.
### Instructions: Place this file and the data files in the same directory.
### 	Create a folder called "figures". Set the R working directory to the main directory.  
### 	At the R command prompt, type
### 	> source("duque_etal_2019_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving figures in the figures directory. If packages do not
###		load properly, install them with install.packages("package_name").
### Data files:
###  duque_etal_2019_data.csv
###   experiment - experiment number (1 or 2)
###   phase - experimental phase: pair formation (pair) or pair maintenance (group)
###   group - squad number for birds
###   pair - pair number for birds
###   condition - hormone manipulation (SAL = saline, MT = mesotocin, OTA = oxytocin antagonist)
###   session - session number for each pair
###   video - unique for each video recording combining experiment, pair, and session
###   coder - ID for video coder to use for reliability check
###   time - number of minutes into the session
###   distance - distance between heads of pair members (cm)
### License: This script is released under the Creative Commons Attribution 4.0 International 
###   Public License (CC BY 4.0). You are free to:
###     Share — copy and redistribute the material in any medium or format
###     Adapt — remix, transform, and build upon the material for any purpose, even commercially. ###   Under the following terms:
###     Attribution — You must give appropriate credit, provide a link to the license, and
###       indicate if changes were made. You may do so in any reasonable manner, but not in any
###       way that suggests the licensor endorses you or your use.
###     No additional restrictions — You may not apply legal terms or technological measures that
###       legally restrict others from doing anything the license permits.
####

####
# Load libraries and functions --------------------------------------------
####
## Load libraries
library(broom)
library(car)
library(here)
library(lme4)
library(rptR)
library(tidyverse)


###
## Convert BIC values to Bayes factor (from Wagenmakers, E.-J. (2007). A practical solution to the pervasive problems of p values. Psychonomic Bulletin & Review, 14(5), 779–804. https://dx.doi.org/10.3758/BF03194105)
###
bic_bf10 <- function(null, alternative) {
  new_bf <- exp((null - alternative) / 2) # convert BICs to Bayes factor
  names(new_bf) <- NULL   # remove BIC label
  return(new_bf)  # return Bayes factor of alternative over null hypothesis
}

## Create themes for plots
# Theme for figures without legends
theme_plots <- function () { 
  theme_bw(base_size=30, base_family="Arial") %+replace% 
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    )
}

# Theme for figures with legends
theme_legend <- function () { 
  theme_bw(base_size=25, base_family="Arial") %+replace% 
    theme(
      panel.grid = element_blank(),  # remove grid lines
      legend.title = element_blank(),  # remove legend title
      legend.key.width = unit(3, "line"),  # increase width of legend lines
      legend.key = element_rect(fill = "transparent", color = NA),  # make legend background transparent
      legend.background = element_rect(fill = "transparent", color = NA)  # make legend background transparent
    )
}

# Creates qqplots using ggplot2 with a linear model
gg_qqplot <- function(LM) # argument: a linear model
{
  y <- quantile(resid(LM), c(0.25, 0.75))  # extract residuals and calculate quantiles of data
  x <- qnorm(c(0.25, 0.75))  # calculate quantiles of normal distribution
  slope <- diff(y) / diff(x)  # calculate slope of expected normal distribution line
  int <- y[1L] - slope * x[1L]  # calculate intercept of expected normal distribution line
  p <- ggplot(LM, aes(sample = .resid)) +
    stat_qq(alpha = 0.5) +  # plot the residuals against standard normal distribution
    geom_abline(slope = slope, intercept = int, color="blue")  # plot expected normal distribution line
  return(p)
}

set.seed(0) # for consistency in model algorithm computations


####
# Import and prepare data to compare conditions --------------------------------------------
####
## Import data
all_raw <- read_csv("duque_etal_2019_data.csv", col_names = TRUE) %>% 
  mutate_if(is.character, as.factor) %>%  # convert all character columns to factors
  mutate(pair = as.factor(pair),  # convert to factor
         condition = relevel(condition, "SAL"))  # make SAL first level

## Pair-formation phase
# Select and clean data
pair_all <- all_raw %>% 
  filter(phase == "pair" & time > 2) %>%  # select pair phase and remove first two minutes from all data to account for behavioral changes due to human experimenter presence
  group_by(experiment, video, coder, group, pair, session, condition) %>% 
  summarize(distance = mean(distance, na.rm = TRUE))  # calculate mean distance over time block

# Remove duplicates
pair_clean <- pair_all[sample(1:nrow(pair_all)),] %>% .[!duplicated(.$video),]  # drop all duplicates by first shuffling duplicated rows, then randomly selecting whichever row occurs first 
# Separate experiments 1 and 2
pair_data_exp1 <- filter(pair_clean , experiment == 1) %>% # filter experiment 1
  group_by(pair, condition) %>%  # group by pair and condition
  summarize(mean_distance = mean(distance, na.rm = TRUE))  # calculate mean distance for each pair and condition
pair_data_exp2 <- filter(pair_clean , experiment == 2) %>% # filter experiment 2
  group_by(pair, condition) %>%  # group by pair and condition
  summarize(mean_distance = mean(distance, na.rm = TRUE))  # calculate mean distance for each pair and condition

## Pair-maintenance (group) phase
# Select and clean data
group_clean <- all_raw %>% 
  filter(phase == "group" & time > 2) %>%  # select group phase and remove first two minutes from all data to account for behavioral changes due to human experimenter presence
  group_by(experiment, group, pair, session, condition) %>% 
  summarize(distance = mean(distance, na.rm = TRUE))  # calculate mean distance over time block

# Separate experiments 1 and 2
group_data_exp1 <- filter(group_clean, experiment == 1) %>% # filter experiment 1
  group_by(pair, condition) %>%  # group by pair and condition
  summarize(mean_distance = mean(distance, na.rm = TRUE))  # calculate mean distance for each pair
group_data_exp2 <- filter(group_clean, experiment == 2) %>% # filter experiment 2
  group_by(pair, condition) %>%  # group by pair and condition
  summarize(mean_distance = mean(distance, na.rm = TRUE))  # calculate mean distance for each pair


####
# Inter-coder reliability --------------------------------------------
####

# Extract videos scored by multiple coders
reliability <- pair_all[which(duplicated(pair_all$video)), ]   # find all videos that are duplicated (with multiple coders)

# Calculate intraclass correlation (ICC): what is the correlation between two randomly selected distances of the same video? 
reliability$video <- factor(reliability$video)  # convert to factor
icc_reliability <- lmer(distance ~ (1 | video), data = reliability)  # conduct linear mixed model with video as random effect 
icc_reliability_CI <- rpt(distance ~ (1 | video), grname = "video", data = reliability, datatype = "Gaussian",  nboot = 1000, npermut = 0)  # calculate ICC and 95% CIs (bootstrap)

icc_value <- as.numeric(icc_reliability_CI$R*100)  # save ICC value
icc_lower <- icc_reliability_CI$CI_emp$`2.5%` * 100  # save ICC lower 95% CI
icc_upper <- icc_reliability_CI$CI_emp$`97.5%` * 100  # save ICC upper 95% CI

# RESULTS: In an empty, random intercept model, 8.74992/(8.74992+0.19519)= 97.82% [0.97, 0.99] of the variation in social distance is accounted for by differences between different videos (i.e., there is excellent agreement in distances quantified from the same video by different coders)


####
# Plot condition effects on distance --------------------------------------------
####

## Produce box plots for both pair and group phase
ggplot(pair_data_exp1, aes(x = condition, y = mean_distance, color = condition)) + 
  # geom_boxplot() +  # plot boxplot
  geom_point(shape = 1, size = 4) +  # plot individual data points
  stat_summary(fun.data = "mean_cl_normal", size = 1.25) +  # plot 95% CIs
  labs(x = "Hormone condition", y = "Pair proximity (cm)") +  # label axes
  theme_plots() # use theme
ggsave(here("figures/box_pair_phase_exp1.png"), width = 7, height = 6)
ggplot(pair_data_exp2, aes(x = condition, y = mean_distance, color = condition)) + 
  geom_boxplot() +  # plot boxplot
  # geom_point(shape = 1, size = 4) +  # plot individual data points
  stat_summary(fun.data = "mean_cl_normal", size = 1.25) +   # plot 95% CIs
  labs(x = "Hormone condition", y = "Pair proximity (cm)") +  # label axes
  scale_colour_discrete(drop = TRUE, limits = levels(pair_data_exp1$condition)) +  # keep three condition colors
  theme_plots() # use theme
ggsave(here("figures/box_pair_phase_exp2.png"), width = 7, height = 6)
ggplot(group_data_exp1, aes(x = condition, y = mean_distance, color = condition)) + 
  # geom_boxplot() +  # plot boxplot
  geom_point(shape = 1, size = 4) +  # plot individual data points
  stat_summary(fun.data = "mean_cl_normal", size = 1.25) +   # plot 95% CIs
  labs(x = "Hormone condition", y = "Pair proximity (cm)") +  # label axes
  theme_plots() # use theme
ggsave(here("figures/box_group_phase_exp1.png"), width = 7, height = 6)
ggplot(group_data_exp2, aes(x = condition, y = mean_distance, color = condition)) + 
  geom_boxplot() +  # plot boxplot
  # geom_point(shape = 1, size = 4) +  # plot individual data points
  stat_summary(fun.data = "mean_cl_normal", size = 1.25) +   # plot 95% CIs
  labs(x = "Hormone condition", y = "Pair proximity (cm)") +  # label axes
  scale_colour_discrete(drop = TRUE, limits = levels(pair_data_exp1$condition)) +  # keep three condition colors
  theme_plots() # use theme
ggsave(here("figures/box_group_phase_exp2.png"), width = 7, height = 6)

####
# Plot session and condition effects on distance --------------------------------------------
####

## Session analysis
# Calculate session means
pair_data_session <- pair_clean  %>%
  group_by(experiment, pair, session, condition) %>%  # for each pair and session
  summarize(mean_distance = mean(distance, na.rm = TRUE))  # calculate mean distance
group_data_session <- group_clean %>%
  group_by(experiment, pair, session, condition) %>%  # for each pair and timeblock
  summarize(mean_distance = mean(distance, na.rm = TRUE))  # calculate mean distance

# Generate plots
experiment_names <- c('1' = 'Experiment 1', '2' = 'Experiment 2')
# Pair-formation phase
ggplot(pair_data_session, aes(x = session, y = mean_distance)) +
  stat_summary(fun.data = "mean_cl_normal", geom = "line", aes(linetype = condition, color = condition), position = position_dodge2(width = 0.5)) +  # plot lines connecting means
  stat_summary(fun.data = "mean_cl_normal", fatten = 6, aes(shape = condition, color = condition), position = position_dodge2(width = 0.35)) +  # plot means as points and CIs
  geom_smooth(method = "lm", col = "black") +
  facet_wrap(~experiment, labeller = as_labeller(experiment_names)) +  # facet by experiment
  labs(x = "Session", y = "Pair proximity (cm)") +  # label axes
  scale_x_continuous(breaks = 1:10) +  # control x-axis breaks
  theme_legend() +  # use theme
  theme(legend.position = c(0.86, 0.85), legend.title=element_blank(), legend.key.width=unit(3, "line"))  # place legend, remove title, and increase legend line length
ggsave(here("figures/session_pair_phase.png"), width = 8, height = 6)

# Pair-maintenance phase
ggplot(group_data_session, aes(x = session, y = mean_distance)) +
  stat_summary(fun.data = "mean_cl_normal", geom = "line", aes(linetype = condition, color = condition), position = position_dodge2(width = 0.35)) +  # plot lines connecting means
  stat_summary(fun.data = "mean_cl_normal", fatten = 6, aes(shape = condition, color = condition), position = position_dodge2(width = 0.35)) +  # plot means as points and CIs
  geom_smooth(method = "lm", col = "black") +
  facet_wrap(~experiment, labeller = as_labeller(experiment_names)) +  # facet by experiment
  labs(x = "Session", y = "Pair proximity (cm)") +  # label axes
  scale_x_continuous(breaks = 1:10) +  # control x-axis breaks
  theme_legend() +  # use theme
  theme(legend.position = c(0.86, 0.85), legend.title=element_blank(), legend.key.width=unit(3, "line"))  # place legend, remove title, and increase legend line length
ggsave(here("figures/session_group_phase.png"), width = 8, height = 6)


####
# Data cleaning for statistical analyses --------------------------------------------
####

# Format all variables as needed for analyses
pair_clean <- pair_clean %>%  
  mutate(session0 = session - 10,  # assign session 10 to 0 for more meaningful intercept interpretation 
    session0_sqrd = -1 * abs(session0)^2) # create quadratic effect of session

group_clean <- group_clean %>%  
  mutate(session0 = session - 10,  # assign session 10 to 0 for more meaningful intercept interpretation 
    session0_sqrd = -1 * abs(session0)^2) # create quadratic effect of session

# Create four datasets: one for each experiment for both pair-formation and pair-maintenance phase
pair_data1 <- filter(pair_clean, experiment == 1)
pair_data2 <- filter(pair_clean, experiment == 2)
group_data1 <- filter(group_clean, experiment == 1)
group_data2 <- filter(group_clean, experiment == 2)


####
# Model selection for four data sets  ----------------------------------------
# First find best-fitting random effect model, then test for fixed effect predictors
# Process is backward by elimintating weakest terms sequentially, starting with full model, until only significant effects remain
# Nested model comparisons (likelihood ratio tests using anova command) are used to select best fitting models
####


### _ Pair Experiment 1 ------------------------------------------------------

## Random effects: selecting empty model
pair1_rand_full <- lmer(distance ~ (1 + session0 | pair) + (1 | group), data = pair_data1) # full random model; trending p=.09 against reduced model 2
pair1_rand2 <- lmer(distance ~ (1 + session0 | pair), data = pair_data1) # drops group ran.intercept; BEST
pair1_rand3 <- lmer(distance ~ (1 | pair), data = pair_data1) # drops random slope

# Likelihood Ratio Tests for nested model comparison
pair1_rand_anova <- anova(pair1_rand_full, pair1_rand2, pair1_rand3) # empty model = pair1_rand2
pair1_rand_anova_tidy <- tidy(pair1_rand_anova)

## Fixed effects (inference does not change whether group random intercept is retained)
pair1_full <- lmer(distance ~ condition * session0 + session0_sqrd + (1 + session0 | pair), data = pair_data1)  # full model
pair1_fixed <- lmer(distance ~ condition + session0 + session0_sqrd + (1 + session0 | pair), data = pair_data1) # drop interaction first, then drop weakest variables
pair1_fixed2 <- lmer(distance ~ condition + session0 + (1 + session0 | pair), data = pair_data1) # drop quadratic session
pair1_fixed3 <- lmer(distance ~ condition + (1 + session0 | pair), data = pair_data1) # drop session (about as good as condition)
pair1_fixed4 <- lmer(distance ~ session0 + (1 + session0 | pair), data = pair_data1) # drop condition (about as good as session)

# Likelihood Ratio Tests for nested model comparison
# No fixed effects significantly change model
pair1_fixed_anova <- anova(pair1_full, pair1_fixed, pair1_fixed2, pair1_fixed3, pair1_fixed4, pair1_rand2) # empty model (pair1_rand2) is equal to all 
pair1_fixed_anova_tidy <- tidy(pair1_fixed_anova)  # create tidy table of model parameters

## Bayes factors for fixed effects
# Extract BICs
pair1_rand_bic <- pair1_fixed_anova_tidy$BIC[1]  # random model
pair1_full_bic <- pair1_fixed_anova_tidy$BIC[6]
pair1_fixed_bic <- pair1_fixed_anova_tidy$BIC[5]
pair1_fixed2_bic <- pair1_fixed_anova_tidy$BIC[4]
pair1_fixed3_bic <- pair1_fixed_anova_tidy$BIC[3]
pair1_fixed4_bic <- pair1_fixed_anova_tidy$BIC[3]

# Convert BICs to BFs compared to null model
pair1_full_bf <- bic_bf10(pair1_rand_bic, pair1_full_bic)
pair1_fixed_bf <- bic_bf10(pair1_rand_bic, pair1_fixed_bic)
pair1_fixed2_bf <- bic_bf10(pair1_rand_bic, pair1_fixed2_bic)
pair1_fixed3_bf <- bic_bf10(pair1_rand_bic, pair1_fixed3_bic)
pair1_fixed4_bf <- bic_bf10(pair1_rand_bic, pair1_fixed4_bic)

## Check assumptions
leveneTest(residuals(pair1_rand2) ~ pair_data1$pair) # assumption of equal variances NOT met; visual suggests slight overdispersion in residuals but not overly problematic; slightly worse at predicting greater distances
plot(pair1_rand2)  # plot residuals vs. predicted values
plot(density(residuals(pair1_rand2)))  # plot distribution of residuals
gg_qqplot(pair1_rand2) # plot qqplot; MOSTLY normally distributed but worse at extremes, slightly more so for large values

## Results
# In experiment 1-pair phase, the best-fitting random effect structure included a random intercept for each unique pair and a random slope over sessions; i.e., allowing pairs to change independently over time (random intercept model for pair with versus without random slope:  χ2(2)=18.53, p<0.001).  However, a random intercept for each group was not warranted (full versus model without group: χ2(1)=2.85, p=0.09).  Although MT-pairs generally perched closer over time than SAL- or OTA- pairs, inclusion of condition, session, their interaction, or quadratic effect of session (time) did not significantly improve an empty model (same random effects with no fixed effects). 


### _ Pair Experiment 2 ------------------------------------------------------

## Random effects: selecting empty model
pair2_rand_full <- lmer(distance ~ (1 + session0 | pair) + (1 | group), data = pair_data2) # fit full model; OVERFIT: group not warranted
pair2_rand2 <- lmer(distance ~ (1 + session0 | pair), data = pair_data2) # drops group ran.int; BEST model
pair2_rand3 <- lmer(distance ~ (1 | pair), data = pair_data2) # drops random slope; worse without

# Likelihood Ratio Tests for nested model comparison
pair2_rand_anova <- anova(pair2_rand_full, pair2_rand2, pair2_rand3)
pair2_rand_anova_tidy <- tidy(pair2_rand_anova)  # create tidy table of model parameters

## Fixed effects (inference does not change whether group random intercept is retained)
pair2_full <- lmer(distance ~ condition * session0 + session0_sqrd + (1 + session0 | pair), data = pair_data2)  # full model
pair2_fixed <- lmer(distance ~ condition + session0 + session0_sqrd + (1 + session0 | pair), data = pair_data2) # drop interaction first, then drop weakest variables
pair2_fixed2 <- lmer(distance ~ session0 + session0_sqrd + (1 + session0 | pair), data = pair_data2) # drop condition; BEST model
pair2_fixed3 <- lmer(distance ~ session0 + (1 + session0 | pair), data = pair_data2) # drop quadratic session

# Likelihood Ratio Tests for nested model comparison
pair2_fixed_anova <- anova(pair2_full, pair2_fixed, pair2_fixed2, pair2_fixed3, pair2_rand2) 
pair2_fixed_anova_tidy <- tidy(pair2_fixed_anova)  # create tidy table of model parameters
pair2_fixed_lmer <- tidy(pair2_fixed2)  # create tidy table of model parameters

## Bayes factors for fixed effects
# Extract BICs
pair2_rand_bic <- pair2_fixed_anova_tidy$BIC[which(pair2_fixed_anova_tidy$term == "pair2_rand2")]  # random model
pair2_full_bic <- pair2_fixed_anova_tidy$BIC[which(pair2_fixed_anova_tidy$term == "pair2_full")]
pair2_fixed_bic <- pair2_fixed_anova_tidy$BIC[which(pair2_fixed_anova_tidy$term == "pair2_fixed")]
pair2_fixed2_bic <- pair2_fixed_anova_tidy$BIC[which(pair2_fixed_anova_tidy$term == "pair2_fixed2")]
pair2_fixed3_bic <- pair2_fixed_anova_tidy$BIC[which(pair2_fixed_anova_tidy$term == "pair2_fixed3")]

# Convert BICs to BFs
pair2_full_bf <- bic_bf10(pair2_rand_bic, pair2_full_bic)
pair2_fixed_bf <- bic_bf10(pair2_rand_bic, pair2_fixed_bic)
pair2_fixed2_bf <- bic_bf10(pair2_rand_bic, pair2_fixed2_bic)
pair2_fixed3_bf <- bic_bf10(pair2_rand_bic, pair2_fixed3_bic)

## Check assumptions
leveneTest(residuals(pair2_fixed2) ~ pair_data2$pair) # assumption of equal variances NOT met; slightly worse at predicting extremes but likely not problematic
plot(pair2_fixed2)   # plot residuals vs. predicted values
plot(density(residuals(pair2_fixed2)))
gg_qqplot(pair2_fixed2) # plot qqplot; MOSTLY normally distributed but worse at extremes, slightly more so for large positive values

## Results
# In experiment 2-pair phase, the best-fitting random effect structure included a random intercept for each unique pair and a random slope over sessions (random intercept model for pair with versus without random slope: χ2(2)=22.31, p<0.001). However, a random intercept for each group was not warranted (overfit full model versus model without group: ; χ2(1)=0.0, p>0.99). Both a linear and quadratic fixed effect of session were warranted (model including linear with versus without quadratic session: χ2(2)=8.16, p=0.004), indicating that pairs perched 0.37 ± 0.12 (standard error) cm closer each subsequent session, but the decrease in distance diminished by 0.04 ± 0.01 cm each session.  That is, though pairs perched more closely over time, the reduction in distance was less pronounced as time progressed. Lastly, inclusion of condition was not warranted (χ2(1)=0.35, p=0.55).  


### _ Group Experiment 1 ------------------------------------------------------

## Random effects: selecting empty model
group1_rand_full <- lmer(distance ~ (1 + session0 | pair) + (1 | group), data = group_data1) # fit full model; OVERFIT: Maximal
group1_rand2 <- lmer(distance ~ (1 | pair) + (1 | group), data = group_data1) # drop random slope 
group1_rand3 <- lmer(distance ~ (1 | pair), data = group_data1) # drop group; Model improves by dropping group, but PAIR ONLY is not different from MAXIMAL: since maximal is overfit, seems only pairs significantly explain where the variability lies
group1_rand4 <- lm(distance ~ 1, data = group_data1) # empty model

# Likelihood Ratio Tests for nested model comparison
group1_rand_anova <- anova(group1_rand_full, group1_rand2, group1_rand3, group1_rand4) # a random intercept for each pair is significantly better than empty, otherwise all other random effects tested were not warranted
group1_rand_anova_tidy <- tidy(group1_rand_anova)  # create tidy table of model parameters
group1_rand_lmer <- tidy(group1_rand3)  # create tidy table of model parameters

## Fixed effects (inference does not change whether group random intercept is retained)
group1_full <- lmer(distance ~ condition * session0 + session0_sqrd + (1 | pair), data = group_data1)  # fit full model
group1_fixed <- lmer(distance ~ condition * session0 + (1 | pair), data = group_data1) # drop slope first
group1_fixed2 <- lmer(distance ~ condition + session0 + (1 | pair), data = group_data1) # drop interaction
group1_fixed3 <- lmer(distance ~ session0 + (1 | pair), data = group_data1) # drop condition (next weakest)

# Likelihood Ratio Tests for nested model comparison
# models with session significantly better than empty (however, keeping condition does not change model); that is, including condition (model 2 vs 3)  does not help or hurt, therefore model without condition chosen for parsimony
group1_fixed_anova <- anova(group1_full, group1_fixed, group1_fixed2, group1_fixed3, group1_rand3)
group1_fixed_anova_tidy <- tidy(group1_fixed_anova)  # create tidy table of model parameters
group1_fixed_lmer <- tidy(group1_fixed3)   # create tidy table of model parameters

## Bayes factors for fixed effects
# Extract BICs
group1_rand_bic <- group1_fixed_anova_tidy$BIC[which(group1_fixed_anova_tidy$term == "group1_rand3")]  # random model
group1_full_bic <- group1_fixed_anova_tidy$BIC[which(group1_fixed_anova_tidy$term == "group1_full")]
group1_fixed_bic <- group1_fixed_anova_tidy$BIC[which(group1_fixed_anova_tidy$term == "group1_fixed")]
group1_fixed2_bic <- group1_fixed_anova_tidy$BIC[which(group1_fixed_anova_tidy$term == "group1_fixed2")]
group1_fixed3_bic <- group1_fixed_anova_tidy$BIC[which(group1_fixed_anova_tidy$term == "group1_fixed3")]

# Convert BICs to BFs
group1_full_bf <- bic_bf10(group1_rand_bic, group1_full_bic)
group1_fixed_bf <- bic_bf10(group1_rand_bic, group1_fixed_bic)
group1_fixed2_bf <- bic_bf10(group1_rand_bic, group1_fixed2_bic)
group1_fixed3_bf <- bic_bf10(group1_rand_bic, group1_fixed3_bic)

## Check assumptions
leveneTest(residuals(group1_fixed3) ~ group_data1$pair) # assumption of equal variances met
plot(group1_fixed3)  # plot residuals vs. predicted values
plot(density(residuals(group1_fixed3)))  # plot distribution of residuals
gg_qqplot(group1_fixed3) # plot qqplot; normally distributed

## Results
# In experiment1-group phase, the best-fitting random effect structure included only a random intercept for each unique pair (against null model with no random effects; x2(1)=5.05, p=0.025). A linear fixed effect of session was warranted (against empty model; x2(1)=6.12, p=0.013), indicating that pairs perched 1.4 ± 0.56 cm closer each subsequent session. However, no other fixed effects tested (condition, session, their interaction, or quadratic effect of session) were warranted. 


### _ Group Experiment 2 ------------------------------------------------------

## Random effects
group2_rand_full <- lmer(distance ~ (1 + session0 | pair) + (1 | group), data = group_data2) # fit full  model; OVERFIT: Maximal
group2_rand2 <- lmer(distance ~ (1 | pair) + (1 | group), data = group_data2) # drop random slope first; BEST
group2_rand3 <- lmer(distance ~ (1 | pair), data = group_data2) # drop group; model is worse without
group2_rand4 <- lmer(distance ~ (1 | group), data = group_data2) # drop pair instead; model is worse without
group2_rand5 <- lm(distance ~ 1, data = group_data2) # drop ALL variables, i.e. a grand mean model

# Likelihood Ratio Tests for nested model comparison
group2_rand_anova <- anova(group2_rand_full, group2_rand2, group2_rand3, group2_rand4, group2_rand5) 
group2_rand_anova_tidy <- tidy(group2_rand_anova)  # create tidy table of model parameters
group2_rand_lmer <- tidy(group2_rand2)  # create tidy table of model parameters

## Fixed effects (inference does not change whether group random intercept is retained)
group2_full <- lmer(distance ~ condition * session0 + session0_sqrd + (1 | pair) + (1 | group), data = group_data2)  # fit full model
group2_fixed <- lmer(distance ~ condition + session0 + session0_sqrd + (1 | pair) + (1 | group), data = group_data2) # drop interaction first
group2_fixed2 <- lmer(distance ~ session0 + session0_sqrd + (1 | pair) + (1 | group), data = group_data2) # drop condition
group2_fixed3 <- lmer(distance ~ session0 + (1 | pair) + (1 | group), data = group_data2) # drop quadratic session

# Likelihood Ratio Tests for nested model comparison
# No significant predictors; all models equivalent therefore empty model chosen for parsimony
group2_fixed_anova <- anova(group2_full, group2_fixed, group2_fixed2, group2_fixed3, group2_rand2) # empty model (group2_rand2 from Line340) is equal to all 
group2_fixed_anova_tidy <- tidy(group2_fixed_anova)  # create tidy table of model parameters

## Bayes factors for fixed effects
# Extract BICs
group2_rand_bic <- group2_fixed_anova_tidy$BIC[which(group2_fixed_anova_tidy$term == "group2_rand2")]  # random model
group2_full_bic <- group2_fixed_anova_tidy$BIC[which(group2_fixed_anova_tidy$term == "group2_full")]
group2_fixed_bic <- group2_fixed_anova_tidy$BIC[which(group2_fixed_anova_tidy$term == "group2_fixed")]
group2_fixed2_bic <- group2_fixed_anova_tidy$BIC[which(group2_fixed_anova_tidy$term == "group2_fixed2")]
group2_fixed3_bic <- group2_fixed_anova_tidy$BIC[which(group2_fixed_anova_tidy$term == "group2_fixed3")]

# Convert BICs to BFs
group2_full_bf <- bic_bf10(group2_rand_bic, group2_full_bic)
group2_fixed_bf <- bic_bf10(group2_rand_bic, group2_fixed_bic)
group2_fixed2_bf <- bic_bf10(group2_rand_bic, group2_fixed2_bic)
group2_fixed3_bf <- bic_bf10(group2_rand_bic, group2_fixed3_bic)

## Check assumptions
leveneTest(residuals(group2_rand2) ~ group_data2$pair) # assumption of equal variances met
plot(group2_rand2)  # plot residuals vs. predicted values
plot(density(residuals(group2_rand2)))  # plot distribution of residuals; slight skew positive tail
gg_qqplot(group2_rand2) # plot qqplot; NOT normally distributed, overshooting extremes

## Results
# In experiment2-group phase, the best-fitting random effect structure included a random intercept for each unique pair and group, but not a random slope (full model with versus without random slope: x2(2)=5.62, p=0.018). Inclusion of condition, session, their interaction, or quadratic effect of session did not significantly improve an empty model (same random effects with no fixed effects). 

# Thus, although a considerable amount of variation is attributable to differences across pairs, and how those pairs change distances across sessions, hormone condition was not warranted in the best-fitting model for either experiment in either phase.
