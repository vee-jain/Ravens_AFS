#' Varalika Jain
#' Project:
#' How common ravens (Corvus corax) exploit anthropogenic food sources (AFSs)
#' through time and space in a semi-transformed, alpine environment

#------------------------------------------------------------------------------
####----(i) LOAD LIBRARIES----####
#' For better understanding of the script, information from the associated 
#' guides or vignettes of the following packages have also been included
library(plyr)
library(dplyr)
library(glmmADMB)
library(car)
library(MuMIn)
library(sjPlot)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(i) READ IN MASTER DATAFRAME----####
#' For details into how 'master.csv' was created, please see the complementary
#' 'master.R' document
master_final <- read.csv("master_final.csv")
head(master_final) 

master_final$id_name <- as.factor(master_final$id_name)
master_final$season <- as.factor(master_final$season)
master_final$year <- as.factor(master_final$year)
master_final$Origin <- as.factor(master_final$Origin)
master_final$Sex <- as.factor(master_final$Sex)
master_final$age_class <- as.factor(master_final$age_class)

####----(1) Summaries----####
#' Summary:
#' The total number of individuals
plyr::count(master_final$id_name)
#' A total of 85 individuals

#' Origin
summary(master_final %>% group_by(id_name, Origin) %>% tally())
#' 29 captive, 56 wild

#' Sex
summary(master_final %>% group_by(id_name, Sex) %>% tally())
#' 51 females, 34 males

#' group adults & sub-adults
summary(master_final %>% group_by(id_name, age_class) %>% tally())
# 12 adults, 66 juveniles, 44 sub-adults
master_final$age_class <- revalue(master_final$age_class, c("sub-adult"="adult"))
summary(master_final %>% group_by(id_name, age_class) %>% tally())
# adult juvenile 
# 53      66  

####----(2) Rescale the predictor variable fixes_by_days----####
class(master_final$fixes_by_days)
master_final$z.fixes_by_days <- scale(master_final$fixes_by_days)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(I) OCCURRENCE DISTRIBUTION----####
####----(1) Prepare dataframe for model----####
#' Average occurrence distribution for each individual-season-year combination 
occ_dist <- master_final %>% dplyr::group_by(id_name, season, year, age_class,
                                             Sex, Origin, z.fixes_by_days, area_95) %>% tally()

plyr::count(occ_dist$id_name)
names(occ_dist)

#' Getting rid of the n column
occ_dist <- occ_dist[,-9]

####----(2) Summaries for each categorical variable----####
summary(occ_dist$age_class)
#' adult juvenile 
#' 224      152 
summary(occ_dist$Origin)
#' captive bred  wild caught 
#' 117          259 
summary(occ_dist$Sex)
#' f   m 
#' 215 161

summary(occ_dist$season)
#' aut spr sum win 
#' 117  65  91 103 

summary(occ_dist$year)
#' 2017 2018 2019 2020 
#' 15  115  197   49

####----(3) The glmm model ----####
#' Log-normal distribution
occ_dist_glmm <- glmmadmb(log(area_95) ~ age_class + season + Sex + Origin + 
                            z.fixes_by_days + year + 
                            (1| id_name),
                          family = "gaussian",
                          data = occ_dist)
summary(occ_dist_glmm)
confint(occ_dist_glmm)
vif(occ_dist_glmm)

####----(4) Dredging for the top model set ----####
#' change the default "na.omit" to "na.fail" to prevent
#' models from being fitted to different datasets in case of missing values
options(na.action = "na.fail") 

occ_dist_set <- dredge(occ_dist_glmm, beta = FALSE, evaluate = TRUE, rank = "AICc")
occ_dist_set

####----(5) Averaging the model with delta AIC < 6 cut-off----####
occ_dist_d6 <- get.models(occ_dist_set, subset = delta < 6)

#' Delta 6 model average
occ_dist_modavg <- model.avg(occ_dist_d6)
summary(occ_dist_modavg)
confint(occ_dist_modavg)
sw(occ_dist_modavg)

#' change back to "na.omit"
options(na.action = "na.omit")

####----(6) Plot effect sizes----####
#' plot model
plot_model(occ_dist_modavg)

####----(7) Set different intercept for season----####
occ_dist$season <- relevel(occ_dist$season, ref = "win")
#' Now run step 3, 4, and 5 again

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(II) AVERAGE MAXIMUM DAILY DISPLACEMENT----####
####----(1) Prepare dataframe for model----####
max_disp <- master_final %>%
  group_by(id_name, date) %>%
  arrange(desc(displacement)) %>%
  slice(1L)
head(max_disp)

#' Average displacement for each individual-season-year combination 
avg_max_disp <- max_disp %>% group_by(id_name, age_class, season, 
                                      Sex, Origin, z.fixes_by_days, year) %>%
  dplyr::summarise(avg_displacement = mean(displacement))

head(avg_max_disp)

#' Summaries same as I.2_

####----(2) The glmm model----####
#' Log-normal distribution
avg_max_disp_glmm <- glmmadmb(log(avg_displacement) ~ age_class + season + Sex + Origin + 
                                z.fixes_by_days + year + 
                                (1| id_name),
                              family = "gaussian",
                              data = avg_max_disp)
summary(avg_max_disp_glmm)
confint(avg_max_disp_glmm)
vif(avg_max_disp_glmm)

####----(3) Dredging for the top model set ----####
#' change "na.omit" to "na.fail"
options(na.action = "na.fail") 

avg_max_disp_set <- dredge(avg_max_disp_glmm, beta = FALSE, evaluate = TRUE, rank = "AICc")
avg_max_disp_set

####----(4) Averaging the model with delta AIC < 6 cut-off----####
avg_max_disp_d6 <- get.models(avg_max_disp_set, subset = delta < 6)

#' Delta 6
avg_max_disp_modavg <- model.avg(avg_max_disp_d6)
summary(avg_max_disp_modavg)
confint(avg_max_disp_modavg)
sw(avg_max_disp_modavg)

#' change back to "na.omit"
options(na.action = "na.omit")

####----(5) Plot effect sizes----####
#' plot model
plot_model(avg_max_disp_modavg)

####----(6) Set different intercept for season----####
avg_max_disp$season <- relevel(avg_max_disp$season, ref = "win")
#' Now run step 2, 3, and 4 again
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(III) NUMBER OF AFSs VISITED----####
####----(1) Prepare dataframe for model----####
#' Group by each individual-season-year combination 
num_afs <- master_final %>% dplyr::group_by(id_name, age_class, season, 
                                            Sex, Origin, z.fixes_by_days, year, 
                                            intersection) %>% tally()

#' Get rid of n column 
names(num_afs)
num_afs <- num_afs[,-9]

#' Count the number of unique intersections between GPS point and radii
num_afs <- num_afs %>% group_by(id_name, age_class, season, 
                                Sex, Origin, z.fixes_by_days, year) %>% 
  tally()

num_afs$n <- (num_afs$n)-1
#' Rename the column
names(num_afs)[names(num_afs) == "n"] <- "number_afs"

####----(2) The glmm model----####
#' Poisson distribution
num_afs_glmm <- glmmadmb(number_afs ~ age_class + season + Sex + Origin + 
                           z.fixes_by_days + year +
                           (1| id_name), 
                         family = "poisson",
                         data = num_afs)
summary(num_afs_glmm)
confint(num_afs_glmm)
vif(num_afs_glmm)

####----(3) Dredging for the top model set----####
#' change "na.omit" to "na.fail"
options(na.action = "na.fail") 

num_afs_set <- dredge(num_afs_glmm, beta = FALSE, evaluate = TRUE, rank = "AICc")
num_afs_set

####----(4) Averaging the model with delta AIC < 6 cut-off----####
num_afs_d6 <- get.models(num_afs_set, subset = delta < 6)

#' Delta 6
num_afs_modavg <- model.avg(num_afs_d6)
summary(num_afs_modavg)
confint(num_afs_modavg)
sw(num_afs_modavg)

#' change back to "na.omit"
options(na.action = "na.omit")

####----(5) Plot effect sizes----####
#' plot model
plot_model(num_afs_modavg)

####----(6) Set different intercept for season----####
num_afs$season <- relevel(num_afs$season, ref = "aut")
#' Now run step 2, 3, and 4 again
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(III) PROBABILITY OF BEING AT ANY AFS----####
####----(1) Prepare dataframe for model----####
#' Group by each individual-season-year combination 
prob_afs <- master_final %>% group_by(id_name, age_class, season, 
                                      Sex, Origin, year, fix_total,
                                      location) %>% tally
#' Rename the column
names(prob_afs)[names(prob_afs) == "n"] <- "number_pts"

#' Updating location column from 'AFS & NA' to 'inside & out'
prob_afs$location <- as.character(prob_afs$location)
prob_afs$location <- revalue(prob_afs$location, c("AFS" = "inside"))
prob_afs$location <- revalue(prob_afs$location, c(" " = "out"))

#' Subset all the inside points from location
#' For this model, we were interested in those individuals that were present at 
#' at least one AFS
prob_afs <- subset(prob_afs, location == "inside")
names(prob_afs)[names(prob_afs) == "number_pts"] <- "inside"

#' Create a separate out column 
prob_afs$out <- (prob_afs$fix_total-prob_afs$inside)

####----(2) Summaries ----####
#' There are 76 individuals

levels(prob_afs$id_name)

summary(prob_afs$age_class)
#' adult juvenile 
#' 220      141 
summary(prob_afs$Origin)
#' captive bred  wild caught 
#' 106          255
summary(prob_afs$Sex)
#' f   m 
#' 206 155

prob_afs$season <- as.factor(prob_afs$season)
summary(prob_afs$season)
#' aut spr sum win 
#' 111  65  82 103

prob_afs$year <- as.factor(prob_afs$year)
summary(prob_afs$year)
#' 2017 2018 2019 2020 
#' 12  111  189   49

####----(3) The glmm model----####
#' Creating an observation level random effect 
prob_afs$obs <- as.factor(seq_len(nrow(prob_afs)))
class(prob_afs$obs)

#' Binomial distribution with observation level random effect (OLRE)
prob_afs_glmm <- glmmadmb(cbind(inside,out) ~ age_class + season +
                            Origin + Sex + year + (1|id_name) + 
                            (1 | obs),
                          data = prob_afs,
                          family = "binomial")
summary(prob_afs_glmm)
confint(prob_afs_glmm)
vif(prob_afs_glmm)

####----(4) Dredging for the top model set----####
#' change "na.omit" to "na.fail"
options(na.action = "na.fail")

prob_afs_set <- dredge(prob_afs_glmm, beta = FALSE, evaluate = TRUE, rank = "AICc")
prob_afs_set

####----(5) Averaging the model with delta AIC < 6 cut-off----####
prob_afs_d6 <- get.models(prob_afs_set, subset = delta < 6)


#' Delta 6
prob_afs_modavg <- model.avg(prob_afs_d6)
summary(prob_afs_modavg)
confint(prob_afs_modavg)
sw(prob_afs_modavg)

#' change back to "na.omit"
options(na.action = "na.omit")

####----(6) Plot effect sizes----####
#' plot model
plot_model(prob_afs_modavg)

####----(7) Set different intercept for season----####
prob_afs$season <- relevel(prob_afs$season, ref = "aut")
#' Now run step 2, 3, and 4 again
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(iii) SUMMARIES---####
#' Function
summary <- function(x, group, z){
  x %>% 
    group_by({{group}}) %>%
    dplyr::summarise(avg=(mean({{z}})), 
                     n=n(),sd=sd({{z}}), se=(sd/sqrt(n)),
                     min=min({{z}}), max=max({{z}}))%>%
    mutate(lower.ci = avg - qt(1 - (0.05 / 2), n - 1) * se,
           upper.ci = avg + qt(1 - (0.05 / 2), n - 1) * se) %>%
    as.data.frame()
}

####----(1) Total number of fixes (whole dataset)----####
totalfixes <- master_final %>% group_by(id_name, fix_total) %>% tally()
totalfixes <- totalfixes %>% group_by(id_name) %>% dplyr::summarise(fixes = sum(fix_total))
sum(totalfixes$fixes)
#' 3620795

summary(x = totalfixes, z = fixes)
#'       avg  n      sd       se min    max lower.ci upper.ci
#'  42597.59 85 69842.3 7575.461  31 252791 27532.95 57662.22

####----(2) Total number of days (whole dataset)----####
totaldays <- master_final %>% group_by(id_name, tracking_days) %>% tally()
totaldays <- totaldays %>% group_by(id_name) %>% dplyr::summarise(days = sum(tracking_days))

summary(x = totaldays, z = days)
#'      avg  n       sd       se min max lower.ci upper.ci
#'  267.7765 85 185.6622 20.13789   3 615 227.7301 307.8229

####----(3) Occurrence distribution----####
#' ha to km
occ_dist$area_95 <- occ_dist$area_95/100

summary(x = occ_dist , z = area_95)
#'        avg  n      sd       se min    max lower.ci upper.ci
#'  130.0419 376 369.9558 19.07901 0.01 3470.11 92.52663 167.5571

summary(occ_dist, age_class, area_95)
#'                  avg  n      sd       se min    max lower.ci upper.ci
#'    adult     171.51036 224 448.9692 29.99801 0.05 3470.11 112.39450 230.62621
#' juvenile  68.93046 152 189.4840 15.36918 0.01 1614.42  38.56405  99.29688

summary(occ_dist, Origin ,area_95)
#'                     avg  n      sd       se min    max lower.ci upper.ci
#' captive bred  39.69376 117 298.2944 27.57733 0.01 3218.53 -14.92661  94.31414
#'  wild caught  170.85552 259 391.8607 24.34904 0.05 3470.11 122.90736 218.80368

summary(occ_dist, Sex ,area_95)
#'         avg  n      sd       se min    max lower.ci upper.ci
#' f 121.4485 215 367.9719 25.09548 0.02 3470.11 71.98249 170.9144
#' m 141.5176 161 373.4288 29.43031 0.01 3218.53 83.39561 199.6395

summary(occ_dist, season ,area_95)
#'           avg  n      sd       se min    max lower.ci upper.ci
#' aut 113.10855 117 202.19470 18.692907 0.05 1182.98  76.084890 150.13220
#' spr  16.96585  65  33.04970  4.099311 0.06  158.33   8.776532  25.15516
#' sum  19.87165  91  38.15291  3.999512 0.01  204.18  11.925920  27.81738
#' win 317.97039 103 632.26258 62.298683 0.11 3470.11 194.401245 441.53953

####----(4) Average maximum daily displacement----####
#' m to km
avg_max_disp$avg_displacement <- avg_max_disp$avg_displacement/1000

summary(x = avg_max_disp , z = avg_displacement)
#'    avg   n       sd       se      min      max lower.ci upper.ci
#'4.83571 376 4.972341 0.2564289 0.02634308 23.79954 4.331491 5.339928

summary(avg_max_disp, age_class, avg_displacement)
#'              avg   n       sd       se      min      max lower.ci upper.ci
#'  adult 5.540434 224 4.793866 0.3203037 0.02634308 23.79954 4.909225 6.171644
#'juvenile 3.797168 152 5.063555 0.4107085 0.05062096 21.43991 2.985691 4.608646


summary(avg_max_disp, Origin ,avg_displacement)
#'                    avg   n       sd       se      min      max lower.ci upper.ci
#' captive bred 1.735190 117 2.028275 0.1875141 0.05062096 17.46578 1.363795 2.106585
#'  wild caught 6.236331 259 5.268568 0.3273729 0.02634308 23.79954 5.591667 6.880994

summary(avg_max_disp, Sex ,avg_displacement)
#'        avg   n       sd       se      min      max lower.ci upper.ci
#' f 4.865475 215 5.159441 0.3518709 0.02634308 23.79954 4.171898 5.559052
#' m 4.795961 161 4.726571 0.3725060 0.05062096 21.43991 4.060298 5.531624


summary(avg_max_disp, season ,avg_displacement)
#'          avg   n       sd       se      min      max lower.ci upper.ci
#' win 6.361086 103 6.166981 0.6076507 0.14757927 23.79954 5.155814 7.566358
#' aut 4.851167 117 4.654847 0.4303408 0.17975614 21.83819 3.998823 5.703511
#' spr 4.094934  65 4.358951 0.5406613 0.20448900 19.38879 3.014839 5.175029
#' sum 3.618436  91 3.730759 0.3910899 0.02634308 16.11753 2.841468 4.395405

####----(5) Number of AFSs visited----####

summary(x = num_afs , z = number_afs)
#' avg   n       sd       se      min  max lower.ci upper.ci
#' 2.611702 376 2.002204 0.1032558   0  12 2.408669 2.814735

summary(num_afs, age_class ,number_afs)
#'          avg   n       sd       se      min  max lower.ci upper.ci
#'    adult 2.946429 224 1.881414 0.1257073   0   9 2.698702 3.194155
#' juvenile 2.118421 152 2.077741 0.1685271   0  12 1.785445 2.451397

summary(num_afs, Origin ,number_afs)
#'                 avg   n       sd       se      min  max lower.ci upper.ci
#' captive bred 1.623932 117 1.165055 0.1077094   0   5  1.41060 1.837264
#'   wild caught 3.057915 259 2.138727 0.1328940   0  12  2.79622 3.319610

summary(num_afs, Sex ,number_afs)
#'     avg   n       sd       se      min  max lower.ci upper.ci
#' f 2.483721 215 1.916106 0.1306774   0   9 2.226141 2.741301
#' m 2.782609 161 2.105634 0.1659472   0  12 2.454879 3.110338

summary(num_afs, season ,number_afs)
#'         avg   n       sd       se      min  max lower.ci upper.ci
#'  aut 2.290598 117 1.650809 0.1526173   0   7 1.988320 2.592876
#'  spr 2.753846  65 2.430476 0.3014634   1  12 2.151604 3.356089
#'  sum 1.978022  91 1.382253 0.1448995   0   6 1.690154 2.265890
#'  win 3.446602 103 2.247998 0.2215018   1   9 3.007254 3.885950

####----(6) Probability of being present at an AFS----####
#' proportion of fixes inside a buffer
prob_afs$p_inside <- (prob_afs$inside/ prob_afs$out)*100

summary(x = prob_afs, z = p_inside)
#'     avg   n       sd       se       min      max lower.ci upper.ci
#' 31.69736 361 63.38842 3.336232 0.6570842 693.5569 25.13641 38.25831

summary(prob_afs, age_class, p_inside)
#'            avg   n       sd       se       min      max lower.ci upper.ci
#' adult 28.47266 220 62.21363 4.194442 0.6570842 693.5569 20.20602 36.73930
#' juvenile 36.72881 141 65.08282 5.480962 1.5360983 557.3723 25.89265 47.56496

summary(prob_afs, Origin, p_inside)
#'                   avg   n       sd       se       min      max lower.ci upper.ci
#' captive bred 31.34683 106 47.49854 4.613466 1.1259676 358.8145 22.19918 40.49448
#'  wild caught 31.84307 255 69.00881 4.321499 0.6570842 693.5569 23.33254 40.35360

summary(prob_afs, Sex, p_inside)
#'        avg   n       sd       se       min      max lower.ci upper.ci
#' f 29.16545 206 57.79590 4.026831 0.6570842 693.5569 21.22613 37.10476
#' m 35.06235 155 70.18816 5.637649 1.1259676 557.3723 23.92525 46.19946

summary(prob_afs, season, p_inside)
#'          avg   n       sd       se       min      max lower.ci upper.ci
#' aut 14.21341 111  10.09418  0.9580973 0.6570842  52.77383 12.31469 16.11214
#' spr 71.05861  65 116.20378 14.4133049 2.8634361 693.55695 42.26472 99.85249
#' sum 39.68355  82  63.82111  7.0478624 1.6528926 358.81455 25.66051 53.70658
#' win 19.34172 103  28.96550  2.8540560 1.2048193 294.71264 13.68072 25.00273

####----(7) Extensively used AFSs----####
popular_afs <- master_final %>% group_by(id_name, intersection) %>% tally()
popular_afs <- popular_afs[,-3]
popular_afs<- popular_afs %>% drop_na("intersection") %>%
  group_by(intersection) %>% tally()
#------------------------------------------------------------------------------

####----END----####