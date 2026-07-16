# Theme: T helper response and exercise - meta-analysis 
# Author: Dr. Florian Javelle 
# Date last update: 05/05/2026

########################################################################### 
#                       PREPARE WORKSPACE                                 #
###########################################################################  

### Clear workspace
rm(list = ls())

### Clear console
cat("\014")

### Install packages 
install.packages("meta")
library(meta)

### Collect required libraries
requiredlibraries <- c("metafor", "tidyr", "dplyr", "data.table", "stringr", "readr",
                       "ggplot2", "reshape2", "cowplot", "gridExtra",
                       "forestplot","metaviz") # funnel plots and Egger's regression

## Loop over required libraries and install if not present (dependencies might have to be installed manually)
for (i in requiredlibraries) {
  if (!(i %in% rownames(installed.packages()))) install.packages(i, dependencies = TRUE)
}

## Loop over required libraries to load
for (i in requiredlibraries) {
  current_lib <- i
  if (!i=="meta"){
    library(current_lib, character.only=TRUE)}
}


########################################################################### 
#                             GATHER DATA Th1                             #
###########################################################################  

### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Th1/pre_post")


### Charge the dataset
### What's the name of the data file to read?
data_file_name <- "Th1_pre_post_r_050.csv"

### For the sensitivity analyses
#data_file_name <- "Th1_pre_post_030.csv"
#data_file_name <- "Th1_pre_post_040.csv"
#data_file_name <- "Th1_pre_post_060.csv"
#data_file_name <- "Th1_pre_post_070.csv"


## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         staining = as.character(staining),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))




######################################################################### 
#                    MAIN ANALYSIS - Th1                               #
#########################################################################

### Meta-analysis from already computed effect sizes 
m1 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m1, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")

# Kostrzewa-Nowak et al. (2019) was detected as very different across meta-analysis and while not influential on the 
# overall effect, its values were unrealistic and brought substantial misleading variance 
fullData <-fullData %>% 
  dplyr::filter(!Study== "Kostrzewa-Nowak et al. (2019)") 
m1 <-metagen(g, seTE=SE_g, data=fullData, studlab=paste(Study), 
               random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
               sm="SMD", fixed=FALSE)

print(m1, overall=TRUE)

#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"

## Creation of functions 
metaoutliers <- function(y, s2, model){
  if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  w <- 1/s2
  y.p <- sum(y*w)/sum(w)
  n <- length(y)
  
  if(missing(model)){
    hetmeasure <- metahet.base(y, s2)
    Ir2 <- hetmeasure$Ir2
    if(Ir2 < 0.3){
      model <- "FE"
      cat("This function uses fixed-effect meta-analysis because Ir2 < 30%.\n")
    }else{
      model <- "RE"
      cat("This function uses random-effects meta-analysis because Ir2 >= 30%.\n")
    }
  }
  
  if(!is.element(model, c("FE", "RE"))) stop("wrong input for the argument model.")
  
  y.p.i <- res <- std.res <- numeric(n)
  if(model == "FE"){
    for(i in 1:n){
      w.temp <- w[-i]
      y.temp <- y[-i]
      y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1/sum(w.temp) + s2[i]
      std.res[i] <- res[i]/sqrt(var.res.i)
    }
  }else{
    for(i in 1:n){
      s2.temp <- s2[-i]
      y.temp <- y[-i]
      tau2.temp <- metahet.base(y.temp, s2.temp)$tau2.DL
      w.temp <- 1/(s2.temp + tau2.temp)
      y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1/sum(w.temp) + s2[i] + tau2.temp
      std.res[i] <- res[i]/sqrt(var.res.i)
    }
  }
  
  outliers <- which(abs(std.res) >= 3)
  if(length(outliers) == 0) outliers <- "All the standardized residuals are smaller than 3"
  
  out <- NULL
  out$model <- model
  out$std.res <- std.res
  out$outliers <- outliers
  
  class(out) <- "metaoutliers"
  return(out)
}

metahet.base <- function(y, s2){
  if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  n <- length(y)
  
  w <- 1/s2
  mu.bar <- sum(w*y)/sum(w)
  
  out <- NULL
  
  out$weighted.mean <- mu.bar
  
  # the conventional methods
  Q <- sum(w*(y - mu.bar)^2)
  H <- sqrt(Q/(n - 1))
  I2 <- (Q - n + 1)/Q
  tau2.DL <- (Q - n + 1)/(sum(w) - sum(w^2)/sum(w))
  tau2.DL <- max(c(0, tau2.DL))
  out$Q <- Q
  out$H <- H
  out$I2 <- I2
  out$tau2.DL <- tau2.DL
  
  # absolute deviation based on weighted mean
  Qr <- sum(sqrt(w)*abs(y - mu.bar))
  Hr <- sqrt((3.14159*Qr^2)/(2*n*(n - 1)))
  Ir2 <- (Qr^2 - 2*n*(n - 1)/3.14159)/(Qr^2)
  tau2.r <- tau2.r.solver(w, Qr)
  out$Qr <- Qr
  out$Hr <- Hr
  out$Ir2 <- Ir2
  out$tau2.r <- tau2.r
  
  # absolute deviation based on weighted median
  
  expit <- function(x) {ifelse(x >= 0, 1/(1 + exp(-x/0.0001)), exp(x/0.0001)/(1 + exp(x/0.0001)))}
  psi <- function(x) {sum(w*(expit(x - y) - 0.5))}
  mu.med <- uniroot(psi, c(min(y) - 0.001, max(y) + 0.001))$root
  out$weighted.median <- mu.med
  Qm <- sum(sqrt(w)*abs(y - mu.med))
  Hm <- sqrt(3.14159/2)*Qm/n
  Im2 <- (Qm^2 - 2*n^2/3.14159)/Qm^2
  tau2.m <- tau2.m.solver(w, Qm)
  out$Qm <- Qm
  out$Hm <- Hm
  out$Im2 <- Im2
  out$tau2.m <- tau2.m
  
  return(out)
}

tau2.r.solver <- function(w, Qr){
  f <- function(tau2){
    out <- sum(sqrt(1 - w/sum(w) + tau2*(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))) - Qr*sqrt(3.14159/2)
    return(out)
  }
  f <- Vectorize(f)
  tau.upp <- Qr*sqrt(3.14159/2)/sum(sqrt(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))
  tau2.upp <- tau.upp^2
  f.low <- f(0)
  f.upp <- f(tau2.upp)
  if(f.low*f.upp > 0){
    tau2.r <- 0
  }else{
    tau2.r <- uniroot(f, interval = c(0, tau2.upp))
    tau2.r <- tau2.r$root
  }
  return(tau2.r)
}

tau2.m.solver <- function(w, Qm){
  f <- function(tau2){
    out <- sum(sqrt(1 + w*tau2)) - Qm*sqrt(3.14159/2)
    return(out)
  }
  f <- Vectorize(f)
  n <- length(w)
  tau2.upp <- sum(1/w)*(Qm^2/n*2/3.14159 - 1)
  tau2.upp <- max(c(tau2.upp, 0.01))
  f.low <- f(0)
  f.upp <- f(tau2.upp)
  if(f.low*f.upp > 0){
    tau2.m <- 0
  }else{
    tau2.m <- uniroot(f, interval = c(0, tau2.upp))
    tau2.m <- tau2.m$root
  }
  return(tau2.m)
}


### Detect outliers
outlier<-metaoutliers(y=m1$TE, s2=m1$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# No outcome have standardized residuals above 3 or below -3. 

## Funnel plot to visualize the results more easily
pdf("Th1_funnel.pdf", width = 10, height = 18)
funnel(m1, random=TRUE, xlim=c(-1.5, 4), ylim=c(0.7, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.99, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)

# One value is <-1 and one is >1 and thus are considered as being influential cases. 
# Nonetheless, those values are not outliers (Kostrzewa-Nowak & Nowak (2020a)a and Tsubakihara et al. (2012)).
# Accordingly, we retain all outcomes and aim to identify the sources of heterogeneity in the moderator analysis.

### Forest plot general
pdf("Th1_no_out_forest.pdf", width = 10, height = 18)

forest(m1, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-2, 5))

dev.off()

####################################### ASYMMETRY ##########################################

### Egger's test (one tailed thus p<.1)
metabias(m1, method="linreg")
#==> Asymmetry


### Fit random-effects model
res <- rma(g, V_g, SE_g, data=fullData, measure="SMD", method="HE")
res

### Trim-and-fill analysis
taf <- trimfill(res)
print(taf)
# No missing studies 

################################################################################################
#                                Moderator analysis - Th1                                      #
################################################################################################

#################################### Categorical moderators ####################################

### RISK of BIAS
## Moderator analysis 
m1mod <- update(m1, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m1mod, overall=FALSE)
# Not significant results  

## Forest plot 
pdf("subgroup_Th1_RoB.pdf", width = 11.25, height = 18)

forest(m1mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

### staining
## Moderator analysis 
m3mod <- update(m1, subgroup=staining, print.subgroup.name=TRUE, fixed=FALSE)
print(m3mod, overall=FALSE)
# Not significant results  

### SEX
## Moderator analysis 
m2mod <- update(m1, subgroup=sex, print.subgroup.name=TRUE, fixed=FALSE)
print(m2mod, overall=FALSE)
# Significant results but 14 vs 1 vs 2

## Forest plot 
pdf("subgroup_Th1_Sex.pdf", width = 11.25, height = 18)

forest(m2mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

### Intensity
## Moderator analysis 
# Let's exclude any studies that do not report 'Intensity' to avoid counting them by default in the analysis.
fullData_bis<-fullData %>% 
  dplyr::filter(!intensity== "") 
m1_c <-metagen(g, seTE=SE_g, data=fullData_bis, studlab=paste(Study), 
               random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
               sm="SMD", fixed=FALSE)

m3mod <- update(m1_c, subgroup=intensity, print.subgroup.name=TRUE, fixed=FALSE)
print(m3mod, overall=FALSE)
# Significant results

## Forest plot 
pdf("subgroup_Th1_intensity.pdf", width = 11.25, height = 18)

forest(m3mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

#################################### Continuous moderators ####################################

### Age 
## We do directly the subgroup_analysis
m1mr <- metareg(m1, age)
summary(m1mr)
# Not significant results

## Print a bubble plot 
bubble(m1mr, lwd=2, col.line="blue", xlab="Age (years)", ylab="Hedges' g")


### BMI 
## We do directly the subgroup_analysis
m2mr <- metareg(m1, bmi)
summary(m2mr)
# Not significant results

## Print a bubble plot 
bubble(m2mr, lwd=2, col.line="blue", xlab="BMI (kg/m2)", ylab="Hedges' g")

### Duration 
## We do directly the subgroup_analysis
m3mr <- metareg(m1, duration)
summary(m3mr)
# Significant results

## Print a bubble plot 
bubble(m3mr, lwd=2, col.line="blue", xlab="duration (min)", ylab="Hedges' g")

### Dose 
## We do directly the subgroup_analysis
m4mr <- metareg(m1, dose)
summary(m4mr)
# Significant results - collinear with duration

## Print a bubble plot 
bubble(m4mr, lwd=2, col.line="blue", xlab="DOSE", ylab="Hedges' g")
# Dose and duration are almost the distributed the same way 

### IFN 
## We do directly the subgroup_analysis
m6mr <- metareg(m1, IFN)
summary(m6mr)
# Not significant results

## Print a bubble plot 
bubble(m6mr, lwd=2, col.line="blue", xlab="[IFN gamma]", ylab="Hedges' g")

### IL_4 
## We do directly the subgroup_analysis
m7mr <- metareg(m1, IL_4)
summary(m7mr)
# Not significant results

## Print a bubble plot 
bubble(m7mr, lwd=2, col.line="blue", xlab="[IL-4]", ylab="Hedges' g")


### IL_6 
## We do directly the subgroup_analysis
m8mr <- metareg(m1, IL_6)
summary(m8mr)
# Not significant results

## Print a bubble plot 
bubble(m8mr, lwd=2, col.line="blue", xlab="[IL-6]", ylab="Hedges' g")

### IL_10 
## We do directly the subgroup_analysis
m9mr <- metareg(m1, IL_10)
summary(m9mr)
# Not significant results

## Print a bubble plot 
bubble(m9mr, lwd=2, col.line="blue", xlab="[IL-10]", ylab="Hedges' g")

### IL_2 
## We do directly the subgroup_analysis
m10mr <- metareg(m1, IL_2)
summary(m10mr)
# Not significant results

## Print a bubble plot 
bubble(m10mr, lwd=2, col.line="blue", xlab="[IL-2]", ylab="Hedges' g")


#################################### MULTIVARIATE ####################################

rma(yi = g, vi=V_g, mods= ~ intensity*duration, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results

rma(yi = g, vi=V_g, mods= ~ IFN*IL_2, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results

#################################### INTERACTION PLOT (HEATMAP) ####################################

## Clean data (remove non-finite values)
fullData_clean <- fullData %>%
  dplyr::filter(is.finite(IFN),
                is.finite(IL_2),
                is.finite(g),
                is.finite(V_g))

## Fit model
model <- metafor::rma(yi = g, vi = V_g,
                      mods = ~ IFN * IL_2,
                      data = fullData_clean,
                      method = "HE",
                      test = "knha")

## Create prediction grid
grid <- expand.grid(
  IFN = seq(min(fullData_clean$IL_10), max(fullData_clean$IFN), length.out = 50),
  IL_2  = seq(min(fullData_clean$IL_2),  max(fullData_clean$IL_2),  length.out = 50)
)

## Predict
X <- model.matrix(~ IFN * IL_2, data = grid)
preds <- predict(model, newmods = X[, -1])

grid$pred <- preds$pred

## Plot heatmap
ggplot2::ggplot(grid, ggplot2::aes(x = IFN, y = IL_2, fill = pred)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::labs(x = "IFN", y = "IL-2", fill = "Hedges' g") +
  ggplot2::theme_minimal()
########################################################################### 
#                 GATHER DATA Th1 - Th1 17 hours post ex                  #
###########################################################################  

### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Th1/pre_post17h")

### Charge the dataset
## What's the name of the data file to read?
data_file_name <- "Th1_pre_post17h_050.csv"

### For the sensitivity analyses
#data_file_name <- "Th1_pre_post17h_030.csv"
#data_file_name <- "Th1_pre_post17h_040.csv"
#data_file_name <- "Th1_pre_post17h_060.csv"
#data_file_name <- "Th1_pre_post17h_070.csv"

## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         met = as.numeric(met),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))



######################################################################### 
#               MAIN ANALYSIS - Th1 17 hours post ex                   #
#########################################################################

### Meta-analysis from already computed effect sizes 
m2 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m2, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")


# Kostrzewa-Nowak et al. (2019) was detected as very different across meta-analysis and while not influential on the 
# overall effect, its values were unrealistic and brought substantial misleading variance 
fullData <-fullData %>% 
  dplyr::filter(!Study== "Kostrzewa-Nowak et al. (2019)") 
m2 <-metagen(g, seTE=SE_g, data=fullData, studlab=paste(Study), 
             random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
             sm="SMD", fixed=FALSE)

print(m2, overall=TRUE)
#################################### OUTLIERS DETECTION ####################################                        

### Detect outliers
outlier<-metaoutliers(y=m2$TE, s2=m2$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# 1 outcome have standardized residuals below -3 (Steensberg et al. (2001)). Before characterizing it as outlier, let's check the funnel plot.

## Funnel plot to visualize the results more easily
pdf("Th1_post17h_funnel.pdf", width = 10, height = 18)
funnel(m2, random=TRUE, xlim=c(-1, 4), ylim=c(0.8, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=1.8184, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)

# One value is <-1 and thus is considered as being an influential case. As it 
# also was an outlier we are going to withdraw it from further analyses 

fullData_bis<-fullData %>% 
  dplyr::filter(!Study== "Steensberg et al. (2001)") 
m2_1 <-metagen(g, seTE=SE_g, data=fullData_bis, studlab=paste(Study), 
               random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
               sm="SMD", fixed=FALSE)
print(m2_1)      

## Funnel plot without the outlier
pdf("Th1_post17h_funnel.pdf", width = 10, height = 18)
funnel(m2_1, random=TRUE, xlim=c(-1, 4), ylim=c(0.8, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=1.84, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()

### Forest plot general
pdf("Th1_post17h_no_out_forest.pdf", width = 10, height = 18)

forest(m2_1, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-2, 5))

dev.off()

####################################### ASYMMETRY ##########################################

### Egger's test (one tailed thus p<.1)
metabias(m2_1, method="linreg")

#==>Asymmetry, thus  trim and fill 

### Fit random-effects model
res <- rma(g, V_g, SE_g, data=fullData_bis, measure="SMD", method="HE")
res

### Trim-and-fill analysis
taf <- trimfill(res)
print(taf)
# No missing studies 

################################################################################################
#                           Moderator analysis - Th1 17h post ex                               #
################################################################################################

#################################### Categorical moderators ####################################
### RISK of BIAS
## Moderator analysis 
m1mod <- update(m2_1, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m1mod, overall=FALSE)
# Not significant results  

## Forest plot 
pdf("subgroup_Th1_17hpost_RoB.pdf", width = 11.25, height = 18)

forest(m1mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

### Intensity
## Moderator analysis 
# All studies have the intensity 

#################################### Continuous moderators ####################################

### Age 
## We do directly the subgroup_analysis
m1mr <- metareg(m2_1, age)
summary(m1mr)
# Not significant 

## Print a bubble plot 
bubble(m1mr, lwd=2, col.line="blue", xlab="Age (years)", ylab="Hedges' g")

### BMI 
## We do directly the subgroup_analysis
m2mr <- metareg(m2_1, bmi)
summary(m2mr)
# Not significant results

## Print a bubble plot 
bubble(m2mr, lwd=2, col.line="blue", xlab="BMI (kg/m2)", ylab="Hedges' g")

### Duration 
## We do directly the subgroup_analysis
m3mr <- metareg(m2_1, duration)
summary(m3mr)
# Not significant results

## Print a bubble plot 
bubble(m3mr, lwd=2, col.line="blue", xlab="duration (min)", ylab="Hedges' g")

### Dose 
## We do directly the subgroup_analysis
m4mr <- metareg(m2_1, dose)
summary(m4mr)
# Not significant results

## Print a bubble plot 
bubble(m4mr, lwd=2, col.line="blue", xlab="DOSE", ylab="Hedges' g")


### IFN 
## We do directly the subgroup_analysis
m6mr <- metareg(m2_1, IFN)
summary(m6mr)
# Not significant results

## Print a bubble plot 
bubble(m6mr, lwd=2, col.line="blue", xlab="[IFN gamma]", ylab="Hedges' g")

### IL_4 
## We do directly the subgroup_analysis
m7mr <- metareg(m2_1, IL_4)
summary(m7mr)
# Not significant results

## Print a bubble plot 
bubble(m7mr, lwd=2, col.line="blue", xlab="[IL-4]", ylab="Hedges' g")


### IL_6 
## We do directly the subgroup_analysis
m8mr <- metareg(m2_1, IL_6)
summary(m8mr)
# Not significant results

## Print a bubble plot 
bubble(m8mr, lwd=2, col.line="blue", xlab="[IL-6]", ylab="Hedges' g")

### IL_10 
## We do directly the subgroup_analysis
m9mr <- metareg(m2_1, IL_10)
summary(m9mr)
# Not significant results

## Print a bubble plot 
bubble(m9mr, lwd=2, col.line="blue", xlab="[IL-10]", ylab="Hedges' g")

### IL_2 
## We do directly the subgroup_analysis
m10mr <- metareg(m2_1, IL_2)
summary(m10mr)
# Not significant results

## Print a bubble plot 
bubble(m10mr, lwd=2, col.line="blue", xlab="[IL-2]", ylab="Hedges' g")


#################################### MULTIVARIATE ####################################

rma(yi = g, vi=V_g, mods= ~ IFN*IL_2, data=fullData_bis, slab = paste(Study), method = "HE", test="knha")
# Not significant results

########################################################################### 
#                             GATHER DATA Th2                             #
###########################################################################  
### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Th2/pre_post")

### Charge the dataset
## What's the name of the data file to read?
data_file_name <- "Th2_pre_post_r_050.csv"

### For the sensitivity analyses
#data_file_name <- "Th2_pre_post_030.csv"
#data_file_name <- "Th2_pre_post_040.csv"
#data_file_name <- "Th2_pre_post_060.csv"
#data_file_name <- "Th2_pre_post_070.csv"

## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         staining = as.character(staining),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))

######################################################################### 
#                    MAIN ANALYSIS - Th2                               #
#########################################################################

### Meta-analysis from already computed effect sizes 
m3 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m3, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")

# Kostrzewa-Nowak et al. (2019) was detected as very different across meta-analysis and while not influential on the 
# overall effect, its values were unrealistic and brought substantial misleading variance 
fullData <-fullData %>% 
  dplyr::filter(!Study== "Kostrzewa-Nowak et al. (2019)") 
m3 <-metagen(g, seTE=SE_g, data=fullData, studlab=paste(Study), 
             random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
             sm="SMD", fixed=FALSE)

print(m3, overall=TRUE)
#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"

### Detect outliers
outlier<-metaoutliers(y=m3$TE, s2=m3$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# One outcome has standardized residuals around -2 (Kostrzewa-Nowak et al. (2019)) . 

## Funnel plot to visualize the results more easily
pdf("Th2_funnel.pdf", width = 10, height = 18)
funnel(m3, random=TRUE, xlim=c(-1, 2), ylim=c(0.7, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.3214, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)
# No influential case 

### Forest plot general
pdf("Th2_no_out_forest.pdf", width = 10, height = 18)

forest(m3, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-4, 3))

dev.off()

####################################### ASYMMETRY ##########################################

### Egger's test (one tailed thus p<.1)
metabias(m3, method="linreg")
#==> No Asymmetry, thus No trim and fill 

### Fit random-effects model
#res <- rma(g, VE_g, SE_g, data=fullData, measure="SMD", method="HE")
#res

### Trim-and-fill analysis
#taf <- trimfill(res)
#print(taf)
# No missing studies 


################################################################################################
#                                Moderator analysis - Th2                                      #
################################################################################################


#################################### Categorical moderators ####################################

### RISK of BIAS
## Moderator analysis 
m1mod <- update(m3, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m1mod, overall=FALSE)
#significant results  

## Forest plot 
pdf("subgroup_Th2_RoB.pdf", width = 11.25, height = 18)

forest(m1mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()



##HEAD TO HEAD comparison 
# Modearate vs Serious 

fullData_bis<-fullData %>% 
  dplyr::filter(!RoB== "Low")

fullData_bis<-fullData_bis %>% 
  dplyr::filter(!RoB== "Critical")

m2_a_ms <-metagen(g, V_g, data=fullData_bis, studlab=paste(Study), 
                  random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
                  sm="SMD", fixed=FALSE)

m2_1mod_ms <- update(m2_a_ms, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m2_1mod_ms, overall=FALSE)
# Significant  

# Moderate vs Low 

fullData_bis<-fullData %>% 
  dplyr::filter(!RoB== "Serious")

fullData_bis<-fullData_bis %>% 
  dplyr::filter(!RoB== "Critical")

m2_a_ml <-metagen(g, V_g, data=fullData_bis, studlab=paste(Study), 
                  random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
                  sm="SMD", fixed=FALSE)

m2_1mod_ml <- update(m2_a_ml, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m2_1mod_ml, overall=FALSE)
# Significant  

# Low vs Serious 

fullData_bis<-fullData %>% 
  dplyr::filter(!RoB== "Critical")

fullData_bis<-fullData_bis %>% 
  dplyr::filter(!RoB== "Moderate")

m2_a_ls <-metagen(g, V_g, data=fullData_bis, studlab=paste(Study), 
                  random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
                  sm="SMD", fixed=FALSE)

m2_1mod_ls <- update(m2_a_ls, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m2_1mod_ls, overall=FALSE)
# Not Significant  




### SEX
## Moderator analysis 
m2mod <- update(m3, subgroup=sex, print.subgroup.name=TRUE, fixed=FALSE)
print(m2mod, overall=FALSE)
# Significant results but 15 vs 1 vs 2

## Forest plot 
pdf("subgroup_Th2_Sex.pdf", width = 11.25, height = 18)

forest(m2mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

### Intensity
## Moderator analysis 
# Let's exclude any studies that do not report 'Intensity' to avoid counting them by default in the analysis.
fullData_bis<-fullData %>% 
  dplyr::filter(!intensity== "") 
m1_c <-metagen(g, seTE=SE_g, data=fullData_bis, studlab=paste(Study), 
               random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
               sm="SMD", fixed=FALSE)

m3mod <- update(m1_c, subgroup=intensity, print.subgroup.name=TRUE, fixed=FALSE)
print(m3mod, overall=FALSE)
# Not significant results

## Forest plot 
pdf("subgroup_Th2_intensity.pdf", width = 11.25, height = 18)

forest(m3mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

### staining
## Moderator analysis 
m4mod <- update(m3, subgroup=staining, print.subgroup.name=TRUE, fixed=FALSE)
print(m4mod, overall=FALSE)
# Significant results 

#################################### Continuous moderators ####################################

### Age 
## We do directly the subgroup_analysis
m1mr <- metareg(m3, age)
summary(m1mr)
#significant results

## Print a bubble plot 
bubble(m1mr, lwd=2, col.line="blue", xlab="Age (years)", ylab="Hedges' g")

### BMI 
## We do directly the subgroup_analysis
m2mr <- metareg(m3, bmi)
summary(m2mr)
# Not significant results

## Print a bubble plot 
bubble(m2mr, lwd=2, col.line="blue", xlab="BMI (kg/m2)", ylab="Hedges' g")

### Duration 
## We do directly the subgroup_analysis
m3mr <- metareg(m3, duration)
summary(m3mr)
# Not significant results

## Print a bubble plot 
bubble(m3mr, lwd=2, col.line="blue", xlab="duration (min)", ylab="Hedges' g")

### Dose 
## We do directly the subgroup_analysis
m4mr <- metareg(m3, dose)
summary(m4mr)
# Not significant results

## Print a bubble plot 
bubble(m4mr, lwd=2, col.line="blue", xlab="DOSE", ylab="Hedges' g")

### IFN 
## We do directly the subgroup_analysis
m6mr <- metareg(m3, IFN)
summary(m6mr)
# Not significant results

## Print a bubble plot 
bubble(m6mr, lwd=2, col.line="blue", xlab="[IFN gamma]", ylab="Hedges' g")

### IL_4 
## We do directly the subgroup_analysis
m7mr <- metareg(m3, IL_4)
summary(m7mr)
# Not significant results

## Print a bubble plot 
bubble(m7mr, lwd=2, col.line="blue", xlab="[IL-4]", ylab="Hedges' g")


### IL_6 
## We do directly the subgroup_analysis
m8mr <- metareg(m3, IL_6)
summary(m8mr)
# Not Significant results

## Print a bubble plot 
bubble(m8mr, lwd=2, col.line="blue", xlab="[IL-6]", ylab="Hedges' g")

### IL_10 
## We do directly the subgroup_analysis
m9mr <- metareg(m3, IL_10)
summary(m9mr)
# Not significant results

## Print a bubble plot 
bubble(m9mr, lwd=2, col.line="blue", xlab="[IL-10]", ylab="Hedges' g")

### IL_2 
## We do directly the subgroup_analysis
m10mr <- metareg(m3, IL_2)
summary(m10mr)
# Not significant results

## Print a bubble plot 
bubble(m10mr, lwd=2, col.line="blue", xlab="[IL-2]", ylab="Hedges' g")


#################################### MULTIVARIATE ####################################

rma(yi = g, vi=V_g, mods= ~ intensity*duration, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results

rma(yi = g, vi=V_g, mods= ~ IL_10*IL_4, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results

########################################################################### 
#                 GATHER DATA Th2  17h post exercise                      #
###########################################################################  

### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Th2/pre_post17h")

### Charge the dataset
## What's the name of the data file to read?
data_file_name <- "Th2_pre_post17h_050.csv"

### For the sensitivity analyses
#data_file_name <- "Th2_pre_post17h_030.csv"
#data_file_name <- "Th2_pre_post17h_040.csv"
#data_file_name <- "Th2_pre_post17h_060.csv"
#data_file_name <- "Th2_pre_post17h_070.csv"

## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         met = as.numeric(met),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))





######################################################################### 
#                 MAIN ANALYSIS - Th2 17h post exercise                 #
#########################################################################

### Meta-analysis from already computed effect sizes 
m4 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m4, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")


# Kostrzewa-Nowak et al. (2019) was detected as very different across meta-analysis and while not influential on the 
# overall effect, its values were unrealistic and brought substantial misleading variance 
fullData <-fullData %>% 
  dplyr::filter(!Study== "Kostrzewa-Nowak et al. (2019)") 
m4 <-metagen(g, seTE=SE_g, data=fullData, studlab=paste(Study), 
             random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
             sm="SMD", fixed=FALSE)

print(m4, overall=TRUE)

#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"


### Detect outliers
outlier<-metaoutliers(y=m4$TE, s2=m4$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# No outlier

## Funnel plot to visualize the results more easily
pdf("Th2_17hpost_funnel.pdf", width = 10, height = 18)
funnel(m4, random=TRUE, xlim=c(-4.5, 3), ylim=c(2, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.2403, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)

# No influential case

### Forest plot general
pdf("Th2_17hpost_no_out_forest.pdf", width = 10, height = 18)

forest(m4, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-5, 2))

dev.off()

####################################### ASYMMETRY ##########################################

### Egger's test (one tailed thus p<.1)
metabias(m4, method="linreg")
#No asymmetry thus no trim-and-fill


################################################################################################
#                      Moderator analysis - Th2 17h post ex                                  #
################################################################################################
#################################### Categorical moderators ####################################

### RISK of BIAS
## Moderator analysis 
m1mod <- update(m4, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m1mod, overall=FALSE)
# Not significant results  

## Forest plot 
pdf("subgroup_Th2_17hpost_RoB.pdf", width = 11.25, height = 18)

forest(m1mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()


### SEX
## Moderator analysis 
# All males

### Intensity
## Moderator analysis 
# Let's exclude any studies that do not report 'Intensity' to avoid counting them by default in the analysis.
fullData_bis<-fullData %>% 
  dplyr::filter(!intensity== "") 
m1_c <-metagen(g, seTE=SE_g, data=fullData_bis, studlab=paste(Study), 
               random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
               sm="SMD", fixed=FALSE)

m3mod <- update(m1_c, subgroup=intensity, print.subgroup.name=TRUE, fixed=FALSE)
print(m3mod, overall=FALSE)
# Not significant but 10 vs 1


#################################### Continuous moderators ####################################

### Age 
## We do directly the subgroup_analysis
m1mr <- metareg(m4, age)
summary(m1mr)
# Not significant results

## Print a bubble plot 
bubble(m1mr, lwd=2, col.line="blue", xlab="Age (years)", ylab="Hedges' g")


### BMI 
## We do directly the subgroup_analysis
m2mr <- metareg(m4, bmi)
summary(m2mr)
# Significant results

## Print a bubble plot 
bubble(m2mr, lwd=2, col.line="blue", xlab="BMI (kg/m2)", ylab="Hedges' g")

### Duration 
## We do directly the subgroup_analysis
m3mr <- metareg(m4, duration)
summary(m3mr)
# Not Significant 

## Print a bubble plot 
bubble(m3mr, lwd=2, col.line="blue", xlab="duration (min)", ylab="Hedges' g")
# Just 2 durations 

### Dose 
## We do directly the subgroup_analysis
m4mr <- metareg(m4, dose)
summary(m4mr)
# Not significant results

## Print a bubble plot 
bubble(m4mr, lwd=2, col.line="blue", xlab="DOSE", ylab="Hedges' g")


### IFN 
## We do directly the subgroup_analysis
m6mr <- metareg(m4, IFN)
summary(m6mr)
# Not significant results

## Print a bubble plot 
bubble(m6mr, lwd=2, col.line="blue", xlab="[IFN gamma]", ylab="Hedges' g")

### IL_4 
## We do directly the subgroup_analysis
m7mr <- metareg(m4, IL_4)
summary(m7mr)
# Not significant results

## Print a bubble plot 
bubble(m7mr, lwd=2, col.line="blue", xlab="[IL-4]", ylab="Hedges' g")


### IL_6 
## We do directly the subgroup_analysis
m8mr <- metareg(m4, IL_6)
summary(m8mr)
# Significant results

## Print a bubble plot 
bubble(m8mr, lwd=2, col.line="blue", xlab="[IL-6]", ylab="Hedges' g")

### IL_10 
## We do directly the subgroup_analysis
m9mr <- metareg(m4, IL_10)
summary(m9mr)
# Not significant results

## Print a bubble plot 
bubble(m9mr, lwd=2, col.line="blue", xlab="[IL-10]", ylab="Hedges' g")

### IL_2 
## We do directly the subgroup_analysis
m10mr <- metareg(m4, IL_2)
summary(m10mr)
# Not significant results

## Print a bubble plot 
bubble(m10mr, lwd=2, col.line="blue", xlab="[IL-2]", ylab="Hedges' g")


#################################### MULTIVARIATE ####################################

rma(yi = g, vi=V_g, mods= ~ intensity*duration, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results


########################################################################### 
#                          GATHER DATA Th1/Th2                            #
###########################################################################  

### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Ratio/pre_post")

### Charge the dataset
## What's the name of the data file to read?
data_file_name <- "Ratio_pre_post_r_050.csv"

### For the sensitivity analyses
#data_file_name <- "Ratio_pre_post_030.csv"
#data_file_name <- "Ratio_pre_post_040.csv"
#data_file_name <- "Ratio_pre_post_060.csv"
#data_file_name <- "Ratio_pre_post_070.csv"

## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         staining = as.character(staining),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))


######################################################################### 
#                    MAIN ANALYSIS - Th1/Th2                               #
#########################################################################

### Meta-analysis from already computed effect sizes 
m5 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m5, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")


# Kostrzewa-Nowak et al. (2019) was detected as very different across meta-analysis and while not influential on the 
# overall effect, its values were unrealistic and brought substantial misleading variance 
fullData <-fullData %>% 
  dplyr::filter(!Study== "Kostrzewa-Nowak et al. (2019)") 
m5 <-metagen(g, seTE=SE_g, data=fullData, studlab=paste(Study), 
             random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
             sm="SMD", fixed=FALSE)

print(m5, overall=TRUE)
#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"

### Detect outliers
outlier<-metaoutliers(y=m5$TE, s2=m5$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# No outlier 

## Funnel plot to visualize the results more easily
pdf("Th1_Th2_funnel.pdf", width = 10, height = 18)
funnel(m5, random=TRUE, xlim=c(-1.5, 2), ylim=c(0.6, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.3544, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.4, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)

# One value >1 and one <-1 and thus are considered as being an influential cases. 
# Nonetheless, those values are not outlier.
# Thus, we keep all the outcomes and try to find out why some of them have better results 
# than others in the moderator analysis. 

### Forest plot general
pdf("Th1_Th2_no_out_forest.pdf", width = 10, height = 18)

forest(m5, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-4, 3))

dev.off()

####################################### ASYMMETRY ##########################################

### Egger's test (one tailed thus p<.1)
metabias(m5, method="linreg")
#==> No Asymmetry, thus No trim and fill 

################################################################################################
#                            Moderator analysis - Th1/Th2                                      #
################################################################################################
#################################### Categorical moderators ####################################

### RISK of BIAS
## Moderator analysis 
m1mod <- update(m5, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m1mod, overall=FALSE)
# Not significant results  

## Forest plot 
pdf("subgroup_Th1_Th2_RoB.pdf", width = 11.25, height = 18)

forest(m1mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()


### SEX
## Moderator analysis 
m2mod <- update(m5, subgroup=sex, print.subgroup.name=TRUE, fixed=FALSE)
print(m2mod, overall=FALSE)
# Significant results but 15 vs 1 vs 2

## Forest plot 
pdf("subgroup_Th1_Th2_Sex.pdf", width = 11.25, height = 18)

forest(m2mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

### Intensity
## Moderator analysis 
# Let's exclude any studies that do not report 'Intensity' to avoid counting them by default in the analysis.
fullData_bis<-fullData %>% 
  dplyr::filter(!intensity== "") 
m1_c <-metagen(g, seTE=SE_g, data=fullData_bis, studlab=paste(Study), 
               random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
               sm="SMD", fixed=FALSE)

m3mod <- update(m1_c, subgroup=intensity, print.subgroup.name=TRUE, fixed=FALSE)
print(m3mod, overall=FALSE)
# Significant

## Forest plot 
pdf("subgroup_Th1_Th2_intensity.pdf", width = 11.25, height = 18)

forest(m3mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()

### staining
## Moderator analysis 
m3mod <- update(m5, subgroup=staining, print.subgroup.name=TRUE, fixed=FALSE)
print(m3mod, overall=FALSE)
# Not Significant results 

#################################### Continuous moderators ####################################

### Age 
## We do directly the subgroup_analysis
m1mr <- metareg(m5, age)
summary(m1mr)
# Not significant results

## Print a bubble plot 
bubble(m1mr, lwd=2, col.line="blue", xlab="Age (years)", ylab="Hedges' g")


### BMI 
## We do directly the subgroup_analysis
m2mr <- metareg(m5, bmi)
summary(m2mr)
# Not significant results

## Print a bubble plot 
bubble(m2mr, lwd=2, col.line="blue", xlab="BMI (kg/m2)", ylab="Hedges' g")

### Duration 
## We do directly the subgroup_analysis
m3mr <- metareg(m5, duration)
summary(m3mr)
# Significant results

## Print a bubble plot 
bubble(m3mr, lwd=2, col.line="blue", xlab="duration (min)", ylab="Hedges' g")

### Dose 
## We do directly the subgroup_analysis
m4mr <- metareg(m5, dose)
summary(m4mr)
# Not significant results

## Print a bubble plot 
bubble(m4mr, lwd=2, col.line="blue", xlab="DOSE", ylab="Hedges' g")


### IFN 
## We do directly the subgroup_analysis
m6mr <- metareg(m5, IFN)
summary(m6mr)
# Not significant results

## Print a bubble plot 
bubble(m6mr, lwd=2, col.line="blue", xlab="[IFN gamma]", ylab="Hedges' g")

### IL_4 
## We do directly the subgroup_analysis
m7mr <- metareg(m5, IL_4)
summary(m7mr)
# Not significant results

## Print a bubble plot 
bubble(m7mr, lwd=2, col.line="blue", xlab="[IL-4]", ylab="Hedges' g")


### IL_6 
## We do directly the subgroup_analysis
m8mr <- metareg(m5, IL_6)
summary(m8mr)
# Not significant results

## Print a bubble plot 
bubble(m8mr, lwd=2, col.line="blue", xlab="[IL-6]", ylab="Hedges' g")

### IL_10 
## We do directly the subgroup_analysis
m9mr <- metareg(m5, IL_10)
summary(m9mr)
# Not significant results

## Print a bubble plot 
bubble(m9mr, lwd=2, col.line="blue", xlab="[IL-10]", ylab="Hedges' g")

### IL_2 
## We do directly the subgroup_analysis
m10mr <- metareg(m5, IL_2)
summary(m10mr)
# Not significant results

## Print a bubble plot 
bubble(m10mr, lwd=2, col.line="blue", xlab="[IL-2]", ylab="Hedges' g")


#################################### MULTIVARIATE ####################################

rma(yi = g, vi=V_g, mods= ~ intensity*duration, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results

rma(yi = g, vi=V_g, mods= ~ IL_10*IL_4, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results

rma(yi = g, vi=V_g, mods= ~ IL_2*IFN, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Not significant results
########################################################################### 
#                 GATHER DATA Th1_Th2  17h post exercise                  #
###########################################################################  
### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Ratio/pre_post17h")

### Charge the dataset
## What's the name of the data file to read?
data_file_name <- "Ratio_pre_post17h_050.csv"

### For the sensitivity analyses
#data_file_name <- "Ratio_pre_post17h_030.csv"
#data_file_name <- "Ratio_pre_post17h_040.csv"
#data_file_name <- "Ratio_pre_post17h_060.csv"
#data_file_name <- "Ratio_pre_post17h_070.csv"

## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         met = as.numeric(met),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))




######################################################################### 
#                 MAIN ANALYSIS - Th1_Th2 17h post exercise             #
#########################################################################

### Meta-analysis from already computed effect sizes 
m6 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m6, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")


# Kostrzewa-Nowak et al. (2019) was detected as very different across meta-analysis and while not influential on the 
# overall effect, its values were unrealistic and brought substantial misleading variance 
fullData <-fullData %>% 
  dplyr::filter(!Study== "Kostrzewa-Nowak et al. (2019)") 
m6 <-metagen(g, seTE=SE_g, data=fullData, studlab=paste(Study), 
             random = TRUE, method.tau = "HE", hakn = TRUE, prediction=TRUE, 
             sm="SMD", fixed=FALSE)

print(m6, overall=TRUE)
#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"

### Detect outliers
outlier<-metaoutliers(y=m6$TE, s2=m6$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# No outlier  

## Funnel plot to visualize the results more easily
pdf("Th1_Th2_17hpost_funnel.pdf", width = 10, height = 18)
funnel(m6, random=TRUE, xlim=c(-3, 6), ylim=c(2.5, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.6433, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)

# No influential case 

### Forest plot general
pdf("Th1_Th2_17hpost_no_out_forest.pdf", width = 10, height = 18)

forest(m6, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-2, 5))

dev.off()

####################################### ASYMMETRY ##########################################

### Egger's test (one tailed thus p<.1)
metabias(m6, method="linreg")
#No asymmetry thus no trim-and-fill analysis 


################################################################################################
#                        Moderator analysis - Th1_Th2 17h post ex                              #
################################################################################################

#################################### Categorical moderators ####################################

### RISK of BIAS
## Moderator analysis 
m1mod <- update(m6, subgroup=RoB, print.subgroup.name=TRUE, fixed=FALSE)
print(m1mod, overall=FALSE)
# Not significant results  

## Forest plot 
pdf("subgroup_Th1_Th2_17hpost_RoB.pdf", width = 11.25, height = 18)

forest(m1mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()


### Intensity
## Moderator analysis 
m2mod <- update(m6, subgroup=intensity, print.subgroup.name=TRUE, fixed=FALSE)
print(m2mod, overall=FALSE)
# Significant results but 10 vs 1

## Forest plot 
pdf("subgroup_Th1_Th2_17hpost_intensity.pdf", width = 11.25, height = 18)

forest(m2mod, sortvar= -TE, layout = "JAMA", rightcols = TRUE, method.tau="HE",weight.study="random", fixed = FALSE,
       xlab = "Hedges' g", xlim=c(-2, 5), print.tau2 = TRUE, spacing= 1.2, ff.random = "bold", 
       colgap = "4mm", prediction=FALSE, subgroup= TRUE, test.subgroup.random = TRUE, col.by="darkblue", print.subgroup.labels=TRUE, 
       test.effect.subgroup.random = TRUE, subgroup.name=NULL, digits.I2=1)

dev.off()
#################################### Continuous moderators ####################################

### Age 
## We do directly the subgroup_analysis
m1mr <- metareg(m6, age)
summary(m1mr)
# Not significant 

## Print a bubble plot 
bubble(m1mr, lwd=2, col.line="blue", xlab="Age (years)", ylab="Hedges' g")


### BMI 
## We do directly the subgroup_analysis
m2mr <- metareg(m6, bmi)
summary(m2mr)
# Not significant results

## Print a bubble plot 
bubble(m2mr, lwd=2, col.line="blue", xlab="BMI (kg/m2)", ylab="Hedges' g")

### Duration 
## We do directly the subgroup_analysis
m3mr <- metareg(m6, duration)
summary(m3mr)
# Not significant results

## Print a bubble plot 
bubble(m3mr, lwd=2, col.line="blue", xlab="duration (min)", ylab="Hedges' g")

### Dose 
## We do directly the subgroup_analysis
m4mr <- metareg(m6, dose)
summary(m4mr)
# Not significant results

## Print a bubble plot 
bubble(m4mr, lwd=2, col.line="blue", xlab="DOSE", ylab="Hedges' g")

### MET 
## We do directly the subgroup_analysis
m5mr <- metareg(m6, met)
summary(m5mr)
# Not significant results

## Print a bubble plot 
bubble(m5mr, lwd=2, col.line="blue", xlab="MET", ylab="Hedges' g")

### IFN 
## We do directly the subgroup_analysis
m6mr <- metareg(m6, IFN)
summary(m6mr)
# Not significant results

## Print a bubble plot 
bubble(m6mr, lwd=2, col.line="blue", xlab="[IFN gamma]", ylab="Hedges' g")

### IL_4 
## We do directly the subgroup_analysis
m7mr <- metareg(m6, IL_4)
summary(m7mr)
# Not significant results

## Print a bubble plot 
bubble(m7mr, lwd=2, col.line="blue", xlab="[IL-4]", ylab="Hedges' g")


### IL_6 
## We do directly the subgroup_analysis
m8mr <- metareg(m6, IL_6)
summary(m8mr)
# Significant results

## Print a bubble plot 
bubble(m8mr, lwd=2, col.line="blue", xlab="[IL-6]", ylab="Hedges' g")

### IL_10 
## We do directly the subgroup_analysis
m9mr <- metareg(m6, IL_10)
summary(m9mr)
# Not significant results

## Print a bubble plot 
bubble(m9mr, lwd=2, col.line="blue", xlab="[IL-10]", ylab="Hedges' g")

### IL_2 
## We do directly the subgroup_analysis
m10mr <- metareg(m6, IL_2)
summary(m10mr)
# Not significant results

## Print a bubble plot 
bubble(m10mr, lwd=2, col.line="blue", xlab="[IL-2]", ylab="Hedges' g")


#################################### MULTIVARIATE ####################################

rma(yi = g, vi=V_g, mods= ~ intensity*duration, data=fullData, slab = paste(Study), method = "HE", test="knha")
# Model significant but not independent predictors


########################################################################### 
#                             GATHER DATA Th1                             #
###########################################################################  

### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Th1/pre_post+2h")


### Charge the dataset
### What's the name of the data file to read?
data_file_name <- "Th1_pre_post+2h_r_050.csv"

### For the sensitivity analyses
#data_file_name <- "Th1_pre_post_030.csv"
#data_file_name <- "Th1_pre_post_040.csv"
#data_file_name <- "Th1_pre_post_060.csv"
#data_file_name <- "Th1_pre_post_070.csv"


## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))




######################################################################### 
#                    MAIN ANALYSIS - Th1                               #
#########################################################################

### Meta-analysis from already computed effect sizes 
m1 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m1, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")



#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"

## Creation of functions 
metaoutliers <- function(y, s2, model){
  if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  w <- 1/s2
  y.p <- sum(y*w)/sum(w)
  n <- length(y)
  
  if(missing(model)){
    hetmeasure <- metahet.base(y, s2)
    Ir2 <- hetmeasure$Ir2
    if(Ir2 < 0.3){
      model <- "FE"
      cat("This function uses fixed-effect meta-analysis because Ir2 < 30%.\n")
    }else{
      model <- "RE"
      cat("This function uses random-effects meta-analysis because Ir2 >= 30%.\n")
    }
  }
  
  if(!is.element(model, c("FE", "RE"))) stop("wrong input for the argument model.")
  
  y.p.i <- res <- std.res <- numeric(n)
  if(model == "FE"){
    for(i in 1:n){
      w.temp <- w[-i]
      y.temp <- y[-i]
      y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1/sum(w.temp) + s2[i]
      std.res[i] <- res[i]/sqrt(var.res.i)
    }
  }else{
    for(i in 1:n){
      s2.temp <- s2[-i]
      y.temp <- y[-i]
      tau2.temp <- metahet.base(y.temp, s2.temp)$tau2.DL
      w.temp <- 1/(s2.temp + tau2.temp)
      y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1/sum(w.temp) + s2[i] + tau2.temp
      std.res[i] <- res[i]/sqrt(var.res.i)
    }
  }
  
  outliers <- which(abs(std.res) >= 3)
  if(length(outliers) == 0) outliers <- "All the standardized residuals are smaller than 3"
  
  out <- NULL
  out$model <- model
  out$std.res <- std.res
  out$outliers <- outliers
  
  class(out) <- "metaoutliers"
  return(out)
}

metahet.base <- function(y, s2){
  if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  n <- length(y)
  
  w <- 1/s2
  mu.bar <- sum(w*y)/sum(w)
  
  out <- NULL
  
  out$weighted.mean <- mu.bar
  
  # the conventional methods
  Q <- sum(w*(y - mu.bar)^2)
  H <- sqrt(Q/(n - 1))
  I2 <- (Q - n + 1)/Q
  tau2.DL <- (Q - n + 1)/(sum(w) - sum(w^2)/sum(w))
  tau2.DL <- max(c(0, tau2.DL))
  out$Q <- Q
  out$H <- H
  out$I2 <- I2
  out$tau2.DL <- tau2.DL
  
  # absolute deviation based on weighted mean
  Qr <- sum(sqrt(w)*abs(y - mu.bar))
  Hr <- sqrt((3.14159*Qr^2)/(2*n*(n - 1)))
  Ir2 <- (Qr^2 - 2*n*(n - 1)/3.14159)/(Qr^2)
  tau2.r <- tau2.r.solver(w, Qr)
  out$Qr <- Qr
  out$Hr <- Hr
  out$Ir2 <- Ir2
  out$tau2.r <- tau2.r
  
  # absolute deviation based on weighted median
  
  expit <- function(x) {ifelse(x >= 0, 1/(1 + exp(-x/0.0001)), exp(x/0.0001)/(1 + exp(x/0.0001)))}
  psi <- function(x) {sum(w*(expit(x - y) - 0.5))}
  mu.med <- uniroot(psi, c(min(y) - 0.001, max(y) + 0.001))$root
  out$weighted.median <- mu.med
  Qm <- sum(sqrt(w)*abs(y - mu.med))
  Hm <- sqrt(3.14159/2)*Qm/n
  Im2 <- (Qm^2 - 2*n^2/3.14159)/Qm^2
  tau2.m <- tau2.m.solver(w, Qm)
  out$Qm <- Qm
  out$Hm <- Hm
  out$Im2 <- Im2
  out$tau2.m <- tau2.m
  
  return(out)
}

tau2.r.solver <- function(w, Qr){
  f <- function(tau2){
    out <- sum(sqrt(1 - w/sum(w) + tau2*(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))) - Qr*sqrt(3.14159/2)
    return(out)
  }
  f <- Vectorize(f)
  tau.upp <- Qr*sqrt(3.14159/2)/sum(sqrt(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))
  tau2.upp <- tau.upp^2
  f.low <- f(0)
  f.upp <- f(tau2.upp)
  if(f.low*f.upp > 0){
    tau2.r <- 0
  }else{
    tau2.r <- uniroot(f, interval = c(0, tau2.upp))
    tau2.r <- tau2.r$root
  }
  return(tau2.r)
}

tau2.m.solver <- function(w, Qm){
  f <- function(tau2){
    out <- sum(sqrt(1 + w*tau2)) - Qm*sqrt(3.14159/2)
    return(out)
  }
  f <- Vectorize(f)
  n <- length(w)
  tau2.upp <- sum(1/w)*(Qm^2/n*2/3.14159 - 1)
  tau2.upp <- max(c(tau2.upp, 0.01))
  f.low <- f(0)
  f.upp <- f(tau2.upp)
  if(f.low*f.upp > 0){
    tau2.m <- 0
  }else{
    tau2.m <- uniroot(f, interval = c(0, tau2.upp))
    tau2.m <- tau2.m$root
  }
  return(tau2.m)
}


### Detect outliers
outlier<-metaoutliers(y=m1$TE, s2=m1$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# No outcome have standardized residuals above 3 or below -3. 

## Funnel plot to visualize the results more easily
pdf("Th1_funnel.pdf", width = 10, height = 18)
funnel(m1, random=TRUE, xlim=c(-1.5, 4), ylim=c(0.7, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.99, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)

# One value is <-1 and one is >1 and thus are considered as being influential cases. 
# Nonetheless, those values are not outliers (Kostrzewa-Nowak & Nowak (2020a)a and Tsubakihara et al. (2012)).
# Accordingly, we retain all outcomes and aim to identify the sources of heterogeneity in the moderator analysis.

### Forest plot general
pdf("Th1_no_out_forest.pdf", width = 10, height = 18)

forest(m1, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-2, 5))

dev.off()

####################################### ASYMMETRY ##########################################

###Number of studies (k=5) too small to test for small study effects (k.min=10)


########################################################################### 
#                             GATHER DATA Th2                             #
###########################################################################  
### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Th2/pre_post+2h")

### Charge the dataset
## What's the name of the data file to read?
data_file_name <- "Th2_pre_post+2h_r_050.csv"

### For the sensitivity analyses
#data_file_name <- "Th2_pre_post_030.csv"
#data_file_name <- "Th2_pre_post_040.csv"
#data_file_name <- "Th2_pre_post_060.csv"
#data_file_name <- "Th2_pre_post_070.csv"

## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))

######################################################################### 
#                    MAIN ANALYSIS - Th2                               #
#########################################################################

### Meta-analysis from already computed effect sizes 
m3 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m3, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")


#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"

### Detect outliers
outlier<-metaoutliers(y=m3$TE, s2=m3$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# One outcome has standardized residuals around -2 (Kostrzewa-Nowak et al. (2019)) . 

## Funnel plot to visualize the results more easily
pdf("Th2_funnel.pdf", width = 10, height = 18)
funnel(m3, random=TRUE, xlim=c(-1, 2), ylim=c(0.7, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.3214, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.8, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)
# No influential case 

### Forest plot general
pdf("Th2_no_out_forest.pdf", width = 10, height = 18)

forest(m3, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-4, 3))

dev.off()

####################################### ASYMMETRY ##########################################


########################################################################### 
#                          GATHER DATA Th1/Th2                            #
###########################################################################  

### Set the working directory
setwd("~/Projekte/Meta Th cells/Data analysis/meta_Th/Ratio/pre_post+2h")

### Charge the dataset
## What's the name of the data file to read?
data_file_name <- "Ratio_pre_post+2h_r_050.csv"

### For the sensitivity analyses
#data_file_name <- "Ratio_pre_post_030.csv"
#data_file_name <- "Ratio_pre_post_040.csv"
#data_file_name <- "Ratio_pre_post_060.csv"
#data_file_name <- "Ratio_pre_post_070.csv"

## Read in csv file of data input
fullData <- read.csv(data_file_name, header=TRUE, sep=";", stringsAsFactors = FALSE)

## if using R 4.0 or later:
fullData <- as.data.frame(fullData) %>% 
  mutate(Study = as.character(Study),
         g = as.numeric(g),
         V_g = as.numeric(V_g), 
         SE_g = as.numeric(SE_g),
         RoB = as.character(RoB),
         intensity = as.character(intensity),
         sex = as.character(sex),
         age = as.numeric(age),
         duration = as.numeric(duration),
         dose = as.numeric(dose),
         IFN = as.numeric(IFN),
         IL_10 = as.numeric(IL_10),
         IL_2 = as.numeric(IL_2),
         IL_4 = as.numeric(IL_4),
         IL_6 = as.numeric(IL_6),
         bmi = as.numeric (bmi))


######################################################################### 
#                    MAIN ANALYSIS - Th1/Th2                               #
#########################################################################

### Meta-analysis from already computed effect sizes 
m5 <-metagen(TE=g, seTE=SE_g,data=fullData, studlab=paste(Study), random = TRUE,
             method.tau = "HE", fixed=FALSE, hakn = TRUE, prediction=TRUE, sm="SMD")

## to see the results shortened 
print(m5, overall=TRUE)

## double-check with the metafor package 
rma(yi = g, sei=SE_g, data=fullData, slab = paste(Study), method = "HE", test="knha")



#################################### OUTLIERS DETECTION ####################################                        
### Rmk
# Calculate the standardized residual for each study in meta-analysis using 
# the methods described in Chapter 12 in Hedges and Olkin (1985) and Viechtbauer 
# and Cheung (2010). A study is considered as an outlier if its standardized 
# residual is greater than 3 in absolute magnitude.If >1.96, it may call for closer inspection.

# WARNING this command requires the creation of the following functions: 
# "metaoutliers", "metahet.base", "tau2.r.solver","tau2.m.solver"

### Detect outliers
outlier<-metaoutliers(y=m5$TE, s2=m5$seTE^2, model="RE")
fullData %>% select (Study, g) %>% mutate(std_res=outlier$std.res) 
stdres<- outlier$std.res 
## To see the results
## Plotted
plot(stdres, type="o", pch=19, xlab="Outcome ID", ylab="standardized residual")

# No outlier 

## Funnel plot to visualize the results more easily
pdf("Th1_Th2_funnel.pdf", width = 10, height = 18)
funnel(m5, random=TRUE, xlim=c(-1.5, 2), ylim=c(0.6, 0),
       xlab="Hedges' g", ylab="Standard error",level =0.95, contour = c(0.9, 0.95, 0.99), 
       shade =c("white", "gray", "darkgray"), col="black", bg="darkgray",
       pch=16,ref=0.3544, lwd=1, cex=1.2) 

legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "p < 0.01"),
       fill=c("white", "gray", "darkgray"), bty="n",
       xjust = -0.4, yjust = 0.60, x.intersp = 1, y.intersp = 1)
dev.off()


## Test for outliers: dfbetas
overall_metafor<-rma.mv(g, V_g, tdist = TRUE, data = fullData)
df_betas<-dfbetas(overall_metafor)
## To see the results
## Printed
fullData %>% select (Study, g) %>% mutate(dfbetas=df_betas)

# One value >1 and one <-1 and thus are considered as being an influential cases. 
# Nonetheless, those values are not outlier.
# Thus, we keep all the outcomes and try to find out why some of them have better results 
# than others in the moderator analysis. 

### Forest plot general
pdf("Th1_Th2_no_out_forest.pdf", width = 10, height = 18)

forest(m5, sortvar= -TE, layout = "JAMA", method.tau="HE",weight.study="random", comb.fixed = FALSE, 
       xlab = "Hedges' g", test.overall.random = TRUE, lwd=1,boxsize=0.8, colgap=unit(4, "mm"),
       digits.mean=2, digits.sd=2, digits.se = 2, digits.tau2 = 2, digits.I2=1, fontsize=6, xlim=c(-4, 3))

dev.off()

####################################### ASYMMETRY ##########################################



