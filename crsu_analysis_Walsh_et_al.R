########################################################################
## Walsh et al. analysis 
## Code author: Enzo Cerullo
## Date: 19/10/2020
########################################################################

## set wd 
setwd("/media/enzo/A05C06A35C0673F6/Users/Enzo/Documents/CRSU")

## packages 
require(gridExtra)
require(egg)
require(plyr)
require(dplyr)
require(readxl)
require(lme4)
require(mada)
require(magic)
require(ggforce)
require(patchwork)
require(rstanarm)
require(bayestestR)
require(lmtest)
require(cmdstanr) # see https://mc-stan.org/users/interfaces/cmdstan fior tutorials on how to install cmdstanr
require(loo)
require(pbmcapply)
require(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

options(scipen = 999)
options(max.print = 1000000000)

filter <- dplyr::filter
mutate <- dplyr::mutate

########################################################################
## data for 2 and 3 way 
data <- read.csv("AllComparativeRevManExport.csv")

# estimates from 3 way studies only
three_way <- filter(data,  Pair.or.triple == "Three way FVR")
write.csv(three_way, "three_way.csv")  # write to csv and read to reset the factor levels (probably a more elegant way to do this)
three_way <- read.csv( "three_way.csv")

#################################
##### select which dataset to use
X <- data.frame(data) ## for three way analysis(for Hoyer & Kuss et al / "cochrane" models)
#X <- data.frame(three_way) ## for full dataset of 5 tests (for NMA Nyaga et al)
#################################
N <- length(X$TP) # num. of studies
X$n1 <- X$TP+X$FN
X$n0 <- X$FP+X$TN
X$true1 <- X$TP
X$true0 <- X$TN 
X$study <- 1:N

X2 <- X %>% arrange(Study.ID)
X2

# Reshape the data from wide to long format 
Y = reshape(X, direction = "long", varying = list( c("n1" , "n0") , c( "true1","true0" ) ) ,
            timevar = "sens" , times = c(1,0) , v.names = c("n","true") ) 
# Sort data by study to cluster the 2 records per study together 
Y = Y[order(Y$id),]
Y$spec<- 1-Y$sens

# dummy variables to add covariates for test type 
# make indicator variables for each test
# this is for use with glmer
Y2 <- mutate(Y, Flu =   ifelse(Test == "Fluoro", 1,0),
                Vis =   ifelse(Test == "Visual", 1,0), 
                Radio = ifelse(Test == "Radio", 1,0),
                Trans = ifelse(Test == "Trans", 1,0), 
                ECM   = ifelse(Test == "ECM",   1,0), 
                seFlu =   Flu*sens, 
                seVis =   Vis*sens, 
                seRadio = Radio*sens,
                seTrans = Trans*sens,
                seECM   = ECM*sens, 
                spFlu =   Flu*spec, 
                spVis =   Vis*spec, 
                spRadio = Radio*spec,
                spTrans = Trans*spec, 
                spECM   = ECM*spec) %>% dplyr::arrange(Study.ID)

Y2

# can also run models with glmer from rstanarm package (at least for non-NMA)
# I didn't use this for the analysis but showing for completeness
mod <- stan_glmer(cbind(true, n-true) ~ 0    
                   #+ sens + spec 
                   + seECM + seFlu +  seRadio +  seTrans + seVis
                   + spECM + spFlu +  spRadio +  spTrans + spVis
                   + (0 + sens + spec | Study.ID:Test), # M3: same variance and corr. for all tests
                   #    + (0 + seECM + spECM  | Study.ID) + (0 + seFlu + spFlu| Study.ID) +  (0 + seRadio + spRadio| Study.ID)  + (0 + seTrans + spTrans| Study.ID) +  (0 + seVis + spVis| Study.ID) ,
                   # + (0 + seECM + spECM + seFlu + spFlu + seRadio + spRadio + seTrans + spTrans + seVis + spVis| Study.ID) ,
                   prior = normal( 0, 2),  
                   warmup =  500, iter = 2000,  adapt_delta = 0.80,
                   data = Y2 , family = binomial)

summary(mod, pars = c("seRadio", "seFlu", "seVis"), digits=5)



X2_threeway <- filter(X2,  Pair.or.triple == "Three way FVR")
#############################################
###################### run models

#### NMA (5 TESTS)
file <- file.path(file = "gs_nma_m2.stan") # nma
mod <- cmdstan_model(file)

# run w/ cmdstanr 
meta_model <- mod$sample(
  data=  list(NS = length(unique(X2$Study.ID)), 
            N = nrow(X2), 
            TP = X2$TP, FN = X2$FN, TN = X2$TN, FP = X2$FP, 
            pos = X2$TP + X2$FN,
            neg = X2$TN + X2$FP,
            Test = as.numeric(as.factor(X2$Test)),
            Study = as.numeric(as.factor(X2$Study.ID))   , 
            total_tests = 5,
            holdout = rep(0,   nrow(X2) )) ,  
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000, 
  adapt_delta = 0.80, 
  max_treedepth = 10,
  refresh = 50,
)

meta_modelr <- rstan::read_stan_csv(meta_model$output_files())

#saveRDS(meta_modelr, file = "nma_m2_fit.rds")



#### DIRECT COMPARISONS (3 TESTS)
file <- file.path(file = "threeway_m3.stan") #  for three-way model
mod <- cmdstan_model(file)

# run w/ cmdstanr 
meta_model <- mod$sample(
  data=  list(NS = length(unique(X2_threeway$Study.ID)), 
              N = nrow(X2_threeway), 
              TP = X2_threeway$TP, FN = X2_threeway$FN, TN = X2_threeway$TN, FP = X2_threeway$FP, 
              pos = X2_threeway$TP + X2_threeway$FN,
              neg = X2_threeway$TN + X2_threeway$FP,
              Test = as.numeric(as.factor(X2_threeway$Test)),
              Study = as.numeric(as.factor(X2_threeway$Study.ID))   , 
              total_tests = 3,
              holdout = rep(0,   nrow(X2_threeway) )) ,  
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000, 
  adapt_delta = 0.80, 
  max_treedepth = 10,
  refresh = 50,
)

# 1 = ecm, 2 = flu, 3 = radio, 4 = trans, 5 = visual

meta_modelr <- rstan::read_stan_csv(meta_model$output_files())

#saveRDS(meta_modelr, file = "threeway_m3_fit.rds")

print(meta_modelr, pars= c( "mu", "Se", "Sp", "Se_pred", "Sp_pred", "diff_Se", "diff_Sp", "ratio_Se", "ratio_Sp", "S",
                            "tau", "Omega", "Sigma", 
                            "sigmasq", "rho", "rho12", "pred"),probs = c(0.025,  0.5, 0.975))


print(meta_modelr, pars= c( "Se", "Sp", "diff_Se", "diff_Sp", "ratio_Se",
                            "ratio_Sp",
                             "Omega", "Sigma"),probs = c(0.025,  0.5, 0.975))

#print(meta_modelr, pars= c( "log_lik2"),probs = c(0.025,  0.5, 0.975))


# 1 = ECM, 2 = FLU, 3 = radio, 4 = trans, 5 = visual

###############
##  output table
out <- summary(meta_modelr,  pars= c( "Se", "Sp", "diff_Se", "diff_Sp", "ratio_Se", "ratio_Sp", "S",
                                      "tau", "Omega", "Sigma") ,probs = c(0.025,  0.5, 0.975))



out_data <- data.frame( Quantile = round((out$summary[,4:6]), 2))
out_data


#write.csv(out_data, file = "table_m3.csv")

View(out_data)
              
################
## density plots
stan_dens(meta_modelr, pars= c( "Se", "Sp", "tau", "Omega"))
################
## loo
loglik <- extract_log_lik(meta_modelr, parameter_name = "log_lik", merge_chains = FALSE)
str(loglik)
length <- length(loglik[1,1,][loglik[1,1,] != "NaN"])
loglik2 <- array(data = loglik[loglik != "NaN"], dim = c(4, 1500, length ))
str(loglik2)
  
r_eff <- relative_eff(exp(loglik2), cores = 2)

## save loo
mod_loo <- cmdstanr::loo(loglik2, r_eff = r_eff, cores = 2)

mod_loo <- rstan::loo(x = meta_modelr,  r_eff = r_eff, cores = 2, moment_match = TRUE)

## output loo
print(mod_loo)


## loo comparison
loo_comparisons <- loo_compare(list(m_loo, ...))
loo_comparisons



#######################################################
##### K fold CV as loo diagnostic is not adequate #####
set.seed(123)

# put data in tibble format
d <-    tibble( TP = X2$TP, FN = X2$FN, TN = X2$TN, FP = X2$FP, 
                pos = X2$TP + X2$FN,
                neg = X2$TN + X2$FP,
                Test = as.numeric(X2$Test),
                Study = as.numeric(X2$Study.ID),
                Study.ID = X2$Study.ID)
d

require(groupdata2)


# Assuming K = 10, the procedure is the following:
# 1) We extract 1/10 of the data and save it in the list G with k = 1;
# 2) We extract 1/9 of the remaining data and save it with k = 2;
# 3) We extract 1/8 of the remaining data and save it with k = 3;
# 4) ...;
# 10) We extract all the data the data and save it with k = 10

K <- 10
dK <- fold(d, k = K, 
           id_col = "Study.ID")
dK

# Create list containing the K datasets
ldata <- plyr::llply(1:K, function(i) {
  list( N = nrow(dK),
        TP = dK$TP, FN = dK$FN, TN = dK$TN, FP = dK$FP, 
        pos = dK$TP + dK$FN,
        neg = dK$TN + dK$FP,
        Test = as.numeric(dK$Test),
        Study = as.numeric(dK$Study),
       NS = length(unique(dK$Study)), 
       total_tests = 5,
       holdout = ifelse(dK$.folds == i, 1, 0))
})

ldata


###################################################################
#### functions to run k-fold cross validation, found on GitHub from:
#### https://github.com/stan-dev/stancon_talks/blob/master/2017/Contributed-Talks/07_nicenboim/kfold.Rmd
###################################################################

# The following function can run all the chains of all the folds of the model in parallel:
stan_kfold <- function(file, list_of_datas, chains, cores,...){
  badRhat <- 1.1
  K <- length(list_of_datas)
  model <- stan_model(file=file)
  # First parallelize all chains:
  sflist <- 
    pbmclapply(1:(K*chains), mc.cores = cores, 
               function(i){
                 # Fold number:
                 k <- round((i+1) / chains)
                 s <- sampling(model, data = list_of_datas[[k]], 
                               chains = 1, chain_id = i,  ...)
                 return(s)
               })
  # Then merge the K * chains to create K stanfits:
  stanfit <- list()
  for(k in 1:K){
    inchains <- (chains*k - 2):(chains*k)
    # Merge `chains` of each fold
    stanfit[[k]] <- sflist2stanfit(sflist[inchains])
  }  
  return(stanfit) 
}

# Wrapper function to extract the log_lik of the held-out data, given a list of stanfits, and a list which indicates with 1 and 0 whether the observation was held out or not:
extract_log_lik_K <- function(list_of_stanfits, list_of_holdout, ...){
  K <- length(list_of_stanfits)
  list_of_log_liks <- plyr::llply(1:K, function(k){
    extract_log_lik(list_of_stanfits[[k]], merge_chains = TRUE , ...)
  })
  # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
  # We define `log_lik_heldout` as a (samples x N_obs) matrix
  # (similar to each log_lik matrix)
  log_lik_heldout <- list_of_log_liks[[1]] * NA
  for(k in 1:K){
    log_lik <- list_of_log_liks[[k]]
    samples <- dim(log_lik)[1] 
    N_obs <- dim(log_lik)[2]
    # This is a matrix with the same size as log_lik_heldout
    # with 1 if the data was held out in the fold k
    heldout <- matrix(rep(list_of_holdout[[k]], each = samples), nrow = samples)
    # Sanity check that the previous log_lik is not being overwritten:
    if(any(!is.na(log_lik_heldout[heldout==1]))){
      warning("Heldout log_lik has been overwritten!!!!")
    }
    # We save here the log_lik of the fold k in the matrix:
    log_lik_heldout[heldout==1] <- log_lik[heldout==1]
  }
  return(log_lik_heldout)
}


kfold <- function(log_lik_heldout)  {
  library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    colLogSumExps(x) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise)
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))

  out <- list(
    pointwise = pointwise,
    elpd_kfold = elpd_kfold,
    se_elpd_kfold = se_elpd_kfold)
  out
 #structure(out, class = "loo")
}


##############################
##### run k-fold #############

# We run all the chains of all the folds of the activation-based model in parallel:
fits <- stan_kfold("gs_nma_m1.stan", list_of_datas = ldata, 
                   chains = 4, cores = 4, seed = 123, iter = 500,
                   warmup=200, control=list(adapt_delta=0.80, 
                                            max_treedepth = 10))

holdout <- lapply(ldata, '[[', "holdout")
# We extract all the held_out log_lik of all the folds
log_lik_ab <- extract_log_lik_K(fits, holdout, "log_lik")
save(log_lik_ab, file = "log_lik_ab.Rda")
str(log_lik_ab)

length <- length(log_lik_ab[1,][log_lik_ab[1,] != "NaN" & !is.na(log_lik_ab[1,]) ])
loglik2 <- array(data = log_lik_ab[log_lik_ab != "NaN" & !is.na(log_lik_ab) ], dim = c(750, length ))
str(loglik2)
sum(is.na(loglik2))

.rs.restartR()

########################################################
## compute elpd_kfold_ic  and model comparison 
# save results - e.g. "kfold_nma_m1" corresponds to kfold from model 1, etc
kfold_nma_m1 <- kfold(loglik2)
kfold_nma_m1
save(kfold_nma_m1, file = "kfold_nma_m1.Rda")

kfold_nma_m1$elpd_kfold*(-2) 
kfold_nma_m2$elpd_kfold*(-2) 

## computer differences in elpd_kfold and se of differences 
kfold_nma_m2$elpd_kfold - kfold_nma_m1$elpd_kfold
sqrt(ncol(loglik2) * var(kfold_nma_m2$pointwise - kfold_nma_m1$pointwise))



# save results - e.g. "kfold_m1" corresponds to kfold from model 1, etc
#kfold_m1 <- kfold(loglik2)
#kfold_m1
#save(kfold_m1, file = "kfold_m1.R")

kfold_m1$elpd_kfold*(-2) 
kfold_m2$elpd_kfold*(-2) 
kfold_m3$elpd_kfold*(-2) 

## computer differences in elpd_kfold and se of differences 
kfold_m3$elpd_kfold - kfold_m1$elpd_kfold
sqrt(ncol(loglik2) * var(kfold_m3$pointwise - kfold_m1$pointwise))

kfold_m3$elpd_kfold - kfold_m2$elpd_kfold
sqrt(ncol(loglik2) * var(kfold_m3$pointwise - kfold_m2$pointwise))


#######################################################################
#######################################################################
## PLOTS
#######################################################################
#######################################################################

# function to generate confidence and prediction regions, hsroc curves 
# from the summary estimates and var-cov parameters extracted from the Stan output data

##############################
####### ROC plot
########################################
# load in the model

# nma (5 tests)

mod <- readRDS("nma_m2_fit.rds")
num_tests <- 5

# direct comparisons (3 tests)
mod <- readRDS("threeway_m3_fit.rds")
num_tests <- 3


###
print(mod, pars= c("Se", "Sp", "lSe" ,"Se_pred", "Sp_pred"),
      probs = c(0.025,0.5, 0.975))


require(bayesSurv)
require(scales)

## credible region
cred_1 <- list()

for (t in 1:num_tests) { 
  cred_1[[t]] <- tibble(y = (rstan::extract(mod, pars = "lSe")$lSe[,t]) , x = (rstan::extract(mod, pars = "lSp")$lSp[,t]))
}


require(data.table)
cred <- rbindlist(cred_1, idcol = TRUE)

# 1 = ECM, 2 = FLU, 3 = radio, 4 = trans, 5 = visual

cred2 <- dplyr::mutate(cred,  Test = factor(cred$.id,   label  = factor(c("ECM", 
                                                                      "Fluorescence",
                                                                      "Imaging", 
                                                                      "Transillumination OCT",
                                                                      "Visual classification"), 
                                                            levels = c("ECM", 
                                                                       "Fluorescence",
                                                                       "Imaging", 
                                                                       "Transillumination OCT",
                                                                       "Visual classification"))) )

cred2 <- dplyr::mutate(cred,  Test = factor(cred$.id,   label  = factor(c("Fluorescence",
                                                                          "Imaging", 
                                                                          "Visual classification"), 
                                                                        levels = c("Fluorescence",
                                                                                   "Imaging", 
                                                                                   "Visual classification"))))


# in inv_probit space
g <- ggplot(data = cred2, aes(x = x, y = y, colour = Test))  + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)
el = pb$data[[1]][c("x","y", "group")]


credible_region <- tibble(x = plogis(el$x), y = plogis(el$y), Test = factor(el$group,    label  = factor(c("ECM", 
                                                                                                         "Fluorescence",
                                                                                                         "Imaging", 
                                                                                                         "Transillumination OCT",
                                                                                                         "Visual classification"), 
                                                                                                       levels = c("ECM", 
                                                                                                                  "Fluorescence",
                                                                                                                  "Imaging", 
                                                                                                                  "Transillumination OCT",
                                                                                                                  "Visual classification"))) )

credible_region <- tibble(x = plogis(el$x), y = plogis(el$y), Test = factor(el$group,    label  = factor(c("Fluorescence",
                                                                                                           "Imaging", 
                                                                                                           "Visual classification"), 
                                                                                                         levels = c("Fluorescence",
                                                                                                                    "Imaging", 
                                                                                                                    "Visual classification"))))

credible_region

g <- ggplot(data = credible_region, aes(x = x, y = y, colour = Test))  + 
  geom_polygon(data = credible_region, aes(x = 1  - x, y = y, colour = Test), alpha=0.05, size=0.4)  + 
  xlim(0,1) + 
  ylim(0,1)

g

####
## prediction region

pred_1 <- list()

for (t in 1:num_tests) { 
  pred_1[[t]] <- tibble(y = (rstan::extract(mod, pars = "lSe_pred")$lSe_pred[,t]), x = (rstan::extract(mod, pars = "lSp_pred")$lSp_pred[,t]))
} 

require(data.table)
pred <- rbindlist(pred_1, idcol = TRUE)

pred2 <- mutate(pred,  Test = factor(.id,  label  = factor(c("ECM", 
                                                             "Fluorescence",
                                                             "Imaging", 
                                                             "Transillumination OCT",
                                                             "Visual classification"), 
                                                           levels = c("ECM", 
                                                                      "Fluorescence",
                                                                      "Imaging", 
                                                                      "Transillumination OCT",
                                                                      "Visual classification"))) )

pred2 <- mutate(pred,  Test = factor(.id,     label  = factor(c("Fluorescence",
                                                                "Imaging", 
                                                                "Visual classification"), 
                                                              levels = c("Fluorescence",
                                                                         "Imaging", 
                                                                         "Visual classification"))))


# in inv_probit space
g <- ggplot(data = pred2, aes(x = x, y = y, colour = Test))  + 
  # geom_point() + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)

el = pb$data[[1]][c("x","y", "group")]


pred_region <- tibble(x = plogis(el$x), y = plogis(el$y), Test = factor(el$group, label  = factor(c("ECM", 
                                                                                                  "Fluorescence",
                                                                                                  "Imaging", 
                                                                                                  "Transillumination OCT",
                                                                                                  "Visual classification"), 
                                                                                                levels = c("ECM", 
                                                                                                           "Fluorescence",
                                                                                                           "Imaging", 
                                                                                                           "Transillumination OCT",
                                                                                                           "Visual classification"))) )

pred_region <- tibble(x = plogis(el$x), y = plogis(el$y), Test = factor(el$group,    label  = factor(c("Fluorescence",
                                                                                                       "Imaging", 
                                                                                                       "Visual classification"), 
                                                                                                     levels = c("Fluorescence",
                                                                                                                "Imaging", 
                                                                                                                "Visual classification"))))
pred_region


g <- ggplot(data = pred_region, aes(x = x, y = y, colour = Test))  + 
  geom_polygon(data = pred_region, aes(x = 1  - x, y = y, colour = Test), alpha=0.05, size=0.4)  + 
  xlim(0,1) + 
  ylim(0,1)

g



## medians
median_sens <- (summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("Se"))$summary[,5])
median_spec <- (summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("Sp"))$summary[,5])

medians <- tibble(median_sens = median_sens, median_spec = median_spec, Test = factor( c(1:num_tests), label  = factor(c("ECM", 
                                                                                                                 "Fluorescence",
                                                                                                                 "Imaging", 
                                                                                                                 "Transillumination OCT",
                                                                                                                 "Visual classification"), 
                                                                                                               levels = c("ECM", 
                                                                                                                          "Fluorescence",
                                                                                                                          "Imaging", 
                                                                                                                          "Transillumination OCT",
                                                                                                                          "Visual classification"))) )

medians <- tibble(median_sens = median_sens, median_spec = median_spec, Test = factor( c(1:num_tests),    label  = factor(c("Fluorescence",
                                                                                                                            "Imaging", 
                                                                                                                            "Visual classification"), 
                                                                                                                          levels = c("Fluorescence",
                                                                                                                                     "Imaging", 
                                                                                                                                     "Visual classification"))))



print(mod, pars= c("Se", "Sp", "lSe" ,"Se_pred", "Sp_pred", "Sigma_bs"),probs = c(0.025,0.5, 0.975))


#############################
## plot


ss<- tibble( 
  Study =as.numeric(as.factor(X2$Study.ID)), 
  TP=X2$TP, FN=X2$FN, FP=X2$FP, TN=X2$TN,
  N=(X2$TP+X2$FN+X2$FP+X2$TN) ,
  Sensitivity= (TP/(TP+FN))    , 
  Specificity= (TN/(TN+FP))  ,
  Test = factor( X2$Test, label  = factor(c("ECM", 
                                                         "Fluorescence",
                                                         "Imaging", 
                                                         "Transillumination OCT",
                                                         "Visual classification"), 
                                                       levels = c("ECM", 
                                                                  "Fluorescence",
                                                                  "Imaging", 
                                                                  "Transillumination OCT",
                                                                  "Visual classification"))), 
  `Two-way or three-way study` = factor(X2$Pair.or.triple, label = c("Three way", "Two way")),
)

ss<- tibble( 
  Study =as.numeric(as.factor(three_way$Study.ID)), 
  TP=three_way$TP, FN=three_way$FN, FP=three_way$FP, TN=three_way$TN,
  N=(three_way$TP+three_way$FN+three_way$FP+three_way$TN) ,
  Sensitivity= (TP/(TP+FN))    , 
  Specificity= (TN/(TN+FP))  ,
  Test = factor( three_way$Test,    label  = factor(c("Fluorescence",
                                               "Imaging", 
                                               "Visual classification"), 
                                             levels = c("Fluorescence",
                                                        "Imaging", 
                                                        "Visual classification"))))
)




ss


print(mod, pars= c("Se", "Sp", "lSe" ,"Se_pred", "Sp_pred", "Sigma_bs", "sigma"),probs = c(0.025,0.5, 0.975))


##########################################################################
# Plot for all 5 tests
##########################################################################
require(RColorBrewer)
g <- ggplot(data = ss, aes(y=Sensitivity, x = 1-Specificity, colour = Test)) + 
  # geom_line(data = roc2,aes(FPR,TPR, colour = Test), size = 1, alpha = 0.9,inherit.aes = F) + 
  geom_point(data = ss, aes(y=Sensitivity, x = 1-Specificity, colour = Test, shape = `Two-way or three-way study`), size = 3, alpha=0.7,inherit.aes = F) + 
  geom_point(data = medians,  aes(y=median_sens, x = 1 - median_spec, colour = Test), size=10, shape = 5,inherit.aes = F)  +      # summary points
  geom_path(data = pred_region, aes(x= 1 - x, y= y, colour = Test), linetype = 2, size = 0.4, inherit.aes = F) +                         # prediction region
  geom_polygon(data = credible_region, aes(x= 1 - x, y= y, colour = Test), alpha=0.05, size=0.4, linetype = 2,inherit.aes = F) + # conf region
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_color_brewer(palette="Set1")
g



tiff("roc_plot_5tests.tif",units = "in", width = 7*1.25, height=5.2*1.25, res=500, compression = "lzw")
g
dev.off()




##########################################################################
# Plots for three-way
########################## PLOT 1 (w/ summaries) ##########################


#############################################
### add in threshold info for fluorescence (for three-way plot only, not NMA )

data_flu <- read_xlsx("comp data with thresholds for fluor (002).xlsx")
#View(data_flu)

data_flu2 <- dplyr::filter(data_flu, `f thresh` %in% abs(c(1:27)) )  %>% 
  dplyr::select(study, tp,fp,fn,tn, `f thresh`) %>% 
  mutate(Study = as.numeric(as.factor(study)), se = tp/(tp+fn), sp = tn/(fp+tn)) %>%
  dplyr::rename(Sensitivity = se, Specificity = sp, TN = tn, FP=fp,FN=fn,TP=tp)
data_flu2
data_flu3 <- (data_flu2[!duplicated(data_flu2[ , c("Study", "TP", "FN","FP","TN" )]),])

ss3 <- left_join( ss,  data_flu3, by = c("Study", "TP", "FN","FP","TN","Sensitivity","Specificity" )) 
ss3

ss_flu <- filter(ss, Test == "Fluorescence")
ss_flu_thresh <- left_join( ss_flu,  data_flu3, by = c( "TP", "FN","FP","TN","Sensitivity","Specificity" ))  %>%
                 dplyr::rename(Study = Study.x) %>%
                 select(-study, -Study.y)

ss_wo_flu  <- filter(ss, Test != "Fluorescence")
ss4 <- rbind(ss_flu_thresh, mutate(ss_wo_flu, `f thresh` = NA))

#View(ss3)
#View(data_flu3)
#View(ss)

#ss2_wo_flu <- filter(ss,  as.integer(Test) != 1 ) %>% dplyr::mutate(`f thresh` = rep(NA, 38))
#ss4 <- rbind(ss3, ss2_wo_flu)

#############################################
### plot

p1 <- ggplot(data = ss4, aes(y=Sensitivity, x = 1-Specificity, colour = Test, group = Test)) + 
  geom_point(data = ss4, aes(y=Sensitivity, x = 1-Specificity, colour = Test, shape = Test), size = 3, alpha=0.7,inherit.aes = F) + 
  geom_text_repel(data = ss4, aes(y=Sensitivity, x = 1-Specificity, colour = Test, label = `f thresh`)) + 
           #  hjust = -0.3, vjust = -0.5, size=4, label.size = 0, nudge_x = 0, show.legend = FALSE, alpha=0)  +
  geom_point(data = medians, aes(y=median_sens, x = 1- median_spec, colour = Test), size=10, shape = 5,inherit.aes = F)  +      # summary points
  geom_path(data = pred_region, aes(x= 1 - x, y= y, colour = Test), linetype = 2,inherit.aes = F) +                         # prediction region
  geom_polygon(data = credible_region, aes(x= 1 - x, y= y, colour = Test), alpha=0.05, size=0.1, linetype = 2,inherit.aes = F) + # conf region
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  +
  theme(legend.position="none") 
p1

########################## PLOT 2 (w/o summaries) ##########################
p2 <- ggplot(data = ss4, aes(y=Sensitivity, x = 1-Specificity, colour = Test, group = Study)) + 
  geom_point(data = ss4, aes(y=Sensitivity, x = 1-Specificity, colour = Test, shape = Test), size = 3, alpha=0.7,inherit.aes = F) + 
  geom_line(data =  ss4, aes(y=Sensitivity, x = 1-Specificity, group = Study), linetype=1, alpha=0.3)  +
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))# +
  #theme(legend.position="none")
 # labs(title= "Study-specific points with study weights, 
#points form the same study are joined by a line")
p2

####### PLOT 1 and PLOT 2 using Patchwork ###########

tiff("roc_plot_3tests.tif",units = "in", width = 11.5, height= 5, res=500, compression = "lzw")

p1 + p2 +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 20)),
   # title = "Direct comparisons of fluorescence, visual and radiograph tests",
  #  tag_levels = c("A")
    )

dev.off()








