library(tidyverse)
library(glmnet)
library(randomForest)
library(grf)
library(hdm)
library(hdm)	
library(randomForest)	
library(xtable)	
library(plm)

set.seed(2018)

#Read dataset
df_a <- read.csv("Downloads/criteo-uplift.csv")
df_treat <- df_a %>%
  filter(treatment == 1) 
df_treat_sampled <- sample_frac(df_treat, 0.002)


df_control <- df_a %>%
  filter(treatment == 0) 
df_control_sampled <- sample_frac(df_control, 0.017)


df <- rbind(df_treat_sampled, df_control_sampled)
df_check <- head(df)
rm(df_treat,df_treat_sampled,df_control,df_control_sampled)

summary(df)
summary(df_mod)
#make Dataset
df <- df %>%
  select(-conversion, -exposure)

df <- df %>%
  plyr::rename(c(treatment = "W",
                 visit= "Y"))
#Naive ATE
naive_ate <- function(dataset, W_var, outcome_var, method = "naive"){
  mean_df <- dataset %>%
    group_by_(W_var) %>%
    summarise_(y = paste("mean(", outcome_var, ")"),
               y_var = paste("var(", outcome_var, ")"),
               count = "n()") %>%
    mutate(y_var_weight = y_var/(count - 1))
  
  E_y0 = mean_df$y[mean_df$W == 0]
  E_y1 = mean_df$y[mean_df$W == 1]
  
  tau_hat <- E_y1 - E_y0
  se_hat <- sum(mean_df$y_var_weight) %>% sqrt()
  
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  
  return(data.frame(Method = method, ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}
result_df <- data.frame()
ate1 <- naive_ate(df , "W", "Y", method = "naive")
result_df <- result_df %>% rbind(ate1)
#Introduce Sampling bias
pt <- .80 # Drop p% of users who satisfy the following condition
pc <- .95

drop_from_treat <-  (df[,"f6"] < 2 | df[,"f0"] > 0.5) 

drop_from_control <-(df[,"f6"] > 2 | df[,"f0"] < 0.5)

drop_treat_idx <- which(df[,"W"] == 1 & drop_from_treat)
drop_control_idx <- which(df[,"W"] == 0 & drop_from_control)

drop_idx <- unique(c(drop_treat_idx[1:round(pt*length(drop_treat_idx))],
                     drop_control_idx[1:round(pc*length(drop_control_idx))]))

print(length(drop_idx))
df_mod <- df[-drop_idx,]
df_mod <- df_mod[sample(nrow(df_mod)),]
rownames(df_mod) <- NULL


#naive_ate with sampling bias
ate2 <- naive_ate(df_mod , "W", "Y", method = "naive_biased")
result_df <- result_df %>% rbind(ate2)

#ggplot
ggplot(result_df, aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#SIngle-eq Lasso
ate_condmean_lasso <- function(dataset, treatment_var, outcome_var) {
  # Covariate names
  regs <- c(covariates, treatment_var)
  
  # glmnet requires inputs as matrices
  x <- as.matrix(dataset[regs])  
  y <- as.matrix(dataset[,outcome_var])
  
  # Set the penalty to betaw to be zero
  pfac <- c(rep(1,length(covariates)), 0) 
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(x, y, 
                     alpha=1, # 
                     penalty.factor=pfac) 
  
  # Automatically performs CV!
  betaw <- coef(model)[treatment_var,]
  return(data.frame(Method = "Single-equation LASSO", ATE = betaw, lower_ci = betaw, upper_ci = betaw))
}

tauhat_lasso_seq <- ate_condmean_lasso(df_mod, treatment_var = "W", outcome_var = "Y")

result_df <- result_df %>% rbind(tauhat_lasso_seq)


#Usual Lasso
ate_lasso <- function(dataset, treatment_var, outcome_var) {
  # Covariate names
  regs <- c(covariates, treatment_var)
  
  # glmnet requires inputs as matrices
  x <- as.matrix(dataset[regs])  
  y <- as.matrix(dataset[,outcome_var])
  
  # Set the penalty to betaw to be zero
  pfac <- c(rep(1,length(covariates)), 1) 
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(x, y, 
                     alpha=1, # 
                     penalty.factor=pfac) 
  
  # Automatically performs CV!
  betaw <- coef(model)[treatment_var,]
  return(data.frame(Method = "Usual LASSO", ATE = betaw, lower_ci = betaw, upper_ci = betaw))
}

tauhat_lasso_all <- ate_lasso(df_mod, treatment_var = "W", outcome_var = "Y")

result_df <- result_df %>% rbind(tauhat_lasso_all)


# Double ML Lasso

DML2.for.PLM <- function(x, d, y, dreg, yreg, nfold=5) {	
  # this implements DML2 algorithm, where there moments are estimated via DML, before constructing	
  # the pooled estimate of theta randomly split data into folds	
  nobs <- nrow(x)	
  foldid <- rep.int(1:nfold,times = ceiling(nobs/nfold))[sample.int(nobs)]	
  I <- split(1:nobs, foldid)	
  # create residualized objects to fill	
  ytil <- dtil <- rep(NA, nobs)	
  # obtain cross-fitted residuals	
  cat("fold: ")	
  for(b in 1:length(I)){	
    dfit <- dreg(x[-I[[b]],], d[-I[[b]]])  #take a fold out	
    yfit <- yreg(x[-I[[b]],], y[-I[[b]]])  # take a folot out	
    dhat <- predict(dfit, x[I[[b]],], type="response")  #predict the fold out	
    yhat <- predict(yfit, x[I[[b]],], type="response")  #predict the fold out	
    dtil[I[[b]]] <- (d[I[[b]]] - dhat) #record residual	
    ytil[I[[b]]] <- (y[I[[b]]] - yhat) #record residial	
    cat(b," ")	
  }	
  rfit <- lm(ytil ~ dtil)               #estimate the main parameter by regressing one residual on the other	
  coef.est <- coef(rfit)[2]             #extract coefficient 	
  se <- sqrt(vcovHC(rfit)[2,2])         #record standard error	
  cat(sprintf("\ncoef (se) = %g (%g)\n", coef.est , se))	
  low_ <- coef.est - se*1.96
  upp_ <- coef.est + se*1.96
  return(data.frame(Method = "Double ML", ATE = coef.est, lower_ci = low_, upper_ci = upp_))
}	

y= as.matrix(df_mod[,14])          #outcome: growth rate	
d= as.matrix(df_mod[,13])          #treatment: initial wealth		
x= as.matrix(df_mod[,-c(13,14)])  #controls: country characteristics	
dreg <- function(x,d){ rlasso(x, d) }  #ML method= rlasso	
yreg <- function(x,y){ rlasso(x, y) }  #ML method = rlasso	
DML2.lasso = DML2.for.PLM(x, d, y, dreg, yreg, nfold=5)	

result_df <- result_df %>% rbind(DML2.lasso)
#ggplot
ggplot(result_df, aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Propensity Lasso
prop_score_lasso <- function(dataset, treatment_var) {
  # glmnet requires inputs as matrices
  x <- as.matrix(dataset[covariates])  
  w <- as.matrix(dataset[,treatment_var])
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(x, w, 
                     alpha=1, 
                     family="binomial") 
  
  # Automatically performs CV
  p <- predict(model, newx=x, type="response")
  return(p)
}

p_lasso <- prop_score_lasso(df_mod, treatment_var = "W")
tauhat_ps_lasso <- prop_score_weight(dataset = df_mod, 
                                     p = p_lasso[,1], 
                                     treatment_var = "W", 
                                     outcome_var = "Y",
                                     method = "Propensity_Weighting_LASSOPS")

result_df <- result_df %>% rbind(tauhat_ps_lasso)

#ggplot
ggplot(result_df, aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
