library(tidyverse)
library(glmnet)
library(hdm)
library(plm)

set.seed(2018)


# download (about 450 MB)
# https://s3.us-east-2.amazonaws.com/criteo-uplift-dataset/
download.file("https://s3.us-east-2.amazonaws.com/criteo-uplift-dataset/criteo-uplift.csv.gz", "criteo-uplift.csv.gz")
# Read dataset
df_a <- read.csv(gzfile("criteo-uplift.csv.gz"))

# stratified sampling by treatment
df_sampled <- bind_rows(
  filter(df_a, treatment==1) %>% sample_frac(.002),
  filter(df_a, treatment==0) %>% sample_frac(.017)
) %>% select(-conversion, -exposure) %>%
  rename(W=treatment, Y=visit)

summary(df_a)
summary(df_sampled)

# Naive ATE
naive_ate <- function(formula=Y~W, data, method="naive"){
  outcome <- as.list(formula)[[2]] %>% as.character
  treatment <- attr(terms(formula), "term.labels")[1]
  print(outcome)
  print(treatment)
  # group mean and variance
  mean <- data.frame(y=data[, outcome],
                     w=data[, treatment]) %>%
    group_by(w) %>%
    summarise(avg=mean(y),
              var=var(y),
              num=n()) %>%
    mutate(var_weight = var/(num - 1))
  # calculate ATE and C.I.
  tau_hat <- mean$avg[mean$w==1] - mean$avg[mean$w==0]
  se_hat <- sum(mean$var_weight) %>% sqrt
  return(data.frame(
    Method = method,
    ATE = tau_hat,
    lower_ci = tau_hat - se_hat * 1.96,
    upper_ci = tau_hat + se_hat * 1.96,
    stringsAsFactors = F
  ))
}

df_result <- data.frame()
df_result <- bind_rows(df_result, naive_ate(Y~W, df_sampled, method = "naive"))

# Introduce Sampling bias
pt <- .80 # Drop p% of users who satisfy the following condition
pc <- .95
drop_idx <- c(which((df_sampled[,"f6"] < 2 | df_sampled[,"f0"] > 0.5) & df_sampled[,"W"] == 1)[1:round(pt*length(drop_treat_idx))],
              which((df_sampled[,"f6"] > 2 | df_sampled[,"f0"] < 0.5) & df_sampled[,"W"] == 0)[1:round(pc*length(drop_control_idx))]
              )
print(length(drop_idx))
df_mod <- df_sampled[-drop_idx, ]
df_mod <- df_mod[sample(nrow(df_mod)),]
rownames(df_mod) <- NULL

# naive_ate with sampling bias
df_result <- bind_rows(df_result, naive_ate(Y~W, df_mod , method = "naive_biased"))


# Single-eq Lasso
ate_condmean_lasso <- function(formula=Y~0+., data, treatment_var="W") {
  
  # glmnet requires inputs as matrices
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data))
  
  # Set the penalty to betaw to be zero
  pfac <- rep(1, NCOL(X)) 
  names(pfac) <- colnames(X)
  pfac[treatment_var] = 0
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(X, y, 
                     alpha=1, # 
                     penalty.factor=pfac) 
  
  # Automatically performs CV!
  betaw <- coef(model)[treatment_var, ]
  return(data.frame(Method = "Single-equation LASSO",
                    ATE = betaw,
                    lower_ci = NA,
                    upper_ci = NA,
                    stringsAsFactors = F))
}

# append the result
df_result <- bind_rows(
  df_result,
  ate_condmean_lasso(Y~0+., df_mod, treatment_var = "W"))


# Usual Lasso
ate_lasso <- function(formula=Y~0+., data, treatment_var="W") {
  # glmnet requires inputs as matrices
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data))
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(X, y, alpha=1) 
  
  # Automatically performs CV!
  betaw <- coef(model)[treatment_var,]
  return(data.frame(Method = "Usual LASSO",
                    ATE = betaw,
                    lower_ci = NA,
                    upper_ci = NA,
                    stringsAsFactors = F))
}


# append the result
df_result <- bind_rows(df_result, ate_lasso(Y~0+., df_mod))

# comparison
ggplot(df_result, aes(x = factor(Method, level=Method), y = ATE, color = Method)) + 
  geom_point() +
  geom_errorbar(aes(ymax = upper_ci, ymin = lower_ci), width = .1) + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

# Double ML Lasso
DML2.for.PLM <- function(y, D, Z, reg_D=hdm::rlasso, reg_Y=hdm::rlasso, n_folds = 5,
                         reg_D_args = NULL, reg_Y_args = NULL) {
  # this implements DML2 algorithm, where there moments are estimated via DML, before constructing	
  # the pooled estimate of theta randomly split data into folds
  n_obs <- NROW(y)
  foldid <- rep.int(1:n_folds, times = ceiling(n_obs/n_folds))[sample.int(n_obs)]
  index_fold <- split(1:n_obs, foldid)
  # initialize y, d
  ytil <- dtil <- rep(NA, n_obs)
  b <- 1  # iteration indicator
  cat("fold: ")
  for(fold in index_fold){
    # fit
    dfit <- do.call(reg_D, append(list(x=as.matrix(Z)[-fold, ],
                                       y=as.matrix(D)[-fold, ]),
                                  reg_D_args)
    )
    yfit <- do.call(reg_Y, append(list(x=as.matrix(Z)[-fold, ],
                                       y=as.matrix(y)[-fold, ]),
                                  reg_Y_args)
    )
    # take residuals
    dtil[fold] <- (as.numeric(D)[fold] - predict(dfit, as.matrix(Z)[fold, ], type="response"))
    ytil[fold] <- (as.numeric(y)[fold] - predict(yfit, as.matrix(Z)[fold, ], type="response"))
    cat(b," ")
    b <- b + 1
  }
  rfit <- lm(ytil ~ dtil) #estimate the main parameter by regressing one residual on the other
  coef.est <- coef(rfit)["dtil"] #extract coefficient
  se <- sqrt(vcovHC(rfit)["dtil", "dtil"])  # record standard error
  cat(sprintf("\ncoef (se) = %g (%g)\n", coef.est , se))
  return(data.frame(Method = "Double ML",
                    ATE = coef.est,
                    lower_ci = coef.est - se * 1.96,
                    upper_ci = coef.est + se * 1.96,
                    stringsAsFactors = F))
}

df_result <- bind_rows(
  df_result,
  DML2.for.PLM(y=df_mod$Y,
               D = df_mod$W,
               Z = dplyr::select(df_mod, starts_with("f")),
               reg_D = rlasso, reg_Y = rlasso)
)
# comparison
ggplot(df_result, aes(x = factor(Method, level=Method), y = ATE, color = Method)) + 
  geom_point() +
  geom_errorbar(aes(ymax = upper_ci, ymin = lower_ci), width = .1) + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))


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

covariates <- dplyr::select(df_mod, starts_with("f")) %>% colnames
p_lasso <- prop_score_lasso(df_mod, treatment_var = "W")

# TODO: ??????????
tauhat_ps_lasso <- prop_score_weight(dataset = df_mod, 
                                     p = p_lasso[,1], 
                                     treatment_var = "W", 
                                     outcome_var = "Y",
                                     method = "Propensity_Weighting_LASSOPS")

df_result <- df_result %>% bind_rows(tauhat_ps_lasso)

#ggplot
ggplot(df_result, aes(x = factor(Method, level=Method), y = ATE, color = Method)) + 
  geom_point() +
  geom_errorbar(aes(ymax = upper_ci, ymin = lower_ci), width = .1) + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))