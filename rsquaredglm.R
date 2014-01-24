#' R-squared and pseudo-rsquared for a list of (generalized) linear (mixed) models
#'
#' This function calls the generic \code{\link{r.squared}} function for each of the
#' models in the list and rbinds the outputs into one data frame
#'
#' @param a list of fitted (generalized) linear (mixed) model objects
#' @return a dataframe with one row per model, and "Class",
#'         "Family", "Marginal", "Conditional" and "AIC" columns
rsquared.glmm=function(modlist) {
  # Iterate over each model in the list
  do.call(rbind,lapply(modlist, r.squared))
}

#' R-squared and pseudo-rsquared for (generalized) linear (mixed) models
#'
#' This generic function calculates the r squared and pseudo r-squared for
#' a variety of(generalized) linear (mixed) model fits.
#' Currently implemented for \code{\link{lm}}, \code{\link{lmerTest::merMod}},
#' and \code{\link{nlme::lme}} objects.
#' Implementing methods usually call \code{\link{.rsquared.glmm}}
#'
#' @param mdl a fitted (generalized) linear (mixed) model object
#' @return Implementing methods usually return a dataframe with "Class",
#'         "Family", "Marginal", "Conditional", and "AIC" columns
r.squared <- function(mdl){
  UseMethod("r.squared")
}

#' Marginal r-squared for lm objects
#'
# This method extracts the variance for fixed and random effects, residuals,
# and the fixed effects for the null model (in the case of Poisson family),
#' and calls \code{\link{.rsquared.glmm}}
#'
#' @param mdl an merMod model (usually fit using \code{\link{lme4::lmer}},
#'        \code{\link{lme4::glmer}}, \code{\link{lmerTest::lmer}},
#'        \code{\link{blme::blmer}}, \code{\link{blme::bglmer}}, etc)
#' @return a dataframe with with "Class" = "lm", "Family" = "gaussian",
#'        "Marginal" = unadjusted r-squared, "Conditional" = NA, and "AIC" columns
r.squared.lm <- function(mdl){
  data.frame(Class=class(mdl), Family="gaussian",
             Marginal=summary(mdl)$r.squared,
             Conditional=NA, AIC=AIC(mdl))
}

#' Marginal and conditional r-squared for merMod objects
#'
#' This method extracts the variance for fixed and random effects, residuals,
#' and the fixed effects for the null model (in the case of Poisson family),
#' and calls \code{\link{.rsquared.glmm}}
#'
#' @param mdl an merMod model (usually fit using \code{\link{lme4::lmer}},
#'        \code{\link{lme4::glmer}}, \code{\link{lmerTest::lmer}},
#'        \code{\link{blme::blmer}}, \code{\link{blme::bglmer}}, etc)
r.squared.merMod <- function(mdl){
  # Get variance of fixed effects by multiplying coefficients by design matrix
  VarF <- var(as.vector(lme4::fixef(mdl) %*% t(mdl@pp$X)))
  # Get variance of random effects by extracting variance components
  VarRand <- colSums(do.call(rbind, lapply(lme4::VarCorr(mdl), function(x) x[1])))
  # Get residual variance
  VarResid <- attr(lme4::VarCorr(mdl),"sc")^2
  if(inherits(mdl, "lmerMod")){
    mdl.aic  <- AIC(update(mdl, REML=F))
    family <- "gaussian"
  }
  else if(inherits(mdl, "glmerMod")){
    family <- summary(mdl)$family
    mdl.aic <- AIC(mdl)
    # Pseudo-r-squared for poisson also requires the fixed effects of the null model
    if(family=="poisson") {
      # Get random effects names to generate null model
      rand.formula <- reformulate(sapply(findbars(formula(mdl)),
                                         function(x) paste0("(", deparse(x), ")")),
                                  response=".")
      # Generate null model (intercept and random effects only, no fixed effects)
      null.mdl <- update(mdl, rand.formula)
      # Get the fixed effects of the null model
      null.fixef <- as.numeric(lme4::fixef(null.mdl))
    }
  }
  # Call the internal function to do the pseudo r-squared calculations
  .rsquared.glmm(VarF, VarRand, VarResid, family = family,
                 mdl.aic = mdl.aic,
                 mdl.class = class(mdl),
                 null.fixef = null.fixef)
}

#' Marginal and conditional r-squared for lme objects
#'
#' This method extracts the variance for fixed and random effects,
#' as well as residuals, and calls \code{\link{.rsquared.glmm}}
#'
#' @param mdl an lme model (usually fit using \code{\link{nlme::lme}})
r.squared.lme <- function(mdl){
  # Get design matrix of fixed effects from model
  Fmat <- model.matrix(eval(mdl$call$fixed)[-2], mdl$data)
  # Get variance of fixed effects by multiplying coefficients by design matrix
  VarF <- var(as.vector(nlme::fixef(mdl) %*% t(Fmat)))
  # Get variance of random effects by extracting variance components
  VarRand <- sum(suppressWarnings(as.numeric(nlme::VarCorr(mdl)
                                          [rownames(nlme::VarCorr(mdl)) != "Residual",
                                           1])), na.rm=T)
  # Get residual variance
  VarResid <- as.numeric(nlme::VarCorr(mdl)[rownames(nlme::VarCorr(mdl))=="Residual", 1])
  # Call the internal function to do the pseudo r-squared calculations
  .rsquared.glmm(VarF, VarRand, VarResid, family = "gaussian",
                 mdl.aic = AIC(update(mdl, method="ML")),
                               mdl.class = class(mdl))
}

#' Marginal and conditional r-squared for glmm given fixed and random variances
#'
#' This function is based on Nakagawa and Schielzeth (2013). It returns the marginal
#' and conditional r-squared, as well as the AIC for each glmm.
#' Users should call the higher-level generic "r.squared", or implement a method for the
#' corresponding class to get varF, varRand and the family from the specific object
#'
#' @param varF Variance of fixed effects
#' @param varRand Variance of random effects
#' @family family family of the glmm (currently works with gaussian, binomial and poisson)
#' @param mdl.aic The model's AIC
#' @param mdl.class The name of the model's class
#' @param null.fixef a numeric vector containing the fixed effects of the null model. This
#'        parameter is only necessary for "poisson" family
#' @return A data frame with "Class", "Family", "Marginal", "Conditional", and "AIC" columns
.rsquared.glmm <- function(varF, varRand, VarResid, family,
                           mdl.aic, mdl.class, null.fixef = NULL){
  Rm <- NULL # Marginal R-squared
  Rc <- NULL # Conditional R-squared
  if(family == "gaussian"){
    # Calculate marginal R-squared (fixed effects/total variance)
    Rm <- varF/(varF+varRand+VarResid)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+VarResid)
  }
  else if(family == "binomial"){
    # Calculate marginal R-squared
    Rm <- varF/(varF+varRand+pi^2/3)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+pi^2/3)
  }
  else if(family == "poisson"){
    # Calculate marginal R-squared
    Rm <- varF/(varF+varRand+log(1+1/exp(null.fixef)))
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+log(1+1/exp(null.fixef)))
  }
  # Bind R^2s into a matrix and return with AIC values
  data.frame(Class=mdl.class, Family = family, Marginal=Rm, Conditional=Rc, AIC=mdl.aic)
}

### Examples ###

# set.seed(9)
# data=data.frame(y=rnorm(100,5,10),y.binom=rbinom(100,1,0.5),
#  y.poisson=rpois(100,5),fixed1=rnorm(100,20,100),
#  fixed2=c("Treatment1","Treatment2"),rand1=LETTERS[1:2],
#  rand1=LETTERS[1:2],
#  rand2=c(rep("W",25),rep("X",25),rep("Y",25),rep("Z",25)))
# 
# library(lme4)
# #Linear model
# mod0=lm(y~fixed1,data)
# #Linear mixed effects model
# mod1=lmer(y~fixed1+(1|rand2/rand1),data)
# mod2=lmer(y~fixed1+fixed2+(1|rand2/rand1),data)
# rsquared.glmm(list(mod0,mod1,mod2))
# #Generalized linear mixed effects model (binomial)
# mod3=glmer(y.binom~fixed1*fixed2+(1|rand2/rand1),family="binomial",data)
# rsquared.glmm(list(mod3))
# #Generalized linear mixed effects model (poisson)
# mod4=glmer(y.poisson~fixed1*fixed2+(1|rand2/rand1),family="poisson",data)
# rsquared.glmm(list(mod4))
# #Get values for all kinds of models
# lmer.models=rsquared.glmm(list(mod0,mod1,mod2,mod3,mod4));lmer.models
# 
# # Load library MuMIn to compare output to function 'r.squaredGLMM'
# library(MuMIn)
# do.call(rbind,lapply(list(mod0,mod1,mod2,mod3,mod4),r.squaredGLMM)) 
# # Error for Poisson model
# 
# # Try with lmerTest package -- output should be the same as above
# library(lmerTest)
# lmerTest.models=rsquared.glmm(list(mod0,mod1,mod2,mod3,mod4)); lmerTest.models
# lmer.models==lmerTest.models #This is generating odd results
# 
# detach(package:lme4,unload=T) #Parts of this package conflict with lme4
# require(nlme)
# mod0=lm(y~fixed1,data)
# mod1=lme(y~fixed1,random=~1|rand2/rand1,data)
# mod2=lme(y~fixed1+fixed2,random=~1|rand2/rand1,data)
# rsquared.glmm(list(mod0,mod1,mod2))
