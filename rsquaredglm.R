# Function rsquared.glmm requires models to be input as a list (can include fixed-
# effects only models,but not a good idea to mix models of class "mer" with models
# of class "lme")

rsquared.glmm=function(modlist) {
  # Iterate over each model in the list
  do.call(rbind,lapply(modlist,function(i) {
    # For models fit using lm
    if(class(i)=="lm") {
      Rsquared.mat=data.frame(Class=class(i),Family="Gaussian",
                              Marginal=summary(i)$r.squared,
                              Conditional=NA,AIC=AIC(i)) }
    # For general linear models fit using lme4, lmerTest, blme...
    else if(inherits(i, "lmerMod")) {
      # Get variance of fixed effects by multiplying coefficients by design matrix
      VarF=var(as.vector(fixef(i) %*% t(i@pp$X)))
      # Get variance of random effects by extracting variance components
      VarRand=colSums(do.call(rbind,lapply(VarCorr(i),function(j) j[1])))
      # Get residual variance
      VarResid=attr(VarCorr(i),"sc")^2
      # Calculate marginal R-squared (fixed effects/total variance)
      Rm=VarF/(VarF+VarRand+VarResid)
      # Calculate conditional R-squared (fixed effects+random effects/total variance)
      Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
      # Bind R^2s into a matrix and return with AIC values
      Rsquared.mat=data.frame(Class=class(i),Family="Gaussian",Marginal=Rm,
                              Conditional=Rc,AIC=AIC(update(i,REML=F))) }
    #For generalized linear models (family=="binomial") fit using lme4, blme...
    else if(inherits(i, "glmerMod") & summary(i)$family=="binomial") {
      # Get variance of fixed effects by multiplying coefficients by design matrix
      VarF=var(as.vector(fixef(i) %*% t(i@pp$X)))
      # Get variance of random effects by extracting variance components
      VarRand=colSums(do.call(rbind,lapply(VarCorr(i),function(j) j[1])))
      # Get residual variance
      VarResid=attr(VarCorr(i),"sc")^2
      # Calculate marginal R-squared
      Rm=VarF/(VarF+VarRand+pi^2/3)
      # Calculate conditional R-squared (fixed effects+random effects/total variance)
      Rc=(VarF+VarRand)/(VarF+VarRand+pi^2/3)
      # Bind R^2s into a matrix and return with AIC values
      Rsquared.mat=data.frame(Class=class(i),Family=summary(i)$family,Marginal=Rm,
                              Conditional=Rc,AIC=AIC(i)) }
    #For generalized linear models (family=="poisson") fit using lme4, blme...
    else if(inherits(i, "glmerMod") & summary(i)$family=="poisson") {
      # Get variance of fixed effects by multiplying coefficients by design matrix
      VarF=var(as.vector(fixef(i) %*% t(i@pp$X)))
      # Get variance of random effects by extracting variance components
      VarRand=colSums(do.call(rbind,lapply(VarCorr(i),function(j) j[1])))
      # Get residual variance
      VarResid=attr(VarCorr(i),"sc")^2
      # Get random effects names to generate null model
      rand.formula=reformulate(sapply(findbars(formula(i)),function(x) 
        paste0("(",deparse(x),")")),response=".")
      # Generate null model (intercept and random effects only, no fixed effects)
      null.mod=update(i,rand.formula)
      # Calculate marginal R-squared
      Rm=VarF/(VarF+VarRand+log(1+1/exp(as.numeric(fixef(null.mod)))))
      # Calculate conditional R-squared (fixed effects+random effects/total variance)
      Rc=(VarF+VarRand)/(VarF+VarRand+log(1+1/exp(as.numeric(fixef(null.mod)))))
      # Bind R^2s into a matrix and return with AIC values
      Rsquared.mat=data.frame(Class=class(i),Family=summary(i)$family,Marginal=Rm,
                              Conditional=Rc,AIC=AIC(i)) }
    # For model fit using nlme
    else if(class(i)=="lme") {
      # Get design matrix of fixed effects from model
      Fmat=model.matrix(eval(i$call$fixed)[-2],i$data)
      # Get variance of fixed effects by multiplying coefficients by design matrix
      VarF=var(as.vector(fixef(i) %*% t(Fmat)))
      # Get variance of random effects by extracting variance components
      VarRand=sum(suppressWarnings(as.numeric(VarCorr(i)
                  [rownames(VarCorr(i))!="Residual",1])),na.rm=T)
      # Get residual variance
      VarResid=as.numeric(VarCorr(i)[rownames(VarCorr(i))=="Residual",1])
      # Calculate marginal R-squared (fixed effects/total variance)
      Rm=VarF/(VarF+VarRand+VarResid)
      # Calculate conditional R-squared (fixed effects+random effects/total variance)
      Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
      # Bind R^2s into a matrix and return with AIC values
      Rsquared.mat=data.frame(Class=class(i),Marginal=Rm,Conditional=Rc,
                              AIC=AIC(update(i,method="ML")))
    } else { print("Function requires models of class lm, lme, mer, or merMod")
} } ) ) }

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
