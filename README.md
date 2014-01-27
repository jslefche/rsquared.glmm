# R-squared for generalized linear mixed-effects models

Created by Jon Lefcheck in Mar. 2013, based on the article by Nakagawa and Schielzeth (2013). 
Original blog post: http://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/

Modified by Juan Sebastian Casallas in Jan. 2014.

## Examples

### Generate a mock dataset

```
set.seed(9)
data <- data.frame(y=rnorm(100, 5, 10), y.binom=rbinom(100, 1, 0.5),
 y.poisson=rpois(100, 5), fixed1=rnorm(100, 20, 100),
 fixed2=c("Treatment1", "Treatment2"),rand1=LETTERS[1:2],
 rand1=LETTERS[1:2],
 rand2=c(rep("W", 25), rep("X", 25), rep("Y", 25), rep("Z", 25)))
```

### lme4
```
library(lme4)
#Linear model
mod0 <- lm(y ~ fixed1, data)
#Linear mixed effects model
mod1 <- lmer(y ~ fixed1 + (1|rand2/rand1), data)
mod2 <- lmer(y ~ fixed1 + fixed2 + (1|rand2/rand1), data)
rsquared.glmm(list(mod0, mod1, mod2))
#Generalized linear mixed effects model (binomial)
mod3 <- glmer(y.binom ~ fixed1*fixed2 + (1|rand2/rand1), family="binomial", data)
mod3.prob <- update(mod3, family = binomial(link = "probit"))
rsquared.glmm(list(mod3, mod3.prob))
#Generalized linear mixed effects model (poisson)
mod4 <- glmer(y.poisson ~ fixed1*fixed2 + (1|rand2/rand1), family="poisson", data)
mod4.sqrt <- update(mod4, family = poisson(link = "sqrt"))
rsquared.glmm(list(mod4, mod4.sqrt))
#Get values for all kinds of models
(lme4.models <- rsquared.glmm(list(mod0, mod1, mod2, mod3, mod3.prob, mod4, mod4.sqrt)))
```
### Compare output to MuMIn::r.squaredGLMM

`MuMIn::r.squaredGLMM` is similar to `r.squared` but cannot calculate r-squared for Poisson models.

```
library(MuMIn)
# Error for Poisson model
mumin.models <- do.call(rbind, lapply(list(mod0, mod1, mod2, mod3, mod3.prob, mod4, mod4.sqrt), r.squaredGLMM))
# Ignoring lm and poisson rows yields same results as ours
all.equal(lme4.models[2:5, 4:5], data.frame(mumin.models)[2:5,], check.attributes = F)
```

### blme

`blme` extends `lme4`, but yields different coefficients for random and fixed effects, which could explain the differences between their conditional r-squared values.

```
library(blme)
#Linear mixed effects model
blme.mod1 <- blmer(y ~ fixed1 + (1|rand2/rand1), data)
blme.mod2 <- blmer(y ~ fixed1 + fixed2 + (1|rand2/rand1), data)
#Generalized linear mixed effects model (binomial)
blme.mod3 <- bglmer(y.binom ~ fixed1*fixed2 + (1|rand2/rand1), family="binomial", data)
blme.mod3.prob <- update(blme.mod3, family = binomial(link = "probit"))
#Generalized linear mixed effects model (poisson)
blme.mod4 <- bglmer(y.poisson ~ fixed1*fixed2 + (1|rand2/rand1), family="poisson", data)
blme.mod4.sqrt <- update(blme.mod4, family = poisson(link = "sqrt"))
#Get values for all kinds of models
(blme.models <- rsquared.glmm(list(mod0, blme.mod1, blme.mod2, blme.mod3, blme.mod3.prob, blme.mod4, blme.mod4.sqrt)))
# blme models yield better conditional r-squared values
all.equal(lme4.models[-(1:2)], blme.models[-(1:2)])
```

### lmerTest

`lmerTest::lmer` extends `lme4::lmer` to allow anova calculations, but their random and mixed effects coefficients are the same.

```
# Try with lmerTest package -- output should be the same as above
library(lmerTest)
#Linear mixed effects model
lmerTest.mod1 <- lmer(y ~ fixed1 + (1|rand2/rand1), data)
lmerTest.mod2 <- lmer(y ~ fixed1 + fixed2 + (1|rand2/rand1), data)
rsquared.glmm(list(mod0, lmerTest.mod1, lmerTest.mod2))
(lmerTest.models <- rsquared.glmm(list(mod0, lmerTest.mod1, lmerTest.mod2)))
# Same results
all.equal(lme4.models[1:3, -(1:3)], lmerTest.models[,-(1:3)])
```

### lme

`nlme::lme` and `lme4::lmer` yield very similar r-squared values.

```
library(nlme)
lme.mod1 <- lme(y~fixed1,random=~1|rand2/rand1,data)
lme.mod2 <- lme(y~fixed1+fixed2,random=~1|rand2/rand1,data)
(lme.models <- rsquared.glmm(list(mod0,lme.mod1,lme.mod2)))
# Results only differ in decimal precision
all.equal(lme4.models[1:3, -(1:3)], lme.models[,-(1:3)], tol = 1e-4)
```