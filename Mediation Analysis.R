rm(list = ls())

install.packages("mediation")
library(mediation)

set.seed(2014)
data("framing", package = "mediation")

?framing

med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)

out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income, data = framing, family = binomial("probit"))

med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", robustSE = TRUE, sims = 1000)

summary(med.out)


#BootStrap

med.out <- mediate(med.fit, out.fit, boot = TRUE, treat = "treat", mediator = "emo", sims = 1000)

summary(med.out)



#Interaction

med.fit <- lm(emo ~ treat + age + educ + gender + income, data=framing)
out.fit <- glm(cong_mesg ~ emo * treat + age + educ + gender + income, data = framing, family = binomial("probit"))

med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", robustSE = TRUE, sims = 100)
summary(med.out)


test.TMint(med.out, conf.level = .95)


#Moderators

#Looking at ACME at different levels of the covariate 

med.fit <- lm(emo ~ treat * age + educ + gender + income, data=framing)
out.fit <- glm(cong_mesg ~ emo + treat * age + emo * age + educ + gender + income, data = framing, family = binomial("probit"))

med.age20 <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", covariates = list(age = 20), sims = 100)
med.age60 <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", covariates = list(age = 60), sims = 100)

summary(med.age20)
summary(med.age60)

#Testing the difference between ACME and ADE at different levels of the covariate

med.init <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", sims=2)
test.modmed(med.init, covariates.1 = list(age = 20), covariates.2 = list(age = 60), sims = 100)


#Having multiple treatments

med.fit <- lm(emo ~ cond + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + cond + age + educ + gender + income, data = framing, family = binomial("probit"))


#Contrasting 2 and 3 treatments
med23.out <- mediate(med.fit, out.fit, treat = "cond", mediator = "emo", control.value = 2, treat.value = 3, sims = 100)
summary(med23.out)

#Contrasting 1 and 4 treatments
med14.out <- mediate(med.fit, out.fit, treat = "cond", mediator = "emo", control.value = 1, treat.value = 4, sims = 100)
summary(med14.out)



#Sensitivity Analysis

med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income, data = framing, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", robustSE = TRUE, sims = 100)

sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "indirect", sims = 100)
summary(sens.out)

sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "direct", sims = 100)
summary(sens.out)

sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "both", sims = 100)
summary(sens.out)


plot(sens.out, sens.par = "rho", main = "Anxiety", ylim = c(-0.2, 0.2))
plot(sens.out, sens.par = "R2", r.type = "total", sign.prod = "positive")

