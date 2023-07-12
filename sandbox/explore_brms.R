library(brms)

fit1 <- brm(formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
            data = kidney, family = lognormal(), prior = c(set_prior("normal(0,5)", class = "b"),
                                                           set_prior("cauchy(0,2)", class = "sd"),
                                                           set_prior("lkj(2)", class = "cor")),
            warmup = 1000, iter = 2000, chains = 4, control = list(adapt_delta = 0.95))
stancode(fit1)
standata(fit1)
summary(fit1, waic = TRUE)
plot(fit1)

hypothesis(fit1, "Intercept- age > 0", class = "sd", group = "patient")
fit2 <- update(fit1, formula. = ~ .- (1 + age|patient) + (1|patient))
LOO(fit1, fit2)
