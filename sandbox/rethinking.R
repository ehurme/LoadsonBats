library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = 'schools.stan', data = schools_dat)


install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
devtools::install_github("rmcelreath/rethinking")
library(rethinking)

f <- alist(
  y ~ dnorm( mu , sigma ),
  mu ~ dnorm( 0 , 10 ),
  sigma ~ dexp( 1 )
)

fit <- quap(
  f ,
  data=list(y=c(-1,1)) ,
  start=list(mu=0,sigma=1)
)

summary(fit)

fit_stan <- ulam( f , data=list(y=c(-1,1)) )
summary(fit_stan)
stancode(fit_stan)

fit_stan <- ulam(
  alist(
    y ~ normal( mu , sigma ),
    mu ~ normal( 0 , 10 ),
    sigma ~ exponential( 1 )
  ), data=list(y=c(-1,1)) )


# prep data
data( UCBadmit )
UCBadmit$male <- as.integer(UCBadmit$applicant.gender=="male")
UCBadmit$dept <- rep( 1:6 , each=2 )
UCBadmit$applicant.gender <- NULL

# varying intercepts model
data("UCBadmit")
m_glmm1 <- ulam(
  alist(
    admit ~ binomial(applications,p),
    logit(p) <- a[dept] + b*male,
    a[dept] ~ normal( abar , sigma ),
    abar ~ normal( 0 , 4 ),
    sigma ~ half_normal(0,1),
    b ~ normal(0,1)
  ), data=UCBadmit )
summary(m_glmm1)

m_glmm2 <- ulam(
  alist(
    admit ~ binomial(applications,p),
    logit(p) <- a[dept] + b[dept]*male,
    c( a , b )[dept] ~ multi_normal( c(abar,bbar) , Rho , sigma ),
    abar ~ normal( 0 , 4 ),
    bbar ~ normal(0,1),
    sigma ~ half_normal(0,1),
    Rho ~ lkjcorr(2)
  ),
  data=UCBadmit )

m_glmm3 <- ulam(
  alist(
    admit ~ binomial(applications,p),
    logit(p) <- v[dept,1] + v[dept,2]*male,
    vector[2]:v[dept] ~ multi_normal( c(abar,bbar) , Rho , sigma ),
    abar ~ normal( 0 , 4 ),
    bbar ~ normal(0,1),
    sigma ~ half_normal(0,1),
    Rho ~ lkjcorr(2)
  ),
  data=UCBadmit )

