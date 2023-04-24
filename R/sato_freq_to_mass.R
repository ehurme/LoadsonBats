duration <- 8 # hours
slope <- 0.01
m1 <- 80
f1 <- 6.8
f2 <- f1+slope*duration

m2 <- (f2/f1)^2*m1
m2


# explore the space of this function
f1 <- seq(6, 8, length = 200)
f2 <- f1+slope*duration
m2 <- (f2/f1)^2*m1
plot(f1, m2)
