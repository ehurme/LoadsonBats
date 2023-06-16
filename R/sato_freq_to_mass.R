duration <- 8 # hours
slope <- 0.001
m1 <- 480
f1 <- 3.2
f2 <- f1+slope*duration

m2 <- (f2/f1)^2*m1
m2


# explore the space of this function
f1 <- seq(6, 8, length = 200)
f2 <- f1+slope*duration
m2 <- (f2/f1)^2*m1
plot(f1, m2)


# look at relationship within M thy
m1 <- 8.28
f1 <- 11.875
f2 <- 12.5# 12.25

m2 <- (f2/f1)^2*m1
m2

