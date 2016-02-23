y <- angrist$totalRepeats
Tobs <- angrist$usesch
z <- angrist$vouch0

# Calculate all the observables
p <- mean(Tobs)
q <- mean(z)
p0 <- mean(Tobs[z == 0])
p1 <- mean(Tobs[z == 1])

p00 <- mean(Tobs == 0 & z == 0)
p01 <- mean(Tobs == 0 & z == 1)
p10 <- mean(Tobs == 1 & z == 0)
p11 <- mean(Tobs == 1 & z == 1)

yb00 <- mean(y[(Tobs == 0) & (z == 0)])
yb01 <- mean(y[(Tobs == 0) & (z == 1)])
yb10 <- mean(y[(Tobs == 1) & (z == 0)])
yb11 <- mean(y[(Tobs == 1) & (z == 1)])

yt00 <- (1 - p0) * yb00
yt01 <- (1 - p1) * yb01
yt10 <- p0 * yb10
yt11 <- p1 * yb11

s2_00 <- var(y[Tobs == 0 & z == 0])
s2_01 <- var(y[Tobs == 0 & z == 1])
s2_10 <- var(y[Tobs == 1 & z == 0])
s2_11 <- var(y[Tobs == 1 & z == 1])

# Calculate the IV coefficients etc
foo <- angrist[,c("strata1", "strata2", "strata3", "strata4",
                  "dmonth1", "dmonth2", "dmonth3", "dmonth4", "dmonth5",
                  "dmonth6", "dmonth7", "dmonth8", "dmonth9", "dmonth10",
                  "dmonth11", "dmonth12")]
bar <- apply(foo, 1, paste, collapse = "")
bar <- strtoi(bar, base = 2)
length(unique(bar))
