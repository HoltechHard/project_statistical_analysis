#################################################################
#            HOMEWORK 02 - CLASSICAL STATISTIC                  #
#                       VARIANT N° 16                           #
#################################################################


# Get input data

x <- scan(file = "data_v16.txt", sep = "")
print("Input data: ")
print(x)

alpha2 <- 0.10
c <- 4.70
d <- 5.12
h <- 0.10
a0 <- 5
sig0 <- 0.30
a1 <- 5.30
sig1 <- 0.30


#                   ********** PROBLEM 01 *************

# 1.1 - sorted vector
xs <- sort(x)
print("Ordened data: ")
print(xs)

# 1.2 - compute empirical cumulative dist. function
# Fn(t) = #{xi <= t} / n
y <- ecdf(x)
plot.ecdf(x, pch = NA, verticals = TRUE, lwd = 4,
            xlab = "x", ylab = "p(x)")

# 1.3 - histogram
m1 <- min(x)
m2 <- max(x)
mm <- m2 - m1

# cut points
brks <- m1 + h * (c(0:ceiling(mm / h)))

#define histogram of probability density

# define histogram
hh <- hist(x, breaks = brks)

# 1.4 - frequency polygon

# coordinates of x and y of the polygon

xx <- c((min(hh$breaks) - (h / 2)), hh$mids, (max(hh$breaks) + (h / 2)))
yy <- c(0, hh$density, 0)

# draw the polygon
polygon(xx, yy, border = "red", lwd = 3)

# draw the function
dd <- points(density(x), type = "l", col = "blue", lwd = 3)


#                   *********** Problem 02 ************

print("\n=============================================")
print("Problem 02: Statistical Estimators")
print("=============================================\n")

# number of samples
n <- length(x)

# mean
m <- mean(x)

# variance
s2 <- mean(x^2) - m^2
s2a <- (n - 1) * var(x) / n

# median
med <- median(xs)

# asimetry
asi <- mean((x - m)^3) / (s2^(3 / 2))

# kurtosis
kur <- mean((x - m)^4) / (s2^2) - 3

# probability
prob <- sum((x >= c) & (x <= d)) / length(x)

estimators <- data.frame(mean = m, variance = s2, median = med,
                asymetry = asi, kurtosis = kur, prob = prob)
print(estimators)


#                      ********* Problem 03/04 ***********

# 4.1 - Define confidence interval for parameter a (mean)
# use distribution t-student
# a = [x* - tal * s/sqrt(n-1) ; x* + tal * s/sqrt(n-1)]

print("=============================================")
print("Problem 04: Confidence Interval")
print("=============================================")

# critical point of t-student dist.
qal <- qt(p = 1 - alpha2 / 2, df = n - 1)
il <- qal * sqrt(s2) / sqrt(n - 1)
ci_a <- c(m - il, m + il)

print(paste("Confidence interval for mean: [", 
            ci_a[1], " , ", ci_a[2], "]", sep = ""))

# 4.2 - Define confidence interval for parameter sigma2 (variance)
# for distribution Chi-quadrad
# sigma2 = [n * s^2/ xi2(alpha/2, n-1) ; n * s^2/ xi2(1 - alpha/2, n-1)]

x1a <- qchisq(alpha2 / 2, n - 1)
x2a <- qchisq(1 - alpha2 / 2, n - 1)
ci_s2 <- n * s2 / c(x2a, x1a)

print(paste("Confidence interval for variance: [",
            ci_s2[1], " , ", ci_s2[2], "]", sep = ""))


#                       *********** Problem 05 ************

# Table of Kolgomorov
#       |     V1     |      V2     |     V3    |          V4              |          V5           |     V6
#   i   |  Fn(xi -)  |   Fn(xi +)  |   F(xi)   |   abs(F(xi) - Fn(xi -))  | abs(F(xi) - Fn(xi+))  | max(V4, V5)
#   1   |  (i-1)/n   |     i/n     |   pnorm   |   abs(F(xi) - (i-1)/n)   |   abs(F(xi) - i/n)    | max(V4, V5)

# Dn = Sup | F(x) - Fn(x) |     => F(x): dist. applied in sample with x, a0, sig0 ;

# critical point: T = sqrt(n) * Dn

print("\n==================================================================")
print("Problem 05: Poof of simple hypothesis with Kolgomorov criterion")
print("==================================================================\n")

xs <- sort(x)

# generate table of Kolgomorov
i <- c(1:n)
v1 <- c(0:(n - 1)) / n
v2 <- v1 + (1 / n)
v3 <- pnorm(xs, a0, sig0)
v4 <- abs(v3 - v1)
v5 <- abs(v3 - v2)
v6 <- pmax(v4, v5)
ks_tab <- data.frame(i = i, lecdf = v1, recdf = v2, hcdf = v3,
                        ldif = v4, rdif = v5, maxdif = v6)
# get the value of Dn
D <- max(v6)

print("Table Kolgomorov criterion")
print(ks_tab)
print(paste("D = ", D))

# take the position
nn <- which(v6 == D)
print(paste("Position: ", nn))

# calculate lambda = sqrt(n) * Dn
D1 <- D * sqrt(n)
print(paste("D1: ", D1))

# Empirical function of probability distribution
plot.ecdf(x, verticals = TRUE, pch = NA, lwd = 5)
# points(xs, v3, "l", col = "blue")

# plot the biggest diference D
xx <- xs[1] - 10 + (xs[n] - xs[1] + 10) * c(0:1000) / 1000
yy <- pnorm(xx, a0, sig0)
points(xx, yy, "l", col = "red", lwd = 5)
# vertical line => difference between the empirical and teorethical function
points(c(xs[nn], xs[nn]), c(v3[nn], v1[nn]), "l", col = "green", lwd = 5)

distf <- function(x){
    pnorm(x, a0, sig0)
}

# function to make Kolgomorov-Smirnov Test
ks <- ks.test(x, distf)
print(ks)

# find the critical point x_alpha

# distribution function of Kolgomorov
pkolm <- function(x, n){
    .Call(stats:::C_pKolmogorov2x, x, n)
}

# quantil function of Kolgomorov
qkolm <- function(q, n){
    uniroot(function(x){pkolm(x,n) - q}, c(0, 1))$root
}

# find the quantil x_alpha in equation: K(x_alpha) = 1 - alpha
xal <- qkolm(1 - alpha2, n)
print(paste("Critical point x_alpha", xal))

# function of Kolmogorov criteria 
# phi(x) = { 0, if sqrt(n) * Dn <= x_alpha
#            1, if sqrt(n) * Dn > x_alpha
phi_kolm <- as.numeric(D1 > xal)

if(phi_kolm == 0){
    print(paste(phi_kolm, " ==> Accept H0, because ", D1, " <= ", xal))
}else{
    print(paste(phi_kolm, " ==> Reject H0, because ", D1, " > ", xal))
}

# calculate the significance level which can't reject H0 => Find p-value ?
# probability of obtaining test results at least as extreme as the result
# under the assumption that the H0 is correct
p_value <- 1 - pkolm(D, n)
print(paste("P-value = ", p_value))



#                       *********** Problem 06 ************

# Table of Xi - Person/Fisher  => intervals of histogram [aj, bj]; F: dist. function of prob. # nolint
#    |   V1   |    V2   |          V3          |          V4           |            V6              |  V7
#  j |   aj   |    bj   | nuj = freq([aj, bj]) |  Pjo = F(bj) - F(aj)  |  (nuj - n*Pjo)/sqrt(n*Pjo) | V6^2

# nuj => frequency of values x belongs to [aj, bj]
# F(aj) => dist. prob. function applied in subsamples [aj, bj] with a0, sig0
# Critical point: Xi2 = sum( ( (nuj - n*Pjo)/sqrt(n*Pjo) )^2)

print("\n===========================================================")
print("Problem 06: Proof of simple hypothesis with Xi2 criterion")
print("=============================================================\n")

# take intervals just with counts >= 5

brk1 <- c(min(x), 4.612, 4.712, 4.812, 4.912, 5.112, 5.212, max(x))
hh1 <- hist(x, breaks = brk1, ylab = "Frequency histogram", freq = TRUE,
                        col = "green", xlim = c(min(x), max(x)))

# redefine the intervals (brk1) of histogram
l_brk <- length(brk1)
brk1[1] <- -Inf
brk1[l_brk] <- Inf
print("Breaks: ")
print(brk1)

# V3: calculate freq(nuj = [aj, bj])
nu <- hh1$counts
print("Count: ")
print(nu)

# define the subintervals [aj, bj]
lw <- brk1[c(1:(l_brk - 1))]
up <- brk1[c(2:l_brk)]

# V4: calculate Pjo = F(bj) - F(aj)
pr <- pnorm(up, a0, sig0) - pnorm(lw, a0, sig0)

# V6: calculate (nu - n*Pjo)/sqrt(n*Pjo)    # nolint
res <- (nu - n * pr) / sqrt(n * pr)

# V7: calculate V6^2
res2 <- res^2

# grades of freedom: #intervals defined
ngr <- length(nu)

# print table of Xi-2 criterion
xi_tab <- data.frame(Group = c(1:ngr), lower = lw, upper = up, count = nu,
                        prob = pr, resid = res, resid2 = res2)
print("Table Xi-2 criterion")
print(xi_tab)

# Calculate Xi2
Xi2 <- sum(res2)
print(paste("Xi^2 = ", Xi2))

# function to make Xi-2 test
xi <- chisq.test(nu, p = pr)
print(xi)

# find quantil x_alpha2, which satisfy equation: Xi_r-1 (x_alpha) = 1 - alpha
x_al2 <- qchisq(p = 1 - alpha2, df = ngr - 1)
print(paste("Critical point x_alpha = ", x_al2))

# function phi for Xi-2 criteria
# phi(x) = { 0, if Xi2 <= x_alpha2
#            1, if Xi2 > x_alpha2
phi_xi <- as.numeric(Xi2 > x_al2)

if(phi_xi == 0){
    print(paste(phi_xi, " ==> Acept H0, because ", Xi2, " <= ", x_al2))
}else{
    print(paste(phi_xi, " ==> Reject H0, because ", Xi2, " > ", x_al2))
}

# calculate the best level of significance which can't reject H0 => p-value = ?
p_value2 <- pchisq(Xi2, ngr - 1, lower.tail = FALSE)
print(paste("p-value for Xi2 = ", p_value2))



#                   ******** Problem 07 *******

print("\n============================================================")
print("Problem 07: Poof of complex hypothesis with Xi2 criterion")
print("============================================================\n")

# calculate pr = Pj = F(bj) - F(aj), where x[1] = m, x[2] = sqrt(s2)
prob.norm <- function(x){
    pnorm(up, x[1], x[2]) - pnorm(lw, x[1], x[2])
}

# function to calculate Xi2 = sum{(nu - n*pj)^2/(n*pj)}
csq.stat <- function(x){
    ex <- n * prob.norm(x)
    return(sum((nu - ex)^2 / (ex)))
}

# apply the critery variant 1 to find optimum xi2*:

# list of Xi2 for each interval
# f: function to be minimized; p: parameters for minimization
qq <- nlm(f = csq.stat, p = c(m, sqrt(s2)))
print(qq)

# optimum xi2* = min{Xi2(theta = m, sqrt(s2))}
Xi2_hat <- qq$min
print(paste("Xi2_hat = ", Xi2_hat))

# quantil value x_alpha3, which satisfy equation: Xi2_r-1-d (x_al3) = 1 - alpha
# d = number of estimated parameters (in our case d = 2, because have m and s2)
d <- 2
x_al3 <- qchisq(p = 1 - alpha2, df = (ngr - 1) - d)
print(paste("Critical point x_alpha = ", x_al3))

# criteria confidence for complex hypothesis
# phi_xihat = { 0, if X2_hat < = x_alpha3
#                1, if x2_hat > x_alpha3
phi_complex <- as.numeric(Xi2_hat > x_al3)

if(phi_complex == 0){
    print(paste(phi_complex, " ==> Accept H0, because ", Xi2_hat, " <= ", x_al3))
}else{
    print(paste(phi_complex, " ==> Reject H0, because ", Xi2_hat, " > ", x_al3))
}

# find the best level of significance for don't reject H0: p-value
# if lower.tail = TRUE => P[X <= xal]; lower.tail = FALSE => P[X > xal]
p_value3 <- pchisq(q = Xi2_hat, df = (ngr - 1) - d, lower.tail = FALSE)
print(paste("p-value for complex hypothesis Xi2 = ", p_value3))



#                  ************ Problem 08 ***********

print("============================================================")
print("Problem 08: Powerfull criteria for proof simple hypotesis")
print("===========================================================")

# proof for H0: a0 = 5, sig0 = 0.3

# quantil C_alpha with a0 and sig0 for H0
c_alpha <- qnorm(1 - alpha2, n * a0, sqrt(n) * sig0)
print(paste("Critical point for H0 c_alpha = ", c_alpha))

sum_x <- sum(x)

# function phi for powerful criterion max. likelihood
phi_pow <- as.numeric(sum_x > c_alpha)

print("Proof H0")
if(phi_pow == 0){
    print(paste(phi_pow, " ==> Accept H0, because ", sum_x, " <= ", c_alpha))
}else{
    print(paste(phi_pow, " ==> Reject H0, because ", sum_x, " > ", c_alpha))
}

# proof for H1: a1 = 5.3, sig1 = 0.3

# quantil C1_alpha with a1 and sig1 for H1
c1_alpha <- qnorm(alpha2, n * a1, sqrt(n) * sig1)
print(paste("Critical point for H1 c1_alpha = ", c1_alpha))

# function phi for powerful criterion max. likelihood
phi1_pow <- as.numeric(sum_x < c1_alpha)

print("Proof H1")
if(phi1_pow == 0){
    print(paste(phi1_pow, " ==> Accept H1, because ", sum_x, " > ", c1_alpha))
}else{
    print(paste(phi1_pow, " ==> Reject H1, because ", sum_x, " < ", c1_alpha))
}
