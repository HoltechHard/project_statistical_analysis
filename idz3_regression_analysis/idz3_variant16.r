#################################################################
#            HOMEWORK 03 - REGRESSION ANALYSIS                  #
#                       VARIANT N° 16                           #
#################################################################


# Get input data

dataset <- read.delim(file = "data_v16_idz3.txt", sep = "", dec = ".")
print("Input data: ")
data <- data.frame(n = dataset["n"], x = dataset["x"], y = dataset["y"])
print(data)

alpha <- 0.20
h <- 2.40

#                   ********** PROBLEM 01 - A *************

# 1.a - formulate linear regressor model y ~ x

print("\n=============================================================")
print("Problem 01 - A: Build lineal model (y = b1 + b2*x)")
print("===============================================================\n")

# number of samples
n <- length(data$y)

# coefficients of independent term
x0 <- array(1, dim = n)

# define extended matrix X and Y
X <- t(matrix(c(x0, data$x), nrow = n, ncol = 2))
Y <- as.matrix(data$y)

# calculate estimators: (X * Xt) * X * Y
S <- X %*% t(X)
S1 <- solve(S)
bhat <- S1 %*% X %*% Y
print("Estimators:")
print(paste("b1 = ", bhat[1], " || b2 = ", bhat[2]))

# 1.a - 1 build graphic results of the experiment
plot(data$x, data$y, main = "Linear regression", 
                    xlab = "X", ylab = "E(y|x)", cex = 1, 
                    sub = "y = 10.787123 - 0.02951252 * x")

# 1.a - 4 - build estimator line
x1 <- c(min(data$x), max(data$x))
y1 <- bhat[1] + bhat[2] * x1
points(x1, y1, "l", col = "blue", lwd = 2)

# 1.a - 5 - build estimator of variance
print("Estimators of error - lineal model")
res <- Y - t(X) %*% as.matrix(bhat)
SS <- sum(res^2)
s2 <- SS / (n - 2)
devs <- sqrt(s2)
print(paste("SSE = ", SS))
print(paste("sigma2_hat = ", s2))
print(paste("stand-dev = ", devs))

# other variant: use libraries

print("Using libraries of R")
q1 <- lsfit(data$x, data$y)
q <- lm(data$y ~ data$x)
qs <- summary(q)
print(qs)
print(paste("coeff = ", q$coefficients))
print(paste("res = ", sum(q$residuals^2) / (n - 2)))


#                   ********** PROBLEM 01 - B *************

# 1.b - build polinomial model with parameters (b1, b2, b3)

print("\n===================================================================")
print("Problem 01 - B: Build polinomial model (y = b1 + b2*x + b3 * x^2)")
print("=====================================================================\n")

# define the input data
n <- length(data$y)
x0 <- array(1, dim = n)
X <- t(matrix(c(x0, data$x, data$x^2), nrow = n, ncol = 3))
Y <- as.matrix(data$y)

# calculate estimators
S <- X %*% t(X)
S1 <- solve(S)
bhat2 <- S1 %*% X  %*% Y
plot(data$x, data$y, main = "Polinomial regressor", 
                xlab = "X", ylab = "E(y|x)", cex = 1, 
                sub = "y = 11.55211564 -0.69029028 * x + 0.084533796 * x^2")
# linear regressor
x1 <- c(min(data$x), max(data$y))
y1 <- bhat[1] + bhat[2] * x1

# polinomial regressor
x2 <- min(data$x) + c(0:1000) * (max(data$x) - min(data$x)) / 1000
y2 <- bhat2[1] + bhat2[2] * x2 + bhat2[3] * x2^2
print("Polinomial estimator")
print(paste("b1 = ", bhat2[1], " || b2 = ", bhat2[2], " || b3 = ", bhat2[3]))

# make graph
points(x1, y1, "l", col = "red", lwd = 2, cex = 2)
points(x2, y2, "l", col = "blue", lwd = 2, cex = 2)
legend("topleft", legend = c("Linear", "Quadratic"),
            col = c("red", "blue"), lwd = 2)

# calculate estimator of variance
print("Estimators of error - polinomial model")
res_pol <- Y - t(X) %*% as.matrix(bhat2)
SS_pol <- sum(res_pol^2)
s2_pol <- SS_pol / (n - 3)
devs_pol <- sqrt(s2_pol)
print(paste("SSE - polinomial = ", SS_pol))
print(paste("sigma2_hat - polinomial = ", s2_pol))
print(paste("stand-dev polinomial = ", devs_pol))


#                   ********** PROBLEM 01 - C *************

print("\n===========================================================")
print("Problem 01 - C: Proof of hypothesis Xi-2 for stand-dev")
print("=============================================================\n")

# build histogram of residuals
h0 <- hist(res_pol)

# build breaks
m1 <- min(res_pol)
m2 <- max(res_pol)
mm <- m2 - m1

# cut points
brk <- m1 + h * (c(0:ceiling(mm / h)))

# define final histogram of frequencies of residuals
hh <- hist(res_pol, breaks = brk, probability = TRUE,
                xlim = c(m1 - 0.5, m2 + 0.5), ylim = c(0, 0.1),
                main = "Histogram of residuals")

# define normal distribution of standard deviation
x3 <- (m1 - 5) + c(0:1000) * (m2 - m1 + 5) / 1000
y3 <- dnorm(x3, 0, devs_pol)
points(x3, y3, "l", col = "red", lwd = 3, cex = 2, xlim = c(m1 - 2, m2 + 2))

# print the breaks and frequencies
print("Breaks and frequencies - histogram of residuals")
print(hh$breaks)
print(hh$counts)

# make proof of confidence of standard-dev for Normal distribution for Xi2

# reformulate breaks: take intervals just with counts >=5

brk1 <- c(min(hh$breaks), -7.0834953, -4.6834953, 0.1165047,
                2.5165047, 4.9165047, 7.3165047, max(hh$breaks))
hh1 <- hist(res_pol, breaks = brk1, ylab = "Freq residuals - XI2", freq = TRUE,
                col = "green")

# prepare proof of hipothesis Xi2

lb <- length(brk1)
brk1[1] <- -Inf
brk1[lb] <- Inf
nu <- hh1$counts

print("new breaks: ")
print(brk1)
print("new frequencies: ")
print(hh1$counts)

# use function Xi2

csq0 <- function(s){
    if(s > 0){
        p <- pnorm(brk1[2:lb], 0, s) - pnorm(brk1[1:lb - 1], 0, s)
        return(sum((nu - n * p)^2 / n / p))
    }else{
        return(Inf)
    }
}

XM <- nlm(csq0, p = devs_pol)
print("some parameters: ")
print(XM)
csq.s <- nlm(csq0, p = devs_pol)$minimum
k <- length(brk1) - 1

# for hypothesis proof: Xi2_k-1-d, where d = 1
print("--- Parameters Xi-2 manually ---")
# best level of significance which can't reject H0
pv <- pchisq(csq.s, df = k - 2, lower.tail = FALSE)
print(paste("p-value = ", pv))
xal <- qchisq(p = 1 - alpha, df = k - 2)
print(paste("critical value x_alpha = ", xal))
print(paste("Xi-2 = ", csq.s))

#function phi for Xi-2 criteria
# phi(x) = { 0, if Xi2 <= x_al
#            1, if Xi2 > x_al
phi_xi <- as.numeric(csq.s > xal)

if(phi_xi == 0){
    print(paste(phi_xi, " ===> Acept H0, because ", csq.s, " <= ", xal))
}else{
    print(paste(phi_xi, " ===> Reject H0, because ", csq.s, " > ", xal))
}

# make proof xi2 using libraries => no es el mismo resultado vale el de arriba
print("--- Parameters Xi-2 with R-library ---")
pr <- pnorm(brk1[2:lb], 0, devs_pol) - pnorm(brk1[1:lb - 1], 0, devs_pol)
xi <- chisq.test(nu, p = pr)
print(xi)

#                   ********** PROBLEM 01 - D *************

print("\n===========================================================")
print("Problem 01 - D: build interval of confidence for b2 and b3")
print("=============================================================\n")

# not change the values of estimators b_hat
C <- diag(c(1, 1, 1))

# pshi = column matrix of estimators
# ph = (b1_hat, b2_hat, b3_hat)
ph <- bhat2

# term S1: (X * XT)^-1
S1 <- solve(S)

# estimator of variance b_pshi
# term V: result of diagonalization: b_pshi = CT * (X * XT)^-1 * C
V <- diag(S1)

xa <- qt(1 - alpha / 2, n - 3)
# standard deviation s2
s1 <- devs_pol
d <- xa * s1 * sqrt(V)

# definition of interval of confidence
# independent coeff: b1, linear: b2; quadratic: b3
CI <- data.frame(parameter = c("b1", "b2", "b3"), lw = ph - d, up = ph + d)
print("Interval of confidence: pshi = pshi_hat +/- xa * s * sqrt(b_pshi)")
print(CI)


#                   ********** PROBLEM 01 - E *************

print("\n===========================================================")
print("Problem 01 - E: Proof of linearity and independence")
print("=============================================================\n")

# proof of linearity
# -------------------
# H0 ==> linear: E(y) = b1 + b2 * x + b3 * x^2
# is linear if H0: b3 = 0 (hypothesis to build for linearity)

# parameter of F
FST1 <- bhat2[3]^2 / V[3] / s2_pol

# critical value x_alpha
xa2 <- qf(1 - alpha, 1, n - 3)

# p-value
pv.f <- pf(FST1, 1, n - 3, lower.tail = FALSE)

print("Proof of linearity")
print("---------------------")
print(paste("F = ", FST1))
print(paste("critical value x_alpha = ", xa2))
print(paste("p-value = ", pv.f))

# criterion F_beta3 and x_alpha

#function phi for F-sendecor
# phi(x) = { 0, if F <= x_al
#            1, if F > x_al
phi_fsendecor <- as.numeric(FST1 > xa2)

if(phi_fsendecor == 0){
    print(paste(phi_fsendecor, " ===> Acept H0, because ", FST1, " <= ", xa2))
}else{
    print(paste(phi_fsendecor, " ===> Reject H0, because ", FST1, " > ", xa2))
}

# with libraries ==> linear dependency
# prepare input variables
x <- data$x
x2 <- data$x^2
y <- data$y

q_lin <- lm(y ~ x + x2)
qs_lin <- summary(q_lin)
qa_lin <- anova(q_lin)
print("linearity with libraries")
print(qa_lin)

# using criterion for proof with p-value and alpha

#function phi2 for F-sendecor
# phi2(x) = { 0, if p-value > alpha
#             1, if p-value < alpha
phi2_fsendecor <- as.numeric(pv.f < alpha)

if(phi_fsendecor == 0){
    print(paste(phi2_fsendecor, " ===> Acept H0, because ",
                pv.f, " > ", alpha))
}else{
    print(paste(phi2_fsendecor, " ===> Reject H0, because ",
                pv.f, " <= ", alpha))
}


# proof of independency
# ----------------------

# H1: independence ==> E(y) = b1 + b2 * x + b3 * x^2
# is independent, if H1: b2 = b3 = 0

C <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 3, ncol = 2)
B <- t(C) %*% S1 %*% C
psi_ind <- t(C) %*% bhat2

# parameter of F
FST2 <- t(psi_ind) %*% solve(B) %*% (psi_ind) / s2_pol / 2

# critical value x_alpha
xa3 <- qf(1 - alpha, 2, n - 3)

# p-value
pv.f2 <- pf(FST2, 2, n - 3, lower.tail = FALSE)

# print results
print("Proof of independency")
print("--------------------------")
print(paste("F2 = ", FST2))
print(paste("critical value x_alpha2 = ", xa3))
print(paste("p-value = ", pv.f2))

# criterion F_beta2_beta3 and x_alpha

#function phi for F-sendecor
# phi(x) = { 0, if F23 <= x_al
#            1, if F23 > x_al
phi_fsendecor2 <- as.numeric(FST2 > xa3)

if(phi_fsendecor2 == 0){
    print(paste(phi_fsendecor2, " ===> Acept H0, because ", FST2, " <= ", xa3))
}else{
    print(paste(phi_fsendecor2, " ===> Reject H0, because ", FST2, " > ", xa3))
}

# with libraries ==> independency
# prepare input variables
q_ind <- lm(y ~ 1)
qaa <- anova(q_ind, q_lin)
print("independence with libraries")
print(qaa)

# using criterion for proof with p-value and alpha

#function phi2 for F-sendecor
# phi2(x) = { 0, if p-value > alpha
#             1, if p-value < alpha
phi2_fsendecor2 <- as.numeric(pv.f2 < alpha)

if(phi2_fsendecor2 == 0){
    print(paste(phi2_fsendecor2, " ===> Acept H0, because ",
                pv.f2, " > ", alpha))
}else{
    print(paste(phi2_fsendecor2, " ===> Reject H0, because ",
                pv.f2, " <= ", alpha))
}


#                   ********** PROBLEM 01 - F *************

print("\n================================================")
print("Problem 01 - F: Use AIC and BIC")
print("==================================================\n")

x <- data$x
y <- data$y

m <- list(formula(y ~ 1), formula(y ~ x), formula(y ~ I(x) + I(x^2)))
aic <- NULL
bic <- NULL

tm <- length(m)

for(i in 1:tm){
    q <- lm(m[[i]])
    aic <- c(aic, AIC(q))
    bic <- c(bic, BIC(q))
}

inf.tab <- data.frame(Model = c("No effect", "lineal", "quadratic"), 
                        AIC = aic, BIC = bic)

# take the minimum ==> the best model in AIC
inf.tab$Model[inf.tab$AIC == min(inf.tab$AIC)]
print(paste("AIC-best: ", inf.tab$Model[inf.tab$AIC == min(inf.tab$AIC)]))
print("AIC-table: ")
print(inf.tab$AIC)

# take the minimum ==> the best model in BIC
inf.tab$Model[inf.tab$BIC == min(inf.tab$BIC)]
print(paste("BIC-best: ", inf.tab$Model[inf.tab$BIC == min(inf.tab$BIC)]))
print("BIC-table: ")
print(inf.tab$BIC)

