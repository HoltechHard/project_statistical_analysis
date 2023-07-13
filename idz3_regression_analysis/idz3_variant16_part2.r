#############################################################
#                   HOMEWORK 03 - ANOVA                     #
#                      VARIANT N° 16                        #
#############################################################


# Get input data

dataset <- read.delim(file = "data_v16_idz3_part2.txt", sep = "", dec = ".")
print("Input data: ")
data <- data.frame(Y = dataset["Y"], A = dataset["A"], B = dataset["B"])
print(data)

# factors: data$A; data$B
# value: data$Y

alpha1 <- 0.01
h <- 0.86


#                   ********** PROBLEM 02 - A *************

# formulate model of 2 factors for ANOVA

print("\n===========================================================")
print("Problem 02 - A: Formulate model of 2 factors for ANOVA")
print("=============================================================\n")

z1 <- data$A
z2 <- data$B
Y <- data$Y

q <- lm(Y ~ as.factor(z1) * as.factor(z2))
qs <- summary(q)
qc <- q$coefficients
print(qc)

print("variance: ")
varx <- sum(q$residuals^2) / (40-20)
print(varx)

print("std-dev: ")
stv <- sqrt(sum(q$residuals^2)/20)
print(stv)


#                   ********** PROBLEM 02 - B *************

# visual representation of additive influence between factors
# x_axis: z2; lines: z1; y_axis: mean for each group eta_ij

print("\n===========================================================")
print("Problem 02 - B: Visual representation of influence of factors")
print("=============================================================\n")

# calculate estimator

data <- data[order(data$B), ]
data <- data[order(data$A), ]
data$ow <- paste0(data$A, ":", data$B)
ow.lev <- levels(as.factor(data$ow))
n.lev <- length(ow.lev)

# matrix X with ones for each corresponding factor

n <- dim(data)[1]

# matrix of coefficients 1 belongs with factors
XM <- matrix(0, nrow = n.lev, ncol = n)
for(i in 1:n.lev){
    XM[i, (data$ow == ow.lev[i])] <- 1
}

# matrix of Y's real values/ observed
YM <- as.matrix(Y)

# model of linear regression: E(Y | (z1, z2)) = XT * eta
# original model: Y = XT * eta
# build estimator: eta_hat = (X * XT)^-1 * X * Y
S <- XM %*% t(XM)
S1 <- solve(S)
eta_hat <- S1 %*% XM %*% YM

## calculate the mean for each group nij

nij_hat <- NULL
nms <- NULL

for(i in unique(data$ow)){
    nij_hat <- c(nij_hat, mean(data$Y[data$ow == i]))
    nms <- c(nms, i)
}

names(nij_hat) <- nms

# result of means nij_hat
print("means for each factor nij_hat:")
print(nij_hat)

# plot the graphic of relationships between factors z1, z2 and Y_mean

# number of factor for each level
nl_z1 <- length(unique(z1))     # z1 = {1, 2, 3, 4, 5}
nl_z2 <- length(unique(z2))     # z2 = {1, 2, 3, 4}

###### Z1 IN X AXIS #####


# possible values for z1 (A)
x <- unique(z1)

# define lines for z2 (B)

# generate z2 = 1
y1_ <- NULL

for(i in unique(z1)){
    y1_ <- c(y1_, nij_hat[paste0(i, ":", 1)])
}

y1_ <- as.numeric(y1_)

# generate z2 = 2
y2_ <- NULL

for(i in unique(z1)){
    y2_ <- c(y2_, nij_hat[paste0(i, ":", 2)])
}

y2_ <- as.numeric(y2_)

# generate z2 = 3
y3_ <- NULL

for(i in unique(z1)){
    y3_ <- c(y3_, nij_hat[paste0(i, ":", 3)])
}

y3_ <- as.numeric(y3_)

# generate z2 = 4

y4_ <- NULL

for(i in unique(z1)){
    y4_ <- c(y4_, nij_hat[paste0(i, ":", 4)])
}

y4_ <- as.numeric(y4_)

# plot graphics

# for y1
plot(x, y1_, "l", col = "red", ylim = c(min(c(y1_, y2_, y3_, y4_)),
                                        max(c(y1_, y2_, y3_, y4_))), lwd = 4)
points(x, y1_, col = "red", pch = 3, cex = 2)

# for y2
points(x, y2_, "l", col = "blue", ylim = c(min(c(y1_, y2_, y3_, y4_)),
                                        max(c(y1_, y2_, y3_, y4_))), lwd = 4)
points(x, y2_, col = "blue", pch = 3, cex = 2)

# for y3
points(x, y3_, "l", col = "green", ylim = c(min(c(y1_, y2_, y3_, y4_)),
                                        max(c(y1_, y2_, y3_, y4_))), lwd = 4)
points(x, y3_, col = "green", pch = 3, cex = 2)

# for y4
points(x, y4_, "l", col = "orange", ylim = c(min(c(y1_, y2_, y3_, y4_)),
                                        max(c(y1_, y2_, y3_, y4_))), lwd = 4)
points(x, y4_, col = "orange", pch = 3, cex = 2)

legend("topleft", legend = c("z2 = 1", "z2 = 2", "z2 = 3", "z2 = 4"),
                col = c("red", "blue", "green", "orange"), lwd = 2)


#                   ********** PROBLEM 02 - C *************

print("\n===========================================================")
print("Problem 02 - C: Proof XI-2 of residuals")
print("=============================================================\n")

# make analysis of residuals
q <- lm(Y ~ as.factor(z1) * as.factor(z2))
qs <- summary(q)
qq <- anova(q)

# calculate the residuals
res <- q$residuals
# calculate the standard-deviation
s1 <- qs$sigma

# build the breaks
m1 <- min(res)
m2 <- max(res)
mm <- m2 - m1

brk <- m1 + h * (c(0:ceiling(mm / h)))

# plot the histogram of residuals
hh <- hist(res, breaks = brk, probability = TRUE,
                xlim = c(m1 - 1, m2 + 1),
                main = "Histogram of residuals")

# plot normal distribution of standard deviation
x3 <- (m1 - 5) + c(0:1000) * (m2 - m1 + 7) / 1000
y3 <- dnorm(x3, 0, s1)
points(x3, y3, "l", col = "red", lwd = 2, xlim = c(m1 - 0.5, m2 + 1))

# take the intervals with count >= 5

brk1 <- c(min(res), -0.455, 0.405, max(res))
hh1 <- hist(res, breaks = brk1, freq = TRUE, col = "green",
                xlim = c(min(res), max(res)))
nu <- hh1$counts

brk1[1] <- -Inf
l.b <- length(brk1)
brk1[l.b] <- Inf

print("Breaks: ")
print(brk1)
print("Frequency: ")
print(nu)

csq0 <- function(s){
    if(s>0){
        p <- pnorm(brk1[2:l.b], 0, s) - pnorm(brk1[1:l.b - 1], 0, s)
        return(sum((nu - n * p)^2 / n / p))
    }else{
        return(Inf)
    }
}

# compute parameter xi2
csq.s <- nlm(csq0, p = s1)$minimum
k <- length(brk1) - 1

# for hypothesis proof: Xi2_k-1-d, where d = 1
print("--- Parameters Xi-2 manually ---")

# best level of significance pv
pv <- pchisq(csq.s, df = k - 2, lower.tail = FALSE)
print(paste("p-value = ", pv))

xal <- qchisq(p = 1 - alpha1, df = k - 2)
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


#                   ********** PROBLEM 02 - D *************

print("\n===========================================================")
print("Problem 02 - D: Proof of simple hypothesis with Xi2 criterion")
print("=============================================================\n")

q <- lm(Y ~ as.factor(z1) * as.factor(z2), data = data)
qa <- anova(q)
res <- q$residuals
SSE <- sum(res^2)
rdf <- n - length(q$coefficients)   # 40 -5 x 4 = 20

# main model: y_hat = mu + alpha_1i + alpha_2j + alpha_12

# linear model for H(12): alpha_12 = 0
q12 <- lm(Y ~ as.factor(z1) + as.factor(z2), data = data)
# linear model for H(1): alpha_1i = 0, alpha_12 = 0
q1 <- lm(Y ~ as.factor(z2), data = data)
# linear model for H(2): alpha_2j = 0, alpha_12 = 0
q2 <- lm(Y ~ as.factor(z1), data = data)
# linear model for H(a): alpha_1i = 0, alpha_2j = 0, alpha_12 = 0
q0 <- lm(Y ~ 1, data = data)

# test of H(12)
q12a <- anova(q12, q)
nms <- names(q12a)
a.o.v <- matrix(unlist(q12a), nrow = 2, ncol = 6)[2, 3:6]
names(a.o.v) <- nms[3:6]
xal <- qf(1 - alpha1, a.o.v[1], rdf)
a.o.v <- c("H_(12)", a.o.v, xal)
AOV_H12 <- as.data.frame(t(a.o.v))
names(AOV_H12) <- c("H", nms[3:6], "x_alpha")
print("proof of H(12): ")
print(AOV_H12)

# test H(1)
q1a <- anova(q1, q)
a.o.v1 <- matrix(unlist(q1a), nrow = 2, ncol = 6)[2, 3:6]
names(a.o.v1) <- nms[3:6]
xal2 <- qf(1 - alpha1, a.o.v1[1], rdf)
a.o.v1 <- c("H_(1)", a.o.v1, xal2)
AOV_H1 <- as.data.frame(t(a.o.v1))
names(AOV_H1) <- c("H", nms[3:6], "x_alpha")
print("proof of H(1): ")
print(AOV_H1)

# test H(2)
q2a <- anova(q2, q)
a.o.v2 <- matrix(unlist(q2a), nrow = 2, ncol = 6)[2, 3:6]
names(a.o.v2) <- nms[3:6]
xal3 <- qf(1 - alpha1, a.o.v2[1], rdf)
a.o.v2 <- c("H_(2)", a.o.v2, xal3)
AOV_H2 <- as.data.frame(t(a.o.v2))
names(AOV_H2) <- c("H", nms[3:6], "x_alpha")
print("proof of H(2): ")
print(AOV_H2)

# test H(120) of no influence of factor 1 and factor 2
q0a <- anova(q0, q)
a.o.v3 <- matrix(unlist(q0a), nrow = 2, ncol = 6)[2, 3:6]
names(a.o.v3) <- nms[3:6]
xal4 <- qf(1 - alpha1, a.o.v3[1], rdf)
a.o.v3 <- c("H_(120)", a.o.v3, xal4)
AOV_H120 <- as.data.frame(t(a.o.v3))
names(AOV_H120) <- c("H", nms[3:6], "x_alpha")
print("proof of H(120): ")
print(AOV_H120)


#                   ********** PROBLEM 02 - E *************

print("\n===========================================================")
print("Problem 02 - E: Select the best model using AIC and BIC")
print("=============================================================\n")

# table of AIC

print("Criteria AIC")
print(paste("y ~ A * B ==> AIC = ", AIC(q)))
print(paste("y ~ A + B ==> AIC = ", AIC(q12)))
print(paste("y ~ B ==> AIC = ", AIC(q1)))
print(paste("y ~ A ==> AIC = ", AIC(q2)))
print(paste("y ~ 1 ==> AIC = ", AIC(q0)))

# table of BIC

print("Criteria BIC")
print(paste("y ~ A * B ==> BIC = ", BIC(q)))
print(paste("y ~ A + B ==> BIC = ", BIC(q12)))
print(paste("y ~ B ==> BIC = ", BIC(q1)))
print(paste("y ~ A ==> BIC = ", BIC(q2)))
print(paste("y ~ 1 ==> BIC = ", BIC(q0)))

