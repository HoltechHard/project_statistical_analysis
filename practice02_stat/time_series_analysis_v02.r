#############################################################
#       TASK 02 - TIME SERIES ANALYSIS (LMM Models)         #
#############################################################
#                       VERSION 2.0

# Author: @HoltechHard
# Topic: Statistical Analysis of dynamics in the changes of
#        observed features belong the time

# LMM MODELS - Linear Mixed Models
# Mixed models has random and fixed effects

# Enunciate:
# ----------------------------
# Cohort Study of Individuals
# ----------------------------
# Variables:
# dependent variable: Y
# times: Visit
# Factors:
# A: constant respect the time
# B: variable respect the time
# Make a cohort study

# add libraries
library("lme4")
library("nlme")
library("ggplot2")
library("stats4")
library("cAIC4")
install.packages("lmerTest")
library("lmerTest")
install.packages("mratios")
library("mratios")
library("RLRsim")

# read the dataset
root <- "practice02_stat/"
dataset <- read.table(paste(root, "IDZ4_Data_30123.csv", sep = ""),
                sep = ",", header = TRUE, as.is = TRUE)

#############################
#          TASK 01          #
#############################

# collect the data which corresponds to variant n° 9
datavar9 <- dataset[dataset$Variant == 9, -1]
print(head(datavar9))
print("Dimensions dataset - Variant 9")
print(dim(datavar9))

# check the levels
lvA <- levels(as.factor(datavar9$A))
lvB <- levels(as.factor(datavar9$B))
lvVisit <- levels(as.factor(datavar9$Visit))

print("Levels for A: ")
print(lvA)      # {0, 1}
print("Levels for B: ")
print(lvB)      # {1, 2}
print("Levels for Visit: ")
print(lvVisit)  # {0, 1, 2, 3, 4}


#############################
#           TASK 02         #
#############################

# Show the graphics of observations without consider time.
# Make ANOVA of 2 factors considering the next structure:
# (A, B) -> Y

# 2.1 - Show the graphics without considering time

# plot relationship between B and Y
plot(datavar9$B, datavar9$Y, "n",
        main = "Graph of observations [not consider time]",
        xlab = "Factor B", ylab = "Feature Y")
legend("topright", title = "Factor A", c("A = 0", "A = 1"),
        fill = c("red", "blue"))

# take values for factor A = 0
ftr <- datavar9$A == 0
points(datavar9$B[ftr], datavar9$Y[ftr], col = "red")

# take values for factor A = 1
ftr <- datavar9$A == 1
points(datavar9$B[ftr], datavar9$Y[ftr], col = "blue")

# 2.2 - Make ANOVA for 2 factors: A, B => Y

# collect the dimension
dim(datavar9)

# collect the names of the columns
names(datavar9)

# take the interaction model
qab <- lm(Y ~ as.factor(A) * as.factor(B), data = datavar9)

# plot model of interaction
qc <- qab$coefficients
# for A = 0
x <- c(0, 1)
y <- qc[1] + qc[3] * x
points((x + 1), y, "l", col = "red", lwd = 3)
# for A = 1
x <- c(0, 1)
y <- qc[1] + qc[2] + qc[3] * x + qc[4] * x
points((x + 1), y, "l", col = "blue", lwd = 3)
dev.off()

# anova for interaction model
q <- anova(qab)
print(q)

# take the additive models
qadd <- lm(Y ~ as.factor(A) + as.factor(B), data = datavar9)
qa <- lm(Y ~ as.factor(A), data = datavar9)
qb <- lm(Y ~ as.factor(B), data = datavar9)
q0 <- lm(Y ~ 1, data = datavar9)

# make 2-factor anova
# main model ==> nij = u + aA_i + aB_j + aAB_ij

# hypothesis for absence interaction ==> H0: aAB_ij = 0
anova(qadd, qab)

# hypothesis for absence effect of factor B ==> HA: aA_i = 0
anova(qa, qab)

# hypothesis for absence effect of factor A ==> HB: aB_j = 0
anova(qb, qab)

# hypothesis of independency ==> H0, HA, HB simultaneously
anova(q0, qab)


#############################
#           TASK 03         #
#############################

# 3.1- Show the dynamics of feature Y in form of trajectories

# Show time effect in feature Y [black lines]
plot(datavar9$Visit, datavar9$Y, "n",
        main = "Dynamics of feature Y in time",
        xlab = "Time [Visit]", ylab = "Feature Y")
for(i in datavar9$ID){
    ftr <- datavar9$ID == i
    d9_idx <- datavar9[ftr, ]
    points(d9_idx$Visit, d9_idx$Y, "l")
}

# check the tendency [red line]
m <- NULL
for(i in as.numeric(lvVisit)){
    m <- c(m, mean(datavar9$Y[datavar9$Visit == i]))
}

points(as.numeric(lvVisit), m, "l", col = "red", lwd = 4)
dev.off()

# 3.2- Evaluate the presence of factor A in the values of Y

# plot effect of each value of A in Y [belong the time]
for(ida in lvA){
    d9a <- datavar9[datavar9$A == ida, ]
    nnm <- paste0()
    pdf(paste0("A = ", ida), width = 9, height = 6)
    plot(d9a$Visit, d9a$Y, "n",
        main = paste0("Trajectories for A = ", ida),
        xlab = "Visit", ylab = "Y")
    for(i in d9a$ID){
        ftr <- d9a$ID == i
        d9ai <- d9a[ftr, ]
        if(ida == 0){
            points(d9ai$Visit, d9ai$Y, "l", col = "blue")
        }else if (ida == 1) {
           points(d9ai$Visit, d9ai$Y, "l", col = "green")
        }
    }

    # plot the tendency line
    m <- NULL
    for(i in as.numeric(lvVisit)){
        m <- c(m, mean(d9a$Y[d9a$Visit == i]))
    }
    points(as.numeric(lvVisit), m, "l", col = "red", lwd = 4)
    dev.off()
}


# 3.3- Evaluate the correlations between Y values for different times
# for each value of A without consider influence of factor B

# function to structure data
f_tab <- function(dd){
  id <- unique(dd$ID)
  nd <- names(dd)
  lvVis <- levels(as.factor(dd$Visit))
  ftb <- matrix(nrow = length(id), ncol = 2 + 2 * length(lvVis))
  nms <- c("ID", "A")
  for (j in lvVis){
    nms <- c(nms, c(paste0("B_V", j), paste0("Y_V", j)))
  }
  colnames(ftb) <- nms
  for (ii in 1:length(id)){
    i <- id[ii]
    dd.i <- dd[dd$ID == i, ]

    if (dim(dd.i)[1] > 0){
      rw <- c(i, dd.i$A[1])
      for (k in lvVis){
        dd.ik <- dd.i[dd.i$Visit == k, ]

        if(dim(dd.ik)[1] > 0){
          rw <- c(rw, c(dd.ik$B, dd.ik$Y))
        } else {
          rw <- c(rw, c(NA, NA))
        }
      }
    }
    ftb[ii, ] <- rw
  }
  return(as.data.frame(ftb))
}

# make new data structure for time series belong columns
dd9 <- f_tab(datavar9)
print(head(dd9))
print("dimensions f-tab: ")
print(dim(dd9))

# correlation without effect of factors
nms <- paste0("Y_V", lvVisit)
# covariation for all data
varr_all <- var(dd9[, names(dd9)%in%nms])
print(varr_all)
# correlation for all data
corr_all <- cor(dd9[, names(dd9)%in%nms])
print(corr_all)

# correlation of Y in each moment respect factor A [not B]
stat_factorA <- list()

for(aix in lvA){
    vnm <- paste0("varr_A", aix)
    vnm1 <- paste0("corr_A", aix)
    d9_aix <- datavar9[datavar9$A == aix, ]
    dd9_aix <- f_tab(d9_aix)
    nms <- paste0("Y_V", lvVisit)
    stat_factorA[[paste0("A=", aix)]] <- list(
        varr = var(dd9_aix[, names(dd9_aix)%in%nms]),
        corr = cor(dd9_aix[, names(dd9_aix)%in%nms])
    )
}

print(stat_factorA)


#############################
#           TASK 04         #
#############################

# check that each different id, defines window of times
# for different values of factors A and B
print(head(datavar9, n = 20))
print(lvA)
print(lvB)

# 4.1- build estimators for mean values for Y
# for different values of factors A and B

# data structure for each value of factor A and B
d9 <- list()

for(va in lvA){
  for(vb in lvB){
    d9[[paste0("a=", va, "_b=", vb)]]  <- 
          datavar9[datavar9$A == va & datavar9$B == vb, ]
  }
}

print(head(d9[["a=0_b=1"]]))

# check each pair of values in each factor A and B
print(names(d9))

# generate the linear mixed effect model for each pair values A, B
m0 <- NULL
for(i in names(d9)){
  m0[[i]] <- lme(Y ~ as.factor(Visit), random = ~1|ID, data = d9[[i]])
}

# check the linear-mixed-models (example for A = 1 and B = 1)
print(m0[["a=1_b=1"]])

# build semivariogram
v0 <- NULL
par(mfrow = c(2, 2))

for(i in names(d9)){
  v0[[i]] <- Variogram(m0[[i]], form =~ Visit | ID)
  plot(v0[[i]]$dist, v0[[i]]$variog, "l", 
        main = paste0("Semivariogram for pair ", i),
        xlab = "distance", ylab = "variogram")
}

# check the table of semivariogram distances and variogram values
print(v0)

# 4.2- evaluate the dependency between correlation of values
# of process centred in time for each pair value of A, B factors

# calculate variance and standard deviation
s2 <- NULL
s1 <- NULL

for(idx in names(d9)){

  # calculate variance
  for(i in lvVisit){
      s2[[idx]] <- c(s2[[idx]], var(d9[[idx]]$Y[d9[[idx]]$Visit == i]))
  }

  # calculate the standard deviation
  s1[[idx]] <- sqrt(s2[[idx]])
}

# define dimensionality
dimvar <- length(lvVisit)

# build correlation and variation values of centred process in time
# get variograms for each pair of factor values A, B
vrg <- NULL
for(idx in names(d9)){
  vrg[[idx]] <- v0[[idx]]$variog
}

# define matrices of covariation and correlation
stat_semivar <- list()

for(idx in names(d9)){
  vrg_aux <- vrg[[idx]]
  s1_aux <- s1[[idx]]
  s2_aux <- s2[[idx]]
  varr <- matrix(nrow = dimvar, ncol = dimvar)
  corr <- matrix(nrow = dimvar, ncol = dimvar)
  
  for(i in 1:dimvar){
    for(j in 1:dimvar){
        if(i != j){
            # formula for correlation
            corr[i, j] <- 1 - vrg_aux[abs(i - j)]
            # formula for variance = si * sj * corr(i, j)
            varr[i, j] <- s1_aux[i] * s1_aux[j] * (1 - vrg_aux[abs(i - j)])
        }else{
            corr[i, j] <- 1
            varr[i, j] <- s2_aux[i]
        }
    }
  }

  stat_semivar[[idx]] <- list(
      corr = corr,
      varr = varr
  )
}

print(stat_semivar)


#############################
#           TASK 05         #
#############################

# Build linear mixed effect model (lme) with simple individual effect
q51 <- lme(fixed = Y ~ as.factor(Visit), random =~ 1|ID,
                data = datavar9, method = "ML")
print(q51)

# Linear mixed model with linear dependency between individual effect and time
q52 <- lme(fixed = Y ~ as.factor(Visit), random =~ (1 + Visit)|ID,
                data = datavar9, method = "ML")
print(q52)

# Using criteria cAIC to check what is the best model
q51_cAIC <- cAIC(q51)
q52_cAIC <- cAIC(q52)

model1 <- "y ~ Visit + (1|ID)"
model2 <- "y ~ Visit + (1+Visit|ID)"

# apply cAIC criteria to establish the best model
tblcAIC <- data.frame(models = c(model1, model2),
                        cAICvalues = c(q51_cAIC$caic, q52_cAIC$caic))

# table of cAIC values
print(tblcAIC)


#############################
#           TASK 06         #
#############################

# Considering the model with simple individual effect, build mixed
# model with dependency of each observed features A and B
# with respect the time of observations

# mixed model: (Y|A, B, Visit)
# model: nijk = u + aA_i + aB_j + aT_k +
#               aA_i * aB_j + aA_i * aT_k + aB_j * aT_k + aABT_ijk
q61 <- lme(fixed = Y ~ as.factor(A) * as.factor(B) * as.factor(Visit),
              random =~ 1|ID, data = datavar9, method = "ML")
print(q61)

# For each value of factor A:
# 1- check the hypothesis of aditivity influence of the factor B
# and the time observation
# 2- check the hypothesis of absence of influence of the factor B
# and the time observation

q6_intBT <- NULL; q6_addBT <- NULL; q6_addB <- NULL; q6_addT <- NULL
anova61 <- NULL; anova62 <- NULL; anova63 <- NULL
pvalues61 <- NULL; pvalues62 <- NULL; pvalues63 <- NULL
alpha <- 0.05

# for each value of factor A
for(idx in lvA){
  # interaction model for factor B and time (for each value of A)
  # model: njk = u + aB_j + aT_k + aBT_jk
  q6_intBT[[paste0("A=", idx)]] <- lme(
                  fixed = Y ~ as.factor(B) * as.factor(Visit),
                  random =~ 1|ID, data = datavar9[datavar9$A == idx, ],
                  method = "ML")
  
  # additive model for factor B and time (for each value of A)
  # model: njk = u + aB_j + aT_k
  q6_addBT[[paste0("A=", idx)]] <- lme(
                  fixed = Y ~ as.factor(B) + as.factor(Visit),
                  random =~ 1|ID, data = datavar9[datavar9$A == idx, ],
                  method = "ML")

  # additive model for time (for each value of A)
  # model: njk = u + aT_k
  q6_addT[[paste0("A=", idx)]] <- lme(
                  fixed = Y ~ as.factor(Visit),
                  random =~ 1|ID, data = datavar9[datavar9$A == idx, ],
                  method = "ML")

  # additive model for factor B (for each value of A)
  # model: njk = u + aB_j
  q6_addB[[paste0("A=", idx)]] <- lme(
                  fixed = Y ~ as.factor(B),
                  random =~ 1|ID, data = datavar9[datavar9$A == idx, ],
                  method = "ML")

  # interaction model: njk = u + aB_j + aT_k + bBT_jk

  # presence of additive influence between factors B and Visit
  # ausence of interaction ==> H0: bBT_jk = 0
  anova61[[paste0("A=", idx)]] <- anova(
            q6_addBT[[paste0("A=", idx)]], q6_intBT[[paste0("A=", idx)]])

  # ausence of effect factor B ==> H0: bBT_jk = 0 and bB_j = 0
  anova62[[paste0("A=", idx)]] <- anova(
            q6_addT[[paste0("A=", idx)]], q6_intBT[[paste0("A=", idx)]])

  # ausence of time observation ==> H0: bBT_jk = 0 and bT_k = 0
  anova63[[paste0("A=", idx)]] <- anova(
            q6_addB[[paste0("A=", idx)]], q6_intBT[[paste0("A=", idx)]])

  # take p-values for each hypothesis
  pvalues61[[paste0("A=", idx)]] <- anova61[[paste0("A=", idx)]]$"p-value"[2]
  pvalues62[[paste0("A=", idx)]] <- anova62[[paste0("A=", idx)]]$"p-value"[2]
  pvalues63[[paste0("A=", idx)]] <- anova63[[paste0("A=", idx)]]$"p-value"[2]
}

# build table of p-values for each hypothesis
tbl_hipothesis <- function(){
    hypoth <- c("H0: bBT_jk = 0", "H0: bBT_jk = 0 and bB_j = 0",
                  "H0: bBT_jk = 0 and bT_k = 0")
    lh <- length(hypoth)
    la <- length(lvA)
    out <- matrix(nrow = lh, ncol = la)

    for(idx in lvA){
      out[1, as.numeric(idx) + 1] <- pvalues61[[paste0("A=", idx)]]
      out[2, as.numeric(idx) + 1] <- pvalues62[[paste0("A=", idx)]]
      out[3, as.numeric(idx) + 1] <- pvalues63[[paste0("A=", idx)]]
    }

    rownames(out) <- hypoth
    colnames(out) <- lvA
    out
}

# check the p-values to make the hypothesis inference
print(tbl_hipothesis())

# check the results and compare with AIC criteria
for(aix in lvA){
  print(cAIC(q6_intBT[[paste0("A=", aix)]])$caic)
  print(cAIC(q6_addBT[[paste0("A=", aix)]])$caic)
  print(cAIC(q6_addB[[paste0("A=", aix)]])$caic)
  print(cAIC(q6_addT[[paste0("A=", aix)]])$caic)
}

# model without factor A
q62 <- lme(fixed = Y ~ as.factor(B) * as.factor(Visit),
                random =~ 1|ID, data = datavar9, method = "ML")
print(q62)

# study the influence of factor A
anov_addA <- anova(q62, q61)
print(anov_addA)
print(paste("p-value = ", anov_addA$"p-value"[2]))

# check cAIC values
q62_aic <- cAIC(q62)$caic
print(q62_aic)
q61_aic <- cAIC(q61)$caic
print(q61_aic)


#############################
#           TASK 08         #
#############################

# calculate Visit^2
datavar9$Visit2 <- datavar9$Visit^2
print(head(datavar9))

# build model mixed-analysis for 2nd polynomial order in time
# model: nijkl = u + aA_i + aB_j + aT_k + aT^2_l +
#                 aAB_ij + aAT_ik + aBT_jk + aAT^2_il + aBT^2_jl +
#                 aABT_ijk + aABT^2_ijl
q8.abt_square <- lme(fixed = Y ~ as.factor(A) * as.factor(B) * Visit +
                  as.factor(A) * as.factor(B) * Visit2,
                  random =~ 1|ID, data = datavar9, method = "ML")
print(q8.abt_square)

# check the significance for each coefficient
anova(q8.abt_square)

# build model mixed-analysis for linear dependency
# model: nijk = u + aA_i + aB_j + aT_j + aAB_ij + aAT_ik + aBT_jk
q8.abt_linear <- lme(fixed = Y ~ as.factor(A) * as.factor(B) * Visit,
                  random =~ 1|ID, data = datavar9, method = "ML")
print(q8.abt_linear)

# proof of hypothesis for linear dependency
# H_T^2: bABT^2_ijl = 0 and bAT^2_il = 0 and bBT^2_jl = 0 and bT^2_l = 0
q8_anova1 <- anova.lme(q8.abt_linear, q8.abt_square)
print(q8_anova1)

# check the p-value for H_T^2
q8_pvalue1 <- q8_anova1$"p-value"[2]
print(q8_pvalue1)

# check the AIC values
q8square_caic <- cAIC(q8.abt_square)$caic
print(q8square_caic)
q8linear_caic <- cAIC(q8.abt_linear)$caic
print(q8linear_caic)


#############################
#           TASK 09         #
#############################

form <- list()
func <- NULL

# build the all possibilities of models

# interaction for 3 variables
func[1] <- "Y ~ A*B*Visit"
form[[1]] <- formula(Y ~ as.factor(A) * as.factor(B) * Visit)

# interaction for 2 variables
func[2] <- "Y ~ A*B + A*Visit + B*Visit"
form[[2]] <- formula(Y ~ as.factor(A) * as.factor(B) + as.factor(A) * Visit +
                as.factor(B) * Visit)
func[3] <- "Y ~ A*B + A*Visit"
form[[3]] <- formula(Y ~ as.factor(A) * as.factor(B) + as.factor(A) * Visit)
func[4] <- "Y ~ A*B + B*Visit"
form[[4]] <- formula(Y ~ as.factor(A) * as.factor(B) + as.factor(B) * Visit)
func[5] <- "Y ~ A*Visit + B*Visit"
form[[5]] <- formula(Y ~ as.factor(A) * Visit + as.factor(B) * Visit)
func[6] <- "Y ~ A*Visit + B"
form[[6]] <- formula(Y ~ as.factor(A) * Visit + as.factor(B))
func[7] <- "Y ~ A + B*Visit"
form[[7]] <- formula(Y ~ as.factor(A) + as.factor(B) * Visit)
func[8] <- "Y ~ A*B + Visit"
form[[8]] <- formula(Y ~ as.factor(A) * as.factor(B) + Visit)
func[9] <- "Y ~ A*Visit"
form[[9]] <- formula(Y ~ as.factor(A) * Visit)
func[10] <- "Y ~ B*Visit"
form[[10]] <- formula(Y ~ as.factor(B) * Visit)

# additive models
func[11] <- "Y ~ A + B + Visit"
form[[11]] <- formula(Y ~ as.factor(A) + as.factor(B) + Visit)
func[12] <- "Y ~ A + Visit"
form[[12]] <- formula(Y ~ as.factor(A) + Visit)
func[13] <- "Y ~ B + Visit"
form[[13]] <- formula(Y ~ as.factor(B) + Visit)
func[14] <- "Y ~ Visit"
form[[14]] <- formula(Y ~ Visit)

# models with Visit^2 for 3 interactions
func[15] <- "Y ~ A*B*Visit + A*B*Visit2"
form[[15]] <- formula(Y ~ as.factor(A) * as.factor(B) * Visit +
                  as.factor(A) * as.factor(B) * Visit2)
func[16] <- "Y ~ A*B*Visit + A*Visit2 + B*Visit2"
form[[16]] <- formula(Y ~ as.factor(A) * as.factor(B) * Visit +
                  as.factor(A) * Visit2 + as.factor(B) * Visit2)
func[17] <- "Y ~ A*B*Visit + A*Visit2"
form[[17]] <- formula(Y ~ as.factor(A) * as.factor(B) * Visit +
                  as.factor(A) * Visit)
func[18] <- "Y ~ A*B*Visit + B*Visit2"
form[[18]] <- formula(Y ~ as.factor(A) * as.factor(B) * Visit +
                  as.factor(B) * Visit2)

# models with Visit^2 for 2 interactions
func[19] <- "Y ~ A*B + A*Visit + A*Visit2 + B*Visit+ B*Visit2"
form[[19]] <- formula(Y ~ as.factor(A) * as.factor(B) +
                  as.factor(A) * Visit + as.factor(A) * Visit2 +
                  as.factor(B) * Visit + as.factor(B) * Visit2)
func[20] <- "Y ~ A*B + A*Visit + A*Visit2"
form[[20]] <- formula(Y ~ as.factor(A) * as.factor(B) + 
                  as.factor(A) * Visit + as.factor(A) * Visit2)
func[21] <- "Y ~ A*B + B*Visit + B*Visit2"
form[[21]] <- formula(Y ~ as.factor(A) * as.factor(B) +
                  as.factor(B) * Visit + as.factor(B) * Visit2)
func[22] <- "Y ~ A*B + Visit + Visit2"
form[[22]] <- formula(Y ~ as.factor(A) * as.factor(B) +
                  Visit + Visit2)
func[23] <- "Y ~ A*Visit + A*Visit2 + B*Visit + B*Visit2"
form[[23]] <- formula(Y ~ as.factor(A) * Visit + as.factor(A) * Visit2 +
                  as.factor(B) * Visit + as.factor(B) * Visit2)
func[24] <- "Y ~ A*Visit + A*Visit2"
form[[24]] <- formula(Y ~ as.factor(A) * Visit + as.factor(A) * Visit2)
func[25] <- "Y ~ B*Visit + B*Visit2"
form[[25]] <- formula(Y ~ as.factor(B) * Visit + as.factor(B) * Visit2)
func[26] <- "Y ~ A*Visit2"
form[[26]] <- formula(Y ~ as.factor(A) * Visit2)
func[27] <- "Y ~ B*Visit2"
form[[27]] <- formula(Y ~ as.factor(B) * Visit2)

# models with Visit^2 for additive
func[28] <- "Y ~ A + B + Visit + Visit2"
form[[28]] <- formula(Y ~ as.factor(A) + as.factor(B) + Visit + Visit2)
func[29] <- "Y ~ A + Visit + Visit2"
form[[29]] <- formula(Y ~ as.factor(A) + Visit + Visit2)
func[30] <- "Y ~ B + Visit + Visit2"
form[[30]] <- formula(Y ~ as.factor(B) + Visit + Visit2)
func[31] <- "Y ~ Visit + Visit2"
form[[31]] <- formula(Y ~ Visit + Visit2)

