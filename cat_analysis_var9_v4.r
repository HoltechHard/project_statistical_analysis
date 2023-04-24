#######################################################
#      TASK 01 - STATISTIC ANALYSIS IN BIG DATA       #
#######################################################
#                      VERSION 4.0
# Comment: Ctrl + enter --> this is a sequence of commands
# to run selected lines of code in R

# Topic: Categorical data analysis
# Enunciate: Variant 09
# The table of data refers to a sample of subjects randomly selected
# for an Italian study on the relation between income and whether
# one posseses a travel credit card.
# At each level of annual income (in mill. lira), the table indicates
# the number of subjects sampled and the number of them possesing at
# least one travel credit card.

# read dataset
data <- read.csv("practice01_stat/Credit.csv", sep = " ")

# structure the data
# variables
# income: salary per person
# numcase: number of total sampled persons
# card: number of persons which have least 1 card
data <- data.frame(data)
colnames(data) <- c("income", "numcase", "card")
print("*** Dataset - Credit Card Analysis ***")
print(data)
print("Dimensionality of dataset: ")
print(dim(data))

######################
#       TASK 01      #
######################

# Divide income levels into 3 groups and graphically
# depict the frequency of having at least one card at
# different income levels
# income < 40  ; 40 <= income < 70  ; income >= 70

# number of samples
samples <- dim(data)[1]

# establish grouping criteria
for(i in 1:samples)
{
    if(data$income[i] < 40)
    {
        data$level[i] <- "low"
    }else if (data$income[i] < 70) {
       data$level[i] <- "medium"
    }else {
       data$level[i] <- "high"
    }
}

# final printed for dataset
# check the new column income levels
print(data)

# redefine 2 columns
for(i in 1:samples){
    data$yescard[i] <- data$card[i]
    data$nocard[i] <- data$numcase[i] - data$card[i]
}

# redefine levels
data$level[which(data$level == "low")] <- 0
data$level[which(data$level == "medium")] <- 1
data$level[which(data$level == "high")] <- 2

print(data)

# build the table of factors

# levels of income
income_lev <- levels(as.factor(data$level))

# levels of credit card possession
possession_lev <- c(0, 1)       # values 0: not |  1: yes

# build table of factors
# factor 1: level of income
# factor 2: level of possession of least 1 credit card

f_table <- function(data){
    a <- levels(as.factor(data$level))
    b <- c(0, 1)
    la <- length(a)
    lb <- length(b)
    out <- matrix(nrow = lb, ncol = la)

    for(i in 1:la){
        ftr <- data[, 4] == a[i]
        for(j in b){
            out[j + 1, i] <- sum(data[ftr, ncol(data) - j])
        }
    }

    rownames(out) <- b
    colnames(out) <- a
    out
}

# final table which grouping levels of income
tablitza <- f_table(data)
print(tablitza)

# check total of values in tablitza
sum(tablitza)

# define structure
tt <- structure(c(39, 7, 24, 14, 6, 10), dim = c(6L))
print(tt)

# barplot: frequencies of least 1 card x income levels
bb <- barplot(tablitza, main = "Frequency of credit card possession",
            xlab = "Income levels", ylab = "Frequency",
            col = c("gray", "black"), ylim = c(0, 42),
            names.arg = c("low", "medium", "high"),
            beside = TRUE)
text(bb, tt + 1, tt, cex = 2)
legend("topright", legend = c("not card", "yes card"),
        col = c("gray", "black"),
                pch = 0, cex = 2, title = "credit cards")


######################
#       TASK 02      #
######################

# Check the dependence between the presence of at least one
# card on the level of income (grouping from point 1)
# by methods of classical categorical (Chi-square)

# print summary-table
print(tablitza)

# horizontal => presence of card (0, 1)
# vertical => income levels (0, 1, 2)

#####  Check the dependency between the categorical variables  #####

# test chi square
# H0: no exists relationships on the categorical variables in the population
# Ha: exists significant association / dependence between categorical variables
# variables are significantly contributing to the
# variability of response variable
# p-value > 0 ==> accept H0
# p-value < 0 ==> accept Ha

# classic chi-square
csq <- chisq.test(tablitza, correct = F)
print(csq)

# p-value
alpha <- 0.05
pv.csq <- csq$p.value
print(paste("p-value = ", pv.csq))
# p-value = 0.00124 < 0.05 ==> reject H0, thus have significant diff.

# define criteria
phi_classic <- as.numeric(pv.csq < alpha)

# define function to get hypothesis inference
hypothesis_proof <- function(phi_xi, pvalue){
    if(phi_xi == 0){
        print(paste(phi_xi, " ===> Accept H0, because ", pvalue, " > ", alpha))
        print("variables don't have significant impact in variability of data")
    }else if(phi_xi == 1){
        print(paste(phi_xi, " ===> Reject H0, because ", pvalue, " <= ", alpha))
        print("variables contribute significantly in variability of data")
    }
}

# inference for classical chi-square test
hypothesis_proof(phi_classic, pv.csq)

# simulation chi-square
csq2 <- chisq.test(tablitza, correct = F, simulate.p.value = T)
print(csq2)

pv.csq2 <- csq2$p.value
print(paste("p-value = ", pv.csq2))

phi_simulated <- as.numeric(pv.csq2 < alpha)

# inference for classical chi-square test
hypothesis_proof(phi_simulated, pv.csq2)


######################
#       TASK 03      #
######################

# Check the dependence of the presence of at least one card on the
# level of income without grouping using glm (model logistic regresion).
# Show graphically the resulting dependence of the probability of
# having at least one card on the level income


# build table with each individual case of 100 samples
glm_table <- function(data){
    out <- data.frame(id = c(1:sum(data$numcase)))
    v1 <- NULL      # column for incomes
    v2 <- NULL      # column for card possession

    for(i in 1:samples)
    {
        lim1 <- data$yescard[i]
        lim2 <- data$nocard[i]

        # fill values for cases that have card = yes
        if(lim1 > 0){
            v1 <- c(v1, array(data$income[i], dim = data$yescard[i]))
            v2 <- c(v2, array(1, dim = data$yescard[i]))
        }

        # fill values for cases that have card = no
        if(lim2 > 0){
            v1 <- c(v1, array(data$income[i], dim = data$nocard[i]))
            v2 <- c(v2, array(0, dim = data$nocard[i]))
        }
    }

    # convert the outputs of each column in numerical data type
    out$v1 <- as.numeric(v1)
    out$v2 <- as.numeric(v2)

    # change the name of the columns in "income" and yescard (1: yes | 0: no)
    names(out)[c(2, 3)] <- names(data)[c(1, 5)]

    return(out)
}

# generate the categorical table with each instance
# for glm analysis
categ_table <- glm_table(data)
print(categ_table)

# create general logistic regression:
# y = f(x), where logit(y) = b0 + b1 * x
# variables
# y: presence or not of least 1 credit card     ==> categorical variable
# x: income (without grouping)                  ==> real variable

# analysis for classic logistic regression
lr_model <- glm(yescard ~ income, data = categ_table, family = "binomial")
lr_ind <- glm(yescard ~ 1, data = categ_table, family = "binomial")

# statistic criteria logarithm relative max-likelihood
llr <- lr_ind$deviance - lr_model$deviance
print(paste("G = ", llr))

# compute the pvalue
llr_pvalue <- pchisq(q = llr, df = 1, lower.tail = FALSE)
print(paste("p-value = ", llr_pvalue))

# apply criteria function and check hypothesis result
phi_llr <- as.numeric(llr_pvalue < alpha)

# inference with relative maximum likelihood and Chi-square
hypothesis_proof(phi_llr, llr_pvalue)

# make anova analysis
# H0: no significant difference between the means of
# the groups (u0 = u1 = ... = uk)
# H1: the variables added contribute significantly to
# the variability of response variable
glm_anova <- anova(lr_ind, lr_model, test = "LRT")
print(glm_anova)

# pvalue = 2.19e-07 *** ==> very high significance

# Conclusion: for that, the impact of categorical variable
# in the logistic regression model is significant


##### variant 2  #####  with probabilities
# analysis which transforms the probability of the response variable
# being a success, resulting in different predicted probabilities for
# the same set of predictor variables.

lr_model2 <- glm(yescard ~ income, data = categ_table,
                    family = binomial("probit"))

lr_ind2 <- glm(yescard ~ 1, data = categ_table,
                    family = binomial("probit"))

glm_anova2 <- anova(lr_ind2, lr_model2, test = "LRT")
print(glm_anova2)

# pvalue = 2.06e-07 *** ==> very high significance


# show graphically the resulting dependence of the
# probability of having at least one card on the
# level income

# check the dependency of logistic regression
print(head(categ_table))

x <- categ_table$income
z <- lr_model$coeff[1] + lr_model$coeff[2] * x
pz <- exp(z) / (1 + exp(z))
print(pz)
plot(x, pz, xlab = "incomes", ylab = "probabilities for card possession",
            main = "Dependency between incomes and prob. to have card", lwd = 3)


######################
#       TASK 04      #
######################

# Check the dependence of the presence of at least one card on the
# level of income with grouping from item (1) using glm.
# Make it using 2 ways: using logistic regression model and
# with poisson model


#################### poisson model ####################
# using frequencies like observations
head(tablitza)

# rebuild data structure to work with frequencies for poisson model

poiss_table <- function(dt){
    at <- rownames(dt)
    bt <- colnames(dt)
    dtpoiss <- data.frame(id = c(1:(length(at) * length(bt))))

    c1 <- NULL
    c2 <- NULL
    c3 <- NULL

    for(i in at){
        for(j in bt){
            c1 <- c(c1, i)
            c2 <- c(c2, j)
            c3 <- c(c3, dt[i, j])
        }
    }

    dtpoiss$c1 <- as.numeric(c1)
    dtpoiss$c2 <- as.numeric(c2)
    dtpoiss$c3 <- as.numeric(c3)
    names(dtpoiss)[c(2, 3, 4)] <- c("possession", "levels", "frequency")

    return(dtpoiss)
}

# new data structure with frequencies
freq_table <- poiss_table(tablitza)
print(freq_table)

# compute the dimension of levels for each categorical variable
nl1 <- length(levels(as.factor(freq_table$possession)))
nl2 <- length(levels(as.factor(freq_table$levels)))

# poisson model based in frequencies
poiss_model <- glm(frequency ~ as.factor(possession) * as.factor(levels),
                    data = freq_table, family = "poisson")

# additive model (with independency)
poiss_ind <- glm(frequency ~ as.factor(possession) + as.factor(levels),
                    data = freq_table, family = "poisson")

# statistic criteria logarithm max-likelihood
llr_poiss <- poiss_ind$deviance - poiss_model$deviance
print(paste("G = ", llr_poiss))

# compute p-value
llrpoiss_pvalue <- pchisq(q = llr_poiss, df = (nl1 - 1) * (nl2 - 1),
                        lower.tail = FALSE)
print(paste("Poisson p-value = ", llrpoiss_pvalue))

# make anova test
poiss_anova <- anova(poiss_ind, poiss_model, test = "LRT")
print(poiss_anova)

# conclusion
phi_poisson <- as.numeric(llr_pvalue < alpha)

# inference for classical chi-square test
hypothesis_proof(phi_poisson, llrpoiss_pvalue)


############## logistic regression ##############
print(freq_table)

group_glm <- function(dtt){
    a <- levels(as.factor(dtt$possession))
    b <- levels(as.factor(dtt$levels))
    la <- length(a)
    lb <- length(b)
    out <- data.frame(id = c(1:sum(freq_table$frequency)))
    v1 <- NULL
    v2 <- NULL

    for(i in 1:la){
        for(j in 1:lb){
            ftr <- freq_table[, 2] == a[i] & freq_table[, 3] == b[j]
            v1 <- c(v1, array(a[i], dim = sum(dtt$frequency[ftr])))
            v2 <- c(v2, array(b[j], dim = sum(dtt$frequency[ftr])))
        }
    }

    out$v1 <- as.numeric(v1)
    out$v2 <- as.numeric(v2)
    names(out)[c(2, 3)] <- names(dtt)[c(2, 3)]

    return(out)
}

# build the new datatable
group_table <- group_glm(freq_table)
print(group_table)

# compute the dimension of levels for each categorical variable
nl3 <- length(levels(as.factor(group_table$levels)))

# build logistic regression model for grouping data
lr_group <- glm(possession ~ as.factor(levels), data = group_table,
                    family = binomial)

# independently model for grouping data
lr_gind <- glm(possession ~ 1, data = group_table, family = binomial)

# statistic criteria logaritm relative max-likelihood
llr2 <- lr_gind$deviance - lr_group$deviance
print(paste("G-logistic = ", llr2))

# calculate p-value
llr2_pvalue <- pchisq(q = llr2, df = nl3 - 1, lower.tail = FALSE)
print(paste("logistic reg. grouping pvalue = ", llr2_pvalue))

# analysis anova
llr2_anova <- anova(lr_gind, lr_group, test = "LRT")
print(llr2_anova)

# define criteria function
phi_llrgroup <- as.numeric(llr2_pvalue < alpha)

# make inference
hypothesis_proof(phi_llrgroup, llr2_pvalue)

