
### Comparison study: LEA::lfmm() and LEA::lfmm2()
### Note: only statistical performances are measured
### lfmm2() is hundred times faster than lfmm()
library(LEA)
library(qvalue)


# Simulation with 100 target loci, with effect sizes ranging between -10 and 10 
# n = 100 individuals
# L = 2000 loci
n = 100 
L = 2000
# Number of hidden factors
K = 3

# measures recorded
rsq <- NULL # Coeff. of determination: squared correlation between X (environment) and U (factors)

pow_lfmm1 = NULL # Estimated power values (lfmm)
fdr_lfmm1 = NULL # Estimated FDR values (lfmm, expected FDR = 5%)

pow_lfmm2 = NULL # Estimated power values (lfmm2)
fdr_lfmm2 = NULL # Estimated FDR values (lfmm2, expected FDR = 5%)

for (i in 1:100){
  # Environmental variable
   x <- as.matrix(rnorm(n)) 
  
  # Create hidden factors and loadings  
    prop <- runif(3, -1.25, 1.25)
    U <- t(tcrossprod(as.matrix(prop), x)) + matrix(rnorm(K*n), ncol = 3)
    sm <- summary(lm(x ~ U))
    rsq[i] <- sm$r.squared #coeff. of determination

    V <- matrix(rnorm(K*L), ncol = 3) #loadings

  # Define effect size of target loci
    B <- rep(0, L) 
    target <- sample(1:L, n) # target loci
    B[target] <- runif(100, -10, +10) # effect sizes

  # Simulate a matrix of haploid genotypes 
  # Simulation performed with the generative model of LFMM
    Y <- tcrossprod(as.matrix(x), B) + tcrossprod(U, V) + 
          matrix(rnorm(n*L, sd = .5), nrow = n)
    Y <- matrix(as.numeric(Y > 0), ncol = L)

  # Fit an LFMM with K = 3 factors 
  # with lfmm2
    mod2 <- lfmm2(input = Y, env = x, K = 3)

  # Compute P-values for all loci
    pv <- lfmm2.test(object = mod2, input = Y, env = x, linear = TRUE)
  # Use FDR control to provide a list of candidate loci
    qv <- qvalue(p=pv$pvalues, fdr = 0.05)$signif
  # Estimate power and FDR   
    pow_lfmm2[i] <- mean((target %in% which(qv))) 
    fdr_lfmm2[i] <- mean(!(which(qv) %in% target))

  # Fit an LFMM with K = 3 factors 
  # with lfmm version 1
    write.lfmm(output.file = "y.lfmm", Y)
    write.env(x, output.file = "x.env")

    project = lfmm( "y.lfmm", 
                "x.env", 
                K = 3, 
                repetitions = 5, 
                project = "new")

  # get adjusted p-values using all runs
    pv_lfmm1 <- lfmm.pvalues(project, K = 3)
    qv <- qvalue(p=pv_lfmm1$pv, fdr = 0.05)$signif
  # Estimate power and FDR   
    pow_lfmm1[i] <- mean((target %in% which(qv))) 
    fdr_lfmm1[i] <- mean(!(which(qv) %in% target))

  # Clean-up lfmm files 
    remove.lfmmProject("y_x.lfmmProject")
    system("rm x*")
    system("rm y*")
}


## Create Figure
par(mfrow = c(1,2))
t.test(fdr_lfmm1, fdr_lfmm2)
boxplot(fdr_lfmm1, fdr_lfmm2, 
        names = c("LFMM1","LFMM2"), 
        col = c("orange", "darkblue"),
        ylab = "False Discovery Rate", las = 1)

plot(rsq, pow_lfmm1, pch = 16, cex = .7, col = "orange", 
     ylab = "Power",
     xlab = "Coefficient of Determination")

mod <- loess(pow_lfmm1~rsq)
points(rsq, fitted(mod), col = "brown", pch = 19 , cex = .6)
points(rsq, pow_lfmm2, pch = 16, cex = .7, col = "lightblue")
mod <- loess(pow_lfmm2~rsq)
points(rsq, fitted(mod), col = "darkblue", pch = 19 , cex = .6)
legend(x = 0.05, y = 0.7, pch = 19, 
       col = c("orange","lightblue"), 
       legend = c("LFMM1","LFMM2"))
