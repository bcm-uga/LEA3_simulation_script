devtools::install_github("bcm-uga/LEA")
library(LEA)


# Simulation study: two-population split and admixture
# Number of individuals in pop samples 
n <-  300 
#pop 1
n1 <- 100 
#pop 2 (admixed pop) 
n2 <- 100 
#pop 3
n3 <-  n - n1 - n2 

# Number of loci in simulated genotypes 
LL = seq(1000, 10000, by = 1000)

# Proportions of missing values in simulated genotypes
pp = seq(0.1, 0.9, by = 0.1)

# Results matrix
res = NULL

for (Le in LL){
  for (prop in pp){
    for (i in 1:20){
      
      L = Le 
      X <- matrix(NA, nrow = n, ncol = L)
      # Drift in F-model
      F = runif(1, 0.02, 0.3) 
      # Admixture coefficient
      alpha <- runif(1,0.05, 0.5) 

      # simulation of a 2-population F-model with admixture
      for (l in 1:L){
          p <- runif(1)
          p1 <- rbeta(1, p*(1-F)/F, (1-p)*(1-F)/F)
          p3 <- rbeta(1, p*(1-F)/F, (1-p)*(1-F)/F)
          p2 <- alpha*p1 + (1-alpha)*p3
          X[1:n1,l] <- rbinom(n1, size = 1, prob = p1)
          X[(n1+1):(n1+n2),l] <- rbinom(n2, size = 1, prob = p2)
          X[(n1+n2+1):n,l] <- rbinom(n3, size = 1, prob = p3)
      }
      
      # Filtering SNPs
      lst.unique <- which( apply(X, 2, FUN = function(x) length(unique(x))) == 1)
      X <- X[,-lst.unique]
      L <- ncol(X)
      
      # Create missing genotypes
      Xo <- X
      X[sample(1:(n*L), round((n*L)*prop))] = 9
      write.lfmm(output.file = "x.lfmm", X)
      
      #Run snmf (10 times)
      project.snmf = snmf("x.lfmm",
                    K = 2, 
                    entropy = TRUE, 
                    repetitions = 10,
                    ploidy = 1,
                    alpha = 100,
                    project = "new")
      
      # get the cross-entropy of the 10 runs for K = 2
      ce <- cross.entropy(project.snmf, K = 2)

      # select the run with the lowest cross-entropy for K = 2
      best <- which.min(ce)
      
      # Impute missing values from the snmf model
      impute(project.snmf, "x.lfmm", method = 'mode', K = 2, run = best)

      # Compare with true genotypes
      # Proportion of correct imputation results:
      p.succ <- mean(Xo[X == 9] == read.lfmm("x.lfmm_imputed.lfmm")[X == 9])
      
      #clean up
      remove.snmfProject("x.snmfProject")
      system(command = "rm x*")

      #Add results to res
      res <-  rbind(res,c(L,
                  prop,
                  F,
                  alpha,
                  p.succ))
    }
  }
}

# Column names
# L: number of loci
# p.miss: prop missing values
# F: drift in F-model
# alpha: admixture coefficient
# p.succ: proportion of correctly reconstructed genotypes

colnames(res) = c("L","p.miss", "F", "alpha", "p.succ")

par(mfrow = c(1,2))
plot(as.factor(res[,"p.miss"]), 
     res[,"p.succ"], 
     xlab = "Rate of missing genotypes",
     ylab = "Accuracy",    
     pch = 19, cex = .5, col = "grey")

plot(res[,"F"], 
     res[,"p.succ"], 
     xlab = "F-value",
     ylab = "Accuracy",    
     pch = 19, cex = .5, col = "grey")

abline(mod, col = "blue", lwd = 3, lty = 2)


#### Imputation on real data 
# Arabidopsis data from Atwell et al. 2010

devtools::install_github("bcm-uga/naturalgwas")
library(naturalgwas)
data(A.thaliana)
Xo <- A.thaliana$genotype
chrpos <- A.thaliana$chrpos 
n = nrow(Xo)
L = ncol(Xo)

res_arabido <- NULL

for (prop in rep(seq(0.2, 0.8,by = 0.1), each = 2)){

  X <-  Xo
  X[sample(1:(n*L), round((n*L)*prop))] <- 9

  write.lfmm(output.file = "x.lfmm", X)

  project.snmf = snmf("x.lfmm",
                    K = 6, 
                    entropy = TRUE, 
                    repetitions = 5,
                    ploidy = 1,
                    alpha = 10,
                    project = "new")

  # get the cross-entropy of the 5 runs for K = 6
  ce = cross.entropy(project.snmf, K = 6)
  
  # select the run with the lowest cross-entropy for K = 6
  best = which.min(ce)
  
  # Imputation based on the snmf model
  impute(project.snmf, "x.lfmm", method = 'mode', K = 6, run = best)

  # Compare with true genotypes
  # Proportion of correct imputation results:
  p.succ = mean(Xo[X == 9] == read.lfmm("x.lfmm_imputed.lfmm")[X == 9])

  # record results and clean-up
  res_arabido <- rbind(res_arabido, c(prop, p.succ))
  remove.snmfProject("x.snmfProject")
  system(command = "rm x*")
}

plot(res_arabido[,1], 
     res_arabido[,2],
     xlab = "Missing",
     ylab = "Correct imputation")

# proportions of correctly reconstructed genotypes
mean(res_arabido[,2])
sd(res_arabido[,2])


x1 = res_arabido[,1] 
x2 = res_arabido[,1]^2
x3 = res_arabido[,1]^3
mod <- lm(res_arabido[,2]  ~  x1  + x3)
summary(mod)
mod


