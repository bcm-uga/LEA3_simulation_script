## LEA
## Genome scan for selection with snmf() and lfmm2() 
library(LEA)

## Read data from Atwell et al 2010 
## Arabidopsis thaliana (Scandinavian accessions)
## Genotype
  X = read.lfmm("sweden/sweden.lfmm")
## Genomic map
  pos <- read.table("sweden/pos.txt", h = TRUE)
  pos <- as.numeric(pos[,1])
  chr <- NULL
  ch = 1; for (i in 2:length(pos)){if (pos[i] < pos[i-1]){ch = ch + 1} 
   ; chr = c(chr, ch)}

  col.chr = NULL
  col.chr[chr %% 2 == 1] <- "darkblue"
  col.chr[chr %% 2 == 0] <- "grey"

## Filter SNPs
  lst.unique <- apply(X, 2, function(x) length(unique(x)))
  X <- X[,lst.unique == 2] 
  
## Genome scan based on PCA
  pc <- prcomp(X, scale = TRUE)
  z <- NULL
  Y <- scale(X)
  for (j in 1:205417){
  z[j] <- summary(lm(Y[,j]~ pc$x[,1]))$coeff[2,3]
  }

  lambda <- median(z^2)/0.456 

  pc.values <- pchisq(z^2/lambda, df = 1, lower.tail = FALSE)

## Calibration check
hist(pc.values, col = "lightblue")

## Genome scan based on snmf
  project.snmf <- snmf("sweden/sweden2.lfmm",
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
  
# Compute p-values for each locus with snmf.pvalues()
  p <-snmf.pvalues(project.snmf, entropy = TRUE, ploidy = 1, K = 2, run = best)
  snmf.values <- p$pvalues

## Calibration check & clean-up
  hist(snmf.values, col = "orange")
  remove.snmfProject("sweden/sweden2.snmfProject")
  

## Genome scan based on lfmm2() used with altitude
## load meta data
  meta_sweden <- read.csv2("sweden/meta_sweden.csv")
  x <- as.numeric(meta_sweden$latitude)

  mod2 <- lfmm2(input = Y, env = x, K = 2)

# Compute P-values
  lfmm.values <- lfmm2.test(object = mod2, 
                   input = Y, 
                   env = x, 
                   linear = TRUE)$pvalues

## Calibration check 
  hist(lfmm.values, col = "palegreen")

## Show positions with -log p > 3
  boo <- -log10(snmf.values) > 3 | -log10(pc.values) > 3 | 
          -log10(lfmm.values) > 3
  
xx <- (1:205417)*765/1000000

## Create Figure
par(mfrow = c(3,1))
par(ps = 18)

plot(xx[boo], -log10(lfmm.values[boo]),
     col = col.chr[lst.unique == 2][boo],
     cex = .6, main = "SNMF",
     pch = 19, xlab = " ", ylab = " ")

plot(xx[boo],-log10(pc.values[boo]),
     col = col.chr[lst.unique == 2][boo],
     cex = .6, main = "PCA",
     pch = 19, xlab = " ", ylab = "-log(Pvalue)")

plot(xx[boo],-log10(lfmm.values[boo]), 
     col = col.chr[lst.unique == 2][boo],
     cex = .6, main = "LFMM",
     pch = 19, xlab = "Position (Megabase)", ylab = " ")
