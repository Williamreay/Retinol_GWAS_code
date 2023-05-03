###############################

## Function for defining credible sets using approximated asymptotic Bayes' factors (Wakefield's method)

## William Reay (2020)

################################

Finemapping_abf <- function(SNP, beta, stderr, W, coverage = 0.95) {
  logsum <- function(x, base = exp(1), na.rm = FALSE) {
    xm <- max(x)
    return(xm + log(sum(exp(x - xm), na.rm = na.rm), base = base))
  }
  
  z <- beta / stderr
  V <- stderr^2
  r <- W^2 / (W^2 + V)
  lbf <- 0.5 * (log(1 - r) + (r * z^2))
  denom <- logsum(lbf)
  PP <- exp(lbf - denom)
  PP_dat <- data.frame(SNP, lbf, PP)
  
  ## Define SNPs in 95% credible set
  
  ordering <- order(PP_dat$PP, decreasing = T)
  idx <- which(cumsum(PP_dat$PP[ordering]) > coverage)[1]
  cs <- as.data.frame(SNP[ordering][1:idx])
  cs <- dplyr::rename(cs, "SNP"="SNP[ordering][1:idx]")
  
  Final_df <- merge(cs, PP_dat, by="SNP")
  return(Final_df)
}
