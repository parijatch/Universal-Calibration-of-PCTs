library(plyr)
library(MASS)
library(data.table)
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(reshape2)

## ---- auto-set working directory to this script's folder ----

this_file <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArg <- "--file="
  match <- grep(fileArg, cmdArgs)
  if (length(match) > 0) {
    return(normalizePath(sub(fileArg, "", cmdArgs[match])))
  }
  
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  
  return(NULL)
}

script_path <- this_file()
if (!is.null(script_path)) {
  script_dir <- dirname(script_path)
  setwd(script_dir)
  message("Working directory set to script folder: ", script_dir)
} else {
  message("Could not determine script path; working directory not changed.
          If using RStudio, change working directory to that of source file.")
}

source("../utils.R")

get.pval <- function(d = 10, n = 1e6, nu = 1, rho=0.5,cor.type = "autoreg") {
  set.seed(1234) 
  cor.mat <- get.cor(d,rho=rho, cor.type = cor.type)
  
  # simulate under null (mu = 0)
  X <- mvrnorm(n = n, mu = rep(0, d), Sigma = cor.mat)
  
  # t distribution denominator sqrt(W/nu), W ~ ChiSq(nu)
  denom <- replicate(n, sqrt(rgamma(n = 1, shape = nu/2, rate = 1/2) / nu))
  X <- X / denom
  
  # one-sided p-values
  pval <- 1 - pt(X, df = nu)
  
  # -----------------------------
  # Save to CSV
  # -----------------------------
  dir.create("pval_data", showWarnings = FALSE)
  
  filename <- sprintf("../pval_data/pvals_d%d_nu%g_%s_rho%g.csv", d, nu, cor.type,rho)
  write.csv(pval, file = filename, row.names = FALSE)
  
  message("Saved p-values to: ", filename)
  
  return(pval)
}

get_or_load_pval <- function(d = 10,
                             n = 1e6,
                             nu = 1,
                             rho=0.5,
                             cor.type = "autoreg") {
  
  # ensure folder exists
  dir.create("../pval_data", showWarnings = FALSE)
  
  # full path to CSV file
  filename <- sprintf("../pval_data/pvals_d%d_nu%g_%s_rho%g.csv", d, nu, cor.type,rho)
  
  # ------------------------------------------------------------
  # CASE 1: file already exists → read it
  # ------------------------------------------------------------
  if (file.exists(filename)) {
    message("File exists: ", filename, " — reading p-values from disk.")
    
    pval_df <- read.csv(filename)
    pval_mat <- as.matrix(pval_df)
    
    # Check dimensions
    n_existing <- nrow(pval_mat)
    d_existing <- ncol(pval_mat)
    
    if (n_existing == n && d_existing == d) {
      message("Using existing p-values (n = ", n_existing,
              ", d = ", d_existing, ").")
      return(pval_mat)
      
    } else {
      message("Existing file dims do not match.",
              "\nExisting: n = ", n_existing, ", d = ", d_existing,
              "\nRequired: n = ", n, ", d = ", d,
              "\nRegenerating with get.pval() ...")
      
      pval <- get.pval(d = d, n = n, nu = nu,rho=rho, cor.type = cor.type)
      return(pval)
    }
  }
  
  # ------------------------------------------------------------
  # CASE 2: file does not exist → generate and save new pvals
  # ------------------------------------------------------------
  message("File not found: ", filename,
          "\nGenerating new p-values with get.pval()...")
  
  pval <- get.pval(d = d, n = n, nu = nu, rho=rho, cor.type = cor.type)
  return(pval)
}




# ----- Calibration Line Plot -----
calibration.lineplot <- function(nu.vec = c(1, 30),
                                 alpha.vec = 10^seq(log10(0.001), log10(0.01), length.out = 100),
                                 d = 10, n = 1e6, rho=0.5,
                                 cor.type = "autoreg") {
  par(mfrow = c(1, length(nu.vec)), mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))
  
  for (nu in nu.vec) {
    message(sprintf("Running nu = %f", nu))
    pval<-get_or_load_pval(d=d,n=n,nu=nu,rho=rho,cor.type=cor.type)
    
    rejections <- t(sapply(alpha.vec, function(alpha) {
      p.pareto <- combine.test(pval, method = "Pareto")
      p.cauchy <- combine.test(pval, method = "Cauchy")
      c(mean(p.pareto < alpha), mean(p.cauchy < alpha))
    }))
    
    ratio <- rejections / alpha.vec
    
    plot(1 / alpha.vec, ratio[,1], type = "l", col = "red", ylim=c(0, max(ratio)),
         xlab=TeX("$1/\\alpha$"), ylab=TeX("Empirical Rejection Rate / $\\alpha$"),
         main=TeX(sprintf("ν = %.2f", nu)))
    lines(1 / alpha.vec, ratio[,2], type = "l", col = "blue")
    abline(h = 1, lty = 2)
    legend("bottomright", legend=c("PCT", "CCT"), col=c("red", "blue"), pch=1)
    
  }
  
  mtext(
    TeX(sprintf("Calibration of PCT & CCT under Multivariate-$t_\\nu(0,\\Sigma)$ with %s $\\Sigma$", cor.type)),
    outer = TRUE, cex = 1.5, line = 1
  )
  
  par(mfrow = c(1, 1))
}


# ----- Calibration Heatmaps -----
calibration.heatmaps <- function(alpha.vec = c( 0.01, 0.005,0.001,0.0005,0.0001),
                                 nu.vec = c(1, 5, 15, 30),
                                 d = 10, n = 1e6, rho=0.5,
                                 cor.type = "autoreg") {
  mat.pareto <- matrix(NA, nrow = length(alpha.vec), ncol = length(nu.vec))
  mat.cauchy <- matrix(NA, nrow = length(alpha.vec), ncol = length(nu.vec))
  
  for (i in seq_along(nu.vec)) {
    nu <- nu.vec[i]
    pval<-get_or_load_pval(d=d,n=n,nu=nu, rho=rho, cor.type=cor.type)
    
    for (j in seq_along(alpha.vec)) {
      alpha <- alpha.vec[j]
      mat.pareto[j, i] <- mean(combine.test(pval, method = "Pareto") < alpha) / alpha
      mat.cauchy[j, i] <- mean(combine.test(pval, method = "Cauchy") < alpha) / alpha
    }
  }
  
  df.pareto <- as.data.frame(as.table(mat.pareto))
  
  df.cauchy <- as.data.frame(as.table(mat.cauchy))

  colnames(df.pareto) <- colnames(df.cauchy) <- c("alpha.idx", "nu.idx", "ratio")
  df.pareto$Method <- "Pareto"
  df.cauchy$Method <- "Cauchy"
  df <- rbind(df.pareto, df.cauchy)
  
  df$alpha <- factor(alpha.vec[df$alpha.idx], levels = rev(alpha.vec))
  df$nu <- factor(nu.vec[df$nu.idx])
  
  gg <- ggplot(df, aes(x = nu, y = alpha, fill = ratio)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 1,limits=c(0,2),
                         name = TeX("Empirical Rejection / $\\alpha$")) +
    facet_wrap(~Method) +
    labs(title = TeX(sprintf("Calibration Heatmaps for Multivariate-$t_\\nu(0,\\Sigma)$ with %s $\\Sigma$", cor.type)),
         x = TeX("$\\nu$"), y = TeX("$\\alpha$")) +
    theme_minimal(base_size = 13) +
    theme(
      panel.spacing = grid::unit(1.75, "lines")
    )
  
  
  print(gg)
}


# Calibration line plots
calibration.lineplot(n=1e6,d=10,nu.vec = c(0.5, 1, 30), rho=0.5,cor.type = "autoreg")

# Calibration heatmaps
calibration.heatmaps(n=1e6,d=10,nu.vec = c(0.1,0.5, 1, 5, 15, 30),rho=0.5,
                     cor.type = "autoreg")
