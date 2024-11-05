#################################################################
# Libraries
#################################################################

# Libraries for spatial analysis
library(spdep)
library(terra)
library(sp)
library(tidyterra)
library(Matrix)
library(readxl)

# Libraries for matrix, distribution and math computation
library(mvtnorm)
library(matrixStats)
library(blockmatrix)
library(invgamma)
library(MASS)

# Libraries for JAGS and bayesian inference
library(rjags)
library(R2jags)
library(coda)

# Library for Figure Combination
library(patchwork)

#Library for plotting
library(ggplot2)

#################################################################
# Auxiliary Functions
#################################################################

histogram_plot <- function(data, x_label) {
  ggplot(data.frame(data), aes(x = data)) +
    geom_histogram(binwidth = diff(range(data)) / 20, fill = "grey", color = "grey20", size = 0.5) +
    theme_minimal() +
    theme(
      plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
      text = element_text(size = 12),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 12),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_line(color = "grey80")
    ) +
    labs(
      x = x_label,
      y = NULL
    )
}


rho_sc_posterior <- function(rho1, rho2, mu1, mu2, tau1, tau2, eta0, eta1, eta2, Adj, Adj2, Minc) {
  n <- nrow(Adj)
  I <- diag(n)
  D <- diag(rowSums(Minc))
  
  rho_sc_values <- numeric(length(rho1))
  
  trace <- function(A) {
    return(sum(diag(A)))
  }
  
  for (i in 1:length(rho1)) {
    mu1_vec <- matrix(mu1[i], nrow = n, ncol = 1)
    mu2_vec <- matrix(mu2[i], nrow = n, ncol = 1)
    
    W <- Adj
    W2 <- Adj2
    E_1 <- solve(tau1[i] * (D - rho1[i] * W)) + (eta0[i] * I + eta1[i] * W + eta2[i] * W2) %*% solve(tau2[i] * (D - rho2[i] * W)) %*% (eta0[i] * I + eta1[i] * W + eta2[i] * W2)
    E_12 <- (eta0[i] * I + eta1[i] * W + eta2[i] * W2) %*% solve(tau2[i] * (D - rho2[i] * W))
    E_2 <- solve(tau2[i] * (D - rho2[i] * W))
    
    J <- matrix(1, nrow = n, ncol = n)
    rho_sc <- trace(J %*% E_12 + J %*% t(E_12)) / (trace(J %*% E_1 + J %*% E_2) + t((mu1_vec - mu2_vec)) %*% J %*% (mu1_vec - mu2_vec))
    
    rho_sc_values[i] <- rho_sc
  }
  
  return(rho_sc_values)
}

#################################################################
# Data Preparation
#################################################################

# Read spatial data: Santiago Metropolitan Region 
rds_santiago <- readRDS("/Users/clementeferrer/Documents/Paper CAR/Datos/gadm41_CHL_3_pk.rds")
rds_santiago <- st_as_sf(rds_santiago) %>% filter(NAME_1 == "Santiago Metropolitan")
ca.neighbors <- poly2nb(rds_santiago)
n <- length(ca.neighbors)
Adj <- sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)

# Create second order neighbors matrix
Adj2 <- Adj %*% Adj
Adj2[Adj > 0] <- 0
diag(Adj2) <- 0 
Adj2 <- (Adj2 > 0) * 1

# Read poverty data
datos_rm <- read_xlsx("/Users/clementeferrer/Documents/Tesis USM/Datos/datos pobreza.xlsx")
Y1 = as.numeric(datos_rm$HT)
Y2 = as.numeric(datos_rm$Casen)
X1 = matrix(1, nrow = nrow(Adj), ncol = 1)
X2 = matrix(1, nrow = nrow(Adj), ncol = 1)

# Model2 : SAE| HT
Y.o2 = c(Y2,Y1)
X.o2 = as.matrix(bdiag(X1, X2))

#################################################################
# JAGS Model Specification
#################################################################

sink("GMCAR.txt")
cat("
    model
    {
    Q1[1:k, 1:k] <- tau1*(D - rho1*Minc)
    Q2[1:k, 1:k] <- tau2*(D - rho2*Minc)
    
    W1[1:k] ~ dmnorm(rep(0, k), Q1)
    W2[1:k] ~ dmnorm(rep(0, k), Q2)
    
    A21 <- eta021 * I + eta121 * Minc + eta221 * Minc2
    
    W[1:k] <- W1
    W[(k+1):(2*k)] <- A21 %*% W1 + W2
    
    for (i in 1:k)
    {
    mu[i] <- X[i,] %*% beta + W[i]
    Y[i] ~ dnorm(mu[i], taue1)
    }
    for (i in (k+1):(2*k))
    {
    mu[i] <- X[i,] %*% beta + W[i]
    Y[i] ~ dnorm(mu[i], taue2)
    }
    
    rho1 ~ dunif(0, 0.999)
    rho2 ~ dunif(0, 0.999)
    tau1 ~ dgamma(0.1, 0.1)
    tau2 ~ dgamma(0.1, 0.1)
    eta021 ~ dnorm(0, 0.01)
    eta121 ~ dnorm(0, 0.01)
    eta221 ~ dnorm(0, 0.01)
    
    taue1 ~ dgamma(0.1, 0.1)
    taue2 ~ dgamma(0.1, 0.1)
    vare1 <- 1/taue1
    vare2 <- 1/taue2
    
    beta[1:nc] ~ dmnorm(rep(0.118,nc), (0.1*I1))
    }
    ", fill = TRUE)
sink()


#################################################################
# Model Data and Initial Values
#################################################################

Minc <- Adj
Minc2 <- Adj2

model.data2 <- list(k = n, nc = ncol(X.o2), I = diag(n), I1 = diag(ncol(X.o2)),
                    Minc = Minc, Minc2 = Minc2,
                    D = diag(rowSums(Minc)), X = X.o2, Y = Y.o2)


model.inits <- rep(list(list(rho1 = 0.5, rho2 = 0.5, tau1 = 1, tau2 = 1, eta021 = 0, 
                             eta121 = 0, eta221 = 0, taue1 = 1, taue2 = 1,
                             beta = rep(0.118, ncol(X.o1)), W1 = rep(0, n), 
                             W2 = rep(0, n))),1)

model.param <- c("beta", "rho1", "rho2", "tau1", "tau2", "eta021", "eta121",
                 "eta221", "vare1", "vare2", "W")

#################################################################
# Model Execution
#################################################################

#Model2 result: Casen | HT
set.seed(123)
est.MCAR2 <- jags(model.data2, model.inits, model.param, "GMCAR.txt",
                  n.chains = 1, n.iter = 30000,n.burnin = 15000, n.thin = 1)
print(est.MCAR2)

#################################################################
# Posterior Analysis
#################################################################

mcmc_samples <- as.mcmc(est.MCAR2)
rho1 <- mcmc_samples[[1]][, "rho1"]
rho2 <- mcmc_samples[[1]][, "rho2"]
mu1 <- mcmc_samples[[1]][, "beta[1]"]
mu2 <- mcmc_samples[[1]][, "beta[2]"]
tau1 <- mcmc_samples[[1]][, "tau1"]
tau2 <- mcmc_samples[[1]][, "tau2"]
eta0 <- mcmc_samples[[1]][, "eta021"]
eta1 <- mcmc_samples[[1]][, "eta121"]
eta2 <- mcmc_samples[[1]][, "eta221"]

rho_sc_values <- rho_sc_posterior(rho1, rho2, mu1, mu2, tau1, tau2, eta0,
                                  eta1, eta2, Adj, Adj2, Minc)

#################################################################
# Plots
#################################################################

rho_sc_plot <- histogram_plot(rho_sc_values, expression(rho[s*","*c]))
rho_sc_plot
plot1 <- histogram_plot(rho1, expression(rho[HT]))
plot2 <- histogram_plot(rho2, expression(rho[SAE]))
plot3 <- histogram_plot(mu1, expression(mu[HT]))
plot4 <- histogram_plot(mu2, expression(mu[SAE]))
plot5 <- histogram_plot(tau1, expression(tau[HT]))
plot6 <- histogram_plot(tau2, expression(tau[SAE]))
plot7 <- histogram_plot(eta0, expression(eta[0]))
plot8 <- histogram_plot(eta1, expression(eta[1]))

combined_plot1 <- (plot1 + plot2) / (plot3 + plot4) / (plot5 + plot6) / (plot7 + plot8)

