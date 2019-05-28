library(zoo)
library(matlib)
library(rio)

Xmatrix <- function(data) {
  count <- as.zoo(c(1:length(data)-1))
  const <- as.zoo(rep(1,length(data)-1))
  coltwo <- as.zoo(count-0.5*length(data))
  coltwo <- window(coltwo,end = (length(data)-1))
  colthree <- window(data,end = (length(data)-1))
  cbind.zoo(const,coltwo, colthree)
}

Ymatrix <- function(data) {
  window(data, start = 2)
}

theta_tau <- function(Xmatrix, Ymatrix) {
  inv(t(Xmatrix) %*% Xmatrix)%*%t(Xmatrix)%*%Ymatrix
}

s_sq_etau <- function(Xmatrix, Ymatrix, data) {
  (length(data)-4)^(-1) %*% t(Ymatrix) %*%
    (
      diag(nrow = nrow(Xmatrix))-Xmatrix %*%
       ((t(Xmatrix) %*% Xmatrix)^(-1)) %*% t(Xmatrix)
      ) %*%
    Ymatrix
}

tau_alphatau <- function(Xmatrix, s_sq_etau_hat, theta) {
  C <- inv(t(Xmatrix) %*% Xmatrix)
  (C[1,1] * s_sq_etau_hat^2)^(-0.5) * theta[1,1]
}

tau_betatau <- function(Xmatrix, s_sq_etau_hat, theta) {
  C <- inv(t(Xmatrix) %*% Xmatrix)
  (C[2,2] * s_sq_etau_hat^2)^(-0.5) * theta[2,1]
}

sigma_sq_zero <- function(data) {
  (length(data)-1)^(-1)*
    sum((window(data,start = 2) - lag(data, i=-1))^2)
}

phi_2 <- function(s_sq_etau,sigma_sq_zero, data) {
  (3*s_sq_etau)^(-1)*
    ((length(data)-1)*sigma_sq_zero-(length(data)-4)*s_sq_etau)
}

y_mean <- function(i, data) {
  ((length(data)-1)^(-1))*sum(lag(data, k = i))
}

phi_3 <- function(s_sq_etau, sigma_sq_zero, data) {
  (2*s_sq_etau)^(-1)*
    (
      (length(data)-1)*
        (
          sigma_sq_zero-(y_mean(i=0,data = data)-y_mean(i=-1,data=data))^2
        )-(length(data)-4)*s_sq_etau
    )
}

setwd("~/Documents/ChStockMkt")
ds <- import("SSEC.csv")
ds <- log(ds$SSEC.Close)
ds <- as.zoo(ds)
X <- Xmatrix(ds)
Y <- Ymatrix(ds)
theta <- theta_tau(X,Y)
s_sq_etau_hat <- s_sq_etau(X,Y,ds)
tau_alphatau_hat <- tau_alphatau(X, s_sq_etau_hat, theta)
tau_betatau_hat <- tau_betatau(X, s_sq_etau_hat, theta)
sigma_sq_zero_hat <- sigma_sq_zero(ds)
phi2 <- phi_2(s_sq_etau_hat, sigma_sq_zero_hat, ds)
phi3 <- phi_3(s_sq_etau_hat, sigma_sq_zero_hat, ds)

ds2 <- import("SZSE.csv")
ds2 <- log(ds2$Open)
ds2 <- as.zoo(ds2)
sz_X <- Xmatrix(ds2)
sz_Y <- Ymatrix(ds2)
sz_theta <- theta_tau(sz_X,sz_Y)
sz_s_sq_etau_hat <- s_sq_etau(sz_X,sz_Y,ds2)
sz_tau_alphatau_hat <- tau_alphatau(sz_X, sz_s_sq_etau_hat, sz_theta)
sz_tau_betatau_hat <- tau_betatau(sz_X, sz_s_sq_etau_hat, sz_theta)
sz_sigma_sq_zero_hat <- sigma_sq_zero(ds2)
sz_phi2 <- phi_2(sz_s_sq_etau_hat, sz_sigma_sq_zero_hat, ds2)
sz_phi3 <- phi_3(sz_s_sq_etau_hat, sz_sigma_sq_zero_hat, ds2)
