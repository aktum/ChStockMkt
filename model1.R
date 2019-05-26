library(zoo)
library(matlib)
library(rio)
setwd("~/Documents/ChStockMkt")
ds <- import("SSEC.csv")
x.data <- as.Date(ds$Date)
ds <- log(ds$SSEC.Close)
ds <- as.zoo(ds)
y_mean <- function(i, data) {
  ((length(data)-1)^(-1))*sum(lag(data, k = i))
}

rho_mu <- function(data) {
  (sum((lag(data, k = -1)-y_mean(i = -1, data = data))^2)^(-1))*
    sum((data-y_mean(i = 0, data = data))*
          (lag(data, k = -1)-y_mean(i = -1,data = data)))
}
rho <- rho_mu(ds)
alpha_mu <- function(rho, data) {
  alpha_mu <- y_mean(i=0,data = data)-rho*y_mean(i=-1, data = data)
}
alpha <- alpha_mu(rho, ds)
s_sq_emu <- function(alpha, rho, data) {
  (length(data)-3)^(-1)*sum((data-alpha-rho*lag(data,k = -1))^2)
}
s_sq_emu_hat <- s_sq_emu(alpha, rho, ds)

s_sq_alphamu <- function(s_sq_emu, data) {
  s_sq_emu*((length(data)-1)^(-1)+(y_mean(i = -1, data))^2*
              (sum((lag(data, k = -1)-
                      y_mean(i = -1, data = data))^2)^(-1)))
}
s_sq_alphamu_hat <- s_sq_alphamu(s_sq_emu_hat, ds)
#the statistic constructed by analogy to the regression "t statistic" for the estimated alpha
tau_alphamu <- function(s_sq_alphamu, alpha_mu) {
  (sqrt(s_sq_alphamu)^(-1))*alpha_mu
}
tau_alphamu_hat <- tau_alphamu(s_sq_alphamu_hat,alpha)
sigma_sq_zero <- function(data) {
  (length(data)-1)^(-1)*
    sum((window(data,end = (length(data)-1))-lag(data, i=-1))^2)
}
sigma_sq_zero_hat <- sigma_sq_zero(ds)
phi_1 <- function(s_sq_emu, sigma_sq_zero) {
  (2*s_sq_emu)^(-1)*
    ((length(data)-1)*sigma_sq_zero-(length(data)-3)*s_sq_emu)
}
phi1 <- phi_1(s_sq_emu_hat, sigma_sq_zero_hat)

ds2 <- import("SZSE.csv")
ds2 <- log(ds2$Open)
ds2 <- as.zoo(ds2)
sz_rho <- rho_mu(ds2)
sz_alpha <- alpha_mu(sz_rho, ds2)
sz_s_sq_emu_hat <- s_sq_emu(sz_alpha, sz_rho, ds2)
sz_s_sq_alphamu_hat <- s_sq_alphamu(sz_s_sq_emu_hat, ds2)
sz_tau_alphamu_hat <- tau_alphamu(sz_s_sq_alphamu_hat,sz_alpha)
sz_sigma_sq_zero_hat <- sigma_sq_zero(ds2)
sz_phi1 <- phi_1(sz_s_sq_emu_hat, sz_sigma_sq_zero_hat)
#Графики для работы
plot.zoo(y = ds, xlab = x.data)
axis(x.data)
