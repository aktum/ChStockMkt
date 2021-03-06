library(rio)
setwd("~/Documents/ChStockMkt")
szse <- import("SZSE.csv")
szse <- log(szse$Open)

colone <- tail(szse,-1)
count <- c(1:(length(szse)-1))
coltwo <- count-0.5*length(szse)
colthree <- head(szse,length(szse)-1)

ds <- as.data.frame(cbind(colone,coltwo,colthree))
#mod_nD_nT = model with no drift and no trend

mod_nD_nT <- lm(colone~0+colthree,ds)
summary(mod_nD_nT)
rss_mod_nD_nT <- sum(mod_nD_nT$residuals^2)
rss_mod_nD_nT <- rss_mod_nD_nT/(length(szse)-4)


mod_D_T <- lm(colone~coltwo+colthree,ds)
summary(mod_D_T)
rss_mod_D_T <- sum(mod_D_T$residuals^2)
rss_mod_D_T_corr <- rss_mod_D_T/(length(szse)-4)

mod_D_nT <- lm(colone~colthree,ds)
summary(mod_D_nT)
rss_mod_D_nT <- sum(mod_D_nT$residuals^2)
rss_mod_D_nT_corr <- rss_mod_D_nT/(length(szse)-4)

#H0 (alpha,beta0,beta1)=(1,0,0)
Phi2 <- (rss_mod_nD_nT-rss_mod_D_T)/3*rss_mod_D_T_corr
#H0 (alpha,beta0,beta1)=(alpha,0,1)
Phi3 <- (rss_mod_D_nT-rss_mod_D_T)/2*rss_mod_D_T_corr
#H0 (alpha,beta1)=(0,1)
Phi1 <- (rss_mod_nD_nT-rss_mod_D_T)/2*rss_mod_D_T_corr
