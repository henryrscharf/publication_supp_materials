library(sp)
library(spdep)
library(spatialreg)

load("../data/cycle_injuries.RData")

highlight <- c(1, 5, 31, 40, 41)
fit <- lm(TotalRateInjury ~ length, data = SRA@data)
y <- fit$residuals
SRA@data$residuals <- y
n <- length(y)
z <- (y - mean(y)) / sqrt(sum(y^2) / n)
SRA_nb <- poly2nb(SRA)
W <- nb2listw(SRA_nb, style = "W")
nsim <- 1e4 # for permutation tests

## moran scatter plot
W_mat <- listw2mat(W)
Wz <- W_mat %*% z
pdf("../fig/moran_scatterplot.pdf")
par(mar = c(4, 4, 0.1, 0.1))
pch <- rep(1, n)
pch[highlight] <- NA
plot(z, Wz, asp = 1, xlab = "z", ylab = "Wz", pch = pch)
points(z[highlight], Wz[highlight], pch = 16)
text(x = z[highlight] - 0.08, y = Wz[highlight], labels = highlight)
abline(0, coef(lm(Wz ~ z))[2], col = "gray")
abline(h = 0, v = 0, lty = 2, col = "gray")
dev.off()

LISA_I <- localmoran(y, listw = W)
SRA@data$LISA_I <- LISA_I[, 1]
SRA@data$LISA_I_p <- LISA_I[, 5]
## double check calculation
range(z * Wz - LISA_I[, 1])

LISA_I_perm <- localmoran_perm(y, listw = W, nsim = nsim)
SRA@data$LISA_I_perm_p <- LISA_I_perm[, 5]
## plot showing close agreement between exact and permutation-based tests
plot(LISA_I[, 5], LISA_I_perm[, 5], xlim = c(0, 0.5), ylim = c(0, 0.5))
cutoff <- 0.05
abline(v = cutoff, h = cutoff)
abline(0, 1, lty = 2, col = "gray")

## there is evidence of global autocorrelation in the residuals
moran.test(y, listw = W)

## fit a new regression model to the data with same predictor and SAR structure for residuals
fit_sar <- errorsarlm(TotalRateInjury ~ length, data = SRA@data, listw = W)
SRA@data$residuals_SAR <- fit_sar$residuals
lm.target <- lm(fit_sar$tary ~ fit_sar$tarX - 1)
Omega <- invIrW(W, fit_sar$lambda)
Omega1 <- tcrossprod(Omega)

## these two methods attempt to adjust for the presence of global autocorrelation
LISA_I_saddle <- localmoran.sad(model = lm.target, nb = SRA_nb, Omega = Omega)
SRA@data$LISA_I_global <- summary(LISA_I_saddle)[, 1]
SRA@data$LISA_I_saddle_p <- summary(LISA_I_saddle)[, 3]

LISA_I_exact <- localmoran.exact.alt(lm.target, nb = SRA_nb, Omega = Omega)
SRA@data$LISA_I_exact_p <- print(LISA_I_exact)[, 3] 

## Two alternative measures of non-stationarity: Geary's c (also a LISA) and Getis-Ord G

## Geary's c
LISA_c <- localC_perm(y, listw = W, nsim = nsim)
SRA@data$LISA_c <- as.numeric(LISA_c)
SRA@data$LISA_c_p <- attr(LISA_c, "pseudo-p")[, 4]

## Getis-Ord G
local_G <- localG_perm(y, listw = nb2listw(SRA_nb), nsim = nsim)
SRA@data$local_G <- as.numeric(local_G)
SRA@data$local_G_p <- pnorm(local_G)

## Summarize
stat_names <- c("LISA_I", "LISA_I_global", "LISA_c", "local_G")
table <- cbind(fit$residuals[highlight], fit_sar$residuals[highlight], 
               SRA@data[highlight, stat_names])
names(table) <- c("residuals", "SAR residuals", "LISA I", "LISA I (SAR adj.)", "LISA c", "local G")
xtable::xtable(table, digits = 3)
pvalue_names <- c('LISA_I_p', 'LISA_I_perm_p', 'LISA_I_saddle_p', 'LISA_I_exact_p', 'LISA_c_p', 'local_G_p')
SRA@data[highlight, pvalue_names]
sapply(pvalue_names, function(test){
  which(SRA@data[, test] < 0.05)
})

pdf("../fig/resid_map.pdf", width = 8.5, height = 7)
spplot(SRA, c("residuals", "residuals_SAR", "TotalRateInjury"), 
       names.attr = c("residuals", "residuals (SAR)", "Injury Rate"),
       col.regions = c(RColorBrewer::brewer.pal(11, "BrBG")), 
       at = seq(-64, 64, l = 11), as.table = F, layout = c(2, 2))
dev.off()

pdf("../fig/map_pvalues.pdf", width = 8.5, height = 5)
spplot(SRA, c('LISA_I_p', 'LISA_I_exact_p', "LISA_c_p", 
              'LISA_I_perm_p', 'LISA_I_saddle_p', "local_G_p"), 
       names.attr = c("Moran I normal", "Moran I exact", "Geary's c perm.", 
                      "Moran I perm.", "Moran I saddle", "Getis-Ord G"),
       col.regions = c("darkred", gray.colors(19, end = 1)), at = seq(0, 1, l= 21),
       sp.layout = list("sp.text", coordinates(SRA)[highlight, ], highlight), as.table = T, layout = c(3, 2))
dev.off()