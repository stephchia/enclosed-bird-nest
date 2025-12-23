# This script performs PCA on climatic variables and pairwise correlations among traits

library(factoextra) # PCA visualization
library(PerformanceAnalytics) # pairwise correlations visualization

# data version
version <- "" # global
version <- "_north" # Northern hemisphere

# import dataset
dt <- read.csv(paste0("data/processed/trait", version, "/sp_traits.csv"), row.names = 1)

## PCA for climatic variables
# (variable explanation: https://chelsa-climate.org/exchelsa-extended-bioclim/)
pca <- prcomp(as.matrix(dt[,c("tmean","dtr","hurs","rsds","vpd","precp","wind")]), scale = T)
# pca <- prcomp(as.matrix(dt[,c("pred","tmean","dtr","hurs","rsds","vpd","precp","wind")]), scale = T)

# Make necessary change of signs of PC1 & PC2 loadings
# (-hurs for PC1 (aridity), +tmean for PC2 (temperature))
if (pca$rotation["precp", 1] > 0) {
  pca$x[, "PC1"] <- -pca$x[,"PC1"]
  pca$rotation[, 1] <- -pca$rotation[,1]
}
if (pca$rotation["tmean", 2] < 0) {
  pca$x[, "PC2"] <- -pca$x[, "PC2"]
  pca$rotation[, 2] <- -pca$rotation[, 2]
}

# PCA summary table (loadings and variance explained) (Table S1)
pca_summary <- rbind(pca$rotation, (pca$sdev)^2 / sum((pca$sdev)^2) * 100)
pca_summary

# Add PC1 & PC2 into the trait dataset and save
dt$PC1 <- pca$x[, "PC1"] # aridity
dt$PC2 <- pca$x[, "PC2"] # temperature

write.csv(dt, paste0("data/processed/trait", version, "/sp_traits_pca.csv"))

# Plot PCA biplot 
rownames(pca$rotation) <- c("Mean tmperature","Diurnal temperature range","Humidity","Solar radiation","Vapor pressure deficit","Precipitation","Wind")
fviz_pca_biplot(pca, repel = TRUE, label = c("var","quali"), axes = c(1, 2),
                col.var = "black", col.ind = "cornflowerblue", alpha.ind = 0.4, labelsize = 4, pointsize = 1, stroke = 0,
                xlab = paste0("PC1 (", round(pca_summary[nrow(pca_summary), 1], 1),"%)"), 
                ylab = paste0("PC2 (", round(pca_summary[nrow(pca_summary), 2], 1),"%)")) +
  theme(panel.grid = element_blank())
ggsave(paste0("output/img/pca", version, ".pdf"), width = 6, height = 4)

# Pairwise correlation between life-history traits 
a <- cor(dt[, c("Enclosed","Dome","Cavity","Ground","Cooperative","Egg","Clutch","PC1","PC2","pred")])
write.csv(a, "pairwise_corr.csv")
cor(dt[, c("tmean","PC2")])

## Pairwise correlation between environmental variables
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x, plot = FALSE)
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
  # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.12), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt, cex = 3) # Resize the text by level of correlation
}

# Creating the scatter plot matrix
pairs(dt[, c("PC1","PC2","pred")],
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      pch = 21, col = alpha("black", .2), 
      labels = c("PC1","PC2","Pred"), main = "",
      row1attop = TRUE, gap = 1, cex.labels = NULL, font.labels = 1)

