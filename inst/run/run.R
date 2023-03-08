# install.packages("devtools", repos = "http://cran.us.r-project.org")
devtools::install_github("saveriofrancini/dataDriven", force = T)
library(dataDriven)
a <- Sys.time()
s <- bootStrap(outDir = getwd())

png("Harmon.png", width = 10000, height = 10000)
spplot(s,
       zcol = "Harmon",
       colorkey = list(width = 30, labels = list(cex = 20)))
dev.off()

png("error.png", width = 10000, height = 10000)
spplot(s,
       zcol = "error",
       colorkey = list(width = 30, labels = list(cex = 20)))
dev.off()

b <- Sys.time()
print(b - a)  