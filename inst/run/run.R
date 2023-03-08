devtools::install_github("saveriofrancini/dataDriven", force = T)
library(dataDriven)

s <- dataDriven::bootStrap(iterations = 100, outDir = getwd())
