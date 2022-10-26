# Fit-hic
library("FitHiC")


fragsfile = "../hicpro/LT1/raw/40000/LT1_40000_abs.bed"
intersfile = "../hicpro/LT1/iced/40000/LT1_40000_iced.matrix"
biasfile = "../hicpro/LT1/iced/40000/LT1_40000_iced.matrix.biases"
outdir <- "../data/FitHiC"
FitHiC(fragsfile, intersfile, outdir, biasfile, libname="LT1",
       distUpThres=2000000, distLowThres=200000, visual=TRUE, useHiCPro=TRUE)
