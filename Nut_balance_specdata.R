#' Soil nutrient mass balance spectral data
#' C,N and Mehlich-3 extractable P,K,S,Ca & Mg, from 60 sentinel sites with HSTXT MIR spectra
#' M. Walsh, Nov. 2015

# install.packages(c("devtools"), dependencies=T)
require(devtools)

# Data setup --------------------------------------------------------------
SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Soils/master/Nut_balance_benchmarks.R"
source_url(SourceURL)

# Load and merge HSTXT MIR spectra
download("https://www.dropbox.com/s/6fvipxqlmg704g3/hstxt_MIR.csv.zip?dl=0", "hstxt_MIR.csv.zip", mode="wb")
unzip("hstxt_MIR.csv.zip", overwrite=T)
mir <- read.table("hstxt_MIR.csv", header=F, sep=",", stringsAsFactors=F)
mir <- as.data.frame(mir)
names(mir) <- c("SSN", paste("m", signif(seq(7497.964, 599.76, length.out=3578), 6),sep=""))
nb60_cal <- merge(nb60_cal, mir, by="SSN")
nb60_val <- merge(nb60_val, mir, by="SSN")

# Write data files --------------------------------------------------------
write.csv(nb60_cal, "nb60_cal.csv", row.names=F)
write.csv(nb60_val, "nb60_val.csv", row.names=F)
