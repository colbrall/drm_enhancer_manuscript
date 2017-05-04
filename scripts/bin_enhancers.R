robust <- read.delim("all_fantom_enhancers.bed",sep = "\t")
for (i in 1:nrow(robust)) {
  if (robust$X.tiss[i] <= 20) {robust$bin[i] <- robust$X.tiss[i]}
  if (robust$X.tiss[i] <= 40 && robust$X.tiss[i] > 20) {robust$bin[i] <- 40}
  if (robust$X.tiss[i] <= 60 && robust$X.tiss[i] > 40) {robust$bin[i] <- 60}
  if (robust$X.tiss[i] <= 80 && robust$X.tiss[i] > 60) {robust$bin[i] <- 80}
  if (robust$X.tiss[i] <= 100 && robust$X.tiss[i] > 80) {robust$bin[i] <- 100}
  if (robust$X.tiss[i] <= 120 && robust$X.tiss[i] > 100) {robust$bin[i] <- 120}
  if (robust$X.tiss[i] <= 140 && robust$X.tiss[i] > 120) {robust$bin[i] <- 140}
  if (robust$X.tiss[i] <= 160 && robust$X.tiss[i] > 140) {robust$bin[i] <- 160}
  if (robust$X.tiss[i] > 160) {robust$bin[i] <- 180} #for sorting; relabelled as >160
}