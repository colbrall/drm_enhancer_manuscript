#colbran_2017_drm_analyses.R
# code used to run analyses and plot things for 
#"Short DNA sequence patterns accurately identify broadly active human enhancers"
library("ggplot2")
library("dplyr")
library("pROC")
library('readr')
library('ppcor')

#-----------------------------------------------#
# Fold Enrichment Analyses
# a through d are DRM counts from drm_finder.py on negative bed files generated using shufflebed
# fant_600 is the same over the entire fantom enhancer dataset, with length set to 600 bp

se <- function(x) sqrt(var(x)/length(x))

# run this for all 4 negative files and for fant_600
a$d_ga <-a$X4/(a$X3 - a$X2)
a$d_ca <-a$X5/(a$X3 - a$X2)
a$d_gc <- a$X6/(a$X3 - a$X2)
a$d_ta <- a$X7/(a$X3 - a$X2)

neg_1 <- mean(fant_600[fant_600$tiss > 45,]$d_ga)/mean(neg$d_ga)
neg_2 <- mean(fant_600[fant_600$tiss > 45,]$d_ga)/mean(a$d_ga)
neg_3 <- mean(fant_600[fant_600$tiss > 45,]$d_ga)/mean(b$d_ga)
neg_4 <- mean(fant_600[fant_600$tiss > 45,]$d_ga)/mean(c$d_ga)
neg_1 <- c(neg_1,mean(fant_600[fant_600$tiss > 45,]$d_ca)/mean(neg$d_ca))
neg_2 <- c(neg_2,mean(fant_600[fant_600$tiss > 45,]$d_ca)/mean(neg$d_ca))
neg_3 <- c(neg_3,mean(fant_600[fant_600$tiss > 45,]$d_ca)/mean(neg$d_ca))
neg_4 <- c(neg_4,mean(fant_600[fant_600$tiss > 45,]$d_ca)/mean(neg$d_ca))
neg_1 <- c(neg_1,mean(fant_600[fant_600$tiss > 45,]$d_gc)/mean(neg$d_gc))
neg_2 <- c(neg_2,mean(fant_600[fant_600$tiss > 45,]$d_gc)/mean(neg$d_gc))
neg_3 <- c(neg_3,mean(fant_600[fant_600$tiss > 45,]$d_gc)/mean(neg$d_gc))
neg_4 <- c(neg_4,mean(fant_600[fant_600$tiss > 45,]$d_gc)/mean(neg$d_gc))
neg_1 <- c(neg_1,mean(fant_600[fant_600$tiss > 45,]$d_ta)/mean(neg$d_ta))
neg_2 <- c(neg_2,mean(fant_600[fant_600$tiss > 45,]$d_ta)/mean(neg$d_ta))
neg_3 <- c(neg_3,mean(fant_600[fant_600$tiss > 45,]$d_ta)/mean(neg$d_ta))
neg_4 <- c(neg_4,mean(fant_600[fant_600$tiss > 45,]$d_ta)/mean(neg$d_ta))

f_enr <- data.frame(DRM = c("GA","CA","GC","TA"),neg_1 = neg_1,neg_2 = neg2,neg_3 = neg_3, neg_4=neg_4)
f_enr$neg_fe <- (f_enr$neg_1 + f_enr$neg_2 + f_enr$neg_3 + f_enr$neg_4)/4
f_enr$neg_log <- (log2(f_enr$neg_1) + log2(f_enr$neg_2) + log2(f_enr$neg_3) + log2(f_enr$neg_4))/4
for (i in 1:nrow(f_enr)) {
  f_enr$neg_se[i] <- se(c(log2(f_enr$neg_1[i]),log2(f_enr$neg_2[i]),log2(f_enr$neg_3[i]),log2(f_enr$neg_4[i])))
}

f_enr$narr_fe <- c(mean(fant_600[fant_600$tiss > 45,]$d_ga)/mean(fant_600[fant_600$tiss == 1,]$d_ga),
                   mean(fant_600[fant_600$tiss > 45,]$d_ca)/mean(fant_600[fant_600$tiss == 1,]$d_ca),
                   mean(fant_600[fant_600$tiss > 45,]$d_gc)/mean(fant_600[fant_600$tiss == 1,]$d_gc),
                   mean(fant_600[fant_600$tiss > 45,]$d_ta)/mean(fant_600[fant_600$tiss == 1,]$d_ta))
f_enr$narr_log <- log2(f_enr$narr_fe)

wilcox.test(fant_600[fant_600$tiss > 45,]$d_ga,a$d_ga)$p.value
wilcox.test(fant_600[fant_600$tiss > 45,]$d_ca,a$d_ca)$p.value
wilcox.test(fant_600[fant_600$tiss > 45,]$d_gc,a$d_gc)$p.value
wilcox.test(fant_600[fant_600$tiss > 45,]$d_ta,a$d_ta)$p.value

wilcox.test(fant_600[fant_600$tiss > 45,]$d_ga,fant_600[fant_600$tiss == 1,]$d_ga)$p.value
wilcox.test(fant_600[fant_600$tiss > 45,]$d_ca,fant_600[fant_600$tiss == 1,]$d_ca)$p.value
wilcox.test(fant_600[fant_600$tiss > 45,]$d_gc,fant_600[fant_600$tiss == 1,]$d_gc)$p.value
wilcox.test(fant_600[fant_600$tiss > 45,]$d_ta,fant_600[fant_600$tiss == 1,]$d_ta)$p.value

ggplot(f_enr,aes(x=DRM,y=neg_log,fill=DRM)) +
  theme(text = element_text(size=18),axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold")) + geom_bar(color="darkgrey",stat="identity")  + scale_fill_manual(values = topo.colors(4,alpha=0.75)) +
  labs(x="DRM",y="Mean log2(FE)") +ylim(-2.5,5.5) + geom_errorbar(aes(ymin=neg_log-neg_se,ymax=neg_log+neg_se,width=0.3),alpha=0.75)
ggplot(f_enr,aes(x=DRM,y=narr_log,fill=DRM)) +
  theme(text = element_text(size=18),axis.text.x = element_text(face = "bold"),axis.text.y = element_text(face = "bold")) + geom_bar(color="darkgrey",stat="identity")  + scale_fill_manual(values = topo.colors(4,alpha=0.75)) +
  labs(x="DRM",y="log2(FE)") +ylim(-2.5,5.5)

#-----------------------------------------------#
# Compare ROC curves
narr_tf <- read_delim("2017-04-24__MKL_1kern_tf_pwms~~all_enh_broad_600~~all_enh_narr_600_cv_test_scores_tf_pwms.txt", "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
narr_6mer <- read_delim("2016-12-08__MKL_1kern_spectrum6~~all_enh_broad_600~~all_enh_narr_600_cv_test_scores_spectrum6.txt", "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
narr_drm <- read_delim("2016-11-30__MKL_1kern_drm_counts~~all_enh_broad_600~~all_enh_narr_600_cv_test_scores_drm_counts.txt", "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
roc.test(roc(narr_tf$X5,narr_tf$X4),roc(narr_6mer$X5,narr_6mer$X4))
roc.test(roc(narr_drm$X5,narr_drm$X4),roc(narr_6mer$X5,narr_6mer$X4))

neg_tf<- read_delim("2017-04-24__MKL_1kern_tf_pwms~~all_enh_broad~~all_enh_broad_rand-neg0_cv_test_scores_tf_pwms.txt", "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
neg_6mer<- read_delim("2016-12-08__MKL_1kern_spectrum6~~all_enh_broad~~all_enh_broad_rand-neg0_cv_test_scores_spectrum6.txt","\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
neg_drm <- read_delim("2016-11-30__MKL_1kern_drm_counts~~all_enh_broad~~all_enh_broad_rand-neg0_cv_test_scores_drm_counts.txt", "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
roc.test(roc(neg_tf$X5,neg_tf$X4),roc(neg_6mer$X5,neg_6mer$X4))
roc.test(roc(neg_drm$X5,neg_drm$X4),roc(neg_6mer$X5,neg_6mer$X4))

gc_tf<- read_delim("2017-04-24__MKL_1kern_tf_pwms~~all_enh_broad~~all_enh_broad_gc-neg0_cv_test_scores_tf_pwms.txt", "\t", escape_double = FALSE, col_names = FALSE,comment = "#", trim_ws = TRUE)
gc_6mer<- read_delim("2016-12-08__MKL_1kern_spectrum6~~all_enh_broad~~all_enh_broad_gc-neg0_cv_test_scores_spectrum6.txt", "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
gc_drm <- read_delim("2016-11-30__MKL_1kern_drm_counts~~all_enh_broad~~all_enh_broad_gc-neg0_cv_test_scores_drm_counts.txt", "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
roc.test(roc(gc_tf$X5,gc_tf$X4),roc(gc_6mer$X5,gc_6mer$X4))
roc.test(roc(gc_drm$X5,gc_drm$X4),roc(gc_6mer$X5,gc_6mer$X4))

roc.test(roc(gc_drm$X5,gc_drm$X4),roc(neg_drm$X5,neg_drm$X4))
roc.test(roc(neg_6mer$X5,neg_6mer$X4),roc(gc_6mer$X5,gc_6mer$X4))
roc.test(roc(gc_6mer$X5,gc_6mer$X4),roc(narr_6mer$X5,narr_6mer$X4))

#-----------------------------------------------#
# TF GC,CpG content, Expression
spec <- read.csv("~/Box Sync/dinucleotide_repeat_motifs/supplementary_files/data/tf_motif_specificity.csv")
spec <- na.omit(spec)
spec$Target.ID <- spec$X.Symbol..Human.
spec$X.Symbol..Human. <- NULL

# Hist of TSPS
for (i in 1:nrow(spec)) {
  if (spec$Tissue.specificity.score[i] < 1) {spec$expression[i] <- "Broad"}
  else {spec$expression[i] <- "Specific"}}
ggplot(spec,aes(Tissue.specificity.score,fill=expression)) + geom_histogram(breaks = seq(0,5,by=0.1),closed="left")+
  theme(text = element_text(size=20))+xlab("TSPS") + ylab("Number of TFs")
 
#JASPAR TFs
# manually fill in some of the NAs that had different syntax/alternate names that I could find online or were in complex in database
# for complexes, chose most specific value
f <- read.delim('~/Box Sync/dinucleotide_repeat_motifs/supplementary_files/data/jaspar_consensus.txt',sep="\t",header=T,stringsAsFactors = F)
t <- left_join(f,spec)
t[is.na(t)]<- 100
for (i in 1:nrow(t)) {
  t$Fantom.ID[i] <- t$Target.ID[i]
  if (t$Target.ID[i] == "BATF::JUN") {t$Tissue.specificity.score[i] <- 1.14
    t$Fantom.ID[i] <- "BATF"} 
  if (t$Target.ID[i] == "EBF1") {t$Tissue.specificity.score[i] <- 0.73
  t$Fantom.ID[i] <- "EBF"} 
  if (t$Target.ID[i] == "HINFP") {t$Tissue.specificity.score[i] <- 0.12
  t$Fantom.ID[i] <- "MIZF"}
  if (t$Target.ID[i] == "JUN(var.2)") {t$Tissue.specificity.score[i] <- 0.45
  t$Fantom.ID[i] <- "JUN"}
  if (t$Target.ID[i] == "JUND(var.2)") {t$Tissue.specificity.score[i] <- 0.18
  t$Fantom.ID[i] <- "JUND"}
  if (t$Target.ID[i] == "MYC::MAX") {t$Tissue.specificity.score[i] <- 0.58
  t$Fantom.ID[i] <- "MYC"}
  if (t$Target.ID[i] == "NFE2::MAF") {t$Tissue.specificity.score[i] <- 3.12
  t$Fantom.ID[i] <- "NFE2"}
  if (t$Target.ID[i] == "NR1H2::RXRA") {t$Tissue.specificity.score[i] <- 0.20
  t$Fantom.ID[i] <- "RXRA"}
  if (t$Target.ID[i] == "Pax6") {t$Tissue.specificity.score[i] <- 2.66
  t$Fantom.ID[i] <- "PAX6"}
  if (t$Target.ID[i] == "RXRA::VDR") {t$Tissue.specificity.score[i] <- 1.08
  t$Fantom.ID[i] <- "VDR"}
  if (t$Target.ID[i] == "SMAD2::SMAD3::SMAD4") {t$Tissue.specificity.score[i] <- 0.20
  t$Fantom.ID[i] <- "SMAD4"}
  if (t$Target.ID[i] == "STAT2::STAT1") {t$Tissue.specificity.score[i] <- 0.33
  t$Fantom.ID[i] <- "STAT1"}
  if (t$Target.ID[i] == "TAL1::GATA1") {t$Tissue.specificity.score[i] <- 3.76
  t$Fantom.ID[i] <- "GATA1"}
  if (t$Target.ID[i] == "TAL1::TCF3") {t$Tissue.specificity.score[i] <- 2.37
  t$Fantom.ID[i] <- "TAL1"}
  if (t$Target.ID[i] == "TP63") {t$Tissue.specificity.score[i] <- 1.99
  t$Fantom.ID[i] <- "TP73L"}
  if (t$Target.ID[i] == "ZEB1") {t$Tissue.specificity.score[i] <- 0.92
  t$Fantom.ID[i] <- "TCF8"}
}
jaspar <- t[-which(t$Tissue.specificity.score == 100),]


# ENCODE TFs
f <- read.table("~/Box Sync/dinucleotide_repeat_motifs/supplementary_files/data/encode_tf_consensus.txt",header=T,sep="\t",stringsAsFactors = F)
for (i in 1:nrow(f)) {f$Target.ID[i] <- strsplit(strsplit(f$Target.ID[i],">")[[1]][2],"_")[[1]][1]}
t <- left_join(f,spec)
t[is.na(t)]<- 100
for (i in 1:nrow(t)) {
  t$Fantom.ID[i] <- t$Target.ID[i]
  if (t$Target.ID[i] == "AHR::ARNT" | t$Target.ID[i] == "AHR::ARNT::HIF1A") {t$Tissue.specificity.score[i] <- 0.96
    t$Fantom.ID[i] <- "AHR"}
  if (t$Target.ID[i] == "BHLHE40") {t$Tissue.specificity.score[i] <- 0.36
    t$Fantom.ID[i] <- "BHLHB2"}
  if (t$Target.ID[i] == "BHLHE41") {t$Tissue.specificity.score[i] <- 0.89
    t$Fantom.ID[i] <- "BHLHB3"}
  if (t$Target.ID[i] == "BPTF") {t$Tissue.specificity.score[i] <- 0.51
    t$Fantom.ID[i] <- "FALZ"}
  if (t$Target.ID[i] == "CLOCK::ARNTL") {t$Tissue.specificity.score[i] <- 0.60
  t$Fantom.ID[i] <- "CLOCK"}
  if (t$Target.ID[i] == "CUX1") {t$Tissue.specificity.score[i] <- 0.19
  t$Fantom.ID[i] <- "CUTL1"}
  if (t$Target.ID[i] == "DDIT3::CEBPA") {t$Tissue.specificity.score[i] <- 1.4
  t$Fantom.ID[i] <- "CEBPA"}
  if (t$Target.ID[i] == "EBF1") {t$Tissue.specificity.score[i] <- 0.73
  t$Fantom.ID[i] <- "EBF"}
  if (t$Target.ID[i] == "GSC2") {t$Tissue.specificity.score[i] <- 3.74
  t$Fantom.ID[i] <- "GSCL"}
  if (t$Target.ID[i] == "HDX") {t$Tissue.specificity.score[i] <- 1
  t$Fantom.ID[i] <- "CXorf43"}
  if (t$Target.ID[i] == "HIF1A::ARNT") {t$Tissue.specificity.score[i] <- 0.39
  t$Fantom.ID[i] <- "HIF1A"}
  if (t$Target.ID[i] == "HINFP") {t$Tissue.specificity.score[i] <- 0.12
  t$Fantom.ID[i] <- "MIZF"}
  if (t$Target.ID[i] == "HLTF") {t$Tissue.specificity.score[i] <- 0.99
  t$Fantom.ID[i] <- "SMARCA3"}
  if (t$Target.ID[i] == "HNF1A" | t$Target.ID[i] == "HNF1") {t$Tissue.specificity.score[i] <- 1.63
  t$Fantom.ID[i] <- "TCF1"}
  if (t$Target.ID[i] == "HNF1B") {t$Tissue.specificity.score[i] <- 2.68
  t$Fantom.ID[i] <- "TCF2"}
  if (t$Target.ID[i] == "HNF4") {t$Tissue.specificity.score[i] <- 3.46
  t$Fantom.ID[i] <- "HNF4A"}
  if (t$Target.ID[i] == "HOMEZ") {t$Tissue.specificity.score[i] <- 0.1
  t$Fantom.ID[i] <- "KIAA1443"}
  if (t$Target.ID[i] == "IKZF1") {t$Tissue.specificity.score[i] <- 1.59
  t$Fantom.ID[i] <- "ZNFN1A1"}
  if (t$Target.ID[i] == "IKZF2") {t$Tissue.specificity.score[i] <- 0.66
  t$Fantom.ID[i] <- "ZNFN1A2"}
  if (t$Target.ID[i] == "IKZF3") {t$Tissue.specificity.score[i] <- 2.49
  t$Fantom.ID[i] <- "ZNFN1A3"}
  if (t$Target.ID[i] == "MEIS1::HOXA9") {t$Tissue.specificity.score[i] <- 1.41
  t$Fantom.ID[i] <- "HOXA9"}
  if (t$Target.ID[i] == "MZF1") {t$Tissue.specificity.score[i] <- 0.11
  t$Fantom.ID[i] <- "ZNF42"}
  if (t$Target.ID[i] == "NFE2L1::MAFG") {t$Tissue.specificity.score[i] <- 0.56
  t$Fantom.ID[i] <- "NFE2L1"}
  if (t$Target.ID[i] == "NKX2-1") {t$Tissue.specificity.score[i] <- 3.35
  t$Fantom.ID[i] <- "TITF1"}
  if (t$Target.ID[i] == "PATZ1") {t$Tissue.specificity.score[i] <- 0.16
  t$Fantom.ID[i] <- "ZNF278"}
  if (t$Target.ID[i] == "PDX1") {t$Tissue.specificity.score[i] <- 4.72
  t$Fantom.ID[i] <- "IPF1"}
  if (t$Target.ID[i] == "RBPJ") {t$Tissue.specificity.score[i] <- 0.3
  t$Fantom.ID[i] <- "RBPSUH"}
  if (t$Target.ID[i] == "RHOXF1") {t$Tissue.specificity.score[i] <- 2.98
  t$Fantom.ID[i] <- "OTEX"}
  if (t$Target.ID[i] == "RHOXF2") {t$Tissue.specificity.score[i] <- 3.63
  t$Fantom.ID[i] <- "PEPP-2"}
  if (t$Target.ID[i] == "TFAP2") {t$Tissue.specificity.score[i] <- 3.07
  t$Fantom.ID[i] <- "TFAP2A"}
  if (t$Target.ID[i] == "TGIF1") {t$Tissue.specificity.score[i] <- 0.34
  t$Fantom.ID[i] <- "TGIF"}
  if (t$Target.ID[i] == "VSX2") {t$Tissue.specificity.score[i] <- 1.38
  t$Fantom.ID[i] <- "CHX10"}
  if (t$Target.ID[i] == "ZBTB14") {t$Tissue.specificity.score[i] <- 0.15
  t$Fantom.ID[i] <- "ZFP161"}
  if (t$Target.ID[i] == "ZBTB18") {t$Tissue.specificity.score[i] <- 1.58
  t$Fantom.ID[i] <- "ZNF238"}
  if (t$Target.ID[i] == "ZEB1") {t$Tissue.specificity.score[i] <- 0.92
  t$Fantom.ID[i] <- "TCF8"}
  if (t$Target.ID[i] == "ZSCAN16") {t$Tissue.specificity.score[i] <- 0.25
  t$Fantom.ID[i] <- "ZNF435"}
}
encode <- t[-which(t$Tissue.specificity.score == 100),]

# CIS-BP TFs
# start by mapping CIS-BP ids to TF names, then match to FANTOM expression
f <- read.table("~/Box Sync/dinucleotide_repeat_motifs/supplementary_files/data/cis_bp_consensus.txt",header=T,sep="\t",stringsAsFactors = F)
ids <- read.table("~/Box Sync/dinucleotide_repeat_motifs/supplementary_files/data/cis_bp_id_map.txt",header=F,sep=" ",stringsAsFactors = F)
for (i in 1:nrow(ids)) {
  ids$name[i] <- strsplit(ids$V3[i],"_")[[1]][1]
  if (startsWith(ids$name[i],"(")) {# remove any parentheses
  ids$name[i] <- substring(ids$name[i],2,nchar(ids$name[i])-1)
  }
}
ids$Target.ID <- ids$V2
ids <- ids[c(4,5)]
f <- left_join(f,ids)
f$Target.ID <- f$name
t <- left_join(f,spec)
t[is.na(t)]<- 100
for (i in 1:nrow(t)){
  t$Fantom.ID[i] <- t$Target.ID[i]
  if (t$Target.ID[i] == "ALX1") {t$Tissue.specificity.score[i] <- 1.91
  t$Fantom.ID[i] <- "CART1"}
  if (t$Target.ID[i] == "BATF3") {t$Tissue.specificity.score[i] <- 0.34
  t$Fantom.ID[i] <- "SNFT"}
  if (t$Target.ID[i] == "BHLHA15") {t$Tissue.specificity.score[i] <- 2.86
  t$Fantom.ID[i] <- "BHLHB8"}
  if (t$Target.ID[i] == "BHLHE40") {t$Tissue.specificity.score[i] <- 0.36
  t$Fantom.ID[i] <- "BHLHB2"}
  if (t$Target.ID[i] == "BHLHE41") {t$Tissue.specificity.score[i] <- 0.89
  t$Fantom.ID[i] <- "BHLHB3"}
  if (t$Target.ID[i] == "BPTF") {t$Tissue.specificity.score[i] <- 0.51
  t$Fantom.ID[i] <- "FALZ"}
  if (t$Target.ID[i] == "CGBP") {t$Tissue.specificity.score[i] <- 0.07
  t$Fantom.ID[i] <- "CXXC1"}
  if (t$Target.ID[i] == "CUX1") {t$Tissue.specificity.score[i] <- 0.19
  t$Fantom.ID[i] <- "CUTL1"}
  if (t$Target.ID[i] == "EBF1") {t$Tissue.specificity.score[i] <- 0.73
  t$Fantom.ID[i] <- "EBF"}
  if (t$Target.ID[i] == "FOXO4") {t$Tissue.specificity.score[i] <- 0.84
  t$Fantom.ID[i] <- "MLLT7"}
  if (t$Target.ID[i] == "GSC2") {t$Tissue.specificity.score[i] <- 3.74
  t$Fantom.ID[i] <- "GSCL"}
  if (t$Target.ID[i] == "HDX") {t$Tissue.specificity.score[i] <- 1
  t$Fantom.ID[i] <- "CXorf43"}
  if (t$Target.ID[i] == "HINFP") {t$Tissue.specificity.score[i] <- 0.12
  t$Fantom.ID[i] <- "MIZF"}
  if (t$Target.ID[i] == "HLTF") {t$Tissue.specificity.score[i] <- 0.99
  t$Fantom.ID[i] <- "SMARCA3"}
  if (t$Target.ID[i] == "HNF1A") {t$Tissue.specificity.score[i] <- 1.63
  t$Fantom.ID[i] <- "TCF1"}
  if (t$Target.ID[i] == "HNF1B") {t$Tissue.specificity.score[i] <- 2.68
  t$Fantom.ID[i] <- "TCF2"}
  if (t$Target.ID[i] == "HOMEZ") {t$Tissue.specificity.score[i] <- 0.1
  t$Fantom.ID[i] <- "KIAA1443"}
  if (t$Target.ID[i] == "IRF9") {t$Tissue.specificity.score[i] <- 0.22
  t$Fantom.ID[i] <- "ISGF3G"}
  if (t$Target.ID[i] == "KDM2B") {t$Tissue.specificity.score[i] <- 0.29
  t$Fantom.ID[i] <- "FBXL10"}
  if (t$Target.ID[i] == "LCOR") {t$Tissue.specificity.score[i] <- 0.34
  t$Fantom.ID[i] <- "MLR2"}
  if (t$Target.ID[i] == "MECOM") {t$Tissue.specificity.score[i] <- 1.25
  t$Fantom.ID[i] <- "EVI1"}
  if (t$Target.ID[i] == "MZF1") {t$Tissue.specificity.score[i] <- 0.11
  t$Fantom.ID[i] <- "ZNF42"}
  if (t$Target.ID[i] == "NKX2-1") {t$Tissue.specificity.score[i] <- 3.35
  t$Fantom.ID[i] <- "TITF1"}
  if (t$Target.ID[i] == "PDX1") {t$Tissue.specificity.score[i] <- 4.72
  t$Fantom.ID[i] <- "IPF1"}
  if (t$Target.ID[i] == "RBPJ") {t$Tissue.specificity.score[i] <- 0.3
  t$Fantom.ID[i] <- "RBPSUH"}
  if (t$Target.ID[i] == "RFX6") {t$Tissue.specificity.score[i] <- 2.32
  t$Fantom.ID[i] <- "RFXDC1"}
  if (t$Target.ID[i] == "RHOXF1") {t$Tissue.specificity.score[i] <- 2.98
  t$Fantom.ID[i] <- "OTEX"}
  if (t$Target.ID[i] == "TGIF1") {t$Tissue.specificity.score[i] <- 0.34
  t$Fantom.ID[i] <- "TGIF"}
  if (t$Target.ID[i] == "TP63") {t$Tissue.specificity.score[i] <- 1.99
  t$Fantom.ID[i] <- "TP73L"}
  if (t$Target.ID[i] == "VSX2") {t$Tissue.specificity.score[i] <- 1.38
  t$Fantom.ID[i] <- "CHX10"}
  if (t$Target.ID[i] == "ZEB1") {t$Tissue.specificity.score[i] <- 0.92
  t$Fantom.ID[i] <- "TCF8"}
  if (t$Target.ID[i] == "ZFHX3") {t$Tissue.specificity.score[i] <- 0.47
  t$Fantom.ID[i] <- "ATBF1"}
  if (t$Target.ID[i] == "ZSCAN16") {t$Tissue.specificity.score[i] <- 0.25
  t$Fantom.ID[i] <- "ZNF435"}
}
cis_bp <- t[-which(t$Tissue.specificity.score == 100),]
cis_bp <- cis_bp[c(1,2,4,5)]

# merge 3 sets together
merged <- full_join(full_join(jaspar,encode),cis_bp) #removes any that are duplicate name/motif pair

# GC, CpG Content of motifs
for (i in 1:nrow(merged)) {
  n <- 0
  cpg <- 0
  mot <- strsplit(as.character(merged$Target.consensus[i]),"")
  for (j in 1:length(mot[[1]])) {
    if (mot[[1]][j] == "C" ) {n <- n+1
      if (j < length(mot[[1]])){if(mot[[1]][j+1] == "G"){cpg <- cpg + 1}}}
    if (mot[[1]][j] == "G" ) {n <- n+1}
  }
  merged$gc[i] <- n/length(mot[[1]])
  merged$cpg[i] <- cpg/length(mot[[1]])
}

# take mean GC/CPG content of motifs for each TF
t <- merged %>% group_by(Target.ID,Tissue.specificity.score) %>% summarise(m_cpg = mean(cpg))
final <- merged %>% group_by(Target.ID,Tissue.specificity.score) %>% summarise(m_gc = mean(gc))
final <- left_join(final,t)
for (i in 1:nrow(final)) {
  if (final$Tissue.specificity.score[i] >= 1) {final$expression[i] <- "specific"}
  if (final$Tissue.specificity.score[i] < 1) {final$expression[i] <- "broad"}
}

# stats, plot
wilcox.test(final[final$expression == "broad",]$m_gc,final[final$expression == "specific",]$m_gc)
ggplot(final,aes(x=expression,y=m_gc,fill=expression)) + geom_violin(color="darkgrey") + scale_fill_manual(values=c("royalblue","lightgreen")) + geom_boxplot(width=0.055,fill="white") +
  xlab("Expression") + ylab("Mean GC Content")+ theme(text = element_text(size=20),axis.title.x= element_text(size = 24),axis.title.y=element_text(size = 24))
wilcox.test(final[final$expression == "broad",]$m_cpg,final[final$expression == "specific",]$m_cpg)
ggplot(final,aes(x=expression,y=m_cpg,fill=expression)) + geom_violin(color="darkgrey") + scale_fill_manual(values=c("royalblue","lightgreen")) + geom_boxplot(width=0.055,fill="white") +
  xlab("Expression") + ylab("Mean CpG Content") + theme(text = element_text(size=20),axis.title.x= element_text(size = 24),axis.title.y=element_text(size = 24))

# regression frameworks
summary(lm(final$Tissue.specificity.score ~ final$m_gc))
spcor(final[c(2,3,4)])

# dummy-coded
final$dummy <- 0
final[final$expression == "broad",]$dummy <- 1
summary(lm(final$dummy ~ final$m_gc))
spcor(final[c(6,3,4)])
