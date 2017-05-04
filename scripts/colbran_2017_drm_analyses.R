#colbran_2017_drm_analyses.R
# code used to run analyses and plot things for 
#"Short DNA sequence patterns accurately identify broadly active human enhancers"
library("ggplot2")
library("dplyr")
library("pROC")
library('readr')
library('ppcor')

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

#-----------------------------------------------#
# TF GC,CpG content, Expression
spec <- read.csv("tf_motif_specificity.csv")
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
f <- read.delim('jaspar_consensus.txt',sep="\t",header=T,stringsAsFactors = F)
t <- left_join(f,spec)
t[is.na(t)]<- 100
for (i in 1:nrow(t)) {
  if (t$Target.ID[i] == "BATF::JUN") {t$Tissue.specificity.score[i] <- 1.14} 
  if (t$Target.ID[i] == "EBF1") {t$Tissue.specificity.score[i] <- 0.73} 
  if (t$Target.ID[i] == "EWSR1-FLI1") {t$Tissue.specificity.score[i] <- 0.81} 
  if (t$Target.ID[i] == "HINFP") {t$Tissue.specificity.score[i] <- 0.12}
  if (t$Target.ID[i] == "JUN(var.2)") {t$Tissue.specificity.score[i] <- 0.45}
  if (t$Target.ID[i] == "JUND(var.2)") {t$Tissue.specificity.score[i] <- 0.18}
  if (t$Target.ID[i] == "MYC::MAX") {t$Tissue.specificity.score[i] <- 0.58}
  if (t$Target.ID[i] == "NFE2::MAF") {t$Tissue.specificity.score[i] <- 3.12}
  if (t$Target.ID[i] == "NR1H2::RXRA") {t$Tissue.specificity.score[i] <- 0.20}
  if (t$Target.ID[i] == "Pax6") {t$Tissue.specificity.score[i] <- 2.66}
  if (t$Target.ID[i] == "RXRA::VDR") {t$Tissue.specificity.score[i] <- 1.08}
  if (t$Target.ID[i] == "SMAD2::SMAD3::SMAD4") {t$Tissue.specificity.score[i] <- 0.20}
  if (t$Target.ID[i] == "STAT2::STAT1") {t$Tissue.specificity.score[i] <- 0.33}
  if (t$Target.ID[i] == "TAL1::GATA1") {t$Tissue.specificity.score[i] <- 3.76}
  if (t$Target.ID[i] == "TAL1::TCF3") {t$Tissue.specificity.score[i] <- 2.37}
  if (t$Target.ID[i] == "TLX1::NFIC") {t$Tissue.specificity.score[i] <- 3.68}
  if (t$Target.ID[i] == "TP63") {t$Tissue.specificity.score[i] <- 1.99}
  if (t$Target.ID[i] == "ZEB1") {t$Tissue.specificity.score[i] <- 0.92}
}
jaspar <- t[-which(t$Tissue.specificity.score == 100),]


# ENCODE TFs
f <- read.table("encode_tf_consensus.txt",header=T,sep="\t",stringsAsFactors = F)
for (i in 1:nrow(f)) {f$Target.ID[i] <- strsplit(strsplit(f$Target.ID[i],">")[[1]][2],"_")[[1]][1]}
t <- left_join(f,spec)
t[is.na(t)]<- 100
for (i in 1:nrow(t)) {
  if (t$Target.ID[i] == "AHR::ARNT" | t$Target.ID[i] == "AHR::ARNT::HIF1A") {t$Tissue.specificity.score[i] <- 0.96}
  if (t$Target.ID[i] == "BHLHE40") {t$Tissue.specificity.score[i] <- 0.36}
  if (t$Target.ID[i] == "BHLHE41") {t$Tissue.specificity.score[i] <- 0.89}
  if (t$Target.ID[i] == "BPTF") {t$Tissue.specificity.score[i] <- 0.51}
  if (t$Target.ID[i] == "CLOCK::ARNTL") {t$Tissue.specificity.score[i] <- 0.60}
  if (t$Target.ID[i] == "CUX1") {t$Tissue.specificity.score[i] <- 0.19}
  if (t$Target.ID[i] == "DDIT3::CEBPA") {t$Tissue.specificity.score[i] <- 1.4}
  if (t$Target.ID[i] == "EBF1") {t$Tissue.specificity.score[i] <- 0.73}
  if (t$Target.ID[i] == "EWSR1::FLI1") {t$Tissue.specificity.score[i] <- 0.85}
  if (t$Target.ID[i] == "GSC2") {t$Tissue.specificity.score[i] <- 3.74}
  if (t$Target.ID[i] == "HDX") {t$Tissue.specificity.score[i] <- 1}
  if (t$Target.ID[i] == "HIF1A::ARNT") {t$Tissue.specificity.score[i] <- 0.39}
  if (t$Target.ID[i] == "HINFP") {t$Tissue.specificity.score[i] <- 0.12}
  if (t$Target.ID[i] == "HLTF") {t$Tissue.specificity.score[i] <- 0.99}
  if (t$Target.ID[i] == "HNF1A" | t$Target.ID[i] == "HNF1") {t$Tissue.specificity.score[i] <- 1.63}
  if (t$Target.ID[i] == "HNF1B") {t$Tissue.specificity.score[i] <- 2.68}
  if (t$Target.ID[i] == "HNF4") {t$Tissue.specificity.score[i] <- 3.46}
  if (t$Target.ID[i] == "HOMEZ") {t$Tissue.specificity.score[i] <- 0.1}
  if (t$Target.ID[i] == "IKZF1") {t$Tissue.specificity.score[i] <- 1.59}
  if (t$Target.ID[i] == "IKZF2") {t$Tissue.specificity.score[i] <- 0.66}
  if (t$Target.ID[i] == "IKZF3") {t$Tissue.specificity.score[i] <- 2.49}
  if (t$Target.ID[i] == "MEIS1::HOXA9") {t$Tissue.specificity.score[i] <- 1.41}
  if (t$Target.ID[i] == "MZF1") {t$Tissue.specificity.score[i] <- 0.11}
  if (t$Target.ID[i] == "NFE2L1::MAFG") {t$Tissue.specificity.score[i] <- 0.56}
  if (t$Target.ID[i] == "NKX2-1") {t$Tissue.specificity.score[i] <- 3.35}
  if (t$Target.ID[i] == "PATZ1") {t$Tissue.specificity.score[i] <- 0.16}
  if (t$Target.ID[i] == "PDX1") {t$Tissue.specificity.score[i] <- 4.72}
  if (t$Target.ID[i] == "RBPJ") {t$Tissue.specificity.score[i] <- 0.3}
  if (t$Target.ID[i] == "RHOXF1") {t$Tissue.specificity.score[i] <- 2.98}
  if (t$Target.ID[i] == "RHOXF2") {t$Tissue.specificity.score[i] <- 3.63}
  if (t$Target.ID[i] == "TFAP2") {t$Tissue.specificity.score[i] <- 3.07}
  if (t$Target.ID[i] == "TGIF1") {t$Tissue.specificity.score[i] <- 0.34}
  if (t$Target.ID[i] == "VSX2") {t$Tissue.specificity.score[i] <- 1.38}
  if (t$Target.ID[i] == "ZBTB14") {t$Tissue.specificity.score[i] <- 0.15}
  if (t$Target.ID[i] == "ZBTB18") {t$Tissue.specificity.score[i] <- 1.58}
  if (t$Target.ID[i] == "ZEB1") {t$Tissue.specificity.score[i] <- 0.92}
  if (t$Target.ID[i] == "ZSCAN16") {t$Tissue.specificity.score[i] <- 0.25}
}
encode <- t[-which(t$Tissue.specificity.score == 100),]

# CIS-BP TFs
# start by mapping CIS-BP ids to TF names, then match to FANTOM expression
f <- read.table("cis_bp_consensus.txt",header=T,sep="\t",stringsAsFactors = F)
ids <- read.table("cis_bp_id_map.txt",header=F,sep=" ",stringsAsFactors = F)
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
  if (t$Target.ID[i] == "ALX1") {t$Tissue.specificity.score[i] <- 1.91}
  if (t$Target.ID[i] == "BATF3") {t$Tissue.specificity.score[i] <- 0.34}
  if (t$Target.ID[i] == "BHLHA15") {t$Tissue.specificity.score[i] <- 2.86}
  if (t$Target.ID[i] == "BHLHE40") {t$Tissue.specificity.score[i] <- 0.36}
  if (t$Target.ID[i] == "BHLHE41") {t$Tissue.specificity.score[i] <- 0.89}
  if (t$Target.ID[i] == "BPTF") {t$Tissue.specificity.score[i] <- 0.51}
  if (t$Target.ID[i] == "CGBP") {t$Tissue.specificity.score[i] <- 0.07}
  if (t$Target.ID[i] == "CUX1") {t$Tissue.specificity.score[i] <- 0.19}
  if (t$Target.ID[i] == "EBF1") {t$Tissue.specificity.score[i] <- 0.73}
  if (t$Target.ID[i] == "FOXO4") {t$Tissue.specificity.score[i] <- 0.84}
  if (t$Target.ID[i] == "GSC2") {t$Tissue.specificity.score[i] <- 3.74}
  if (t$Target.ID[i] == "HDX") {t$Tissue.specificity.score[i] <- 1}
  if (t$Target.ID[i] == "HINFP") {t$Tissue.specificity.score[i] <- 0.12}
  if (t$Target.ID[i] == "HLTF") {t$Tissue.specificity.score[i] <- 0.99}
  if (t$Target.ID[i] == "HNF1A") {t$Tissue.specificity.score[i] <- 1.63}
  if (t$Target.ID[i] == "HNF1B") {t$Tissue.specificity.score[i] <- 2.68}
  if (t$Target.ID[i] == "HOMEZ") {t$Tissue.specificity.score[i] <- 0.1}
  if (t$Target.ID[i] == "IRF9") {t$Tissue.specificity.score[i] <- 0.22}
  if (t$Target.ID[i] == "KDM2B") {t$Tissue.specificity.score[i] <- 0.29}
  if (t$Target.ID[i] == "LCOR") {t$Tissue.specificity.score[i] <- 0.34}
  if (t$Target.ID[i] == "MECOM") {t$Tissue.specificity.score[i] <- 1.25}
  if (t$Target.ID[i] == "MZF1") {t$Tissue.specificity.score[i] <- 0.11}
  if (t$Target.ID[i] == "NKX2-1") {t$Tissue.specificity.score[i] <- 3.35}
  if (t$Target.ID[i] == "PDX1") {t$Tissue.specificity.score[i] <- 4.72}
  if (t$Target.ID[i] == "RBPJ") {t$Tissue.specificity.score[i] <- 0.3}
  if (t$Target.ID[i] == "RFX6") {t$Tissue.specificity.score[i] <- 2.32}
  if (t$Target.ID[i] == "RHOXF1") {t$Tissue.specificity.score[i] <- 2.98}
  if (t$Target.ID[i] == "TGIF1") {t$Tissue.specificity.score[i] <- 0.34}
  if (t$Target.ID[i] == "TP63") {t$Tissue.specificity.score[i] <- 1.99}
  if (t$Target.ID[i] == "VSX2") {t$Tissue.specificity.score[i] <- 1.38}
  if (t$Target.ID[i] == "ZEB1") {t$Tissue.specificity.score[i] <- 0.92}
  if (t$Target.ID[i] == "ZFHX3") {t$Tissue.specificity.score[i] <- 0.47}
  if (t$Target.ID[i] == "ZSCAN16") {t$Tissue.specificity.score[i] <- 0.25}
}
cis_bp <- t[-which(t$Tissue.specificity.score == 100),]
cis_bp <- cis_bp[c(1,2,4)]

# merge 3 sets together
merged <- full_join(full_join(jaspar,encode),cis_bp) #removes any that are duplicate name/motif

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

final$dummy <- 0
final[final$expression == "broad",]$dummy <- 1
summary(lm(final$dummy ~ final$m_gc))
spcor(final[c(6,3,4)])

