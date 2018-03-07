### Guide to Scripts and data for "Short DNA sequence patterns accurately identify broadly active human enhancers"
##### Laura Colbran April 21, 2017
###### https://doi.org/10.1186/s12864-017-3934-9
**scripts/**  
kmer_enrichment.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;calculates fold enrichment of kmers between 2 fasta files  
drm_finder.py  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;counts DRMs in a fasta file  
bin_enhancers.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;bins enhancers by activity  
colbran_2017_drm_analyses.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;code for DRM fold enrichments, TF mapping, pROC  

**data/**  
all_fantom_enhancers.bed  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Broad enhancers = all with #tiss >45  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Context Specific = random subset of those with #tiss = 1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;when directly comparing the two, I set all lengths in both to 600 bp  
dhs_broad.bed  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;DHS-defined broad enhancers  
hist_broad.bed  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Histone-mark-defined broad enhancers  
broad_fantom_kmer_enrichments.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;summary file of fold enrichments calculated using kmer_enrichment.py  
broad_fantom_kmer_weights.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;summary file of weights given to all 6-mers by all SVM classifiers  
cis_bp_consensus.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;CIS-BP IDs, consensus motifs  
cis_bp_id_map.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;map of CIS-BP IDs to TF names  
encode_tf_consensus.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ENCODE IDs, consensus motifs  
jaspar_consensus.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;JASPAR IDs, consensus motifs  
tf_motif_specificity.csv  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;FANTOM TSPS scores, IDs  
supp_table_1.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;contains the TF mapping- names, motifs, expression  

**classifiers/**  
contains output and scripts from all SVM classifiers

**Files relevant for main text figures:**

**Figure 1**- Fold Enrichments of DRM counts from all_fantom_enhancers.bed  
       _See colbran_2017_drm_analyses.R for code_  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**A**- FE = mean(broad enhancer count)/mean(negative count)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Obtained FEs from 4 random genomic background sets  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*see Methods*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Took mean and SD of log2(FE)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**B**- FE = mean(broad enhancer count)/mean(context-spec enhancer count)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For both, significance as determined by Wilcoxon Rank Sum on full count distributions.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;plotted using ggplot2 geom_bar  

**Figure 2**- DRM ROC curve, counts  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**A**- scripts and curves in classifiers/2016-11-30*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;features = count/bp for all 4 DRMs  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**B**- ggplot geom_violin, scaled by area, of counts in all_fantom_enhancers.bed  

**Figure 3**- 6-mer and TF ROC curves  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**A**- scripts, curves in classifiers/2016-12-08*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**B**- scripts, curves in classifiers/2017-04-24*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;features = count/bp for all motifs in CIS-BP, JASPAR, and ENCODE databases  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_see colbran_2017_drm_analyses.R for stats code_  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**C**- data in broad_fantom_kmer_enrichments.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;as in 1A, 'Negatives' and 'GC-matched' are mean FE over 4 random sets  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ggplot2 geom_violin, scaled by area

**Figure 4**- GC content/Enrichment  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**A,B,D**- data in broad_fantom_kmer_enrichments.txt  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;binned 6-mers by  GC content  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ggplot2 geom_violin, scaled by area  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;used Spearman's rho to find correlation between GC content and enrichment  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**C**- data in all_fantom_enhancers.bed  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;binned using bin_enhancers.R  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;plotted GC content vs bin

**Figure 5**- TF motifs  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**A**- GC content  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**B**- CpG Content  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;_see colbran_2017_drm_analyses.R for code_
