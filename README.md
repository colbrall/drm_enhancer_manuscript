# Guide to Scripts and data for "Short DNA sequence patterns accurately identify broadly active human enhancers"
# Laura Colbran April 21, 2017
#

scripts/
        kmer_enrichment.py
                calculates fold enrichment of kmers between 2 fasta files
        drm_finder.py
                counts DRMs in a fasta file
        bin_enhancers.R
                bins enhancers by activity
        colbran_2017_drm_analyses.R
                code for TF mapping, pROC

data/
        all_fantom_enhancers.bed
                Broad enhancers = all with #tiss >45
                Context Specific = random subset of those with  #tiss = 1
                when directly comparing the two, I set all lengths in both to 600 bp (see Methods)
        dhs_broad.bed
                DHS-defined broad enhancers
        hist_broad.bed
                Histone-mark-defined broad enhancers
        kmer_enrichments.txt
                summary file of fold enrichments calculated using kmer_enrichment.py 
        cis_bp_consensus.txt
                CIS-BP IDs, consensus motifs
        cis_bp_id_map.txt
                map of CIS-BP IDs to TF names
        encode_tf_consensus.txt
                ENCODE IDs, consensus motifs
        jaspar_consensus.txt
                JASPAR IDs, consensus motifs
        tf_motif_specificity.csv
                FANTOM TSPS scores, IDs
        supp_table_1.txt
                contains the TF mapping- names, motifs, expression
        
classifiers/
        contains output and scripts from all SVM classifiers


Files relevant for main text figures:
Figure 1- Fold Enrichments of DRM counts from all_fantom_enhancers.bed 
       A- FE = mean(broad enhancer count)/mean(negative count)
          Obtained FEs from 4 random genomic background sets (see Methods)
          Took mean and SD of log2(FE)
       B- FE = mean(broad enhancer count)/mean(context-spec enhancer count)
       For both, significance as determined by Wilcoxon Rank Sum on full count distributions.
       plotted using ggplot2 geom_bar
       
Figure 2- DRM ROC curve, counts
       A- scripts and curves in classifiers/2016-11-30*
          features = count/bp for all 4 DRMs
       B- ggplot geom_violin, scaled by area, of counts in all_fantom_enhancers.bed

Figure 3- 6-mer and TF ROC curves
       A- scripts, curves in classifiers/2016-12-08*
       B- scripts, curves in classifiers/2017-04-24*
          features = count/bp for all motifs in CIS-BP, JASPAR, and ENCODE databases 
          see R script for stats code
       C- data in broad_fantom_kmer_enrichments.txt
          as in 1A, 'Negatives' and 'GC-matched' are mean FE over 4 random sets
          ggplot2 geom_violin, scaled by area

Figure 4- GC content/Enrichment
       A,B,D- data in broad_fantom_kmer_enrichments.txt
              binned 6-mers by  GC content
              ggplot2 geom_violin, scaled by area
              used Spearman's rho to find correlation between GC content and enrichment
       C- data in all_fantom_enhancers.bed
          binned using bin_enhancers.R
          plotted GC content vs bin 
