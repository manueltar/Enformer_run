
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

interpret_enformer = function(option_list)
{

  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")

  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform K562_features and HL60_features ----
  
  K562_features = unlist(strsplit(opt$K562_features, split=","))
  
  cat("K562_features_\n")
  cat(sprintf(as.character(K562_features)))
  cat("\n")
  
  HL60_features = unlist(strsplit(opt$HL60_features, split=","))
  
  cat("HL60_features_\n")
  cat(sprintf(as.character(HL60_features)))
  cat("\n")
  
  #### READ Table_S6 ----
  
  Table_S6<-as.data.frame(fread(file=opt$Table_S6, sep="\t", header=T), stringsAsFactors=F)
  
  
  cat("Table_S6_0\n")
  cat(str(Table_S6))
  cat("\n")
  cat(str(unique(Table_S6$VAR_38)))
  cat("\n")
 
  Table_S6$chr<-gsub("_.+$","",Table_S6$VAR_38)
  Table_S6$chr<-gsub("^chr","",Table_S6$chr)
  Table_S6$pos38<-gsub("^[^_]+_","",Table_S6$VAR_38)
  Table_S6$pos38<-as.integer(gsub("_.+$","",Table_S6$pos38))
  Table_S6$ref<-gsub("^[^_]+_[^_]+_","",Table_S6$VAR_38)
  Table_S6$ref<-gsub("_.+$","",Table_S6$ref)
  Table_S6$alt<-gsub("^[^_]+_[^_]+_[^_]+_","",Table_S6$VAR_38)
 

  cat("Table_S6_1\n")
  cat(str(Table_S6))
  cat("\n")
  cat(str(unique(Table_S6$chr)))
  cat("\n")
  cat(str(unique(Table_S6$pos38)))
  cat("\n")
  cat(str(unique(Table_S6$ref)))
  cat("\n")
  cat(str(unique(Table_S6$alt)))
  cat("\n")
  
  #### READ Enformer_result ----
  
  Enformer_result<-as.data.frame(fread(file=opt$Enformer_result, sep="\t", header=T), stringsAsFactors=F)
  
  Enformer_result$VAR_38<-paste(Enformer_result$chr,Enformer_result$pos,Enformer_result$ref,Enformer_result$alt,sep='_')
  
  cat("Enformer_result_0\n")
  cat(str(Enformer_result))
  cat("\n")
  
  #indx.const<-c(which(colnames(Enformer_result) == 'chrom'),which(colnames(Enformer_result) == 'pos'),which(colnames(Enformer_result) == 'id'),which(colnames(Enformer_result) == 'ref'),which(colnames(Enformer_result) == 'alt'))
  indx.const<-c(which(colnames(Enformer_result) == 'VAR_38'))
  
  indx.col<-grep(paste(K562_features, collapse="|"), colnames(Enformer_result))
  
  Enformer_result_subset_K562_matrix<-as.matrix(Enformer_result[,c(indx.col)])
  row.names(Enformer_result_subset_K562_matrix)<-Enformer_result$VAR_38
  
  cat("Enformer_result_subset_K562_matrix_0\n")
  cat(str(Enformer_result_subset_K562_matrix))
  cat("\n")
  
  mean_matrix_K562<-apply(Enformer_result_subset_K562_matrix,1,mean)
  
  cat("mean_matrix_K562_0\n")
  cat(str(mean_matrix_K562))
  cat("\n")
  
  mean_df_K562<-as.data.frame(mean_matrix_K562)
  mean_df_K562$VAR_38<-row.names(mean_df_K562)
  row.names(mean_df_K562)<-NULL
  colnames(mean_df_K562)[which(colnames(mean_df_K562) == 'mean_matrix_K562')]<-'mean_DNASE_K562'
  
  cat("mean_df_K562_0\n")
  cat(str(mean_df_K562))
  cat("\n")
  
  ###
  
  indx.col<-grep(paste(HL60_features, collapse="|"), colnames(Enformer_result))
  
  Enformer_result_subset_HL60_matrix<-as.matrix(Enformer_result[,c(indx.col)])
  row.names(Enformer_result_subset_HL60_matrix)<-Enformer_result$VAR_38
  
  cat("Enformer_result_subset_HL60_matrix_0\n")
  cat(str(Enformer_result_subset_HL60_matrix))
  cat("\n")
  
  mean_matrix_HL60<-apply(Enformer_result_subset_HL60_matrix,1,mean)
  
  cat("mean_matrix_HL60_0\n")
  cat(str(mean_matrix_HL60))
  cat("\n")
  
  mean_df_HL60<-as.data.frame(mean_matrix_HL60)
  mean_df_HL60$VAR_38<-row.names(mean_df_HL60)
  row.names(mean_df_HL60)<-NULL
  colnames(mean_df_HL60)[which(colnames(mean_df_HL60) == 'mean_matrix_HL60')]<-'mean_DNASE_HL60'
  
  cat("mean_df_HL60_0\n")
  cat(str(mean_df_HL60))
  cat("\n")
  
  ### Merges
  
  mean_df_K562<-merge(mean_df_K562,
                      mean_df_HL60,
                      by='VAR_38')
  
  cat("mean_df_K562_1\n")
  cat(str(mean_df_K562))
  cat("\n")
  cat(str(unique(mean_df_K562$VAR_38)))
  cat("\n")
  
 
 
 Table_S6<-merge(mean_df_K562,
                 Table_S6,
                 by='VAR_38')
 
 cat("Table_S6_0\n")
 cat(str(Table_S6))
 cat("\n")
 cat(str(unique(Table_S6$VAR_38)))
 cat("\n")
  
 #### SAVE ----

  setwd(out)

  
  write.table(Table_S6, file=paste("Enformer_summary_",'K562','_','HL60',".tsv",sep=''),sep="\t", quote=F,row.names = F)
  
  

}

TF_motif_annotation = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform context_initiator_TFs_subset ----
  
  context_initiator_TFs_subset = unlist(strsplit(opt$context_initiator_TFs_subset, split=","))
  
  cat("context_initiator_TFs_subset_0\n")
  cat(sprintf(as.character(context_initiator_TFs_subset)))
  cat("\n")
  
  #### READ and transform context_only_TFs_subset ----
  
  context_only_TFs_subset = unlist(strsplit(opt$context_only_TFs_subset, split=","))
  
  cat("context_only_TFs_subset_0\n")
  cat(sprintf(as.character(context_only_TFs_subset)))
  cat("\n")
  
  #### READ and transform version_string ----
  
  version_string = unlist(strsplit(opt$version_string, split=","))
  
  cat("version_string_0\n")
  cat(sprintf(as.character(version_string)))
  cat("\n")
  
  
  #### READ and transform Intersect_SNP ----
  
  Intersect_SNP = opt$Intersect_SNP
  
  cat("Intersect_SNP_0\n")
  cat(sprintf(as.character(Intersect_SNP)))
  cat("\n")
  
  
  #### READ TF_motifs_annotated ----
  
  TF_motifs_annotated<-as.data.frame(fread(file=opt$TF_motifs_annotated, sep="\t", header=T), stringsAsFactors=F)
  
  
  cat("TF_motifs_annotated_0\n")
  cat(str(TF_motifs_annotated))
  cat("\n")
  cat(str(unique(TF_motifs_annotated$rs)))
  cat("\n")
  
  ###### comparison ----
  
  
  df_comparison<-unique(TF_motifs_annotated[which(TF_motifs_annotated$version%in%version_string &
                                                TF_motifs_annotated$Intersect_SNP%in%Intersect_SNP),])
  
  
  cat("df_comparison_0\n")
  cat(str(df_comparison))
  cat("\n")
  cat(str(unique(df_comparison$rs)))
  cat("\n")
  
  
  
  
  
  indx.context_initiator<-grep(paste(context_initiator_TFs_subset,collapse="|"),df_comparison$Motif_ID)
  
  indx.context_only<-grep(paste(context_only_TFs_subset,collapse="|"),df_comparison$Motif_ID)
  
  
  df_context_initiator<-df_comparison[indx.context_initiator,]
  df_context_initiator$TF_class<-'context_initiator'
  
  cat("df_context_initiator_0\n")
  cat(str(df_context_initiator))
  cat("\n")
  cat(str(unique(df_context_initiator$rs)))
  cat("\n")
  
  df_context_only<-df_comparison[indx.context_only,]
  df_context_only$TF_class<-'context_only'
  
  cat("df_context_only_0\n")
  cat(str(df_context_only))
  cat("\n")
  cat(str(unique(df_context_only$rs)))
  cat("\n")
  
  
  indx.ambivalent<-which(df_context_initiator$rs%in%df_context_only$rs)
  
  cat("indx.ambivalent_0\n")
  cat(str(indx.ambivalent))
  cat("\n")
  
  cat(sprintf(as.character(unique(df_context_initiator$rs[indx.ambivalent]))))
  cat("\n")
  
  df_context_initiator<-df_context_initiator[-indx.ambivalent,]
  
  cat("df_context_initiator_1\n")
  cat(str(df_context_initiator))
  cat("\n")
  cat(str(unique(df_context_initiator$rs)))
  cat("\n")
  
  
  
  
  ALL_TF<-rbind(df_context_initiator,df_context_only)
  
  cat("ALL_TF_0\n")
  cat(str(ALL_TF))
  cat("\n")
  
  ALL_TF.dt<-data.table(ALL_TF, key=c('rs'))
  
  
  ALL_TF_collapsed<-as.data.frame(ALL_TF.dt[,.(TF_class=paste(unique(TF_class), collapse='|')), by=key(ALL_TF.dt)])
  
  cat("ALL_TF_collapsed_0\n")
  cat(str(ALL_TF_collapsed))
  cat("\n")
  cat(str(unique(ALL_TF_collapsed$rs)))
  cat("\n")
  
  
 
  #### READ Enformer_results ----
  
  Enformer_results<-as.data.frame(fread(file=opt$Enformer_results, sep="\t", header=T), stringsAsFactors=F)
  
  
  cat("Enformer_results_0\n")
  cat(str(Enformer_results))
  cat("\n")
  
 
  ##### Merge, create factors and save -----
  
  
  Enformer_results<-merge(Enformer_results,
                          ALL_TF_collapsed,
                          by='rs',
                          all.x=T)
  
  Enformer_results$TF_class[is.na(Enformer_results$TF_class)]<-'NOT_ANNOTATED'
  
  Enformer_results$TF_class<-factor(Enformer_results$TF_class,
                                    levels=c('NOT_ANNOTATED','context_initiator','context_only'),
                                    ordered=T)
  
  cat("Enformer_results_1\n")
  cat(str(Enformer_results))
  cat("\n")
  cat(sprintf(as.character(names(summary(Enformer_results$TF_class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Enformer_results$TF_class))))
  cat("\n")
  
  #### SAVE ----
  
  setwd(out)
  
  
  write.table(Enformer_results, file=paste("Enformer_summary_",'K562','_','HL60','_','TF_class',".tsv",sep=''),sep="\t", quote=F,row.names = F)
  
  saveRDS(Enformer_results, file=paste("Enformer_summary_",'K562','_','HL60','_','TF_class',".rds", sep=''))
  
 
  
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--Table_S6"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Enformer_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_features"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--HL60_features"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--context_initiator_TFs_subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--context_only_TFs_subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--version_string"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Intersect_SNP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_motifs_annotated"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Enformer_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  interpret_enformer(opt)
  TF_motif_annotation(opt)
  
  
}


###########################################################################

system.time( main() )