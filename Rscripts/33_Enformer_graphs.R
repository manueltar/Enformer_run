
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

merge_with_MPRA = function(option_list)
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
  
 
  #### READ MPRA_results ----
  
  MPRA_results<-as.data.frame(fread(file=opt$MPRA_results, sep="\t", header=T), stringsAsFactors=F)
  
  
  cat("MPRA_results_0\n")
  cat(str(MPRA_results))
  cat("\n")
  cat(str(unique(MPRA_results$VAR)))
  cat("\n")
  
  MPRA_results$carried_to_VAR<-paste('chr',MPRA_results$carried_variants,sep='')
  
  
  MPRA_results_subset<-MPRA_results[which(MPRA_results$carried_to_VAR == MPRA_results$VAR),]
  
  cat("MPRA_results_subset_0\n")
  cat(str(MPRA_results_subset))
  cat("\n")
  cat(str(unique(MPRA_results_subset$VAR)))
  cat("\n")
  
  MPRA_results_subset<-MPRA_results_subset[which(MPRA_results_subset$TILE == 'TILE_3'),]
  
  cat("MPRA_results_subset_0\n")
  cat(str(MPRA_results_subset))
  cat("\n")
  cat(str(unique(MPRA_results_subset$VAR)))
  cat("\n")
  
  indx.int<-c(which(colnames(MPRA_results_subset) == 'VAR'),
              which(colnames(MPRA_results_subset) == 'Cell_Type'),
              which(colnames(MPRA_results_subset) == 'TILE'),
              which(colnames(MPRA_results_subset) == 'ASE'),
              which(colnames(MPRA_results_subset) == 'Per_tile_experimental_class'))
  
 
  MPRA_results_subset_2<-unique(MPRA_results_subset[,indx.int])
  
  cat("MPRA_results_subset_2_0\n")
  cat(str(MPRA_results_subset_2))
  cat("\n")
  cat(str(unique(MPRA_results_subset_2$VAR)))
  cat("\n")
  
  #### READ Enformer_results ----
  
  Enformer_results<-readRDS(file=opt$Enformer_results)
  
  Enformer_results$chr<-paste('chr',Enformer_results$chr,sep='')
  

  cat("Enformer_results_0\n")
  cat(str(Enformer_results))
  cat("\n")
  
  
  #### K562 ----
  
  MPRA_results_subset_2_K562<-MPRA_results_subset_2[which(MPRA_results_subset_2$Cell_Type == 'K562'),]
  
  cat("MPRA_results_subset_2_K562_0\n")
  cat(str(MPRA_results_subset_2_K562))
  cat("\n")
  
  
  
  ### Z -score ASE
  
  mean_K562<-mean(MPRA_results_subset_2_K562$ASE)
  sd_K562<-sd(MPRA_results_subset_2_K562$ASE)
  
  cat("mean_K562 & sd_K562\n")
  cat(sprintf(as.character(mean_K562)))
  cat("\n")
  cat(sprintf(as.character(sd_K562)))
  cat("\n")
  
  MPRA_results_subset_2_K562$ASE_Z_score<-((MPRA_results_subset_2_K562$ASE-mean_K562)/sd_K562)
  
  cat("MPRA_results_subset_2_K562_1\n")
  cat(str(MPRA_results_subset_2_K562))
  cat("\n")
  
  MPRA_results_subset_2_K562<-merge(MPRA_results_subset_2_K562,
                                  Enformer_results, by=c('VAR'))
  
  cat("MPRA_results_subset_2_K562_1\n")
  cat(str(MPRA_results_subset_2_K562))
  cat("\n")
  
  indx.int<-c(which(colnames(MPRA_results_subset_2_K562) == 'VAR'),
              which(colnames(MPRA_results_subset_2_K562) == 'Cell_Type'),
              which(colnames(MPRA_results_subset_2_K562) == 'TILE'),
              which(colnames(MPRA_results_subset_2_K562) == 'ASE_Z_score'),
              which(colnames(MPRA_results_subset_2_K562) == 'mean_DNASE_K562'),
              which(colnames(MPRA_results_subset_2_K562) == 'TF_class'),
              which(colnames(MPRA_results_subset_2_K562) == 'Per_tile_experimental_class'))
  
  
  K562_df<-unique(MPRA_results_subset_2_K562[,indx.int])
  
 
  
 #### SAVE ----

  setwd(out)

  
  write.table(K562_df, file=paste("Enformer_comparison_",'K562',".tsv",sep=''),sep="\t", quote=F,row.names = F)
  saveRDS(K562_df, file=paste("Enformer_comparison_",'K562',".rds",sep=''))
  
  
  #### HL60 ----
  
  MPRA_results_subset_2_HL60<-MPRA_results_subset_2[which(MPRA_results_subset_2$Cell_Type == 'HL60'),]
  
  cat("MPRA_results_subset_2_HL60_0\n")
  cat(str(MPRA_results_subset_2_HL60))
  cat("\n")
  
  
  
  ### Z -score ASE
  
  mean_HL60<-mean(MPRA_results_subset_2_HL60$ASE)
  sd_HL60<-sd(MPRA_results_subset_2_HL60$ASE)
  
  cat("mean_HL60 & sd_HL60\n")
  cat(sprintf(as.character(mean_HL60)))
  cat("\n")
  cat(sprintf(as.character(sd_HL60)))
  cat("\n")
  
  MPRA_results_subset_2_HL60$ASE_Z_score<-((MPRA_results_subset_2_HL60$ASE-mean_HL60)/sd_HL60)
  
  cat("MPRA_results_subset_2_HL60_1\n")
  cat(str(MPRA_results_subset_2_HL60))
  cat("\n")
  
  MPRA_results_subset_2_HL60<-merge(MPRA_results_subset_2_HL60,
                                    Enformer_results, by=c('VAR'))
  
  cat("MPRA_results_subset_2_HL60_1\n")
  cat(str(MPRA_results_subset_2_HL60))
  cat("\n")
  
  indx.int<-c(which(colnames(MPRA_results_subset_2_HL60) == 'VAR'),
              which(colnames(MPRA_results_subset_2_HL60) == 'Cell_Type'),
              which(colnames(MPRA_results_subset_2_HL60) == 'TILE'),
              which(colnames(MPRA_results_subset_2_HL60) == 'ASE_Z_score'),
              which(colnames(MPRA_results_subset_2_HL60) == 'mean_DNASE_HL60'),
              which(colnames(MPRA_results_subset_2_HL60) == 'TF_class'),
              which(colnames(MPRA_results_subset_2_HL60) == 'Per_tile_experimental_class'))
  
  
  HL60_df<-unique(MPRA_results_subset_2_HL60[,indx.int])
  
  
  
  #### SAVE ----
  
  setwd(out)
  
  
  write.table(HL60_df, file=paste("Enformer_comparison_",'HL60',".tsv",sep=''),sep="\t", quote=F,row.names = F)
  saveRDS(HL60_df, file=paste("Enformer_comparison_",'HL60',".rds",sep=''))

}


corr.graphs = function(option_list)
{
  suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
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
  
  setwd(out)
  
  K562_df<-readRDS(file=paste("Enformer_comparison_",'K562',".rds",sep=''))
  
  cat("K562_df_0\n")
  cat(str(K562_df))
  cat("\n")
  cat(str(unique(K562_df$VAR)))
  cat("\n")
  
  K562_df.dt<-data.table(K562_df, key="TF_class")
  
  K562_cor_df<-as.data.frame(K562_df.dt[,.(Pearson_cor=cor(ASE_Z_score,mean_DNASE_K562, method="pearson"), n=.N), by=key(K562_df.dt)], stringsAsFactors=F)
  
  K562_cor_df$Cell_type<-"K562"
  
  K562_cor_df$MOCK_COORD<-1
  
  cat("K562_cor_df_0\n")
  cat(str(K562_cor_df))
  cat("\n")
  
  breaks.Rank<-seq(-1,1,by=0.5)
  labels.Rank<-as.character(breaks.Rank)
  
  cat("labels.Rank_0\n")
  cat(str(labels.Rank))
  cat("\n")
  
  # K562_cor_df %>%
  #   mutate(myaxis = paste0(K562_cor_df$TF_class, "\n", K562_cor_df$n,)) %>%
  #   mutate(myaxis=fct_reorder(myaxis,as.numeric(K562_cor_df$TF_class))) %>%

  
  ggheatmap <-ggplot(data=K562_cor_df,
                      aes(x=TF_class, y=MOCK_COORD, fill = Pearson_cor))+
                geom_tile(color = "white")+
                scale_fill_gradient2(name=paste("Pearson","Correlation","Coefficient", sep="\n"),
                                     low = "red", high = "green",mid="white",midpoint=0,
                                     na.value = NA,
                                     breaks=breaks.Rank,
                                     labels=labels.Rank,
                                     limits=c(breaks.Rank[1],
                                              breaks.Rank[length(breaks.Rank)]))+
                theme_minimal()+ # minimal theme
                scale_x_discrete(name=element_blank(), drop =F) +
                scale_y_discrete(name=NULL, drop=F)+
                theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1, size = 14, family="sans"))+
                coord_fixed()
  
  
  breaks.ASE_Z_score<-seq(-2,+2,by=0.5)
  labels.ASE_Z_score<-as.character(breaks.ASE_Z_score)
  
  cat("labels.ASE_Z_score_0\n")
  cat(str(labels.ASE_Z_score))
  cat("\n")
  
  breaks.mean_DNASE_K562<-seq(-50,+50,by=10)
  labels.mean_DNASE_K562<-as.character(breaks.mean_DNASE_K562)
  
  cat("labels.mean_DNASE_K562_0\n")
  cat(str(labels.mean_DNASE_K562))
  cat("\n")
  
  
  scatter_plot<-ggplot(data=K562_df,
                  aes(x=ASE_Z_score, 
                      y=mean_DNASE_K562,
                      color=TF_class))+
    geom_point(size=4)+
    theme_bw()+
    scale_x_continuous(name="ASE Z-score K562",
                       breaks=breaks.ASE_Z_score,labels=labels.ASE_Z_score, limits=c(breaks.ASE_Z_score[1],breaks.ASE_Z_score[length(breaks.ASE_Z_score)]))+
    scale_y_continuous(name="Enformer DNASE K562",
                       breaks=breaks.mean_DNASE_K562,labels=labels.mean_DNASE_K562, limits=c(breaks.mean_DNASE_K562[1],breaks.mean_DNASE_K562[length(breaks.mean_DNASE_K562)]))+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
          axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_color_manual(values=c('gray',"palevioletred2",'#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',"black",
                                "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F)+
    geom_hline(yintercept=0, color="black", linetype='dashed',size=1.5) +
    geom_vline(xintercept=0, color="black", linetype='dashed',size=1.5) +
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=16))+
    ggeasy::easy_center_title()
  
  graph_DEF<-plot_grid(scatter_plot,ggheatmap,
                       ncol = 1,
                       nrow=2,
                       rel_heights=c(1,1))
  
  svgname<-paste("Graph_corr_","K562",".svg", sep='')
  
  setwd(out)
  
  ggsave(svgname, plot= graph_DEF,
         device="svg",
         height=10, width=12)
  
 ############## 
  
  HL60_df<-readRDS(file=paste("Enformer_comparison_",'HL60',".rds",sep=''))
  
  cat("HL60_df_0\n")
  cat(str(HL60_df))
  cat("\n")
  cat(str(unique(HL60_df$VAR)))
  cat("\n")
  
  
  HL60_df.dt<-data.table(HL60_df, key="TF_class")
  
  HL60_cor_df<-as.data.frame(HL60_df.dt[,.(Pearson_cor=cor(ASE_Z_score,mean_DNASE_HL60, method="pearson"), n=.N), by=key(HL60_df.dt)], stringsAsFactors=F)
  
  HL60_cor_df$Cell_type<-"HL60"
  
  HL60_cor_df$MOCK_COORD<-1
  
  cat("HL60_cor_df_0\n")
  cat(str(HL60_cor_df))
  cat("\n")
  
  breaks.Rank<-seq(-1,1,by=0.5)
  labels.Rank<-as.character(breaks.Rank)
  
  cat("labels.Rank_0\n")
  cat(str(labels.Rank))
  cat("\n")
  
  # HL60_cor_df %>%
  #   mutate(myaxis = paste0(HL60_cor_df$TF_class, "\n", HL60_cor_df$n,)) %>%
  #   mutate(myaxis=fct_reorder(myaxis,as.numeric(HL60_cor_df$TF_class))) %>%
  
  
  ggheatmap <-ggplot(data=HL60_cor_df,
                     aes(x=TF_class, y=MOCK_COORD, fill = Pearson_cor))+
    geom_tile(color = "white")+
    scale_fill_gradient2(name=paste("Pearson","Correlation","Coefficient", sep="\n"),
                         low = "red", high = "green",mid="white",midpoint=0,
                         na.value = NA,
                         breaks=breaks.Rank,
                         labels=labels.Rank,
                         limits=c(breaks.Rank[1],
                                  breaks.Rank[length(breaks.Rank)]))+
    theme_minimal()+ # minimal theme
    scale_x_discrete(name=element_blank(), drop =F) +
    scale_y_discrete(name=NULL, drop=F)+
    theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1, size = 14, family="sans"))+
    coord_fixed()
  
  
  breaks.ASE_Z_score<-seq(-2,+2,by=0.5)
  labels.ASE_Z_score<-as.character(breaks.ASE_Z_score)
  
  cat("labels.ASE_Z_score_0\n")
  cat(str(labels.ASE_Z_score))
  cat("\n")
  
  breaks.mean_DNASE_HL60<-seq(-50,+50,by=10)
  labels.mean_DNASE_HL60<-as.character(breaks.mean_DNASE_HL60)
  
  cat("labels.mean_DNASE_HL60_0\n")
  cat(str(labels.mean_DNASE_HL60))
  cat("\n")
  
  
  scatter_plot<-ggplot(data=HL60_df,
                       aes(x=ASE_Z_score, 
                           y=mean_DNASE_HL60,
                           color=TF_class))+
    geom_point(size=4)+
    theme_bw()+
    scale_x_continuous(name="ASE Z-score HL60",
                       breaks=breaks.ASE_Z_score,labels=labels.ASE_Z_score, limits=c(breaks.ASE_Z_score[1],breaks.ASE_Z_score[length(breaks.ASE_Z_score)]))+
    scale_y_continuous(name="Enformer DNASE HL60",
                       breaks=breaks.mean_DNASE_HL60,labels=labels.mean_DNASE_HL60, limits=c(breaks.mean_DNASE_HL60[1],breaks.mean_DNASE_HL60[length(breaks.mean_DNASE_HL60)]))+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
          axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_color_manual(values=c('gray',"palevioletred2",'#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',"black",
                                "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F)+
    geom_hline(yintercept=0, color="black", linetype='dashed',size=1.5) +
    geom_vline(xintercept=0, color="black", linetype='dashed',size=1.5) +
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=16))+
    ggeasy::easy_center_title()
  
  graph_DEF<-plot_grid(scatter_plot,ggheatmap,
                       ncol = 1,
                       nrow=2,
                       rel_heights=c(1,1))
  
  svgname<-paste("Graph_corr_","HL60",".svg", sep='')
  
  setwd(out)
  
  ggsave(svgname, plot= graph_DEF,
         device="svg",
         height=10, width=12)
  
  
  #### SAVE ----
  
  DEF<-rbind(K562_cor_df,HL60_cor_df)
  
  
  setwd(out)
  
  write.table(DEF, file=paste("Enformer_DNASE_vs_ASE_comparison_",'correlation_values',".tsv",sep=''),sep="\t", quote=F,row.names = F)
  
  
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
    make_option(c("--MPRA_results"), type="character", default=NULL, 
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
  
  merge_with_MPRA(opt)
  corr.graphs(opt)
  
  
}


###########################################################################

system.time( main() )