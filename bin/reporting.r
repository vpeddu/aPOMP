library(ggplot2)
library(reshape2)
library(viridis)
library(wesanderson)
library(ggvegan)

setwd('/Users/vikas/Documents/UCSC/lab/kim_lab/evmeta/bin/')

reporting_stats = read.csv('/Users/vikas/Documents/UCSC/lab/kim_lab/evmeta/bin/zymo_reporting_stats.csv')
#reporting_stats = reporting_stats[-c(2),]
#reporting_stats$name = c('aPOMP','Kraken2','Megan-LR')
colnames(reporting_stats) = c('Name','Precision','Recall','F1')

reporting_stats_long = melt(reporting_stats)
reporting_stats_plot = ggplot(reporting_stats_long, aes(fill=variable, y=value, x=Name)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete=TRUE) +
  scale_y_continuous(expand = c(NA, 1), limits = c(0,1)) +
  theme_classic() + 
  #xlab('Pipeline') + 
  #ylab('Percent') + 
  theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position="bottom")
reporting_stats_plot

ggsave(plot = reporting_stats_plot, file = 'reporting_stats_plot.pdf', height = 3, width = 3)

rarefaction_file = as.data.frame(t(read.csv('/Users/vikas/Documents/UCSC/lab/kim_lab/evmeta/bin/raraefaction_stats.csv')))
rarefactions_long <- data.frame(matrix(ncol=3, dimnames=list(NULL, c('name','taxid','count'))))


for (i in 1:nrow(rarefaction_file)){ 
  name = rownames(rarefaction_file)[i]
  taxids = rarefaction_file$V1[i]
  split_taxids = strsplit(split = ",",noquote(gsub("\\[|\\]", "", taxids)))[[1]]

  counts = rarefaction_file$V2[i]
  split_counts = strsplit(split = ",",noquote(gsub("\\[|\\]", "", counts)))[[1]]
  
  tmp = data.frame(rep(name,length(split_counts)) , split_taxids, split_counts)
  colnames(tmp)<-c('name','taxid','count')
  rarefactions_long<-rbind(rarefactions_long,tmp)
  }




library(biomformat)

ass = rich_sparse_file = system.file("extdata", "/Users/vikas/Downloads/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_sup.biom'", 
                                     package = "biomformat")

x = read_biom('/Users/vikas/Downloads/Zymo-GridION-EVEN-BB-SN_Guppy_6.0.1_sup.biom')

