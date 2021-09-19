packages = c("readxl", "stringr", "reshape2", "ggplot2", "ggprism", 
             "magrittr", "rstatix", "ggpubr", "optparse", "data.table")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "https://cloud.r-project.org/")
      library(x, character.only = TRUE)
    }
  }
)

package.check

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar=".xlsx"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name for stats", metavar=".csv"),
  make_option(c("-c", "--control"), type="character", default="N2a Control",
              help="name of control sample to be used", metavar="Control"),
  make_option(c("-p", "--pdf"), type="character", default=NULL,
              help="output file name for PDF", metavar=".pdf"),
  make_option(c("-w", "--width"), type="integer", default=10,
              help="width of PDF", metavar="integer"),
  make_option(c("-t", "--height"), type="integer", default=5,
              help="height of PDF", metavar="integer"),
  make_option(c("-g", "--gene"), type="character", default="Gapdh",
              help="gene to be used for comparison")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

dat<-read_excel(opt$file, skip = 6)
dat<-dat[1:which(is.na(dat))[1]-1,]
dat<-dat[,1:3]
dat[dat == "Undetermined"] <- NA
genes<-unique(dat$Target)
list<-lapply(genes, function(x){
  tmp<-dat[dat$Target == x,]
})

names(list)<-genes

list<-lapply(genes, function(x){
  ifelse(anyNA(list[[x]]$Cт, recursive = T), return(), return(list[[x]]))
})

dat<-do.call(rbind, list)
genes<-unique(dat$Target)
dat$Cт<-as.numeric(dat$Cт)
list<-lapply(genes, function(x){
  tmp<-dat[dat$Target == x,]
})

names(list)<-genes

list<-lapply(genes, function(x){
  tmp2<-list[[x]]$Cт - list[[opt$gene]]$Cт
  list[[x]]$delta_ct<-tmp2
  list[[x]]$rel_qty<-2^-list[[x]]$delta_ct
  return(list[[x]])
})

dat<-do.call(rbind, list)
dat$rel_gapdh<-dat$rel_qty/dat[dat$Target %in% opt$gene,]$rel_qty
dat<-reshape2::melt(dat)
dat<-dat[dat$variable == "rel_gapdh",]
dat$`Sample Name`
dat<-dat[dat$Target != "Pantr2",]
dat$`Sample Name` %<>%
  gsub(" 1", "", .) %>% gsub(" 2", "", .) %>% gsub(" 3", "", .)
colnames(dat)<-c("sample", "target", "variable", "value")
dat<-dat[,-3]
dat<-dat[dat$target != opt$gene,]
dat$target<-as.factor(dat$target)
dodge<-position_dodge(width = 0.8)

pwc <- group_by(dat, target) %>%
  pairwise_t_test(value ~ sample, ref.group = opt$control)
pwc <- pwc %>% add_xy_position(x = "target")
stat.test <- pwc %>% add_xy_position(x = "target")

pdf(file = opt$pdf, width=opt$width , height=opt$height )
ggbarplot(dat, x = "target", y = "value", 
          add = c("mean_sd", "point"),
          color = "black", palette = "jco",
          fill = "sample",
          position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype = 2, size = 1.0) + 
  ylab(paste0("Expression relative to ", opt$gene)) +
  xlab("Gene Name") +
  labs(fill = "Sample Name") +
  theme_prism() +
  stat_pvalue_manual(stat.test, label="p.signif", hide.ns = T)
dev.off()

fwrite(stat.test, file = opt$out)

