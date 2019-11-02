library(readr)
library(ggplot2)

source("similarity_scores.R")
source("coexpression_comparison.R")

set.seed(13456)

file_name <- c("All_Tissue_Site_Details.combined.reads.gct")
# read in raw expression matrix
gct_file <- read.table(file_name, skip = 2, header = TRUE,check.names=FALSE, 
                       row.names=1)

gct_file[is.na(gct_file)] <- 0
 
#remove lowly expressed genes
#gct_file <- gct_file[rowMeans(gct_file) > ncol(gct_file),]


#reduce to genes of interest and samples from tissues of interest
GTEx_sampleinfo <- read_delim("GTEx_sampleinfo.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
GTEx_sampleinfo <- as.data.frame(GTEx_sampleinfo)


#extract set of genes to construct networks with
heart_genes <- read_delim("heart_genes.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
brain_genes <- read_delim("brain_genes.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
skeletalmuscle_genes <- read_delim("skeletalmuscle_genes.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)

# remove any gene shared across these tables
remove <- unique(c(intersect(heart_genes$`Gene name`, skeletalmuscle_genes$`Gene name`),
                   intersect(heart_genes$`Gene name`, brain_genes$`Gene name`),
                   intersect(brain_genes$`Gene name`, skeletalmuscle_genes$`Gene name`)))

heart_genes <- heart_genes[!(heart_genes$`Gene name` %in% remove),]
brain_genes <- brain_genes[!(brain_genes$`Gene name` %in% remove),]
skeletalmuscle_genes <- skeletalmuscle_genes[!(skeletalmuscle_genes$`Gene name` %in% remove),]

# select 100 random genes that are all in the gene expression data
gene_list <- list(heart_genes$`Gene stable ID`, brain_genes$`Gene stable ID`, 
                  skeletalmuscle_genes$`Gene stable ID`)

l <- list()
i <- 1
for(g in gene_list){
  genes <- g[g %in% gsub("\\..*","",rownames(gct_file))]
  select <- sample(1:length(genes), 100, replace = FALSE)
  l[[i]] <- genes[select]
  i <- i + 1
}


#extract profiles for each tissue of interest seperately
heart <- GTEx_sampleinfo[GTEx_sampleinfo$histological_type=="Heart",]$Sample_Name
brain <- GTEx_sampleinfo[GTEx_sampleinfo$histological_type=="Brain",]$Sample_Name
muscle <- GTEx_sampleinfo[GTEx_sampleinfo$histological_type=="Muscle",]$Sample_Name


#for each tissue type
i <- 1
t <- c("heart", "brain", "muscle")
profiles <- list()

for (tissue in list(heart, brain, muscle)){
  #extract expression table
  tissue_file <- gct_file[,colnames(gct_file) %in% tissue]
  
  #get only the genes of interest
  tissue_file <- tissue_file[gsub("\\..*","",rownames(gct_file)) %in% l[[i]],]
  print(dim(tissue_file))
  profiles[[paste(t[i])]] <- tissue_file
  
  #graph distribution of biological variation in each tissue type
  cq <- calculate_quartile_dev(tissue_file)
  dat <- as.data.frame(list("cq"=cq))
  
  pdf(paste(t[i],"_quartile_dev.pdf",sep=""))
  g = ggplot(dat, aes(x=cq)) +
    geom_histogram(binwidth=.5, colour="black", fill="white") +
    theme(panel.background=element_blank()) +
    theme(axis.line = element_line(colour = "black", 
                      size = 1, linetype = "solid")) +
    geom_vline(aes(xintercept=mean(cq, na.rm=T)),   # Ignore NA values for mean
               color="red", linetype="dashed", size=1) +
    labs(title="Coefficient of Quartile Deviation",
         x ="Quartile Deviation", y = "Frequency")
  print(g)
  dev.off()

  cv <- coefficient_variation(tissue_file)
  dat <- as.data.frame(list("cv"=cv))
  
  pdf(paste(t[i],"_cov.pdf",sep=""))
  g = ggplot(dat, aes(x=cv)) +
    geom_histogram(binwidth=.5, colour="black", fill="white") +
    theme(panel.background=element_blank()) +
    theme(axis.line = element_line(colour = "black",
                      size = 1, linetype = "solid")) +
    geom_vline(aes(xintercept=mean(cv, na.rm=T)),   # Ignore NA values for mean
               color="red", linetype="dashed", size=1) +
    labs(title="Coefficient of Variation",
         x ="Variation", y = "Frequency")
  
  print(g)
  dev.off()
  i <- i + 1 
}


tissues=t
sizes=c(100)
sim_measure="spearman"
rep=TRUE

  similarity <- list()
  results <- c()
  i <- 1
  for (p in profiles){
    for(s in sizes){
      construct_networks(samples=p, tissue=tissues[i], N=10, s_size=s,
                         sim_measure=sim_measure, rep=rep)
      
    }
    i <- i+1
  }

