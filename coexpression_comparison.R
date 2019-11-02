library(multivariance)
library(readr)
library(ggplot2)
library(vegan)
library(stats)

################################################################################
# This script can be used to compare different co-expression networks.
# Experiments may be conducted for different sample subsets.
# Capable of constructing networks using Pearson Correlation, Spearman
# Correlation, and Mutual Information.

# Different similarity measures can also be applied to compare the networks.
# Currently capable of measuring absolute difference between co-expression
# networks as well as fuzzy set similarity from CoexSim.
# Also contains methods to analyze RNA-seq data to construct networks, but 
# it is also possible to use with microarray data.

# Author: Katie Ovens
# Date: July 2019
################################################################################

#'  Similarity measure between two graphs.
#'
#'  Calculate the difference in edge weights between two graphs. Each graph
#'  should contain the same nodes and edges.
#'  
#'  @param graph1 vector representing first network containing edge 
#'                weights.
#'  @param graph2 vector representing second network containing edge 
#'                weights.
#'  @param norm normalization factor depending on the range of scores
#'
#'  @return The normalized sum of the absolute differences between 
#'          corresponding edges in the graphs.
#'
edge_similarity <- function(graph1, graph2, norm=2){
  n_edges <- length(graph1)
  #check that the dataframes match for first 2 columns
  if ( n_edges != length(graph2) ){
    #stop, the graphs do not contain the same nodes to compare 
    stop("The graphs do not contain the same nodes to compare.")
  }
  #extract weight columns for each graph and get the absolute difference
  abs_diff <- abs(graph1 - as.numeric(graph2))
  #sum all together and normalize
  total <- sum(abs_diff)/(norm*n_edges)
  return(1-total)
}


################################################################################
#'  Assign rankings to items in a list
#'
#'  Calculate the rank of each item in a list. There is an option to
#'  bin the weights.
#'  
#'  @param w List of sorted weights that need to be ranked.
#'  @param epsilon Numeric value that cannot be exceeded between items in the 
#'                 list in order for them to be ranked the same.
#'                 Default is that each gene will get a different rank unless
#'                 the value is exactly equal.
#'  
#'  @return A vector of ranks where each corresponds to an item in the list.
#'
ranking <- function(w, epsilon=0){
  rank <- 1
  pos <- 1
  first <- 1
  ranks <- c()
  
  # go through weights
  while(pos <= length(w)){
    # when the difference with the first exceeds epsilon
    if(abs(w[first] - w[pos]) > epsilon){
      # set item to first and increase rank
      first <- pos
      rank <- rank + 1
      ranks <- c(ranks, rank)
    }
    else{
      #keep the rank the same
      ranks <- c(ranks, rank)
    }
    pos <- pos + 1
  }
  return(ranks)
}

################################################################################
#'  Perform Kendall test on a table of values
#'
#'  Calculate the rank of each item in columns of a table while allowing for 
#'  ties.
#'  
#'  @param values List of sorted weights that need to be ranked.
#'  @param eps Numeric value that cannot be exceeded between items in the 
#'                 list in order for them to be ranked the same.
#'                 Default is that each gene will get a different rank unless
#'                 the value is exactly equal.
#'  @param n.permutation number of permutations to perform Kendal test
#'  
#'  @return Result of kendall concordance coefficient test 
#'          (one for each column of values) across values of epsilon.
#'

kendall <- function(values, eps=c(0.00, 0.01, 0.05, 0.10), 
                    n.permutations=1000){
  rankings <- values
  kendall.result <- list()
  #for each epsilon
  for (e in eps){
    #for each edge, we have a correlation for each group of samples N
    #for each column
    for(c in 1:ncol(values)){
      #sort entire dataframe by column
      values <- values[order(values[,c]),]
      
      #pass the column to ranking
      r <- ranking(values[,c], e)
      
      #replace cor values by ranking
      rankings[,c] <- r
    }
    #perform kendall on the dataframe of ranks and store result for each value 
    #epsilon
    kendall.result[[paste(e)]] <- kendall.global(values, nperm=n.permutations,
                                                 mult = "BH")
                                     
  }
  return(kendall.result)
}

################################################################################

#'  Write correlation networks to file (to avoid memory errors).
#'
#'  Method to calculate correlation between genes from gene expression data.
#'  
#'  @param gct_file Table of samples/gene expression matrix 
#'                 (columns are samples, rows are genes).
#'  @param out_name The name of the file to write the correlations to.
#'  @param cor_method The correlation method to use. Default is "pearson".
#'
calculate_correlation <- function(gct_file, out_name, cor_method="pearson"){
  write("to,from,cor", file=out_name, append=TRUE)
  for(i in 1:nrow(gct_file)){
    for(j in i:nrow(gct_file)){
      # do not require edge between the same gene (uninformative)
      if (i == j){
        c <- 0
      }
      else{
        if (cor_method == "multi"){
            c <- multicorrelation(cbind(as.numeric(gct_file[i,]), as.numeric(gct_file[j,])),type="m.multi.2",multicorrelation.type = "normalized")
        }
        else{
            c <- cor(x=as.numeric(gct_file[i,]), y=as.numeric(gct_file[j,]),method=cor_method,  use="complete.obs")
        }
        if(is.na(c)){
            c <- 1
        }
        write(paste(rownames(gct_file)[i],rownames(gct_file)[j],c,sep=","), file=out_name, append=TRUE)
      }
    }
  }
}


################################################################################

#'  Generate N networks from random combinations of samples of size s_size.
#'
#'  Method to generate non-overlapping subsets of gene expression samples and
#'  apply a similarity measure to construct networks.
#'  
#'  @param samples Table of samples/gene expression matrix 
#'                 (columns are samples, rows are genes).
#'  @param N The number of networks to construct.
#'  @param s_size The number of samples that should be used to construct each
#'                of the N networks.
#'  @param sim_measure Measure to use to construct the networks.
#'  @param rep Boolean to allow or avoid repetition of samples.
#'
#'
construct_networks <- function(samples, tissue, N, s_size, 
                               sim_measure="pearson", rep=FALSE){
  # number of samples required in total
  num_samples <- N * s_size
  selected <- 0
  if(num_samples > ncol(samples)){
    # randomly select this number of samples without replacement
    selected <- sample(1:ncol(samples), num_samples, replace=rep)
  }
  else{
    stop("Not enough samples in the expression matrix.")
  }
  
  #get a list of sample sets
  tables <- split(selected, ceiling(seq_along(selected)/s_size))
  edge_list <- list()
  n <- 1
  for (net in names(tables)){
    if (sim_measure == "pearson"){
      #construct network using pearson correlation
      calculate_correlation(samples[,tables[[net]]], 
                            out_name = paste("network_",n,"_", s_size,"_",tissue,".csv", sep=""),
                            "pearson")
      
    }
    else if (sim_measure == "spearman"){
      #construct network using spearman correlation
      calculate_correlation(samples[,tables[[net]]], 
                            out_name = paste("network_",n,"_", s_size, "_", tissue, ".csv", sep=""),
                            "spearman")
    }
    else if (sim_measure == "multi"){
      #construct network using multivariance
      calculate_correlation(samples[,tables[[net]]],
                            out_name = paste("network_",n,"_", s_size, "_", tissue, ".csv", sep=""),
                            "multi")

    }
    
  n <- n + 1
  }
  
}

################################################################################

#'  Calculate the average similarity score for a list of networks.
#'  
#'  @param networks A list containing all networks to compare based on a 
#'                  similarity measure.
#'  @param n1 The name of the first node column (string).
#'  @param n2 The name of the second node column (string).
#'
#'  @return A average similarity score for all networks compared.
#'
compare_networks <- function(networks,n1, n2){
  sim_total <- 0
  unique_comparisons <- combn(x = names(networks), m = 2)
  for (c in 1:ncol(unique_comparisons)){
    sim_total = sim_total + edge_similarity(networks[[unique_comparisons[,c][1]]],
                                            networks[[unique_comparisons[,c][2]]],
                                            n1, n2)
  }
  
  #average scores
  sim_total = sim_total/ncol(unique_comparisons)
  return(sim_total)
}

################################################################################

#'  Generate N networks from random combinations of samples of size s_size.
#'
#'  Method to generate non-overlapping subsets of gene expression samples and
#'  apply a similarity measure to construct networks.
#'  
#'  @param profiles List of expression profiles to analyze.
#'  @param N The number of networks to construct for each sample size.
#'  @param sizes A vector containing different numbers of samples that should 
#'               be used to construct N networks.
#'  @param sim_measure Measure to use to construct the networks.
#'  @param rep Boolean to allow or avoid repetition of samples.
#'
#'  @return A list of similarity scores for each sample size in sizes.
#'
main_experiment <- function(profiles, tissues, N=5, sizes=seq(3,10,1),
                               sim_measure="pearson", rep=FALSE){
  similarity <- list()
  i <- 1
  results <- c()
  for (p in profiles){
    for(s in sizes){
      nets <- construct_networks(samples=p, tissue=tissues[i], N=N, s_size=s, 
                                 sim_measure=sim_measure, rep=rep)
      #calculate similarity and store in list
      similarity[[s]] <- compare_networks(nets, "to", "from")
    }
    
    #convert list to data frame and cbind to previous results
    df <- do.call(rbind.data.frame, similarity)
    df$replicate <- rep(i,nrow(df))
    if (length(results)==0){
      results <- df
    }
    else{
      results <- rbind(results, df)
    }
    i <- i+1
  }
  colnames(results) <- c("sample_size", "similarity","replicate")
  #make a graph using ggplot2 of the change in agreement between scores
  #Do the graphs stabilize over time or do they remain very different?
  g<-ggplot(results, aes(x=sample_size, y=similarity, group=replicate)) +
    geom_line(aes(color=replicate))+
    geom_point
  
  pdf(paste("results_",sim_measure, ".pdf", sep=""))
  print(g)
  dev.off()
  
}


##############################################################################
get_similarity <- function(tissue, N=10, sample_sizes = seq(3,50,1), 
                           method="pearson"){
  require("vegan") || stop("Package vegan is not available!")
  require("PMCMR") || stop("Package PMCMR is not available!")
  require("RVAideMemoire") || stop("Package RVAideMemoire is not available!")
  similarity <- list()
  values <- c()
  for(t in tissue){
    for(s in sample_sizes){
      sim_total <- c()
      for(n in 1:N){
        name <- paste(method,"/network_",n,"_", s,"_",t,".csv",sep="")
        values <- read.csv(name, header = TRUE, sep = ",",
                           colClasses=c("character","character","numeric"))
        
        sim_total = cbind(sim_total, values[,3])
      }
      #perform friedman test and kendall concordance test
      #kendall concordance coefficient test
      kendall.result <- kendall(sim_total, eps=c(0.00), 
                                n.permutations=1000)
      
      #Friedman rank sum test
      friedman.result <- friedman.test(as.matrix(sim_total))
      similarity[[t]][[paste(s)]] <- list(friedman.result, kendall.result)
    }
  }

  return(similarity)
}

##############################################################################
visualize_kendall <- function(results, tissue, e="0", title){
  #Plot Kendall W
  df <- c()
  for (t in tissue){
    sizes <- as.numeric(names(results[[t]]))
    res <- c()
    for(size in names(results[[t]])){
      res <- c(res, results[[t]][[paste(size)]][[2]][[e]]$Concordance_analysis["W",])
    }
    if(length(df)==0){
      df <- cbind(sizes, res, rep(t, length(sizes)))
    }
    else{
      df <- rbind(df, cbind(sizes, res, rep(t, length(sizes))))
    }
    
  }
  df <- as.data.frame(df)
  colnames(df) <- c("sample_size", "W", "tissue")
  df$sample_size <- as.numeric(as.character(df$sample_size))
  df$W <- as.numeric(as.character(df$W))
  
  plt <- ggplot(data=df, aes(x=sample_size, y=W, group=tissue)) +
    geom_line(aes(color=tissue)) +
    geom_point(aes(color=tissue)) +
    theme(plot.title=element_text(hjust=0.5, face = "bold", size=20))+
    theme(panel.background=element_blank(),
          axis.line = element_line(color="gray15"),
          axis.text.x = element_text(size=18, color="gray15"),
          axis.text.y = element_text(size=18, color="gray15"),
          axis.title = element_text(size=20, face = "plain", color="gray15"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="gray15"),
          legend.key=element_blank(),
          legend.title.align=0.5,
          legend.position = c(0.85, 0.22),
          legend.text=element_text(size=18),
          legend.title=element_text(size=20)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("K e n d a l l   W") +
    xlab("Sample Size") +
    scale_y_continuous(limits=c(0, 1))+
    xlim(0,50)
  
  pdf(paste(title, ".pdf", sep=""), width=8, height=5.5)
  print(plt)
  dev.off()
}
