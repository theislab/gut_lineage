library(rhdf5)
library(Matrix)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(RColorBrewer)

#load batch corrected data (Combat batch correction)
f_path <- '~/Desktop/tmp_folder/data/gut_AB_AL_control_mutants_anno.h5ad' 
dset <- h5read(f_path, '/',compoundAsDataFrame=FALSE)
figure_path <- '~/Desktop/tmp_folder/figures/'
file_path <- '~/Desktop/tmp_folder/table/'
today <- format(Sys.Date(), '%y%m%d')
#get genes/cell ID
barcodes <- unlist(dset$obs$index)
genes <- unlist(dset$var$index)
#for the record: As of 16th October 2018, I have looked up, if the 'raw data' are Combat corrected or not.
#They are not. They are log-transformed and that's it. 
data_raw <- dset$raw.X$data
index_raw <- dset$raw.X$indices
ptr_raw <- dset$raw.X$indptr
sparse_mat <- sparseMatrix(p = as.numeric(ptr_raw), 
                           x=as.numeric(data_raw), 
                           i = as.numeric(index_raw)+1)
genes_raw <- dset$raw.var$index
n_genes <- dset$obs$n_genes
#set levels
sample_levels <- c("CD_1" ,"CD_2", "CD_3", "Control_1",  "Control_2" , "Control_6" ,          
                   "Control_3_FVR" ,  "Control_4_FVR", "Control_5_FVR","Control_7_FVR_only",
                   "Mutant_1" ,     "Mutant_2",   "Mutant_3_FVR",    "Mutant_4_FVR"
                   )
new_sample_levels <- c("FVF" ,"FVF", "FVF", "Control_1",  "Control_2" , "Control_6" ,          
                       "Control_3_FVR" ,  "Control_4_FVR", "Control_5_FVR","Control_7_FVR_only",
                       "Mutant_1" ,     "Mutant_2",   "Mutant_3_FVR",    "Mutant_4_FVR"
)

cell_type_levels <- c("ISC", "Enterocyte", "Enterocyte progenitor",  
                          "Goblet progenitor" , "early Goblet"  ,   "Goblet cell"  ,  
      
                          "Paneth primed ISC", "Paneth progenitor" ,  "Paneth 1" ,   "Paneth 2",     
                          "Lgr5+ EEC" ,  "Sox4+ early EE progenitor" , 
                          "Ngn3 progenitor", "Isl1/Arx progenitor",  
                          "Pax4 progenitor" ,"Ghrl progenitor",  
                          "EC",  "EC-Reg4",  "SAKD", "SIA", "SIK", "SIL-P", "SILA", "SIN", 
                          "Tuft progenitor" ,  "Tuft 1" , "Tuft 2")

new_cell_type_levels <- c("ISC", "Enterocyte", "Enterocyte progenitor",  
                          "Goblet progenitor" , "early Goblet"  ,   "Goblet cell"  ,  
                       
                          "Paneth primed ISC", "Paneth progenitor" ,  "Paneth 1" ,   "Paneth 2",     
                          "Lgr5+ EEC" ,  "Sox4+ early EE progenitor" , 
                          "Ngn3 progenitor", "Isl1-Arx progenitor",  
                          "Pax4 progenitor" ,"Ghrl progenitor",  
                          "EC",  "EC-Reg4",  "SAKD", "SIA", "SIK", "SIL-P", "SILA", "SIN", 
                          "Tuft progenitor" ,  "Tuft 1" , "Tuft 2")

proma_cell_type_levels <- c("ISC", "Enterocyte", "Enterocyte progenitor",  
                            "Goblet progenitor" , "early Goblet"  ,   "Goblet cell"  ,  
                         
                            "Paneth primed ISC", "Paneth progenitor" ,  "Paneth cell" ,   "Paneth cell",     
                            "EE progenitor" ,  "Sox4+ early EE progenitor" , 
                            "Ngn3 progenitor",  rep("EEC",11), 
                            "Tuft progenitor" ,"Tuft cell" , "Tuft cell")
#get attributes
cellData <- data.frame(
  sample=factor(unlist(dset$uns$sample_categories)[unlist(dset$obs$sample)+1], 
                levels=sample_levels),
  n_genes=n_genes,
  genetics = factor(unlist(dset$uns$genetics_categories)[unlist(dset$obs$genetics)+1]), 
  cell_type = factor(unlist(dset$uns$cell_type_test_categories)[unlist(dset$obs$cell_type_test+1)], 
                     cell_type_levels),
  progenitor_type = factor(unlist(dset$uns$cell_type_test_categories)[unlist(dset$obs$cell_type_test+1)], 
                           cell_type_levels)
) 

levels(cellData$cell_type) <- new_cell_type_levels
levels(cellData$progenitor_type) <- proma_cell_type_levels
###############################
#remove FVF samples from test #
###############################
FVF_idx <- !(cellData$sample %in% c('CD_1', 'CD_2', 'CD_3'))
cellData <- cellData[FVF_idx,]
#levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
sparse_mat <- sparse_mat[,FVF_idx]

#prepare pairwise tests mutant vs. control
for (type in levels(cellData$cell_type)){
  print(paste0('Testing ', type))
  testIdx <- cellData$cell_type == type
  sizes <- table(cellData$genetics[testIdx])
  if (any(sizes==0)){
    print('One test group is empty. Skip group.')
  }else{
    
    data_test <- as.matrix(sparse_mat[,cellData$cell_type %in% type])
    rownames(data_test) <- genes_raw
    batch_test <- cellData$sample[cellData$cell_type %in% type]
    batch_test <- droplevels(batch_test)
    genetics_test <- cellData$genetics[cellData$cell_type %in% type]
    n_genes_test <- cellData$n_genes[cellData$cell_type %in% type]
    #?limma
    #prepare model matrix
    dge <- edgeR::DGEList(data_test,sample=batch_test,group=genetics_test)
    design <- model.matrix(~genetics_test+n_genes_test)
    colnames(design) <- gsub("genetics_test|batch_test", "", colnames(design))
    #contrasts <- makeContrasts(mutant-control, levels = design)
    y <- new("EList")
    y$E <- dge
    rm(dge)
    
    fit <- limma::lmFit(y, design = design)
    rm(y)
    fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
    #fit2 <- contrasts.fit(fit, contrasts)
    rest <- limma::topTable(fit,coef='mutant',
                            number = Inf,adjust.method='BH')
    
    subset_rest <- rest[abs(rest$logFC)>0.05 & rest$adj.P.Val<0.05,]
    
    write.csv(x = subset_rest, 
              file=paste0(file_path, 'limma_',today, '_', type,'_mutant_vs_control.csv'))
    write.csv(x = rest, 
              file=paste0(file_path,'limma_',today, '_', type,'_mutant_vs_control_all.csv'))
    
  }
}

#################################
# test progenitor groups (coarse grained)
#################################
#prepare pairwise tests mutant vs. control
for (type in levels(cellData$progenitor_type)){
  print(paste0('Testing ', type))
  testIdx <- cellData$progenitor_type == type
  sizes <- table(cellData$genetics[testIdx])
  if (any(sizes==0)){
    print('One test group is empty. Skip group.')
  }else{
    
    data_test <- as.matrix(sparse_mat[,cellData$progenitor_type %in% type])
    rownames(data_test) <- genes_raw
    batch_test <- cellData$sample[cellData$progenitor_type %in% type]
    batch_test <- droplevels(batch_test)
    genetics_test <- cellData$genetics[cellData$progenitor_type %in% type]
    n_genes_test <- cellData$n_genes[cellData$progenitor_type %in% type]
    #?limma
    #prepare model matrix
    dge <- edgeR::DGEList(data_test,sample=batch_test,group=genetics_test)
    design <- model.matrix(~genetics_test+n_genes_test)
    colnames(design) <- gsub("genetics_test|batch_test", "", colnames(design))
    #contrasts <- makeContrasts(mutant-control, levels = design)
    y <- new("EList")
    y$E <- dge
    rm(dge)
    
    fit <- limma::lmFit(y, design = design)
    rm(y)
    fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
    #fit2 <- contrasts.fit(fit, contrasts)
    rest <- limma::topTable(fit,coef='mutant',
                            number = Inf,adjust.method='BH')
    
    subset_rest <- rest[abs(rest$logFC)>0.05 & rest$adj.P.Val<0.05,]
    
    write.csv(x = subset_rest, 
              file=paste0(file_path,'limma_',today, '_', type,'_mutant_vs_control_progen.csv'))
    write.csv(x = rest, 
              file=paste0(file_path,'limma_',today, '_', type,'_mutant_vs_control_progen_all.csv'))
    
  }
}



#################################
# remove newer samples for test #
#################################

# #remove FVF samples from test
# FVFplus_idx <- !(cellData$sample %in% c('CD_1', 'CD_2', 'CD_3',  "Control_6" , 
#"Control_5_FVR","Control_7_FVR_only"))
# cellData <- cellData[FVFplus_idx,]
# cellData$sample <- droplevels(cellData$sample)
# #levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
# sparse_mat <- sparse_mat[,FVFplus_idx]
# 
# #prepare pairwise tests mutant vs. control
# for (type in levels(cellData$cell_type)){
#   print(paste0('Testing ', type))
#   testIdx <- cellData$cell_type == type
#   sizes <- table(cellData$genetics[testIdx])
#   if (any(sizes==0)){
#     print('One test group is empty. Skip group.')
#   }else{
#     
#     data_test <- as.matrix(sparse_mat[,cellData$cell_type %in% type])
#     rownames(data_test) <- genes_raw
#     keep_genes <- rowMeans(data_test)>0.05
#     print(paste0('Keep ', sum(keep_genes), ' genes.'))
#     data_test <- data_test[keep_genes,]  
#     batch_test <- cellData$sample[cellData$cell_type %in% type]
#     
#     genetics_test <- cellData$genetics[cellData$cell_type %in% type]
#     
#     #?limma
#     #prepare model matrix
#     dge <- edgeR::DGEList(data_test,sample=batch_test,group=genetics_test)
#     design <- model.matrix(~0+genetics_test+batch_test)
#     colnames(design) <- gsub('genetics_test','', colnames(design))
#     colnames(design) <- gsub('batch_test','', colnames(design))
#     y <- new("EList")
#     y$E <- dge
#     rm(dge)
#     
#     contr.matrix <- makeContrasts(
#       MutvsCont = mutant - control, 
#       levels = colnames(design))
#     
#     fit <- limma::lmFit(y, design = design)
#     vfit <- contrasts.fit(fit, contrasts=contr.matrix)
#     rm(y)
#     # fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
#     efit <- eBayes(vfit)
#     #plotSA(efit, main="Final model: Mean−variance trend")
#     dt <- decideTests(efit)
#     print(summary(dt))
#     
#     rest <- limma::topTable(efit,coef=1,
#                             number = Inf,adjust.method='BH')
#     subset_rest <- rest[abs(rest$logFC)>0.05 & rest$adj.P.Val<0.05,]
#     write.csv(x = subset_rest, 
#               file=paste0(file_path, 'limma_',type,'_mutant_vs_control_paired.csv'))
#     
#   }
# }
# 
# for (type in levels(cellData$cell_type)){
#   print(paste0('Testing ', type))
#   testIdx <- cellData$cell_type == type
#   sizes <- table(cellData$genetics[testIdx])
#   if (any(sizes==0)){
#     print('One test group is empty. Skip group.')
#   }else{
#     
#     data_test <- as.matrix(sparse_mat[,cellData$cell_type %in% type])
#     rownames(data_test) <- genes_raw
#     batch_test <- cellData$sample[cellData$cell_type %in% type]
#     
#     genetics_test <- cellData$genetics[cellData$cell_type %in% type]
#     
#     # create scatter plot for gene means
#     mean_control <- rowMeans(data_test[,genetics_test=='control'])
#     mean_mutant <- rowMeans(data_test[,genetics_test=='mutant'])
#     mean_df <- data.frame(mean_control= mean_control, 
#                           mean_mutant = mean_mutant)
#     g <- ggplot(mean_df, aes(mean_control, mean_mutant)) + labs(title=type) +
#       geom_point(alpha=0.5) + geom_abline(slope = 1, intercept = 0) + theme_classic()
#     ggsave(filename = paste0(figure_path,'mean_', type, 'mutant_control.pdf'), 
#            plot = g, width = 6, height=6)
#     #write.csv(x = subset_rest, 
#     #          file=paste0(file_path, 'limma_',type,'_mutant_vs_control_paired.csv'))
#     
#   }
# }


########################
# create volcano plots #
########################
library(ggplot2)
library(RColorBrewer)
library(dplyr)

files <- list.files(path = file_path, pattern = 'paired.csv')
files <- files[grepl('limma_', files)]
for (filing in files){
gsub('limma_', '',gsub('.csv', '', filing)) -> title
limma_res <- read.csv(paste0(file_path, filing), row.names = 1)
data_res <- limma_res %>% select(logFC,AveExpr, adj.P.Val)
data_res$adj.P.Val <- -log10(data_res$adj.P.Val)

g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
  labs(x='log Fold Change (mutant/WT)', y='-log10(adjusted p)', title=title) + 
  theme_classic()
ggsave(filename = paste0(figure_path,'limma_', title, '.pdf'), plot = g0, width = 6, height=6)

g1 <- ggplot(data_res, aes(logFC, AveExpr)) + geom_point(alpha=0.7) + 
  labs(x='log Fold Change (mutant/WT)', y='Average Expression', title=title) + 
  theme_classic()
ggsave(filename = paste0(figure_path,'limma_lfc_mean_', title, '.pdf'), plot = g1, width = 6, height=6)
}

###################################
# get TFs from mutant vs. control (progenitors) #
###################################
library(biomartr)
files <- list.files(path = file_path, pattern = '_mutant_vs_control_progen.csv')
files <- files[grepl(today, files)]
files <- files[grepl('(progeni|ISC)', files)]


gene_list <- list()
threshold_score <- 1e-5

for (file in files){
  tmp <- read.csv(paste0(file_path, file))
  ix_score <- tmp$adj.P.Val< threshold_score & abs(tmp$logFC)>0.05
  gene_list <- rbind(gene_list, tmp[ix_score,])
}

gene_list <- gene_list[!duplicated(gene_list$X),]
ix <- 0.05

#mart = useMart('ensembl', dataset='mmusculus_gene_ensembl')
GO_tbl <- getGO(organism = "Mus musculus", 
                genes    = gene_list$X[abs(gene_list$logFC)>ix],
                filters  = 'mgi_symbol')

#GO-term for TF activity is GO:0003700;
TFs <- GO_tbl[GO_tbl$goslim_goa_accession %in% c('GO:0003700'),]
a_ix <- match(gene_list$X, TFs$mgi_symbol)
TFs_sort <- gene_list$X[!is.na(a_ix)]
test <- levels(TFs_sort)[as.numeric(TFs_sort)]
write.csv(test[order(test)], 
          file=paste0(file_path, 'TFs_progenitors_mutant_control_2021.csv'), row.names = FALSE)


