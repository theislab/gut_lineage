library(rhdf5)
library(Matrix)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

f_path <- '~/Desktop/tmp_folder/data/gut_AB_AL_log_cor_control_anno.h5ad'
dset <- h5read(f_path, '/',compoundAsDataFrame=FALSE)
figure_path <- '~/Desktop/tmp_folder/figures/'
file_path <- '~/Desktop/tmp_folder/table/'
today <- format(Sys.Date(), '%y%m%d')
#get genes/cell ID
barcodes <- unlist(dset$obs$index)
genes <- unlist(dset$var$index)

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
                   "Control_3_FVR" ,  "Control_4_FVR", "Control_5_FVR","Control_7_FVR_only")
cell_type_levels <- c("ISC", "Enterocyte", "Enterocyte progenitor",  
                      "Goblet progenitor" , "early Goblet"  ,   "Goblet cell"  ,  
                      "Paneth primed ISC", "Paneth progenitor" ,  "Paneth 1" ,   "Paneth 2",     
                      "Lgr5+ EEC" ,  "Sox4+ early EE progenitor" , 
                      "Ngn3 progenitor", "Isl1/Arx progenitor",  
                      "Pax4 progenitor" ,"Ghrl progenitor",  
                      "EC",  "EC-Reg4",  "SAKD", "SIA", "SIK", "SIL-P", "SILA", "SIN", 
                      "Tuft progenitor" ,  "Tuft 1" , "Tuft 2")

progen_test_levels <- c("ISC", "Enterocyte", "Enterocyte progenitor",  
                        "Goblet progenitor" ,  "Goblet progenitor"  ,   "Goblet cell"  ,  
                        "Paneth primed ISC", "Paneth progenitor" ,  "Paneth 1" ,   "Paneth 2",     
                        "Lgr5+ EEC" ,  "Sox4+ early EE progenitor" , 
                        "Ngn3 progenitor", "Isl1/Arx progenitor",  
                        "Pax4 progenitor" ,"Ghrl progenitor",  
                        "EC",  "EC-Reg4",  "SAKD", "SIA", "SIK", "SIL-P", "SILA", "SIN", 
                        "Tuft progenitor" ,  "Tuft 1" , "Tuft 2")

#get attributes
cellData <- data.frame(
  sample=factor(unlist(dset$uns$sample_categories)[unlist(dset$obs$sample)+1], 
                levels=sample_levels),
  n_genes = n_genes,
  cell_type = factor(unlist(dset$uns$refined_clustering_categories)[unlist(dset$obs$refined_clustering+1)], 
                     cell_type_levels),
  progen_test_type = factor(unlist(dset$uns$refined_clustering_categories)[unlist(dset$obs$refined_clustering+1)], 
                     cell_type_levels)
) 

levels(cellData$progen_test_type) <- progen_test_levels
#create boxplot 
ggplot(cellData, aes(progen_test_type, n_genes)) +geom_boxplot() +theme_bw()
###############################
#remove FVF samples from test #
###############################
FVF_idx <- !(cellData$sample %in% c('CD_1', 'CD_2', 'CD_3'))
cellData <- cellData[FVF_idx,]
#levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
sparse_mat <- sparse_mat[,FVF_idx]

#prepare pairwise tests
celltypes_OI <- c('ISC','Paneth primed ISC', "Enterocyte progenitor", "Goblet progenitor",
                  "Paneth progenitor", "Sox4+ early EE progenitor", 
                  "Tuft progenitor","Ngn3 progenitor")
celltypes_OI_PC <- c('ISC','Paneth primed ISC','Paneth progenitor', 'Paneth 1', 'Paneth 2')
pairs1 <-  t(combn(celltypes_OI, 2))
pairs2 <-  t(combn(celltypes_OI_PC, 2))


for (i in 1:(dim(pairs1)[1])){
  cell_types <-pairs1[i,]
  print(cell_types)
  
  data_test <- as.matrix(sparse_mat[,cellData$progen_test_type %in% cell_types])
  rownames(data_test) <- genes_raw
  
  keep_genes <- rowMeans(data_test)>0.05
  print(paste0('Keep ', sum(keep_genes), ' genes.'))
  data_test <- data_test[keep_genes,]  
  
  batch_test <- cellData$sample[cellData$progen_test_type %in% cell_types]
  batch_test <- droplevels(batch_test)
  ct_test <- cellData$progen_test_type[cellData$progen_test_type %in% cell_types]
  ct_test <- droplevels(ct_test)
  #cell_types_test <- levels(ct_test)
  #levels(ct_test) <- c('ct1', 'ct2')
  #n_genes_test <- cellData$n_genes[cellData$progen_test_type %in% cell_types]
  #?limma
  #prepare model matrix
  
  dge <- edgeR::DGEList(data_test,sample=batch_test,group=ct_test)
  design <- model.matrix(~ct_test+batch_test)
  colnames(design) <- gsub('ct_test|batch_test','', colnames(design))
  #contrasts <- makeContrasts(ct1-ct2, levels = design)
  y <- new("EList")
  y$E <- dge
  rm(dge)
  
  fit <- limma::lmFit(y, design = design)
  rm(y)
  fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
  #fit2 <- contrasts.fit(fit, contrasts)
  rest <- limma::topTable(fit,coef= 2,
                          number = Inf,adjust.method='BH')
  
  subset_rest <- rest[abs(rest$logFC)>0.05 & rest$adj.P.Val<0.05,]
  title <- paste0(colnames(fit$coefficients)[2],'_vs_', 
                  cell_types[!cell_types==colnames(fit$coefficients)[2]], 
                  '_controls')
  write.csv(x = subset_rest, 
            file=paste0(file_path, 'limma_',today, '_', title,'.csv'))
  
 #create volcano plot
  data_res <- subset_rest %>% select(logFC, adj.P.Val)
  data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
  
  g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
    labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
    theme_classic()
  ggsave(filename = paste0(figure_path,'limma_', today, '_', title, '.pdf'), 
         plot = g0, width = 6, height=6)
}

#second test pair (Paneth lineage)
for (i in 1:(dim(pairs2)[1])){
  cell_types <-pairs2[i,]
  print(cell_types)
  
  
  data_test <- as.matrix(sparse_mat[,cellData$progen_test_type %in% cell_types])
  rownames(data_test) <- genes_raw
  
  keep_genes <- rowMeans(data_test)>0.05
  print(paste0('Keep ', sum(keep_genes), ' genes.'))
  data_test <- data_test[keep_genes,]  
  
  batch_test <- cellData$sample[cellData$progen_test_type %in% cell_types]
  batch_test <- droplevels(batch_test)
  ct_test <- cellData$progen_test_type[cellData$progen_test_type %in% cell_types]
  ct_test <- droplevels(ct_test)
  #cell_types_test <- levels(ct_test)
  #levels(ct_test) <- c('ct1', 'ct2')
  n_genes_test <- cellData$n_genes[cellData$progen_test_type %in% cell_types]
  #?limma
  #prepare model matrix
  
  dge <- edgeR::DGEList(data_test,sample=batch_test,group=ct_test)
  design <- model.matrix(~ct_test+batch_test)
  colnames(design) <- gsub('ct_test|batch_test','', colnames(design))
  #contrasts <- makeContrasts(ct1-ct2, levels = design)
  y <- new("EList")
  y$E <- dge
  rm(dge)
  
  fit <- limma::lmFit(y, design = design)
  rm(y)
  fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
  #fit2 <- contrasts.fit(fit, contrasts)
  rest <- limma::topTable(fit,coef= 2,
                          number = Inf,adjust.method='BH')
  
  subset_rest <- rest[abs(rest$logFC)>0.05 & rest$adj.P.Val<0.05,]
  title <- paste0(colnames(fit$coefficients)[2],'_vs_', 
                  cell_types[!cell_types==colnames(fit$coefficients)[2]], 
                  '_controls')
  write.csv(x = subset_rest, 
            file=paste0(file_path, 'limma_',today, '_', title,'_PC.csv'))
  
  #create volcano plot
  data_res <- subset_rest %>% select(logFC, adj.P.Val)
  data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
  
  g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
    labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
    theme_classic()
  ggsave(filename = paste0(figure_path,  'limma_',today,'_', title, '_PC.pdf'), plot = g0, width = 6, height=6)
}


#####################
# paired sample fit #
#####################
# FVFplus_idx <- !(cellData$sample %in% c('CD_1', 'CD_2', 'CD_3', 
#                                         'Control_5_FVR', 'Control_6', 'Control_7_FVR_only'))
# cellData <- cellData[FVFplus_idx,]
# #levels(cellData$sample) <- dropEmptyLevels(cellData$sample)
# sparse_mat <- sparse_mat[,FVFplus_idx]
# 
# #prepare pairwise tests
# celltypes_OI <- c('ISC', 'Paneth primed ISC',"Enterocyte progenitor", "Goblet progenitor",
#                   "Paneth progenitor", "Sox4+ early EE progenitor", 
#                   "Tuft progenitor","Ngn3 progenitor")
# pairs <-  t(combn(celltypes_OI, 2))
# #cell_types <- c('ISC', 'Sox4+ early EE progenitor')
# #cell_types <- c('ISC', 'Ngn3 pro')
# #cell_types <- c('ISC', 'Lgr5 pos EEC')
# #cell_types <- c('Sox4+ early EE progenitor', 'Lgr5 pos EEC')
# #cell_types <- c('ISC', 'Paneth primed ISC')
# #cell_types <- c('Paneth primed ISC', 'Paneth pro')
# #cell_types <- c('Paneth primed ISC', 'Paneth mature')
# #cell_types <- c('ISC', 'Paneth mature')
# #cell_types <- c('Lyz2 pos', 'Paneth mature')
# 
# for (i in 1:(dim(pairs)[1]-1)){
#   cell_types <-pairs[i,]
#   print(cell_types)
#   
#   
#   data_test <- as.matrix(sparse_mat[,cellData$cell_type %in% cell_types])
#   rownames(data_test) <- genes_raw
#   batch_test <- cellData$sample[cellData$cell_type %in% cell_types]
#   ct_test <- cellData$cell_type[cellData$cell_type %in% cell_types]
#   ct_test <- droplevels(ct_test)
#   
#   #?limma
#   #prepare model matrix
#   dge <- edgeR::DGEList(data_test,sample=batch_test,group=ct_test)
#   design <- model.matrix(~ct_test+batch_test)
#   y <- new("EList")
#   y$E <- dge
#   rm(dge)
#   
#   fit <- limma::lmFit(y, design = design)
#   rm(y)
#   fit <-  limma::eBayes(fit, trend = TRUE, robust = TRUE)
#   
#   rest <- limma::topTable(fit,coef=paste0('ct_test',cell_types[2]),
#                           number = Inf,adjust.method='BH')
#   
#   subset_rest <- rest[abs(rest$logFC)>0.05 & rest$adj.P.Val<0.05,]
#   title <- paste0(cell_types[2],'_vs_', cell_types[1], '_controls_paired')
#   write.csv(x = subset_rest, 
#             file=paste0(file_path,'limma_',today, '_', title,'.csv'))
#   
#   #create volcano plot
#   data_res <- subset_rest %>% select(logFC, adj.P.Val)
#   data_res$adj.P.Val <- -log10(data_res$adj.P.Val)
#   
#   g0 <- ggplot(data_res, aes(logFC, adj.P.Val)) + geom_point(alpha=0.7) + 
#     labs(x='log Fold Change', y='-log10(adjusted p)', title=title) + 
#     theme_classic()
#   ggsave(filename = paste0(figure_path,'limma_',today, '_', title, '.pdf'), plot = g0, width = 6, height=6)
# }

##################################
# rank plots for all comparisons #
##################################
# files <- list.files(path = file_path, pattern = 'controls.csv')
# files <- files[grepl('181017', files)]
# files <- files[!grepl('mutant', files)]

#####
# create ranking plots for progenitor vs ISC 
#####
# for (file in files){
#   titling <- paste(sapply(strsplit(file, '_', fixed=TRUE), '[', 3:5), collapse=' ') 
#   tmp <- read.csv(paste0(file_path, file))
#   shorten <- 1:min(30, dim(tmp)[1])
#   tmp <- tmp[order(abs(tmp$logFC), decreasing = TRUE),][shorten,]
#   tmp$color <- c( 'cornflowerblue','black')[1+as.numeric(tmp$logFC>0)]
#   
#   p0 <- ggplot(tmp, aes(x=shorten, y=abs(logFC), color=color)) + geom_blank() + 
#     labs(x='rank', y='log fold change') + 
#     ggtitle(titling) +  
#     expand_limits(y=c(0,5)) +
#     annotate('text', shorten, abs(tmp$logFC), label=tmp$X, color=tmp$color, angle = 90) + 
#     theme_bw() + theme(text = element_text(size=15))
#   ggsave(plot = p0, filename = paste0(figure_path,  'rank_plot_', today,'_', titling, '.pdf'), width=5, height=5)
# }
# 
# for (file in files){
#   titling <- paste(sapply(strsplit(file, '_', fixed=TRUE), '[', rev(3:5)), collapse=' ') 
#   tmp <- read.csv(paste0(file_path, file))
#   shorten <- 1:min(30, dim(tmp)[1])
#   tmp <- tmp[order(abs(tmp$logFC), decreasing = TRUE),][shorten,]
#   tmp$color <- c( 'cornflowerblue','black')[1+as.numeric(tmp$logFC<0)]
#   
#   p0 <- ggplot(tmp, aes(x=shorten, y=abs(logFC), color=color)) + geom_blank() + 
#     labs(x='rank', y='log fold change') + 
#     ggtitle(titling) +  
#     expand_limits(y=c(0,5)) +
#     annotate('text', shorten, abs(tmp$logFC), label=tmp$X, color=tmp$color, angle = 90) + 
#     theme_bw() + theme(text = element_text(size=15))
#   ggsave(plot = p0, filename = paste0(figure_path,  'rank_plot_swap_', today,'_', titling, '.pdf'), width=5, height=5)
# }
# 
# for (file in files){
#   titling <- paste(sapply(strsplit(file, '_', fixed=TRUE), '[', 3:5), collapse=' ') 
#   tmp <- read.csv(paste0(file_path, file))
#   shorten <- 1:min(30, dim(tmp)[1])
#   tmp <- tmp[order(tmp$logFC, decreasing = TRUE),][shorten,]
#   tmp$color <- c( 'cornflowerblue','black')[1+as.numeric(tmp$logFC>0)]
#   
#   p0 <- ggplot(tmp, aes(x=shorten, y=logFC, color=color)) + geom_blank() + 
#     labs(x='rank', y='log fold change') + 
#     ggtitle(titling) +  
#     expand_limits(y=c(-0.2,5)) +
#     annotate('text', shorten, abs(tmp$logFC), label=tmp$X, color=tmp$color, angle = 90) + 
#     theme_bw() + theme(text = element_text(size=15))
#   ggsave(plot = p0, filename = paste0(figure_path, 'rank_plot_only_pos_',today,'_', titling, '.pdf'), width=5, height=5)
# }
# 
# for (file in files){
#   titling <- paste(sapply(strsplit(file, '_', fixed=TRUE), '[', rev(3:5)), collapse=' ') 
#   tmp <- read.csv(paste0(file_path, file))
#   shorten <- 1:min(30, dim(tmp)[1])
#   tmp <- tmp[order(tmp$logFC, decreasing = FALSE),][shorten,]
#   tmp$color <- c( 'cornflowerblue','black')[1+as.numeric(tmp$logFC<0)]
#   
#   p0 <- ggplot(tmp, aes(x=shorten, y=abs(logFC), color=color)) + geom_blank() + 
#     labs(x='rank', y='log fold change') + 
#     ggtitle(titling) +  
#     expand_limits(y=c(-0.2,5)) +
#     annotate('text', shorten, abs(tmp$logFC), label=tmp$X, color=tmp$color, angle = 90) + 
#     theme_bw() + theme(text = element_text(size=15))
#   ggsave(plot = p0, filename = paste0(figure_path, 'rank_plot_only_pos_swap_',today,'_', titling, '.pdf'), width=5, height=5)
# }

###################################
# get TFs from ISC vs. progenitor #
###################################

#######
# get complete results to select for TFs
#######
library(biomartr)
files <- list.files(path = file_path, pattern = 'controls.csv')
files <- files[grepl('210209', files)]
filesPC_GC <- files[grep('Paneth progenitor_vs_Goblet pro', files)] 
filesISC <- files[grep('_ISC_(c|v)', files)]
files <- c( filesPC_GC, filesISC)


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
write.csv(test[order(test)], file=paste0(file_path, 'TFs_progenitors_2021.csv'), row.names = FALSE)

#get receptor expression 
#   GO:0007165 signal transduction
# comment: I checked for signalling receptors and related GO terms, but those genes are not differential
Recs <- GO_tbl[GO_tbl$goslim_goa_accession %in% c('GO:0007165'),]
a_ix <- match(gene_list$X, Recs$mgi_symbol)
Recs_sort <- gene_list$X[!is.na(a_ix)]
test <- levels(Recs_sort)[as.numeric(Recs_sort)]
write.csv(test, file=paste0(file_path, 'receptors_progenitors.csv'), row.names = FALSE)

