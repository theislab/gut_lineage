library(rhdf5)
library(Matrix)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

today <- format(Sys.Date(), '%y%m%d')
######
# limma results
######
#load differential expression results for mutant vs control
file_path <- '~/Desktop/tmp_folder/table/'
files <- list.files(path = file_path, pattern = 'mutant_vs_control_progen.csv')
files <- files[grepl('limma_', files)]
files <- files[grepl('210209', files)]
files <- files[grepl('progeni|ISC', files)] #select only progenitor and stem cell populations, not mature cells
limma_results <- data.frame(X='', logFC=numeric(1),AveExpr=0, adj.P.Val=0, cell_type='')
for (filing in files){
  cells_tested <- sapply(strsplit(filing, '_'),'[',3)
  limma_res <- read.csv(paste0(file_path, filing))
  data_res <- limma_res %>% select(X, logFC,AveExpr, adj.P.Val)
  data_res$cell_type <- cells_tested
  limma_results <- rbind(limma_results, data_res)
}  
limma_results <- limma_results[-1,]
write.csv(x = limma_results, file = paste0(file_path,  'limma_results_',today,'_mutant_vs_control_progenitors.csv'))

# load comparison between progenitors and mature cell types
files <- list.files(path = file_path, pattern = '_controls.csv')
files <- files[grepl('210209', files)]
limma_results2 <- data.frame(X='', logFC=numeric(1),AveExpr=0, 
                             adj.P.Val=0, cell_type1='', cell_type2='')
for (filing in files){
  cells_tested1 <- sapply(strsplit(filing, '_'),'[',3)
  cells_tested2 <- sapply(strsplit(filing, '_'),'[',5)
  limma_res <- read.csv(paste0(file_path, filing))
  data_res <- limma_res %>% select(X, logFC,AveExpr, adj.P.Val)
  data_res$cell_type1 <- cells_tested1
  data_res$cell_type2 <- cells_tested2
  limma_results2 <- rbind(limma_results2, data_res)
}  
limma_results2 <- limma_results2[-1,]
write.csv(x = limma_results2, file = paste0(file_path,  'limma_results_',today,'_controls_progen.csv'))

a <- unique(limma_results2$X[abs(limma_results2$logFC)>0.1])
write.csv(a[order(a)], file = paste0(file_path, 'limma_results_', today, '_differential_genes_lfc05.csv'), row.names = FALSE)

# load comparison between ISCs and maturing Paneth cells
files <- list.files(path = file_path, pattern = '_controls_PC.csv')
files <- files[grepl('210209', files)]
limma_results3 <- data.frame(X='', logFC=numeric(1),AveExpr=0, adj.P.Val=0, 
                             cell_type1='', cell_type2='')
for (filing in files){
  cells_tested1 <- sapply(strsplit(filing, '_'),'[',3)
  cells_tested2 <- sapply(strsplit(filing, '_'),'[',5)
  limma_res <- read.csv(paste0(file_path, filing))
  data_res <- limma_res %>% select(X, logFC,AveExpr, adj.P.Val)
  data_res$cell_type1 <- cells_tested1
  data_res$cell_type2 <- cells_tested2
  limma_results3<- rbind(limma_results3, data_res)
}  
limma_results3 <- limma_results3[-1,]
write.csv(x = limma_results3, file = paste0(file_path, 'limma_results_',today,'_controls_PC_maturation.csv'), row.names = FALSE)

#######
#check Wnt genes for significance
#######
canonical_Wnt_genes <- c('Dvl2', 'Axin2', 'Ascl2', 'Lgr5', 'Myc', 'Hopx')
some_more_wnt_genes <- c('Wnt4', 'Wnt5a','Wnt5b', 'Wnt9a','Wnt9b', 'Wnt6',
                         'Fzd2', 'Fzd3',   'Fzd6','Fzd7',
                         'Dvl1','Dvl2','Dvl3','Invs',
                         'Vangl1', 'Vangl2', 'Prickle1',
                         'Celsr1', 'Celsr2', 'Celsr3',
                         'Scrib', 'Fuz', 'Intu',
                         'Ror2', 'Ryk', 'Ptk7',
                         'Smurf1', 'Smurf2', 'Mapk8',
                         'Cdc42', 'Rhoa','Exoc3', 'Exoc4', 'Exoc5',  'Jun')
#check receptors for significance
receptors <- c(
  'Fgfr1','Fgfr2','Fgfr3','Fgfr4', #FGF
  'Egfr', 'Erbb2','Erbb3','Lrig1', #EGF 
  'Bmpr1a', 'Bmpr2','Id1','Id2','Id3', #BMP
  'Notch1','Notch2','Notch3','Jag1','Dll1','Dll4','Hes1','Lfng', #Notch 
  'Sfrp5','Sfrp1','Fzd1','Fzd2','Fzd3','Fzd6','Fzd7','Fzd8','Lrp5','Lrp6', #Wnt 
  'Smo', 'Shh', 'Ihh', #Hedgehog
  'Yap1','Tead2', 'Tead3', #Hippo
  'Fzd3', 'Fzd6', 'Ror2', 'Ptk7', 'Celsr1', 'Vangl1', 'Vangl2', 'Prickle1', 'Jun', #Wnt/PCP
  'Ephb2', 'Ephb3', 'Efnb1', #Ephrin signalling
  'Itgb1', 'Itga1','Itga2', 'Itga3', 'Itga5', 'Itga6', 'Itga9' #integrin signalling
  
)

####
# check some more genes (including TFs) for significance
####
more_interesting_genes =c('Cdkn1a', "Lbh", 'Klf15', "Ier2", 'Klf3', 'Atoh1', "Spdef", "Btg2", "Insm1", "Neurog3", "Sox4" )
#differential TFs mutant vs control
tf_mutant_control <- read.csv(paste0(file_path,'TFs_progenitors_mutant_control_Oct.csv'))
#differential TFs ISC vs progenitor or Paneth progenitor vs Goblet progenitor
tf_progen <- read.csv(paste0(file_path,'TFs_progenitors_2021.csv'))
tfs_differential = intersect(unlist(tf_mutant_control), unlist(tf_progen))
tfs_diffset = setdiff(unlist(tf_mutant_control), unlist(tf_progen))
write.csv(tfs_differential, 
          file=paste0(file_path, 'TFs_diff_progenitors_mutant_control_Oct.csv'), 
          row.names = FALSE)
progenitors <- c('ISC','Paneth primed ISC', 'Paneth progenitor',  "Sox4+ early EE progenitor" , 
                 "Ngn3 progenitor", 'EEC', 'Goblet progenitor', 'Tuft progenitor')

write.csv(x = limma_results[limma_results$X %in% canonical_Wnt_genes & limma_results$cell_type %in% progenitors,],
          row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_mutant_control_Wnt.csv'))
write.csv(x = limma_results[limma_results$X %in% some_more_wnt_genes & limma_results$cell_type %in% progenitors,],
          row.names =FALSE,
          file= paste0(file_path,'limma_results_',today,'_mutant_control_Wnt-PCP.csv'))
write.csv(x = limma_results[limma_results$X %in% more_interesting_genes & limma_results$cell_type %in% progenitors,],
          row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_mutant_control_genes_of_interest.csv'))
write.csv(x = limma_results[limma_results$X %in% tfs_differential & limma_results$cell_type %in% progenitors,],
          row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_mutant_control_TFs_intersect.csv'))
write.csv(x = limma_results[limma_results$X %in% tfs_diffset & limma_results$cell_type %in% progenitors,],
          row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_mutant_control_TFs_diffset.csv'))
write.csv(x = limma_results[limma_results$X %in% receptors& limma_results$cell_type %in% progenitors,],
          row.names =FALSE,
          file= paste0(file_path,'limma_results_',today,'_mutant_control_receptors.csv'))
#progenitor tests
write.csv(x = limma_results2[limma_results2$X %in% canonical_Wnt_genes ,],row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_controls_Wnt.csv'))
write.csv(x = limma_results2[limma_results2$X %in% some_more_wnt_genes,],row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_controls_Wnt-PCP.csv'))
write.csv(x = limma_results2[limma_results2$X %in% more_interesting_genes ,],row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_controls_genes_of_interest.csv'))
write.csv(x = limma_results2[limma_results2$X %in% tfs_differential ,],row.names =FALSE,
          file= paste0(file_path, 'limma_results_',today,'_controls_TFs_intersect.csv'))
write.csv(x = limma_results2[limma_results2$X %in% receptors,],
          row.names =FALSE,
          file= paste0(file_path,'limma_results_',today,'_controls_receptors.csv'))
#########
#load data for summary plots
########
f_path <- '~/Desktop/tmp_folder/data/gut_AB_AL_log_cor_control_anno.h5ad'
dset <- h5read(f_path, '/',compoundAsDataFrame=FALSE)
figure_path <- '~/Desktop/tmp_folder/figures/'

f_path_mut <- '~/Desktop/tmp_folder/data/gut_AB_AL_control_mutants_anno.h5ad' 
dset_mut <- h5read(f_path_mut, '/', compoundAsDataFrame = FALSE)

#get genes/cell ID
barcodes <- unlist(dset$obs$index)
genes <- unlist(dset$var$index)



sample_levels <- c( "Control_1",  "Control_2" , "Control_6" ,          
                   "Control_3_FVR" ,  "Control_4_FVR", "Control_5_FVR","Control_7_FVR_only","CD_1" ,"CD_2", "CD_3")
new_sample_levels <- c("Control_1",  "Control_2" , "Control_6" ,          
                       "Control_3_FVR" ,  "Control_4_FVR", "Control_5_FVR","Control_7_FVR_only",'FVF', 'FVF', 'FVF')
cell_type_levels <- c("ISC",  "Enterocyte progenitor",  "Enterocyte",
                      "Goblet progenitor" , "early Goblet"  ,   "Goblet cell"  ,  
      
                      "Paneth primed ISC", "Paneth progenitor" ,  "Paneth 1" ,   "Paneth 2",     
                      "Lgr5+ EEC" ,   "Sox4+ early EE progenitor" , 
                      "Ngn3 progenitor", "Isl1/Arx progenitor",  
                      "Pax4 progenitor" ,"Ghrl progenitor",  
                      "EC",  "EC-Reg4",  "SAKD", "SIA", "SIK", "SIL-P", "SILA", "SIN", 
                      "Tuft progenitor"  , "Tuft 1" , "Tuft 2")
split_Gpc_levels <- c("ISC", "Enterocyte progenitor",   "Enterocyte",
                      "Goblet progenitor" , "early Goblet"  ,   "Goblet cell"  ,  
                      "ISC", "Paneth progenitor" ,  "Paneth cell" ,   "Paneth cell",     
                      "EE progenitor" ,   "EE progenitor" , 
                      "EE progenitor", "EEC",  
                      "EEC" ,"EEC",  
                      "EEC",  "EEC",  "EEC", "EEC", "EEC",  "EEC",  "EEC", "EEC", 
                      "Tuft progenitor"  , "Tuft cell" , "Tuft cell")

proma_ct_mutant_levels <- c("ISC","Enterocyte progenitor",   "Enterocyte", 
                            "Goblet progenitor" ,  "early Goblet"   ,   "Goblet cell"  ,  
                            "ISC", "Paneth progenitor" ,  "Paneth cell" ,   "Paneth cell",     
                            "EE progenitor" ,   "EE progenitor" , 
                            "EE progenitor", "EEC",  
                            "EEC" ,"EEC",  
                            "EEC",  "EEC",  "EEC", "EEC", "EEC", "EEC", "EEC", "EEC", 
                            "Tuft progenitor"  , "Tuft cell" , "Tuft cell")

major_cell_type_levels = c( "ISC"  , "Enterocyte", "Goblet cell" ,  "Paneth cell" ,"EEC" ,       
                            "Tuft cell")
proma_cell_type_levels = c( "ISC"  ,  "Enterocyte progenitor" ,"Enterocyte", 
                            "Goblet progenitor", "early Goblet" , "Goblet cell" , 
                            "Paneth progenitor" , "Paneth cell" ,
                            "EE progenitor"   ,"EEC" ,        
                            "Tuft progenitor", "Tuft cell")

genetics <- c('control', 'mutant')

colors_proma <-c('#A30059', #ISC
                 '#FF4A46', #Enterocyte progenitor
                 '#5A0007', #Enterocyte
                 '#006FA6', #Goblet progenitor
                 '#0000A6', #Goblet
                 '#fdb462', #Paneth progenitor
                 '#ffd92f', #Paneth
                 '#63FFAC', #EEprogenitor
                 '#0A8201', #EEC
                 '#9083CF', #Tuft progenitor
                 '#304130' #Tuft
                 )

colors_split_Gpc <-c('#A30059', #ISC
                 '#FF4A46', #Enterocyte progenitor
                 '#5A0007', #Enterocyte
                 '#006FA6', #Goblet progenitor
                 '#CB63CC', #early Goblet
                 '#0000A6', #Goblet
                 '#fdb462', #Paneth progenitor
                 '#ffd92f', #Paneth
                 '#63FFAC', #EEprogenitor
                 '#0A8201', #EEC
                 '#9083CF', #Tuft progenitor
                 '#304130' #Tuft
)

#get cell attributes for Mutants and Controls
cellData_mut <- data.frame(mouse_line = factor(unlist(dset_mut$uns$genetics_categories)[unlist(dset_mut$obs$genetics)+1]),
                           sample =  factor(unlist(dset_mut$uns$sample_categories)[unlist(dset_mut$obs$sample)+1]),
                           cc_score = factor(unlist(dset_mut$uns$phase_categories)[unlist(dset_mut$obs$phase)+1], 
                                         levels=c('G1', 'G2M', 'S')),
                      # cell_type = factor(unlist(dset_mut$uns$cell_type_test_categories)[unlist(dset_mut$obs$cell_type_test)+1],
                      #                    levels = cell_type_levels),
                       proma_cell_type = factor(unlist(dset_mut$uns$proma_cell_type_categories)[unlist(dset_mut$obs$proma_cell_type)+1],
                                                levels = proma_cell_type_levels)
                       )
#levels(cellData_mut$proma_cell_type) <- proma_ct_mutant_levels

#count cells per cell type
cellData_mut_new <- cellData_mut
cellData_mut.m <- melt(table(cellData_mut_new), 
                   id.vars=c('mouse_line', 'sample', 'cc_score', 'proma_cell_type')) 
cellData_mut.m <- cellData_mut.m[cellData_mut.m$value>0,]

#select only non enriched samples
cellData_mut.m <- cellData_mut.m[cellData_mut.m$sample %in% c('Control_1', 'Control_2', 'Mutant_1', 'Mutant_2'),]

#get cell attributes for Controls only
cellData <- data.frame(
  sample=factor(unlist(dset$uns$sample_categories)[unlist(dset$obs$sample)+1], 
                levels=sample_levels),
  genetics = unlist(dset$uns$genetics_categories)[unlist(dset$obs$genetics)+1],
  cell_type = factor(unlist(dset$uns$refined_clustering_categories)[unlist(dset$obs$refined_clustering+1)], 
                     cell_type_levels),
  cc_score = factor(unlist(dset$uns$phase_categories)[unlist(dset$obs$phase)+1], 
                    levels=c('G1', 'G2M', 'S')),
  major_cell_type = factor(unlist(dset$uns$major_cell_type_categories)[unlist(dset$obs$major_cell_type+1)], 
                           major_cell_type_levels),
  proma_cell_type = factor(unlist(dset$uns$proma_cell_type_categories)[unlist(dset$obs$proma_cell_type+1)], 
                           proma_cell_type_levels),
  split_Gpc_type = factor(unlist(dset$uns$refined_clustering_categories)[unlist(dset$obs$refined_clustering+1)], 
                          cell_type_levels)
                ) 
#get sample sizes
levels(cellData$split_Gpc_type) <- split_Gpc_levels
sample_size <- table(cellData$sample)


cellData_new <- cellData
levels(cellData_new$sample) <- new_sample_levels
cellData.m <- melt(table(cellData_new), 
                   id.vars=c('cell_type', 'major_cell_type', 'proma_cell_type', 'split_Gpc_type')) 
cellData.m <- cellData.m[cellData.m$value>0,]

CD <- cellData.m %>%
  select(sample, cell_type,major_cell_type, proma_cell_type,value) %>%
  group_by(cell_type) %>% 
  summarise(total_count= sum(value))
CD_major <- cellData.m %>%
  select(sample, cell_type,major_cell_type, proma_cell_type,value) %>%
  group_by(major_cell_type) %>% 
  summarise(total_count= sum(value))
CD_proma<- cellData.m %>%
  select(sample, cell_type,major_cell_type, proma_cell_type,value) %>%
  group_by(proma_cell_type) %>% 
  summarise(total_count= sum(value))

CD_split_Gpc<- cellData.m %>%
  select(sample, cell_type,split_Gpc_type,value) %>%
  group_by(split_Gpc_type) %>% 
  summarise(total_count= sum(value))

CD_sample <- cellData.m %>%
  select(sample, cell_type,major_cell_type, proma_cell_type,value) %>%
  group_by(sample, cell_type) %>% 
  summarise(total_count= sum(value))

CD_sample_tot <- cellData.m %>%
  select(sample, cell_type,major_cell_type, proma_cell_type,value) %>%
  group_by(sample) %>% 
  summarise(total_count= sum(value))

CD_sample_tot2 <- data.frame(sample=names(sample_size), value=sample_size)

CD_sample_M <- cellData.m %>%
  select(sample, cell_type,major_cell_type, proma_cell_type,value) %>%
  group_by(sample, major_cell_type) %>% 
  summarise(total_count= sum(value))
CD_sample_pro <- cellData.m %>%
  select(sample, cell_type,major_cell_type, proma_cell_type, value) %>%
  group_by(sample, proma_cell_type) %>% 
  summarise(total_count= sum(value))
CD_sample_pro_Gpc <- cellData.m %>%
  select(sample, cell_type,major_cell_type, split_Gpc_type, value) %>%
  group_by(sample, split_Gpc_type) %>% 
  summarise(total_count= sum(value))


################
# create plots #
################

g06 <- ggplot(cellData_mut.m[cellData_mut.m$proma_cell_type %in% c('ISC', 'Paneth progenitor', 
                                                                   'Paneth cell','EE progenitor',
                                                                   'EEC'),], 
             aes(mouse_line, value, fill=cc_score))+ 
  scale_fill_manual(values =c('#1f77b4', '#ff7f0e', '#2ca02c'),
                    name='Cell cycle\nphase') +
  geom_bar(position="fill", stat='identity')  +  
  labs(x='Mouse line', y='Proliferation (%)') +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0, 25,50 ,75, 100)
  ) + facet_grid(~proma_cell_type) + 
  #scale_x_discrete(breaks=c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only','FVF'),
  #                 labels=c('Control 6', 'Control 5 FVR', 'Control 7 FVR only','FVF')) +
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g06
ggsave(filename = paste0(figure_path,today, '_cellcycle_per_cell_type_mutants.pdf'), plot = g06, width = 10, height=7) 

#Ext Data Fig. 8c
g06 <- ggplot(cellData_mut.m[cellData_mut.m$proma_cell_type %in% c("ISC"  ,  
                                                                   "Enterocyte progenitor" ,"Enterocyte",  
                                                                   "Goblet progenitor",'early Goblet', "Goblet cell" , 
                                                                   "Paneth progenitor" , "Paneth cell" ,
                                                                   "EE progenitor"   ,"EEC" ,        
                                                                   "Tuft progenitor", "Tuft cell"),], 
              aes(mouse_line, value, fill=cc_score))+ 
  scale_fill_manual(values =c('#1f77b4', '#ff7f0e', '#2ca02c'),
                    name='Cell cycle\nphase') +
  geom_bar(position="fill", stat='identity')  +  
  labs(x='Mouse line', y='Proliferation (%)') +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0, 25,50 ,75, 100)
  ) + facet_grid(~proma_cell_type) + 
  #scale_x_discrete(breaks=c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only','FVF'),
  #                 labels=c('Control 6', 'Control 5 FVR', 'Control 7 FVR only','FVF')) +
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g06
ggsave(filename = paste0(figure_path,today, '_cellcycle_per_cell_type_mutants.pdf'), plot = g06, width = 15, height=7) 

#Ext Data Fig. 4 
g07 <- ggplot(cellData.m[cellData.m$proma_cell_type %in% c("ISC"  ,  
                                                           "Enterocyte progenitor" ,"Enterocyte",  
                                                           "Goblet progenitor",'early Goblet', "Goblet cell" , 
                                                           "Paneth progenitor" , "Paneth cell" ,
                                                           "EE progenitor"   ,"EEC" ,        
                                                           "Tuft progenitor", "Tuft cell"),], 
              aes(proma_cell_type, value, fill=cc_score))+ 
  scale_fill_manual(values =c('#1f77b4', '#ff7f0e', '#2ca02c'),
                    name='Cell cycle\nphase') +
  geom_bar(position="fill", stat='identity')  +  
  labs(x='Cell type', y='Proliferation (%)') +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0, 25,50 ,75, 100)
  ) +
  #scale_x_discrete(breaks=c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only','FVF'),
  #                 labels=c('Control 6', 'Control 5 FVR', 'Control 7 FVR only','FVF')) +
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g07
ggsave(filename = paste0(figure_path, today, '_cellcycle_per_cell_type_progen.pdf'), plot = g07, width = 7, height=7) 


g01 <- ggplot(CD_sample_tot, aes(sample, total_count)) +# facet_wrap(~genetics)+
  geom_col(position="dodge") + #scale_fill_brewer(type='qual', palette = 'Paired', guide=FALSE) + 
  lims(y=c(0,20000))+
  #scale_y_sqrt(breaks=c(50, 500, 1000, 5000, 10000, 15000)) + 
  scale_x_discrete(breaks= CD_sample_tot$sample, 
                   labels=c('Control 1', 'Control 2', 'Control 6', 'Control 3 (FVR)',
                            'Control 4 (FVR)', 'Control 5 (50% FVR)', 'Control 7 (90% FVR)',
                            'FVF (3 samples)')
                   ) +
  theme_classic() + labs(title='Cell number per sample', x='Sample', y='Size') +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g01
ggsave(filename = paste0(figure_path, 'counts_per_sample_all_controls.pdf'), plot = g01, width = 5, height=7)  

g02 <- ggplot(CD_sample_tot2, aes(sample, value.Freq)) +# facet_wrap(~genetics)+
  geom_col(position="dodge") + #scale_fill_brewer(type='qual', palette = 'Paired', guide=FALSE) + 
  lims(y=c(0,20000))+
  #scale_y_sqrt(breaks=c(50, 500, 1000, 5000, 10000, 15000)) + 
 scale_x_discrete(breaks= levels(CD_sample_tot2$sample), 
                   labels=c('FVF 1', 'FVF 2', 'FVF 3',
                     'Control 1', 'Control 2', 'Control 6', 'Control 3 (FVR)',
                            'Control 4 (FVR)', 'Control 5 (50% FVR)', 'Control 7 (90% FVR)'
                            )
  ) +
  theme_classic() + labs(title='Cell number per sample', x='Sample', y='Size') +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g02
ggsave(filename = paste0(figure_path, 'counts_per_sample_all_controls2.pdf'), plot = g02, width = 5, height=7)  



g0 <- ggplot(CD, aes(cell_type, total_count)) +# facet_wrap(~genetics)+
  geom_col(position="dodge") + #scale_fill_brewer(type='qual', palette = 'Paired', guide=FALSE) + 
  scale_y_sqrt(breaks=c(50, 500, 1000, 5000, 10000, 15000)) + 
  scale_x_discrete(breaks= cell_type_levels) +
  theme_classic() + labs(title='Cell number per cell type', x='Cell type', y='Population size') +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g0
ggsave(filename = paste0(figure_path, 'counts_per_refined_cell_type_all_controls.pdf'), plot = g0, width = 10, height=7)  

g1 <- ggplot(CD_major, aes(major_cell_type, total_count)) +# facet_wrap(~genetics)+
  geom_col(position="dodge") + #scale_fill_brewer(type='qual', palette = 'Paired', guide=FALSE) + 
  scale_y_sqrt(breaks=c(1000,2000,3000,4000, 5000, 10000, 20000)) + 
  theme_classic() + labs(title='Cell number per cell type', x='Cell type', y='Population size') +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g1
ggsave(filename = paste0(figure_path, today, '_counts_per_major_cell_type_all_controls.pdf'), plot = g1, width = 4, height=7) 

g2 <- ggplot(CD_proma, aes(proma_cell_type, total_count)) +# facet_wrap(~genetics)+
  geom_col(position="dodge") + #scale_fill_brewer(type='qual', palette = 'Paired', guide=FALSE) + 
  scale_y_sqrt(breaks=c(500, 1000,2000,3000,4000, 5000, 10000, 20000)) +
  theme_classic() + labs(title='Cell number per cell type', x='Cell type', y='Population size') +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g2
ggsave(filename = paste0(figure_path, today, '_counts_per_main_progen_cell_type_all_controls.pdf'), plot = g2, width = 6, height=7) 

#Fig 3e 
g21 <- ggplot(CD_split_Gpc, aes(split_Gpc_type, total_count, fill=split_Gpc_type)) +# facet_wrap(~genetics)+
  geom_col(position="dodge") + scale_fill_manual(values=colors_split_Gpc, guide=FALSE)+#scale_fill_brewer(type='qual', palette = 'Paired', guide=FALSE) + 
  scale_y_sqrt(breaks=c(500, 1000,2000,3000,4000, 5000, 10000, 15000)) + #guide_legend() + 
  theme_classic() + labs(title='Cell number per cell type', x='Cell type', y='Population size') +
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g21
ggsave(filename = paste0(figure_path, today, '_counts_per_cell_type_all_controls.pdf'), plot = g21, width = 6, height=7) 


g2 <- ggplot(CD_sample[CD_sample$sample %in% c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only'),], 
       aes(cell_type, total_count, fill=sample)) +
  scale_fill_manual(breaks=c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only'), values= c('grey', 'firebrick1', 'cornflowerblue') ) +
  # facet_wrap(~genetics)+
  geom_bar(position="fill", stat='identity') +   
  labs(x='Cell type', y='Population size (%)') +
 scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                    labels=c(0, 25,50 ,75, 100)
                    ) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g2
ggsave(filename = paste0(figure_path, 'counts_per_refined_cell_type_c5-7.pdf'), plot = g2, width = 10, height=7)  

g2 <- ggplot(CD_sample_pro[CD_sample_pro$sample %in% c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only'),], 
             aes(proma_cell_type, total_count, fill=sample)) +
  scale_fill_manual(breaks=c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only'), 
                    values= c('grey', 'firebrick1', 'cornflowerblue') ) +
  # facet_wrap(~genetics)+
  geom_bar(position="fill", stat='identity') +   
  labs(x='Cell type', y='Population size (%)') +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0, 25,50 ,75, 100)
  ) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g2
ggsave(filename = paste0(figure_path, 'counts_per_main_progen_cell_type_c5-7.pdf'), plot = g2, width = 10, height=7) 



#Fig 3c
g6 <- ggplot(CD_sample_pro_Gpc[CD_sample_pro_Gpc$sample %in% c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only','FVF'),], 
             aes(sample, total_count, fill=split_Gpc_type)) +
  scale_fill_manual(breaks=levels(cellData$split_Gpc_type), 
                    values=colors_split_Gpc, 
                    name='Cell type') +
  # facet_wrap(~genetics)+
  geom_bar(position="fill", stat='identity') +   
  labs(x='Sample', y='Cell type frequency (%)') +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0, 25,50 ,75, 100)
  ) + 
  scale_x_discrete(breaks=c('Control_6', 'Control_5_FVR', 'Control_7_FVR_only','FVF'),
                   labels=c('no enrichment', 'FVR enrichment (50%)', 'FVR enrichment (90%)','FVF enrichment (50%)')) +
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.y  = element_text(size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=15))
g6
ggsave(filename = paste0(figure_path, today, '_counts_per_cell_type_split_Gpc_FVF_c5-7_stacked.pdf'), plot = g6, width = 6, height=10) 


#############################
# load Wnt/PCP score        #
#############################

wnt_data <- read.csv('~/Documents/Cooperations/Böttcher_IDR/10X_data/analysis/notebooks/table/Wnt_scale_per_cell_type.csv')
wnt_data$refined_clustering <- factor(wnt_data$refined_clustering, levels = cell_type_levels)
wnt_data2 <- read.csv('~/Documents/Cooperations/Böttcher_IDR/10X_data/analysis/notebooks/table/Wnt_scale_per_promain_cell_type.csv')
wnt_data2$proma_cell_type <- factor(wnt_data2$proma_cell_type, levels = proma_cell_type_levels)
wnt_data <- melt(wnt_data, id.vars='refined_clustering')
wnt_size_new <- c(sum(wnt_data$value[wnt_data$variable=='none']),
                  sum(wnt_data$value[wnt_data$variable=='Wnt.PCP']))
wnt_data$value_norm <- wnt_data$value/rep(wnt_size_new, each=dim(wnt_data)[1]/2)
wnt_size_new2 <- c(sum(wnt_data2['none']),
                  sum(wnt_data2['Wnt.PCP']))
wnt_data2 <- melt(wnt_data2, id.vars='proma_cell_type')
wnt_data2$value_norm <- wnt_data2$value/rep(wnt_size_new2, each=dim(wnt_data2)[1]/2)

wnt_data3 <- data.frame(type=rep('overall',2), 
                        variable=c('none', 'Wnt.PCP'),
                        value=wnt_size_new
                        )
wnt_frac <- wnt_data %>% select(refined_clustering, variable, value) %>%
        group_by(refined_clustering, variable) %>% summarise(total=sum(value))
wnt_frac2 <- wnt_data %>% select(refined_clustering, variable, value) %>%
  group_by(refined_clustering) %>% summarise(total=sum(value))
wnt_frac <- wnt_frac[wnt_frac$variable=='Wnt.PCP',]
wnt_frac$total = wnt_frac$total/wnt_frac2$total

wnt_fracM <- wnt_data2 %>% select(proma_cell_type, variable, value) %>%
  group_by(proma_cell_type, variable) %>% summarise(total=sum(value))
wnt_fracM2 <- wnt_data2 %>% select(proma_cell_type, variable, value) %>%
  group_by(proma_cell_type) %>% summarise(total=sum(value))
wnt_fracM <- wnt_fracM[wnt_fracM$variable=='Wnt.PCP',]
wnt_fracM$total = wnt_fracM$total/wnt_fracM2$total

g81 <- ggplot(wnt_frac, 
             aes(x=refined_clustering, y=total, fill=variable)) +
  scale_fill_manual( values=c('#b82531'),  
                     name='Wnt/PCP score') +
  guides(fill=FALSE) +
  # facet_wrap(~genetics)+
  geom_bar(position="dodge",stat="identity") +   
  labs(x='Cell type', y='Fraction Wnt/PCP\n positive cells (%)') +
  scale_y_sqrt(breaks=c(0,0.01,0.1),
               labels=c(0,1, 10)
  ) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g81
ggsave(filename = paste0(figure_path, 'Wnt_score_per_cell_type_frac.pdf'), plot = g81, width = 10, height=7)

g82 <- ggplot(wnt_fracM, 
              aes(x=proma_cell_type, y=total, fill=variable)) +
  scale_fill_manual( values=c('#b82531'),  
                     name='Wnt/PCP score') +
  guides(fill=FALSE) +
  # facet_wrap(~genetics)+
  geom_bar(position="dodge",stat="identity") +   
  labs(x='Cell type', y='Fraction Wnt/PCP\n positive cells (%)') +
  scale_y_sqrt(breaks=c(0,0.01,0.1),
               labels=c(0,1, 10)
  ) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g82
ggsave(filename = paste0(figure_path, 'Wnt_score_per_progen_cell_type_frac.pdf'), plot = g82, width = 10, height=7)


g8 <- ggplot(wnt_data, 
             aes(x=refined_clustering, y=value_norm, fill=variable)) +
  scale_fill_manual( values=c('#aaaaaa', '#b82531'), 
                     name='Wnt/PCP score',
                     breaks=c('none', 'Wnt.PCP'),
                     labels=c('negative', 'positive')
  ) +
  # facet_wrap(~genetics)+
  geom_bar(position="dodge",stat="identity") +   
  labs(x='Cell type', y='Population size (%)') +
  scale_y_sqrt(breaks=c(0,0.01,0.1,0.5,0.75,1),
                     labels=c(0,1, 10,50 ,75, 100)
  ) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g8
ggsave(filename = paste0(figure_path, 'Wnt_score_per_cell_type_dodge.pdf'), plot = g8, width = 10, height=7)

g9 <- ggplot(wnt_data2, 
             aes(x=proma_cell_type, y=value_norm, fill=variable)) +
  scale_fill_manual( values=c('#aaaaaa', '#b82531'), 
                     name='Wnt/PCP score',
                     breaks=c('none', 'Wnt.PCP'),
                     labels=c('negative', 'positive')
  ) +
  # facet_wrap(~genetics)+
  geom_bar(position="dodge",stat="identity") +   
  labs(x='Cell type', y='Population size (%)') +
  scale_y_sqrt(breaks=c(0,0.01,0.1,0.5,0.75,1),
                     labels=c(0,1,10,50 ,75, 100)
  ) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g9
ggsave(filename = paste0(figure_path, 'Wnt_score_per_progen_cell_type_dodge.pdf'), plot = g9, width = 10, height=7)

g10 <- ggplot(wnt_data3, 
             aes(x=type, y=value, fill=variable)) +
  scale_fill_manual( values=c('#aaaaaa', '#b82531'), 
                     name='Wnt/PCP score',
                     breaks=c('none', 'Wnt.PCP'),
                     labels=c('negative', 'positive')
  ) +
  scale_x_discrete(breaks = 'overall', labels='') +
  # facet_wrap(~genetics)+
  geom_bar(position="dodge",stat="identity") +   
  labs(x='Wnt/PCP signaling', y='Population size (total)') +
  scale_y_sqrt(breaks=c(1000, 10000, 60000)
  ) + 
  theme_classic() + 
  theme(axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        axis.text.x  = element_text(angle=45, hjust = 1, size=10))
g10
ggsave(filename = paste0(figure_path, 'Wnt_score_overall.pdf'), plot = g10, width = 10, height=7)

