library(reshape2)         # for melting data
library(vegan)            # for mantel test
library(ggplot2)
library(phangorn)
library(ggforce)
library(tanggle)
library(RSplitsTree)
library(igraph)
library(ggtree)
library(ape)
library(tidytree)
library(RRtools)
library(ade4)
library(adegenet)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(diveRsity)
library(dplyr)
library(geosphere)
# library(ggfortify)
# library(ggmap)
library(ggrepel)
library(ggpubr)
library(ggthemes)
# library(heatmaply)
library(lattice)
library(openxlsx)
library(ozmaps)
library(RColorBrewer)
library(SNPRelate)
library(stringr)
library(tidyr)
library(fastDiversity)
library(reshape)

source('https://github.com/eilishmcmaster/SoS_functions/blob/33bd7065baa91b7e2a1800b6481d63573fb38d28/dart2svdquartets.r?raw=TRUE')
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")

topskip   <- 6
nmetavar  <- 18
RandRbase <- "" #main directory 
missingness <- 0.3
species   <- "RV" # folder name
dataset   <- "DRanu22-7417" # report name

basedir <- ""

d1        <- new.read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE)

d2 <- d1
meta      <- read.meta.data.full.analyses.df(d2, basedir, species, dataset)

d2        <- dart.meta.data.merge(d2, meta) 

#keep <- d2$sample_names[!is.na(d2$meta$analyses[,'keep'])]

#d2 <- remove.by.list(d2, keep)

d3        <- remove.poor.quality.snps(d2, min_repro=0.96, max_missing=0.3)%>% 
  remove.fixed.snps()
d4        <- sample.one.snp.per.locus.random(d3, seed=12345) 

d5 <- remove.by.missingness(d4, 0.5)

dms2 <- d5
dms <- dms2

m2 <- dms$meta$analyses %>% as.data.frame()
#na_samples <- d2$sample_names[is.na(d2$meta$analyses[,'keep'])]

# Convert to a data frame
na_samples_df <- data.frame(Sample_Names = na_samples)

# Save the list of NA samples to an Excel file
write.xlsx(na_samples_df, file = "na_samples.xlsx", rowNames = FALSE)

#### colours ####

sp_colours <-   named_list_maker(dms$meta$analyses[,'sp'], 'Spectral',11)
site_colours <-   named_list_maker(dms$meta$site, 'Paired',11)
region_colours <-   named_list_maker(dms$meta$analyses[,'region'], 'Paired',8)

sp_shapes <- (1:length(unique(dms$meta$analyses[,'sp']))) 
names(sp_shapes) <- unique(dms$meta$analyses[,'sp'])

region_shapes <- (1:length(unique(dms$meta$analyses[,'region']))) 
names(region_shapes) <- unique(dms$meta$analyses[,'region'])

#### kiship ####

kin <- individual_kinship_by_pop(dms2, RandRbase, species, dataset, dms$meta$analyses[,"region"], maf=0.1, mis=0.3, as_bigmat=TRUE)

#kin2 <- as.data.frame(kin) %>%mutate_all(~replace(.,.<0.35355339, 0)) #removes all of the pairwise connections that are k<0.35355339
#kin3 <- kin

kin2 <- as.data.frame(kin) %>%mutate_all(~replace(.,.<0.45, 0)) #removes all of the pairwise connections that are k<0.35355339
kin3 <- kin2

kin3<- as.data.frame(kin3) %>%mutate_all(~replace(.,.>0, 1))
network <- graph_from_adjacency_matrix(as.matrix(kin3), mode="undirected", diag=F,weighted=T) #makes the network based on k>0.45
plot(network)

ceb <- cluster_fast_greedy(network) # make the bubbles for the network plot
ceb_net_plot <- plot(ceb, network, vertex.label.color="transparent", vertex.size=2, edge.width=0.4) #make the network plot
ceb_net_plot
# 
clones <-as.data.frame(cbind(genet=ceb$membership, sample=ceb$names)) #get the clones from the network as a df
clones_out <- merge(clones, m2[,c("sample","lat","long","site","sp","region")], by="sample") #add some metadata
clones_out_filtered <- clones_out %>%
  group_by(genet) %>%
  filter(n() > 1) %>%
  ungroup()
# 
write.xlsx(clones_out_filtered, paste0(species,"/outputs/",species,"_coanalysis_coanalysis_clones_out.xlsx"), rowNames=FALSE, colNames=TRUE)
# 

clones_to_keep <- clones_out %>%
  group_by(genet) %>%
  slice(1) %>%  # This keeps the first clone for each genet group
  ungroup() %>%
  pull(sample)  # Extract the 'sample' column (or whichever column holds your sample names)

#Identify the clones to remove (those not in the 'clones_to_keep' list)
clones_to_remove <- clones_out %>%
  filter(!sample %in% clones_to_keep) %>%  # Remove the clones that are kept
  pull(sample)  # Extract the 'sample' column (or whichever column holds your sample names)

print(clones_to_remove)

# Save the list of clones to remove to a CSV file
write.csv(data.frame(Clones_to_Remove = clones_to_remove), "clones_to_remove.csv", row.names = FALSE)

# #### families ####
# 
# #families
# fam_network <- graph_from_adjacency_matrix(as.matrix(kin), mode="undirected", diag=F,weighted=T) #makes the network based on k>0.45
# plot(fam_network)
# 
# fam_ceb <- cluster_fast_greedy(fam_network) # make the bubbles for the network plot
# fam_ceb_net_plot <- plot(fam_ceb, fam_network, vertex.label.color="transparent", vertex.size=2, edge.width=0.4) #make the network plot
# fam_ceb_net_plot
# 
# families <-as.data.frame(cbind(family=fam_ceb$membership, sample=fam_ceb$names)) #get the clones from the network as a df
# families <- merge(families, m2[,c("sample","lat","long","site","sp","region")], by="sample") #add some metadata
# families <- families[order(as.numeric(families$family)),] #order the table by genet
# 
# write.xlsx(families, paste0(species,"/outputs/",species,"_coanalysis_coanalysis_families.xlsx"), rowNames=FALSE, colNames=TRUE)
# 
# families
# 
# # Remove duplicate clones from working data
# #close clones with highest quality data (least missingess)
# missingness_gt <- as.data.frame(rowSums(is.na(dms[["gt"]]))/ncol(dms[["gt"]])) #get missingness of all samples
# colnames(missingness_gt) <- "missingness"
# clone_missing <- merge(clones_out, missingness_gt, by.x="sample", by.y=0) # merge the clone data
# clone_missing <- clone_missing[order(clone_missing$missingness),] #order by missingness low to high
# 
# unique_genets <- distinct(clone_missing, genet, .keep_all=TRUE) #keeps the top result for each genet (lowest missingness)
# 
# #make a list removing duplicate clones
# non_clones <- dms$sample_names[which(!dms$sample_names %in% clones_out$sample)]
# approved_clones <- unique_genets$sample
# 
# clones_dereplicated <- c(non_clones, approved_clones) # list of sample to keep
# 
# mx <- as.data.frame(meta$analyses)
# mx$keep <- NA
# mx$keep[mx$sample %in% clones_dereplicated] <- 'y'
# 
# write.xlsx(mx, paste0(species,"/meta/",species,"_",dataset,"_meta2.xlsx"), rowNames=FALSE, colNames=TRUE)
# 
# # remove clones from dms
# dms_no_clones <- remove.by.list(dms, clones_dereplicated)
# dms <- dms_no_clones
dms_maf2 <- remove.by.maf(dms, 0.02)

#### pca ####
gen_d5 <- new("genlight", dms_maf2$gt) #convert df to genlight object for glPca function
gen_pca <- glPca(gen_d5, parallel=TRUE, nf=6) #do pca -- this method somehow allows the input to hav1 NAs

g_pca_df <- gen_pca[["scores"]] #extract PCs
g_pca_df2 <- merge(g_pca_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=FALSE) # some in DArT are not in meta?

pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:6]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)") #create names for axes

pca_plot1 <- ggplot(g_pca_df2, aes(x=PC1, y=PC2, colour=region, shape=region))+ xlab(pcnames[1])+ylab(pcnames[2])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=region_colours)+
  scale_shape_manual(values=region_shapes)
# geom_text_repel(mapping=aes(label=sp_label), color="black", min.segment.length = 0, size=2)

pca_plot1

pca_plot2 <- ggplot(g_pca_df2, aes(x=PC3, y=PC4, colour=region, shape=region))+ xlab(pcnames[3])+ylab(pcnames[4])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=region_colours)+
  scale_shape_manual(values=region_shapes)

pca_plot3 <- ggplot(g_pca_df2, aes(x=PC5, y=PC6, colour=region, shape=region))+ xlab(pcnames[5])+ylab(pcnames[6])+
  geom_point(size=2)+
  theme_few()+geom_vline(xintercept = 0, alpha=0.2)+geom_hline(yintercept = 0, alpha=0.2)+
  labs(colour="", shape="")+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right",
        legend.text = element_text(face="italic"),
        axis.title = element_text(size=10), axis.text = element_text(size=8))+
  guides(colour = guide_legend(title.position = "top"))+
  scale_colour_manual(values=region_colours)+
  scale_shape_manual(values=region_shapes)

all3_pca_plots <- ggarrange(pca_plot1, pca_plot2, pca_plot3,# labels=c("A","B","C"),
                            common.legend = TRUE, ncol=3, legend = "right")
all3_pca_plots

ggsave(paste0(species,"/outputs/",species,"_coanalysis_pca_maf2.pdf"),
       all3_pca_plots, width = 35, height = 10, units = "cm", dpi=600)


dist_matrix <- as.matrix(dist(dms$gt, diag=TRUE))
Heatmap(dist_matrix)

# Heatmap(dms$gt)



#### splitstree ####

splitstree(dist(dms$gt, method = "euclidean"), paste0(species,"/outputs/",species,"_coanalysis_splits.nex"))


# Read network data from Nexus file
Nnet <- phangorn::read.nexus.networx(paste0(species,"/outputs/",species,"_coanalysis_splits.nex"))

x <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2], 
                sample=rep(NA, nrow(Nnet$.plot$vertices)))

x[Nnet$translate$node,"sample"] <- Nnet$translate$label
x <- merge(x, m2, by="sample", all.x=TRUE, all.y=FALSE)
# x$svdq_pop[!is.na(x$sample)&is.na(x$svdq_pop)] <- "ungrouped"

net_x_axis <- max(x$x)-min(x$x)
net_y_axis <- max(x$y)-min(x$y)

hull <- x %>% group_by(region) %>% 
  slice(chull(x, y))


x_filtered <- x[!is.na(x$sp),]


splitstree_plot_svdq <- ggplot(Nnet, aes(x = x, y = y)) +
  geom_shape(data = hull, alpha = 0.4, expand = 0.02, radius = 0.02,
             aes(fill = region, color = "transparent")) +
  geom_point(data = x_filtered, shape=19, color="white", size=3) +
  geom_splitnet(layout = "slanted", size = 0.2) +
  geom_point(data = x_filtered, aes(color = sp, shape = sp)) +
  scale_shape_manual(values = sp_shapes, na.translate = FALSE)+
  scale_fill_manual(values = region_colours, na.translate = FALSE) +
  scale_colour_manual(values = sp_colours, na.translate = FALSE) +
  theme_void() +
  expand_limits(x = c(min(x_filtered$x) - 0.2 * net_x_axis, max(x_filtered$x) + 0.2 * net_x_axis),
                y = c(min(x_filtered$y) - 0.1 * net_y_axis, max(x_filtered$y) + 0.1 * net_y_axis)) +
  theme(
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.text = element_text(face = "italic"),  # Italicize legend labels
    legend.key.size = unit(0.75, 'lines')
  ) +
  coord_fixed() +
  labs(color = "Species", shape = "Species", fill = "") +
  guides(colour = guide_legend(title.position = "top", nrow = 4, override.aes = list(fill = NA, linetype = 0)),
         fill = "none",
         shape = guide_legend(title.position = "top", nrow = 4))+
  geom_tiplab2(size=1.5, hjust=-0.2)


splitstree_plot_svdq 
# Save the plot
ggsave(paste0(species,"/outputs/",species,"_coanalysis_splitstree.pdf"),
       splitstree_plot_svdq , width = 25, height = 25, units = "cm", dpi = 600)



#### stats ####

fst_stats <- faststats(dms$gt, 
                       genetic_group_variable = dms$meta$analyses[,"region"],
                       site_variable = dms$meta$analyses[,"site"], 
                       max_missingness=0.3, 
                       maf = 0.05)
fst_stats <- fst_stats %>% arrange(genetic_group)
fst_stats


fst_stats2 <- faststats(dms$gt, 
                        genetic_group_variable = dms$meta$analyses[,"sp"],
                        site_variable = dms$meta$analyses[,"region"], 
                        max_missingness=0.3,
                        minimum_n = 2)
fst_stats2 <- fst_stats2 %>% arrange(genetic_group)
fst_stats2 

hos <- (rowSums(dms$gt==1, na.rm=TRUE)/ rowSums(!is.na(dms$gt)))
hos2 <- cbind(hos, dms$meta$analyses)
ggplot(hos2, aes(x = as.numeric(hos))) +
  geom_histogram( fill = "skyblue", color = "black") +  # Adjust binwidth as needed
  facet_wrap(~ region) +                                              # Facet by the 'sp' column
  theme_few()



### fst ####

sppop_freq <- as.data.frame(table(dms$meta$site))
not_n1_sites <- as.vector(sppop_freq[sppop_freq$Freq<2,1]) #remove sps where n<=1
not_n1_sites <- c(not_n1_sites, NA) 
not_n1_samples <- dms$sample_names[which(!(dms$meta$site %in% not_n1_sites))]
# not_n1_samples <- not_n1_samples[not_n1_samples!='NSW1191727']

fst_dms <- remove.by.list(dms, not_n1_samples)

# fst_dms <- dms
length(fst_dms$sample_names)
length(fst_dms$locus_names)

gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pFst      <- population.pw.Fst(fst_dms, fst_dms$meta$site, RandRbase,species,dataset, maf_val=0.05, miss_val=0.3) #calculates genetic distance
pS        <- population.pw.spatial.dist(fst_dms,  fst_dms$meta$site) #calculates geographic distance between populations

####plot IBD plot

# Make self comparisons NA
diag(pFst$Fst) <- NA
diag(pS$S) <- NA

#Mantel test 
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 10000, na.rm = TRUE) #mantel test, finds if matrices are signficantly similar
man

# mantel plot
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <-Fst_sig$Geo_dist/1000 

# # adding metadata for sites
Fst_sig2 <- merge(Fst_sig, distinct(m2[,c("site","region")]), by.x="Var1", by.y="site", all.y=FALSE)
Fst_sig2 <- merge(Fst_sig2, distinct(m2[,c("site","region")]), by.x="Var2", by.y="site", all.y=FALSE)
Fst_sig2$same_group <- ifelse(Fst_sig2$region.x == Fst_sig2$region.y, "Within", "Between")

fstp1 <- ggplot(Fst_sig2, aes(x= Geo_dist2, y=Fst, colour=same_group))+geom_point(size=1, alpha=0.3)+
  labs(x="Distance (km)", y="FST", colour="Comparison")+
  # facet_zoom(x=Geo_dist2<25, zoom.size=1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
fstp1

ggsave(paste0(species,"/outputs/",species,"_manning_fst.pdf"),
       fstp1, width = 15, height = 15, units = "cm", dpi=600)

paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif)

# Make heatmaps
# geo dist
geo_d <-pS$S #this is a square matrix
mat <- geo_d/1000 # convert to km 

#FST
mat2 <-pFst$Fst
diag(mat2) <- NA

order_hm <- Heatmap(mat2,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE)
od <- colnames(mat2)[column_order(order_hm)]

mat = mat[od, od]
mat2 = mat2[od, od]

agg <- unique(fst_dms$meta$analyses[, c("site", "region","sp")]) %>% as.data.frame() # create aggregated df of pop_largeecies and site
mat2 <- merge(mat2, agg, by.x=0, by.y="site", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat2) <- mat2$Row.names

mat2$Row.names <- NULL
mat2 <- mat2[match(colnames(mat2)[1:nrow(mat2)],rownames(mat2)),]

row_ann <- rowAnnotation(region = mat2$region,
                         col=list(region=region_colours),
                         na_col="white",
                         annotation_legend_param = list(labels_gp=gpar(fontface="italic",fontsize=8),
                                                        title_gp=gpar(fontsize=10)),
                         annotation_name_gp = gpar(fontsize = 0),
                         annotation_name_side="top")

row_ann2 <- rowAnnotation(sp = mat2$sp,
                          col=list(sp=sp_colours),
                          na_col="white",
                          annotation_legend_param = list(labels_gp=gpar(fontface="italic",fontsize=8),
                                                         title_gp=gpar(fontsize=10)),
                          annotation_name_gp = gpar(fontsize = 0),
                          annotation_name_side="top")

bottom_ann <- HeatmapAnnotation(region = mat2$region, col = list(region = region_colours),
                                annotation_name_gp = gpar(fontsize = 0),
                                annotation_legend_param = list(labels_gp=gpar(fontface="italic", fontsize=8),
                                                               title_gp=gpar(fontsize=10)),
                                annotation_name_side="left",
                                na_col = "white")

# specify fst heatmap colours 
gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))


#specify geo heatmap colours
palette <-  colorRamp2(c(0, max(mat, na.rm=TRUE)), c("white", "#80B1D3"))

geo <- Heatmap(mat,rect_gp = gpar(type = "none"),
               width = nrow(mat)*unit(6, "mm"),
               height = nrow(mat)*unit(6, "mm"),
               col=palette,na_col="white",
               bottom_annotation = bottom_ann,
               row_names_gp = gpar(fontsize = 8, fontface="italic"),
               column_names_gp = gpar(fontsize = 8),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               name="Distance (km)",
               heatmap_legend_param = list(title_gp = gpar(fontsize = 10),
                                           labels_gp = gpar(fontsize = 8)),
               # cluster_rows = TRUE, 
               # cluster_columns = TRUE,
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(i >= j) {
                   grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                   grid.text(sprintf("%.f", mat[,1:nrow(mat)][i, j]), x, y, gp = gpar(fontsize = 6))
                 }
               }
)

# make fst heatmap
gene <- Heatmap(as.matrix(mat2[,1:nrow(mat2)]), rect_gp = gpar(type = "none"),
                width = nrow(mat2)*unit(6, "mm"),
                height = nrow(mat2)*unit(6, "mm"),
                right_annotation = c(row_ann, row_ann2),
                col=gene_col,na_col="grey",
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 0),
                border_gp = gpar(col = "black", lty = 1),
                name="FST",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10),
                                            labels_gp = gpar(fontsize = 8)),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    grid.text(sprintf("%.2f", mat2[,1:nrow(mat2)][i, j]), x, y, gp = gpar(fontsize = 6))
                  }
                })

gene_width <- nrow(mat2)*unit(6, "mm")

draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE)

# Set the file name and parameters
filename <- paste0(species,"/outputs/",species,"_fst_plot.pdf")
width <- ncol(mat2) * 0.25 +2
height <- ncol(mat2) * 0.25 +0.5
dpi <- 300
units <- "cm"

# Set up the PNG device
pdf(filename, width = width, height = height)

# Draw the plot
draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE)

# Turn off the PNG device
dev.off()

##### migrate ####

make_genepop_file <- function(dms, maf,missing, group, grouping){
  
  dmsx <- dms %>%
    remove.poor.quality.snps(., min_repro=0.96,max_missing=missing) %>%
    remove.by.maf(., maf)
  
  ds <- dmsx$gt
  
  print(paste("Loci:", ncol(ds)))
  print(paste("Samples:",nrow(ds)))
  
  
  # if(ncol(ds) >=50){# make into genepop format
  old <- c("0","1","2", NA)
  new <- c("0101","0102","0202","0000")
  ds[ds %in% old] <- new[match(ds, old, nomatch = 0000)]
  
  # populations
  pops <- unique(grouping) # get population names
  
  # write genepop file
  gf <- paste0(species, "/popgen/genepop_",group,".gen") # make genepop file path
  
  cat(paste0("genepop file: ",species, " with MAF ", paste0(maf)), # first line of genepop file
      file=gf,sep="\n")
  cat(colnames(ds),file=gf,sep="\n", append=TRUE) # one loci name per line
  
  remove <- c() # vector for populations excluded from the analysis
  
  for (i in 1:length(pops)){ #loop for making the population groups
    if (length(which(grouping %in% pops[i]))<=1){ # find if the population is n=1
      cat("Removing population ", pops[i], " due to n=1")
      remove <- c(remove, pops[i]) # add the pop name to remove vector
    }else{
      cat("pop",file=gf,sep="\n", append=TRUE) # add the data to the genepop file
      df <- ds[which(grouping %in% pops[i]),]
      for (j in 1:nrow(df)){
        cat(c(paste0(pops[i],","),df[j,], "\n"),file=gf,sep="\t", append=TRUE)
      }
    }
    
  } #end of pops loop
  return(gf)
}

#Migration was estimated using the divMigrate method (19) with Jost’s D metric of differentiation (47), as implemented in the R package diveRsity (46). This approach uses allele frequency differences between population pairs to estimate rates of migration in each direction; note that these rates are relative to other population pairs in the same data set and cannot be compared across data sets.

c <- make_genepop_file(dms, maf=0.05,missing=0.2, species, dms$meta$analyses[,"site"])

v <- diveRsity::divMigrate(infile=c, outfile=NULL, stat="gst",plot_network=TRUE, filter_threshold = 0, boots=1000, para=TRUE)

# Save a single object to a file
saveRDS(v, paste0(species,"/outputs/plots/gst_m_10000_reps.rds"))
# Restore it under a different name
v <- readRDS(paste0(species,"/outputs/plots/gst_m_10000_reps.rds"))

# d_mig <- v$gRelMig #v$dRelMig # all migration values
d_mig <- v$gRelMigSig #v$dRelMigSig # significant only

gpop <- read.genepop(c)

# from is rows, to is columns
colnames(d_mig) <- unique(gpop@pop) 
rownames(d_mig) <- unique(gpop@pop) 

# Filter by significance -- significance is if there is a significant difference in the directions
# Test for overlap of the estimated 95% confidence intervals. Where there is no overlap, the directional gene flow components are said to be significantly different (asymmetric).
# mig_sig <- v$gRelMigSig #v$dRelMigSig
# colnames(mig_sig) <- unique(gpop@pop)
# rownames(mig_sig) <- unique(gpop@pop)

qgraph::qgraph(d_mig,legend = TRUE, edge.labels = TRUE, curve = 2.5, mar = c(2, 2, 5, 5))

long_mig <- melt(d_mig)
long_mig2 <- long_mig[(long_mig$value>0.2),] # only keep connections where m>0.2
# long_mig2 <- long_mig2[(long_mig2$Var1!=long_mig2$Var2), ] # remove self connections


meta_agg <- m2 %>%
  group_by(site, sp, region) %>%
  dplyr::reframe(lat = mean(as.numeric(lat), na.rm=TRUE),
                 long = mean(as.numeric(long),na.rm=TRUE))

long_mig2 <- merge(long_mig2, meta_agg, by.x = "Var1", by.y = "site", all.y = FALSE)
long_mig2 <- merge(long_mig2, meta_agg, by.x = "Var2", by.y = "site", all.y = FALSE)
long_mig2 <- distinct(long_mig2)  

colnames(long_mig2)[1:3] <- c("from", "to", "m")

# Plot 2: long_mig
long_mig2 <- long_mig2[order(long_mig2$m), ]


oz_states <- ozmaps::ozmap_states

aus_map_point <- ggplot() +
  geom_sf(data = oz_states, color="black", size=2) 

gst_no_map <- ggplot() + coord_cartesian() + coord_fixed() +
  theme_few() +
  geom_curve(data = long_mig2[long_mig2$m>0.5,], #
             aes(x = long.x, y = lat.x,
                 xend = long.y, yend = lat.y, colour = m), 
             size = 0.5, na.rm = TRUE, curvature = 0.3, 
             arrow = arrow(angle = 20, ends = "first", type = "open", length = unit(2, "mm"))) +
  scale_color_gradient(low = "white", high = "red") + # midpoint = 0.5, mid = "blue",
  geom_point(data = meta_agg, mapping = aes(x = long, y = lat), colour = "black") +
  labs(x = "Longitude", y = "Latitude", colour = "Migration (m)") + guides(size = "none", alpha = "none") +
  ggrepel::geom_label_repel(data = meta_agg, aes(x = long, y = lat, label = site),
                            min.segment.length = 0.25, color = "black", fill = "white", size = 3,
                            segment.colour = "white", alpha = 0.9, label.size = 0, nudge_y = 0.003) +
  theme(legend.position = "bottom") 


gst_no_map




