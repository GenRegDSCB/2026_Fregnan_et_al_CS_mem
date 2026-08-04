library("DESeq2")
library("dplyr")
library("RColorBrewer")
library("ggfortify")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

## Study 1 PCA -------
cdata <- "0_Input_FPKM.tsv"
pdata <- "0_Input_metadata.tsv"

countData <- read.delim(cdata, row.names=1, check.names=FALSE)
sdata <- read.delim(pdata, row.names=1, check.names=FALSE)

sdata <- dplyr::filter(sdata, Study %in% c("S1", "Both"))
sdata <- dplyr::filter(sdata, Tissue %in% c("CN", "MEM"))
sdata <- dplyr::filter(sdata, Class %in% c("PENIS_CTRL", "PENIS_CS_FLAT", "PENIS_CS_GRAT", "PENIS_CS_SCA", "PENIS_DEN"))

countData <- countData[, row.names(sdata)]
countData <- countData[rowSums(countData > 0) > 0, ]

cdatat <- data.frame(t(countData), check.names = F)

cdata.pca<-prcomp(log(cdatat+0.1,2), center=T, scale=T)

dat <- data.frame(cbind(sdata, cdatat, ID=row.names(sdata)), check.names=F)
colorh<-colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(12)

class.c<-colorRampPalette(brewer.pal(12,"Paired"))(24)

class <- c("CN_CTRL"=class.c[2], "CN_CS_MEM_FLAT"=class.c[10],
           "CN_CS_MEM_GRAT"=class.c[13], "CN_CS_MEM_SCA"=class.c[17],
           "PENIS_CTRL"=class.c[2], "PENIS_CS_FLAT"=class.c[10],
           "PENIS_CS_SCA"=class.c[17], "PENIS_CS_GRAT"=class.c[13],
           "PENIS_DEN"=class.c[22])

dat$Class <- factor(dat$Class, levels=c("CN_CTRL", "CN_CS_MEM_FLAT", "CN_CS_MEM_GRAT", "CN_CS_MEM_SCA"))

autoplot(cdata.pca, data = dat, colour = "Class", size=2, 
         alpha=2, label=F, label.size=1, shape="Tissue") + 
  theme_bw()+
  scale_color_manual(values = class)+
  theme(legend.position="right")

## Sample correlation
library("corrplot")

a <- cor(countData, method="spearman")

colcor <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(256))

corrplot(a, 
         type="upper", 
         tl.cex=0.6, 
         tl.col="black", 
         col=colcor, 
         order = "hclust",
         hclust.method="ward.D2",
         addCoef.col="white",
         number.cex=0.7,)

## Study 2 PCA -------
cdata <- "0_Input_FPKM.tsv"
pdata <- "0_Input_metadata.tsv"

countData <- read.delim(cdata, row.names=1, check.names=FALSE)
sdata <- read.delim(pdata, row.names=1, check.names=FALSE)

sdata <- dplyr::filter(sdata, Study %in% c("S2", "Both"))

countData <- countData[, row.names(sdata)]
countData <- countData[rowSums(countData > 0) > 0, ]

cdatat <- data.frame(t(countData), check.names = F)

cdata.pca<-prcomp(log(cdatat+0.1,2), center=T, scale=T)

dat <- data.frame(cbind(sdata, cdatat, ID=row.names(sdata)), check.names=F)
colorh<-colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(12)

autoplot(cdata.pca, data = dat, colour = "Class", size=2,
         alpha=2, label=F, label.size=1, shape="Tissue") + 
  theme_bw()+
  #scale_color_manual(values = colorh)+
  theme(legend.position="right")

## Study 1 Upset plot -------
a <- read.delim("~/Desktop/a", check.names=F)
b <- list()

for(i in 1:ncol(a)){
  b[[i]] <- a[,i][is.na(a[,i])==F]
  names(b)[i] <- names(a)[i]
}

UpSetR::upset(UpSetR::fromList(b), nsets = 10, number.angles = 0, order.by = "freq", set_size.numbers_size = 15, text.scale = c(1.3, 1.3, 1.3, 1.3, 1.2, 1.3), sets.bar.color = "black", att.color = "black", main.bar.color = "black")

### Heatmap DE ----
library("ComplexHeatmap")
library("circlize")

DE <- read.delim("0_Input_DE_membrane.tsv", row.names=1, check.names=F)

colorh<-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"PuOr"))(256))

ch<-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(11))

DE <- DE[unique(genes$Gene),]

cellh <- c("APOC1+ FB"=ch[1], "CC EC"=ch[2], "CC SMC"=ch[3], "COMP+ FB"=ch[4],"MAC"=ch[5], "PC"=ch[6], "PI16+ FB"=ch[7], "SWC"=ch[8],"T"=ch[9], "Vessel EC"=ch[10], "Vessel SMC"=ch[11])
  
row_ha = rowAnnotation(Cell=DE$cluster,
                              col = list(Cell = cellh),
                              border = TRUE,
                              simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize =5),
                              annotation_legend_param = list(border="black"))

row_ha = rowAnnotation(Class=genes$Class,
                       border = TRUE,
                       simple_anno_size = unit(0.3, "cm"),
                       annotation_name_gp = gpar(fontsize =5),
                       annotation_legend_param = list(border="black"))

h1 <- Heatmap(DE[,1:3],
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 7),
        column_title_gp = gpar(fontsize = 8),
        #left_annotation = row_ha,
        cluster_columns = F,
        cluster_rows = T,
        name = "log2FC",
        border = "black",
        column_names_rot = 45,
        col=colorRamp2(breaks=c(-4,-2,0,2,4), 
                       colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
        heatmap_legend_param = list(border = "black"),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        clustering_distance_columns =  "euclidean", 
        clustering_distance_rows =  "euclidean")

### Heatmap DE Z-score ----
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

# Table containing the count data
cdata <- "0_Input_FPKM.tsv"
genes <- "0_Input_Selected_gene_class.tsv"

# Table containing the covariates
pdata <- "0_Input_metadata.tsv"

#### Read the input files
countData <- read.delim(cdata, row.names=1, check.names=FALSE)
sdata <- read.delim(pdata, row.names=1, check.names=FALSE)
genes <- read.delim(genes, check.names=FALSE)

DE <- read.delim("0_Input_DE_membrane.tsv", row.names=1, check.names=F)

sdata <- dplyr::filter(sdata, Study %in% c("S1", "Both"))
sdata <- dplyr::filter(sdata, Tissue %in% c("MEM"))

countData <- countData[row.names(DE), row.names(sdata)]
countData <- countData[unique(genes$Gene),]

colorh<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(256))

class.c<-colorRampPalette(brewer.pal(12,"Paired"))(24)

class <- c("CN_CTRL"=class.c[2], "CN_CS_MEM_FLAT"=class.c[10],
           "CN_CS_MEM_GRAT"=class.c[13], "CN_CS_MEM_SCA"=class.c[17],
           "PENIS_CTRL"=class.c[2], "PENIS_CS_FLAT"=class.c[10],
           "PENIS_CS_SCA"=class.c[17], "PENIS_CS_GRAT"=class.c[13],
           "PENIS_DEN"=class.c[22])

expression <- scale(t(log10(countData+1)))

column_ha = HeatmapAnnotation(Class=sdata$Class,
                              col = list(Class = class),
                              border = TRUE,
                              simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize =5),
                              annotation_legend_param = list(border="black"))

a <- read.delim("0_Input_DE_membrane_target.tsv")

mark_at = which(colnames(expression) %in% a$ID)

ra = rowAnnotation(foo = anno_mark(at = mark_at, labels = a$ID, labels_gp = gpar(fontsize=3), link_width=unit(5, "mm")))

h2 <- Heatmap(t(expression),
              top_annotation = column_ha,
              #right_annotation = ra,
              row_names_gp = gpar(fontsize = 4),
              column_names_gp = gpar(fontsize = 6),
              cluster_columns = T,
              cluster_rows = T,
              name = "Z-score",
              border = "black",
              column_names_rot = 45,
              col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[255])),
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              clustering_distance_columns = "euclidean",
              clustering_distance_rows = "euclidean")

h1 + h2

#### Boxplot miRNA
library("ggpubr")
dat <- data.frame(miR=t(countData["hsa-miR-452-5p",]), Class=phenoData$Class, Time=phenoData$Time, PSQI=phenoData$PSQI_cat)

names(dat)[1] <- "miR"

ggplot(dat, aes(x=Time, y=log10(miR), col=Class))+
  geom_boxplot()+
  geom_jitter(width=0.15)+
  theme_bw()+
  facet_grid(cols = vars(Class), rows = vars(PSQI))+
  geom_pwc()

### Scatterplot
library("ggpubr")
a <- read.delim("~/Desktop/a")

colorh<-colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(12)

ggplot(a, aes(x=EV, y=MIV))+
  geom_point(aes(col=Coherent))+
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_cor(method="spearman")+
  labs(x="log2FC(EV)", y="log2FC(MIV)")+
  scale_color_manual(values = c(colorh[1], colorh[5], colorh[2], colorh[6]))+
  ggtitle(label = "7d vs 3d")
  
### Heatmap dotplot enrichment
a <- read.delim("~/Desktop/a")

a$Description <- factor(a$Description, levels=rev(unique(a$Description)))

ggplot(a, aes(y=Description, col=-Log.q.value., size=Enrichment, x=GeneList))+
  geom_point()+
  theme_bw()+
  ggtitle(label = "7d vs 3d")

### Volcano plot ----
a <- read.delim("~/Desktop/a")

pal2 <- colorRampPalette(brewer.pal(9, "RdBu"))(256)

ggplot(a, aes(x=log2FC, y=-log10(pval), fill=log2FC, label = Symbol)) +
  geom_point(shape=21, alpha=0.7, col="black") +
  theme_bw() +
  scale_fill_gradientn(colours=rev(pal2), limits=c(-11,11)) +
  labs(x="log2FC", y="-log10(adj. p-value)", col="log2FC") +
  geom_hline(yintercept = 1.3, linetype = "dashed", col="red")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", col="red")+
  xlim(c(-11,11))

### Booble plot ----
a <- read.delim("0_Input_Enrichment_terms.tsv")

pal2 <- colorRampPalette(brewer.pal(9, "YlOrRd"))(256)

a$Class <- factor(a$Class, levels=c("Axon/Synaptic", "Vascular", "Inflammatory", "Metabolic/Ischemic"))

a$Description <- factor(a$Description, levels=rev(unique(a$Description)))

ggplot(a, aes(x=Comparison, y=Description, fill=log10(Enrichment))) +
  geom_point(shape=21, alpha=0.7, col="black", aes(size=Gene)) +
  theme_bw() +
  scale_fill_gradientn(colours=pal2) +
  facet_grid(rows=vars(Class), space = "free", scale="free")+
  labs(x="Comparison", y="Term", size="N° DE genes")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  geom_text(aes(label=-LogP), size=2)


## Heatmap selected genes
genes <- read.delim("0_Input_Selected_gene_class.tsv")

## Alluvial plot
library("tidyverse")
library("ggalluvial")

gene_ontology <- read.delim("~/Desktop/a")
sorted_genes <- read.delim("~/Desktop/b")

# Match the ontology annotation to the sorted gene list
gene_positions <- sorted_genes %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(gene_order = row_number(), gene_y = rev(gene_order))

plot_data <- gene_positions %>%
  select(Gene, gene_order, gene_y) %>%
  inner_join(gene_ontology, by = "Gene") %>%
  distinct(Gene, Class, .keep_all = TRUE)

ontology_positions <- plot_data %>%
  count(Class, name = "n_genes") %>%
  arrange(desc(n_genes), Class) %>%
  mutate(
    ontology_y = seq(
      from = nrow(gene_positions),
      to = 1,
      length.out = n()))

plot_data <- plot_data %>%
  left_join(ontology_positions %>% select(Class, ontology_y), 
            by = "Class")

plot_data$Class <- factor(plot_data$Class, levels=c("Axon/Synaptic", "Vascular", "Inflammatory", "Metabolic/Ischemic"))

ggplot() +
  geom_curve(data = plot_data,
    aes(x = 0, y = gene_y, xend = 1, 
        yend = ontology_y, colour = Class),
    curvature = 0.15, alpha = 0.55, linewidth = 0.4) +
  geom_point(data = gene_positions, aes(x = 0, y = gene_y),
    size = 1.8) +
  geom_text(data = gene_positions,aes(x = -0.03, y = gene_y,
      label = Gene), hjust = 1, size = 2) +
  geom_point(data = ontology_positions, aes(x = 1, y = ontology_y),size = 2) +
  geom_text(data = ontology_positions, aes(x = 1.03, y = ontology_y, label = paste0(Class, " (", n_genes, ")")), hjust = 0, size = 3) +
  scale_x_continuous(
    limits = c(-0.3, 1.35), breaks = c(0, 1),labels = c("Sorted genes", "Ontology")) +
  
  labs(x = NULL, y = NULL,colour = "Ontology") +
  theme_void() +
  theme(axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "none", plot.margin = margin(10, 80, 10, 80))
