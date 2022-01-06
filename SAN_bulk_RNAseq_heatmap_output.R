# load the data
#load the phenotype table
pdata=read.table("phenotype_SAN_drugs.txt",sep = "",row.names = 1)
# load the expression data
edata=read.table("expression_SAN_drugs.txt",sep = "",row.names = 1)

# table for factor/character variables
pdata$group=as.factor(pdata$group)

table(pdata$group)



# DEG analysis use EdgeR package
# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(edata)

y<-DGEList(counts=edata,genes=genelist)

# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]


# Normalization
y <- calcNormFactors(y, method="TMM")

y$samples$group <- pdata$group

group<-pdata$group

design <- model.matrix(~0+group)

colnames(design) <- levels(group)

design


# add gene symbol
library(org.Hs.eg.db) # options(connectionObserver = NULL)
y$genes$symbol <- mapIds(org.Hs.eg.db, rownames(y),
                         keytype="ENSEMBL", column="SYMBOL")
head(y$genes)

y <- y[!is.na(y$genes$symbol), ]


# heatmap of selected pacemaker genes
SAN_marker<-read.table("pacemaker_markers.txt",sep = "",header = TRUE)
SAN_marker<-SAN_marker$Gene

logCPM <- cpm(y, prior.count=2, log=TRUE)

rownames(logCPM) <- y$genes$symbol

colnames(logCPM) <- paste(y$samples$group, 1:3, sep="-")

logCPM<-logCPM[,4:9]

index<-which(rownames(logCPM) %in% SAN_marker)

logCPM <- logCPM[index,]

coolmap(logCPM,margin=c(10,6))


# heatmap of inflammation genes
inflammation<-read.csv("inflammatory_response.csv")

inflammation_gene<-toupper(inflammation$inflammatory_response)


logCPM <- cpm(y, prior.count=2, log=TRUE)

rownames(logCPM) <- y$genes$symbol

colnames(logCPM) <- paste(y$samples$group, 1:3, sep="-")

logCPM<-logCPM[,4:9]

index2<-which(rownames(logCPM) %in% inflammation_gene)

logCPM <- logCPM[index2,]

coolmap(logCPM,margin=c(10,6))

# heatmap of chemokine genes
chemokine<-read.csv("chemokine.csv")

chemokine_gene<-toupper(chemokine$chemokine)


logCPM <- cpm(y, prior.count=2, log=TRUE)

rownames(logCPM) <- y$genes$symbol

colnames(logCPM) <- paste(y$samples$group, 1:3, sep="-")

logCPM<-logCPM[,4:9]

index3<-which(rownames(logCPM) %in% chemokine_gene)

logCPM <- logCPM[index3,]

coolmap(logCPM,margin=c(10,6))


# heatmap of immune response genes
immune<-read.csv("immune_response.csv")

immune_response_gene<-toupper(immune$immune_response)

logCPM <- cpm(y, prior.count=2, log=TRUE)

rownames(logCPM) <- y$genes$symbol

colnames(logCPM) <- paste(y$samples$group, 1:3, sep="-")

logCPM<-logCPM[,4:9]

index4<-which(rownames(logCPM) %in% immune_response_gene)

logCPM <- logCPM[index4,]

coolmap(logCPM,margin=c(10,6))



# heatmap of ferroptosis genes
ferroptosis<-read.csv(file = "ferroptosis.csv",header = TRUE)

ferroptosis_gene<-ferroptosis$ferroptosis[1:10]

logCPM <- cpm(y, prior.count=2, log=TRUE)

rownames(logCPM) <- y$genes$symbol

colnames(logCPM) <- paste(y$samples$group, 1:3, sep="-")

logCPM<-logCPM[,4:9]

index5<-which(rownames(logCPM) %in% ferroptosis_gene)

logCPM <- logCPM[index5,]

coolmap(logCPM,margin=c(10,6))





# heatmap use DESeq2 package
# perform DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(edata, pdata, design = ~ group)
dds <- DESeq(dds)

# transform raw counts into normalized values
# rlog transformed and variance stabilization
rld <- rlogTransformation(dds,blind = T)
vsd <- varianceStabilizingTransformation(dds,blind = T)

# heatmap of top 100 expressed genes
library("RColorBrewer")
library("gplots")
# 100 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:100]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="100 Top Expressed Genes Heatmap")

# heatmap of selected genes
table_human<-read.table("table_human_index.csv",header = TRUE,sep = ",",row.names = 1)

gene_matrix<-assay(vsd)

id<-match(rownames(gene_matrix),table_human$ensembl_gene_id)

rownames(gene_matrix)<-table_human$external_gene_name[id]

pacemaker_markers<-read.table("pacemaker_markers.txt",header = TRUE,sep = "",row.names = NULL)

index<-which(rownames(gene_matrix) %in% pacemaker_markers$Gene)

gene_matrix<-gene_matrix[index,]

heatmap.2(gene_matrix, col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6,main="pacemaker marker Genes Heatmap")


