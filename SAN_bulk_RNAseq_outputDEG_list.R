# load the data
#load the phenotype table
pdata=read.table("phenotype_SAN_new.txt",sep = "",row.names = 1)
# load the expression data
edata=read.table("expression_SAN_new.txt",sep = "",row.names = 1)

# table for factor/character variables
pdata$group=as.factor(pdata$group)

table(pdata$group)


# DEG analysis use EdgeR package
# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(edata)

y<-DGEList(counts=edata,genes=genelist)

# add transcript length
library(biomaRt)
human<-useMart(dataset="hsapiens_gene_ensembl",biomart='ensembl')
Attributes<-listAttributes(human)
table_human<-getBM(attributes = c('ensembl_gene_id','transcript_length','external_gene_name','entrezgene_id'),mart = human)

id<-match(rownames(y$genes),table_human$ensembl_gene_id)

y$genes$Length<-table_human$transcript_length[id]

RPKM<-rpkm(y)

RPKM<-as.data.frame(RPKM)

write.csv(RPKM,file = "SAN_SARSCOV2_mock_RPKM.csv")


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

# the DEA result for all the genes
# dea <- lrt$table

y <- estimateDisp(y, design, robust = TRUE)

fit<-glmQLFit(y,design,robust = TRUE)

IN_vs_MOCK<-makeContrasts(SARSCOV2-mock,levels = design)

res<-glmQLFTest(fit,contrast = IN_vs_MOCK)

toptag <- topTags(res, n = nrow(y$genes), p.value = 1)

dea <- toptag$table 

dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

write.csv(dea,file = "SAN_SARS_mock_DEG_edgeR.csv")



# Make a basic volcano plot

with(dea, plot(logFC, -log10(PValue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(dea, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(dea, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(dea, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)

id2<-match(dea$genes,table_human$ensembl_gene_id)
dea$names<-table_human$external_gene_name[id2]

names<-subset(dea, FDR<.05 & abs(logFC)>1)$names

with(subset(dea, FDR<.05 & abs(logFC)>1), textxy(logFC, -log10(PValue), labs=names, cex=.5))

legend("topleft",legend=c("FDR<0.05","Foldchange>2","FDR<0.05 & Foldchange>2"),pch=15,col=c("red","orange","green"),horiz = TRUE,text.width = 1,bty = "n")

