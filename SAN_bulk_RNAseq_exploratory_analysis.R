# load the data
#load the phenotype table
pdata=read.table("phenotype_SAN_drugs.txt",sep = "",row.names = 1)
# load the expression data
edata=read.table("expression_SAN_drugs.txt",sep = "",row.names = 1)

# pdata=read.table("phenotype_SAN_drugs_with_viral.txt",sep = "",row.names = 1)
# edata=read.table("expression_SAN_drugs_with_viral.txt",sep = "",row.names = 1)

# table for factor/character variables
pdata$group=as.factor(pdata$group)

table(pdata$group)


# remove low expression data
edata = edata[rowMeans(edata) > 10, ]


# look at overall distributions
boxplot(edata[,1])

boxplot(log2(edata[,1]+1))

boxplot(log2(edata+1),col=2,range=0)

# transform the edata
edata=log2(edata+1)


# get the exact principal components use prcomp

pc=prcomp(edata)

pch<-c(15,16,17,18)
colors<-rep(c("purple","blue","red","green"),3)
group<-pdata$group
plot(pc$rotation[,1],pc$rotation[,2],pch=pch[group],col=colors[group],ylab="PC2 (5.23%)",xlab="PC1 (93.80%)")
legend("topright",legend=levels(group),pch=pch,col=colors)


# get the exact principal components use svd function
edata_centered = t(t(edata) - colMeans(edata))

svd = svd(edata_centered)

# look at the percent variance explained
plot(svd$d,ylab = "Singular value",col=2)

plot(svd$d^2/sum(svd$d^2),ylab = "Percent Variance Explained",col=2)


plot(svd$v[,1],svd$v[,2],pch=19,ylab="2nd PC",
     xlab="1st PC",col=as.numeric(pdata$group))


# clustering
# sample distance plot
library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(edata))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(edata))
colnames(sampleDistMatrix) <- paste(colnames(edata))
colors2 <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sample_distance_plot<-
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors2)


# tree clustering
hclust1=hclust(sampleDists)
plot(hclust1)
plot(hclust1,hang=0.5)


# 3D PCA
library(rgl)
plot3d(pc$rotation[,1:3],col=colors[group])

plot3d(pc$rotation[,1:3],col = colors[group],pch=pch[group],zlab = "PC3 (0.46%)", ylab = "PC2 (5.23%)", xlab = "PC1 (93.80%)")
text3d(pc$rotation[,1],pc$rotation[,2],pc$rotation[,3],texts = c(rownames(pc$rotation)))
legend3d("topright",legend=levels(group),pch=pch,col=colors)

plot3d(pc$rotation[,1:3],col =rep(c("purple","blue","red","green"),each=3),size=10,pch=pch[group],zlab = "PC3 (0.46%)", ylab = "PC2 (5.23%)", xlab = "PC1 (93.80%)")


rgl.snapshot('3D_PCA_im_de_IN.png', fmt = 'png')

rgl.postscript('3D_PCA_im_de_IN.pdf', fmt = 'pdf')


dir.create("animation_merge")
for (i in 1:360) {
  view3d(userMatrix=rotationMatrix(2*pi * i/360, 0, 1, 0))
  rgl.snapshot(filename=paste("animation_merge/frame-",
                              sprintf("%03d", i), ".png", sep=""))}


