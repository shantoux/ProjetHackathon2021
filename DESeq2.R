library(DESeq2)
library(ggplot2)

#Récupération des comptages
files <- list.files(path="counts_files/", pattern="*0.counts$", full.names = T)
gene_names <- read.table(files[1], sep="", header=F)[-1,1]# get gene names
df    <- as.data.frame(do.call(cbind,lapply(files,function(fn)
  read.table(fn,header=F,sep="", check.names = F)[,7])))
colnames(df) <- df[1,] #get sample names
df = df[-1,] #delete from data frame
rownames(df) = gene_names
df[]<- lapply(df, as.integer)

info <- data.frame(strain=c("mutated","mutated","mutated","wildtype","wildtype","wildtype","wildtype","mutated"),
                    row.names=c("SRR628582.bam","SRR628583.bam","SRR628584.bam","SRR628585.bam","SRR628586.bam","SRR628587.bam","SRR628588.bam","SRR628589.bam"))


#retrait de gènes non exprimés
df = df[rowSums(df)>0,]
#matrix de comptage
matrix_df=as.matrix(df)
#construction de l'objet
dds <- DESeqDataSetFromMatrix(matrix_df, colData = info, ~ strain)

#PCA après avoir normalisé les données de comptage
vsd <- varianceStabilizingTransformation(dds)
dir.create("plots")#nouveau dossier afin de stocker les plots
png("plots/rplot_PCA.png", width = 700, height = 700)
plotPCA(vsd, intgroup="strain") + geom_text(aes(label=name),vjust=1) #plotPCA + nom des échantillons
dev.off()

#analyse de l'expression différentielle des gènes
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),] #Résultats de l'analyse ordonnées par p ajustée croissante

diff_expressed_genes = resOrdered[1:10,] #10 gènes différentiellement exprimés avec les plus petites p ajustées

#Nombre de gènes différentiellement exprimés
table(res$padj<0.05)

png("plots/rplot_plotCounts.png", width = 700, height = 700) #plotCounts du gène avec la p-ajustée la plus faible                                           
plotCounts(dds, gene=which.min(res$padj), intgroup="strain", main = "plotCounts du gène PHLDB2", col = 'blue')
dev.off()

png("plots/rplot_MA.png", width = 700, height = 700) #MA-plot 
plotMA(res, ylim=c(-5,5))
dev.off()
                                            
save.image("image.RData")



             
