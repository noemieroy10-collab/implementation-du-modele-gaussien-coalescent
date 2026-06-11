#JEU DE DONNÉES POISSONS


library(ape)
library(geiger)
library(ggplot2)
library(arbutus)
library(tidyr)
library(phylolm)


# arbre phylogénétique :

arbre_poissons = read.tree("fishes/data/recodedTreeNamed.tre")
plot(arbre_poissons)


arbre_poissons$tip.label

data_poissons = read.csv("fishes/data/master_fpkm.csv",
                         header = TRUE,
                         row.names = 1,)
head(data_poissons)
dim(data_poissons) #16740 gènes , 118 individus
colnames(data_poissons)

#moyenne ou pas ?

poissons20 = data.frame(
  GholNS = rowMeans(data_poissons[, 1:6]),
  GsexNS = rowMeans(data_poissons[, 7:12]),
  PbimNS = rowMeans(data_poissons[, 13:18]),
  XhelNS = rowMeans(data_poissons[, 19:24]),
  LperNS = rowMeans(data_poissons[, 25:30]),
  PlatNS = rowMeans(data_poissons[, 31:36]),
  PlimNS = rowMeans(data_poissons[, 37:42]),
  PmexPuyNS = rowMeans(data_poissons[, 43:48]),
  PmexPichNS = rowMeans(data_poissons[, 49:54]),
  PmexTacNS = rowMeans(data_poissons[, 55:60]),
  PmexPuyS = rowMeans(data_poissons[, 61:65]),
  PmexTacS = rowMeans(data_poissons[, 66:70]),
  PmexPichS = rowMeans(data_poissons[, 71:76]),
  PlatS = rowMeans(data_poissons[, 77:82]),
  LsulS = rowMeans(data_poissons[, 83:88]),
  GsexS = rowMeans(data_poissons[, 89:94]),
  GeurS = rowMeans(data_poissons[, 95:100]),
  GholS = rowMeans(data_poissons[, 101:106]),
  PbimS = rowMeans(data_poissons[, 107:112]),
  XhelS = rowMeans(data_poissons[, 113:118])
)

head(poissons20)



fit_gene_modeles_avecGCbis <- function(x, arbre) {
  
  trait <- as.numeric(x)
  names(trait) <- names(x)
  
  fit_BM = try(fitContinuous(arbre, trait, model = "BM"), silent = TRUE)
  fit_OU = try(fitContinuous(arbre, trait, model = "OU"), silent = TRUE)
  fit_EB = try(fitContinuous(arbre, trait, model = "EB"), silent = TRUE)
  fit_GC = try(phylolm(trait ~ 1, phy = arbre, model = "GC"), silent = TRUE)
  
  aic_BM = fit_BM$opt$aic 
  aic_OU = fit_OU$opt$aic 
  aic_EB = fit_EB$opt$aic 
  aic_GC = AIC(fit_GC) 
  
  aics = c(BM = aic_BM, OU = aic_OU, EB = aic_EB, GC = aic_GC)
  
  cat("\nAIC:\n")
  print(aics)
  
  names(which.min(aics))
}

modeles100pois <- apply(poissons20[1:100, ], 1, function(x) { #sur les 100 premiers genes
  fit_gene_modeles_avecGCbis(x, arbre_poissons)
})


#saveRDS(modelesp100GCpois, file = "modeles100GCpois.rds")


meilleurmodele <- as.data.frame(table(modeles100pois))

colnames(meilleurmodele) <- c("modele", "nombre_de_genes")

meilleurmodele$modele <- factor(
  meilleurmodelepGCbis$modele,
  levels = c("BM", "OU", "EB", "GC")
)

ggplot(meilleurmodele,
       aes(x = modele, y = nombre_de_genes)) +
  geom_bar(stat = "identity", fill = "#E5B80B") +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")





