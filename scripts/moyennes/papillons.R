#JEU DE DONNEES PAPILLONS

library(ape)
library(geiger)
library(ggplot2)
library(arbutus)
library(tidyr)
library(phylolm)


# arbre phylogénétique :

arbre_papillons = read.nexus("Heliconius_Butterflies/Data/Heliconiini.multiple_uniform_constraints.MCC.tree")
plot(arbre_papillons)

library(ape)




# garder uniquement certaines espèces
sous_arbre <- keep.tip(arbre_papillons, c("Heliconius_charithonia_8830_P", "Heliconius_sara_8862_P", "Heliconius_erato_erato_NCS2556_FG","Heliconius_doris_02_1939_Pe","Heliconius_melpomene_rosina_546_P"))

plot(sous_arbre)


sous_arbre$tip.label
sous_arbre$tip.label <- c("Heliconius_charithonia", "Heliconius_doris", "Heliconius_rato","Heliconius_melpomene","Heliconius_sara")
plot(sous_arbre)

data_papillons = read.csv("Heliconius_Butterflies/Data/expression_matrix.csv",
                         header = TRUE,
                         row.names = 1,)

head(data_papillons)
dim(data_papillons) #2393 gènes 49 individus
colnames(data_papillons)


papillons5= data.frame(
  Heliconius_charithonia = rowMeans(data_papillons[, 1:12]),
  Heliconius_doris = rowMeans(data_papillons[, 13:24]),
  Heliconius_rato = rowMeans(data_papillons[, 25:30]),
  Heliconius_melpomene = rowMeans(data_papillons[, 31:38]),
  Heliconius_sara = rowMeans(data_papillons[, 39:49])
)

head(papillons5)



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

modeles100pap <- apply(papillons5[1:100, ], 1, function(x) {
  fit_gene_modeles_avecGCbis(x, sous_arbre)
})

modeles100pap


meilleurmodele <- as.data.frame(table(modeles100pap))

colnames(meilleurmodele) <- c("modele", "nombre_de_genes")

meilleurmodele$modele <- factor(
  meilleurmodele$modele,
  levels = c("BM", "OU", "EB", "GC")
)

ggplot(meilleurmodele,
       aes(x = modele, y = nombre_de_genes)) +
  geom_bar(stat = "identity", fill = "#E5B80B") +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")


