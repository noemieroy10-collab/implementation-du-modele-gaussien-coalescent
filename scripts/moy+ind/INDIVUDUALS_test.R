#ECREVISSES INDIVIDUALS

library(ape)
library(geiger)
library(ggplot2)
library(phylolm)
library(tidyverse)

# arbre phylogénétique :

arbre_ecrevisses = read.tree("interspecific_rnaseq/data/crayfish.nodelabels.tre")
plot(arbre_ecrevisses)
arbre_ecrevisses$tip.label

data_ecrevisses = read.table("interspecific_rnaseq/data/orthogroups.TMM.EXPR.matrix",
                             header = TRUE,
                             row.names = 1,)
head(data_ecrevisses)
dim(data_ecrevisses) #3560 gènes , 34 individus


#on transpose le tableau car phylolm veut les individus en ligne et gènes en colonnes
data_ecrevissesT <- as.data.frame(t(data_ecrevisses))

#On crée la colonne espèce en extrayant le début du nom de l'individu
# (CCRYP_KC8782 devient CCRYP)
data_ecrevissesT$espece = sub("_.*", "", rownames(data_ecrevissesT))

#phylolm a besoin d'une colonne qui correspond
# aux tip.labels de l'arbre pour faire le lien avec les individus.
data_ecrevissesT $tip_label <- data_ecrevissesT $espece

head(data_ecrevissesT[, c("espece", "tip_label")]) # Pour vérifier


#data_ecrevissesT$individu <- rownames(data_ecrevissesT)



# %>% permet de repndre le resultat de la ligne précédente et l'envoie dans la fonction qui suit

stats_ecrevissesT <- data_ecrevissesT%>%# nouveau tableau avec moyenne et ecart type pour inclure la variance entre individus
  pivot_longer(cols = starts_with("OG"),names_to = "gene",values_to = "expression") %>%
  group_by(gene, espece) %>%
  summarise(Moyenne = mean(expression), #moyenne pour chaque trait par espece
    N = n(), #compter le nbre de ligne par espèce
    SD = sd(expression), #calcul ecart type entre les individus
    SE = if_else(is.na(SD) | N <= 1, 0, SD / sqrt(N)), #SE
    .groups = "drop"
  )
 
head(stats_ecrevissesT) # il y a NA quand il y a un seul ind pour 1 espece
stats_ecrevissesT$SE

#le script résume les individus par (moyenne,variance)

fit_gene_modeles_individus <- function(nom_gene, df_statistiques, arbre) {
  
  
  #on recrée un tableau à chaque itération pour chaque gène
  #on récupère la moyenne pour chaque espèce et son SE pour prendre en compte la variabilité entre inds
  donnees_gene =df_statistiques[df_statistiques$gene == nom_gene, ] #data frame classique
  donnees_gene = as.data.frame(donnees_gene)
  rownames(donnees_gene) = donnees_gene$espece
  
  #alignement strict sur l'ordre des feuilles de l'arbre
  df_local = donnees_gene[arbre$tip.label, ]
  
  # Création de la colonne de variance exigée par phylolm (se2)
  df_local$input_error = df_local$SE^2
  df_local$input_error[df_local$input_error == 0] <- 1e-5
  #on peut avoir un SE = 0 si un seul individu pour une espèce mais
  #si zero phylolm pas content donc on met valeur faible pour ne pas trop modifier
  

  fit_BM = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "BM", 
                        input_error = df_local$input_error, REML = TRUE), silent = TRUE)
  
  fit_OU = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "OUfixedRoot", 
                        input_error =df_local$input_error, REML = TRUE), silent = TRUE)
  
  fit_EB = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "EB", 
                        input_error = df_local$input_error, REML = TRUE), silent = TRUE)
  
  fit_GC = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "GC", 
                        input_error = df_local$input_error, REML = TRUE), silent = TRUE)
  

  aic_BM = if(inherits(fit_BM, "try-error")) NA else AIC(fit_BM)
  aic_OU = if(inherits(fit_OU, "try-error")) NA else AIC(fit_OU)
  aic_EB = if(inherits(fit_EB, "try-error")) NA else AIC(fit_EB)
  aic_GC = if(inherits(fit_GC, "try-error")) NA else AIC(fit_GC)
  
  aics =c(BM = aic_BM, OU = aic_OU, EB = aic_EB, GC = aic_GC)
  
  cat("\n", nom_gene, ":\n")
  print(aics) 
  
  if(all(is.na(aics))) return(NA)
  names(which.min(aics))
}


modeles30_ind = sapply(unique(stats_ecrevissesT$gene)[1:30], function(g) {
  fit_gene_modeles_individus(g, stats_ecrevissesT, arbre_ecrevisses)
})



# Voir le résultat final
table(modeles30_ind)



# 1. On transforme le vecteur en data.frame de comptage
df_visu= as.data.frame(table(modeles30_ind))
colnames(df_visu) = c("modele", "nombre_de_genes")

# 2. On force les niveaux du facteur pour inclure tous les modèles, même ceux à 0
df_visu$modele =factor(
  df_visu$modele,
  levels = c("BM", "OU", "EB", "GC")
)

# 3. Ton code ggplot mis à jour avec le bon tableau
ggplot(df_visu, aes(x = modele, y = nombre_de_genes)) +
  geom_bar(stat = "identity", fill = "#E5B80B") +
  scale_x_discrete(drop = FALSE) +  # Garde les barres à 0 visibles
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")








#POISSONS INDIVIDUALS

library(ape)
library(geiger)
library(ggplot2)
library(phylolm)
library(tidyverse)


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


#on transpose le tableau car phylolm veut les individus en ligne et gènes en colonnes
data_poissonsT <- as.data.frame(t(data_poissons))

#On crée la colonne espèce en extrayant le début du nom de l'individu
# (CCRYP_KC8782 devient CCRYP)
data_poissonsT$espece = sub("_.*", "", rownames(data_poissonsT))

#unique(data_poissonsT$espece)
#arbre_poissons$tip.label
#setdiff(unique(data_poissonsT$espece), arbre_poissons$tip.label)

#phylolm a besoin d'une colonne qui correspond
# aux tip.labels de l'arbre pour faire le lien avec les individus.
data_poissonsT$tip_label <- data_poissonsT$espece

head(data_poissonsT[, c("espece", "tip_label")]) # Pour vérifier


#data_ecrevissesT$individu <- rownames(data_ecrevissesT)



# %>% permet de repndre le resultat de la ligne précédente et l'envoie dans la fonction qui suit

stats_poissonsT= data_poissonsT%>%# nouveau tableau avec moyenne et ecart type pour inclure la variance entre individus
  pivot_longer(cols = starts_with("gene"),names_to = "gene",values_to = "expression") %>%
  group_by(gene, espece) %>%
  summarise(Moyenne = mean(expression), #moyenne pour chaque trait par espece
            N = n(), #compter le nbre de ligne par espèce
            SD = sd(expression), #calcul ecart type entre les individus
            SE = if_else(is.na(SD) | N <= 1, 0, SD / sqrt(N)), #SE
            .groups = "drop"
  )

head(stats_poissonsT) # il y a NA quand il y a un seul ind pour 1 espece
stats_ecrevissesT$SE

# le script résume les individus par (moyenne,variance)

fit_gene_modeles_individus <- function(nom_gene, df_statistiques, arbre) {
  
  
  #on recrée un tableau à chaque itération pour chaque gène
  #on récupère la moyenne pour chaque espèce et son SE pour prendre en compte la variabilité entre inds
  donnees_gene =df_statistiques[df_statistiques$gene == nom_gene, ] #data frame classique
  donnees_gene = as.data.frame(donnees_gene)
  rownames(donnees_gene) = donnees_gene$espece
  
  #alignement strict sur l'ordre des feuilles de l'arbre
  df_local = donnees_gene[arbre$tip.label, ]
  
  # Création de la colonne de variance exigée par phylolm (se2)
  df_local$input_error = df_local$SE^2
  df_local$input_error[df_local$input_error == 0] <- 1e-5
  #on peut avoir un SE = 0 si un seul individu pour une espèce mais
  #si zero phylolm pas content donc on met valeur faible pour ne pas trop modifier
  
  
  fit_BM = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "BM", 
                       input_error = df_local$input_error, REML = TRUE), silent = TRUE)
  
  fit_OU = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "OUfixedRoot", 
                       input_error =df_local$input_error, REML = TRUE), silent = TRUE)
  
  fit_EB = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "EB", 
                       input_error = df_local$input_error, REML = TRUE), silent = TRUE)
  
  fit_GC = try(phylolm(Moyenne ~ 1, phy = arbre, data = df_local, model = "GC", 
                       input_error = df_local$input_error, REML = TRUE), silent = TRUE)
  
  
  aic_BM = if(inherits(fit_BM, "try-error")) NA else AIC(fit_BM)
  aic_OU = if(inherits(fit_OU, "try-error")) NA else AIC(fit_OU)
  aic_EB = if(inherits(fit_EB, "try-error")) NA else AIC(fit_EB)
  aic_GC = if(inherits(fit_GC, "try-error")) NA else AIC(fit_GC)
  
  aics =c(BM = aic_BM, OU = aic_OU, EB = aic_EB, GC = aic_GC)
  
  cat("\n", nom_gene, ":\n")
  print(aics) 
  
  if(all(is.na(aics))) return(NA)
  names(which.min(aics))
}


modeles30pois_ind = sapply(unique(stats_poissonsT$gene)[1:30], function(g) {
  fit_gene_modeles_individus(g, stats_poissonsT, arbre_poissons)
})



# Voir le résultat final
table(modeles30pois_ind)



# 1. On transforme le vecteur en data.frame de comptage
df_visu= as.data.frame(table(modeles30pois_ind))
colnames(df_visu) = c("modele", "nombre_de_genes")

# 2. On force les niveaux du facteur pour inclure tous les modèles, même ceux à 0
df_visu$modele =factor(
  df_visu$modele,
  levels = c("BM", "OU", "EB", "GC")
)

# 3. Ton code ggplot mis à jour avec le bon tableau
ggplot(df_visu, aes(x = modele, y = nombre_de_genes)) +
  geom_bar(stat = "identity", fill = "#E5B80B") +
  scale_x_discrete(drop = FALSE) +  # Garde les barres à 0 visibles
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")
