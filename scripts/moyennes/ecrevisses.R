#reproduction des graphs de l'article sur les modeles avec les données brutes
# du jeu de données écrevisses

library(ape)
library(geiger)
library(ggplot2)


# arbre phylogénétique :

arbre_ecrevisses = read.tree("interspecific_rnaseq/data/crayfish.nodelabels.tre")
plot(arbre_ecrevisses)
arbre_ecrevisses$tip.label

data_ecrevisses = read.table("interspecific_rnaseq/data/orthogroups.TMM.EXPR.matrix",
                        header = TRUE,
                        row.names = 1,)

head(data_ecrevisses)
dim(data_ecrevisses) #3560 gènes , 34 individus

# étant donné que l'arbre ne contient que 14 espèces et qu'on en 34 dans les données, 
# je fais la moyenne des individus par espèce pour avoir 14 colonnes
#ecrevisses14 <- ecrevisses[, c("CCRYP_KC8782","CDUBI_KC8779", "CGRAY_KC8769","CHAMU_KC8582","CNERT_KC8482","CRUST_KC8774","CSETO_KC8659","CTENE_KC8583","OAUST_KC8577","OINCO_KC8597","PFALL_KC8497","PHORS_KC8784","PLUCI_KC8540","PPALL_KC8786")]

ecrevisses14 = data.frame(
  groupe1 = rowMeans(data_ecrevisses[, c("CCRYP_KC8782","CCRYP_KC8782")]),
  groupe2 = rowMeans(data_ecrevisses[, c("CDUBI_KC8779", "CDUBI_KC8780", "CDUBI_KC8781")]),
  groupe3 = rowMeans(data_ecrevisses[, c("CGRAY_KC8769", "CGRAY_KC8770","CGRAY_KC8771")]),
  groupe4 = rowMeans(data_ecrevisses[, c("CHAMU_KC8582", "CHAMU_KC8601","CHAMU_KC8602")]),
  groupe5 = rowMeans(data_ecrevisses[, c("CNERT_KC8482", "CNERT_KC8483","CNERT_KC8484")]),
  groupe6 = rowMeans(data_ecrevisses[, c("CRUST_KC8774", "CRUST_KC8775")]),
  groupe7 = rowMeans(data_ecrevisses[, c("CSETO_KC8659", "CSETO_KC8660","CSETO_KC8661")]),
  groupe8 = rowMeans(data_ecrevisses[, c("CTENE_KC8583", "CTENE_KC8585","CTENE_KC8587")]),
  groupe9 = rowMeans(data_ecrevisses[, c("OAUST_KC8577", "OAUST_KC8579","OAUST_KC8580")]),
  groupe10 = rowMeans(data_ecrevisses[, c("OINCO_KC8597", "OINCO_KC8598")]),
  groupe11 = rowMeans(data_ecrevisses[, c("PFALL_KC8497", "PFALL_KC8498","PFALL_KC8499")]),
  groupe12 = rowMeans(data_ecrevisses[, c("PHORS_KC8784", "PHORS_KC8785")]),
  groupe13 = rowMeans(data_ecrevisses[, c("PLUCI_KC8540", "PLUCI_KC8783","CTENE_KC8587")]),
  groupe14 = rowMeans(data_ecrevisses[, c("PPALL_KC8786","PPALL_KC8786")])
)

dim(ecrevisses14)

colnames(ecrevisses14) = c("CCRYP", "CDUBI","CGRAY","CHAMU","CNERT","CRUST","CSETO","CTENE","OAUST","OINCO","PFALL","PHORS","PLUCI","PPALL")
head(ecrevisses14)



fit_gene_models = function(x, arbre) { #fonction pour tester chaque modele au gène
  
  trait = as.numeric(x) #transfo en valeur numérique
  names(trait)= names(x)
  
  #fitcontinous est une fonction de geiger qui permet d'ajuster un modele sur les données
  fit_BM = try(fitContinuous(arbre, trait, model = "BM"), silent = TRUE) #coller le modele BM
  fit_OU= try(fitContinuous(arbre, trait, model = "OU"), silent = TRUE) # modele OU
  fit_EB = try(fitContinuous(arbre, trait, model = "EB"), silent = TRUE) #modele EB
  
  aics <- c( #récupération des aic
    BM = fit_BM$opt$aic,
    OU = fit_OU$opt$aic,
    EB = fit_EB$opt$aic
  )
  names(which.min(aics)) # chosiir celui qui minimise l'AIC
  
}


modeles30 = apply(ecrevisses14[1:30, ], 1, function(x) { # test sur les 30 premiers gènes déjà
  fit_gene_models(x, arbre_ecrevisses)
})

modeles30

meilleurmodele = as.data.frame(table(modeles30))
colnames(meilleurmodele) = c("modele", "nombre_de_genes")
meilleurmodele

ggplot(meilleurmodele, aes(x = modele, y = nombre_de_genes)) + # distribution des modèles pour les 30 gènes
  geom_bar(stat = "identity", fill = "#E5B80B") +
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")




#fig2b
# pour la figure 2b, ils utilisent le package arbutus pour faire les 5 tests
#mais pas de maj pour dernière version de R
# demain il faut essayer avec 4.2.3 pour voir si ça marche
library(geiger)
library(arbutus) # librairie arbutus pour pouvoir faire les 5 tests statstiques de l'article

##### test pour 1 gène
trait <- as.numeric(ecrevisses14[2, ])
names(trait) <- colnames(ecrevisses14)

trait <- trait[arbre$tip.label]

fit <- fitContinuous(arbre, trait, model = modeles30[2])
arb <- arbutus(fit)

arb #ok fonctionne
#####

library(ape)
library(geiger)
library(arbutus)

pvalues = data.frame() # on crée un tableau vide pour y mettre les pvalues des différents tests

for(i in 1:30){ #boucle pour les 30 premiers gènes
  
  trait = as.numeric(ecrevisses14[i, ]) # on convertit les valeurs du tableau en valeurs numérique
  names(trait)= colnames(ecrevisses14) # on associe chaque valeur à une espèce
  trait= trait[arbre_ecrevisses$tip.label] #meme noms pour les traits que celui de l'arbre
  
  modele = as.character(modeles30[i]) #on passe les modeles en carac
  
  fit = try( # on ajuste le modele evol 
    fitContinuous(arbre_ecrevisses, trait, model = modele),
    silent = TRUE
  )
  if(inherits(fit, "try-error")) next # si jamais ça échoue, on ignore et on passe au gène suivant
  
  arb= try( # calcul pvals
    arbutus(fit),
    silent = TRUE 
  )
  if(inherits(arb, "try-error")) next #idem 
  
  p= arb$p.values #on récupère les pvalues
  
  pvalues <- rbind( # et on les ajoute à notre data frame resultats
    pvalues,
    data.frame(
      gene = i,
      modele = modele,
      m.sig = p["m.sig"],
      c.var = p["c.var"],
      s.var = p["s.var"],
      s.asr = p["s.asr"],
      s.hgt = p["s.hgt"],
      d.cdf = p["d.cdf"]
    )
  )
}

pvalues

library(tidyr)

tableau = resultats |> pivot_longer( 
    cols = c(c.var, s.var, s.asr, s.hgt, d.cdf), # pas le m.sig dans l'article
    names_to = "test",
    values_to = "pvalue"
  )

ggplot(tableau, aes(x = pvalue)) +
  geom_histogram(bins = 30, fill = "#0097a7") + #intervalles = 30 
  facet_wrap(~test, ncol=3) + #sépare le graphique par type de test
  theme_classic() + #style graphique
  xlab("p-values") +
  ylab("nombre de gènes") +
  geom_vline(xintercept = 0.05, linetype = "dashed") # on ajoute la limite


# H0 : le modele peut expliquer les données 
# donc si pval > 0.05,
# on ne rejette pas H0 et le modele représente assez bien les données



# le modele coalescent considere le fait que les gènes n'évoluent pas forcemeent
# de la mm facon que les especes, donc leur arbre est différent pour chaque trait et 
# que chaque trait dépend de l'expression de plusieurs loci, ayant eux meme un arbre différent
# donc un trait depend de la somme d'effet de plusieurs loci et quand L tend vers l'infini,
# trait tend vers une loi normale et c'est parce que ça suit une loi normale qu'on peut estimer mu 
# et la matrice de cova meme si on connait pas l'expression de chaque loci 





#AJOUTER LE MODELE GC
#figure 2a
# dans l'article ils disent que GC peut être apparenté à un BM dont
# la taille des branches a été modifiée

#L’expression génique évolue comme un BM, mais sur un arbre modifié pour tenir 
#compte de la discordance gène/arbre d’espèces liée au coalescent.

transform_tree_coal= function(tree, lambda = 2) { # fonction pour transformer longueur branches (controle intensité avec lambda)
  tree2 = tree #copie de l'arbre pour ne pas le modifier directmeent
  
  profondeur = node.depth.edgelength(tree)# profondeur de chaque noeud depuis la racine
  #ex : racine = 0, feuilles = valeurs les plus grandes
  
  
  edges= tree$edge #on récupère la matrice de connexion entre parent/enfant
  edge.length = tree$edge.length #longueur de chaque branche
  new_lengths = numeric(length(edge.length)) #prepa nouveau vect
  
  for(i in seq_along(edge.length)){ # on boucle sur toutes les branches
    parent <- edges[i,1] #pour la branche i, parent = noeud de depart
    enfant  <- edges[i,2] #pour la branche i, enfant = noeuf d'arrivée
    
    l_u <- edge.length[i] #longueur originale de la branche i
    
    
    t_pa <- profondeur[parent] # temps jusqu'au parent, plus t_pa est grand plus on est éloigné dans le temps
    
    new_lengths[i] <- # transformation de l'article(proposition 11)
      l_u + # la nouvelle branche = longueur originale 
      (lambda - 1) * #+effet coalescent
      (1 - exp(-l_u)) * # + dépendence à la longueurde la branche
      exp(-t_pa) # + dépendance au temps : effet décroissant avec le temps (si proche racine = t petit donc effet fort et inversement)
  }
  
  tree2$edge.length <- new_lengths # on remplace les anciennes longueurs par les nouvelles
  tree2 # on a notre nouvel arbre
}

library(ape)

arbre_gc <- transform_tree_coal(arbre, lambda = 2) 

plot(arbre_gc)


library(phylolm)

trait = as.numeric(ecrevisses14[1, ])
names(trait) = colnames(ecrevisses14) 
tabl = data.frame(species = names(trait),trait = trait)
rownames(tabl) <- df$species #nouveau tableau

head(tabl)


# cconstruction d'une cova phylog,ajustement modele, estimation de moyenne evolutive et variance evol
fit_gc = phylolm( #pour faire le fit avec GC, ressemble à un BM avec des longueurs de branches différentes
  trait ~ 1, #trait = μ + erreur_phylogénétique, 1 représente intercept seul
  data = tabl,
  phy = arbre_gc,
  model = "BM"
)

fit_gc

fit_gene_modeles_avecGC <- function(x, arbre) { # meme fonction de précedemment mais en ajoutant le fit de GC
  
  arbre_gc <- transform_tree_coal(arbre, lambda = 2)
  
  trait <- as.numeric(x)
  names(trait) <- names(x)
  
  fit_BM <- try(fitContinuous(arbre, trait, model = "BM"), silent = TRUE) #coller le modele BM
  fit_OU <- try(fitContinuous(arbre, trait, model = "OU"), silent = TRUE) # modele OU
  fit_GC <- phylolm(trait ~ 1,data = df,phy = arbre_gc,model = "BM") #GC
  fit_EB <- try(fitContinuous(arbre, trait, model = "EB"), silent = TRUE) #EB
  
  aics <- c(
    BM = fit_BM$opt$aic,
    OU = fit_OU$opt$aic,
    EB = fit_EB$opt$aic,
    GC = AIC(fit_GC)
  )
  
  cat("\nAIC:\n")
  print(aics)
  
  names(which.min(aics)) 
}

modeles30GC <- apply(ecrevisses14[1:30, ], 1, function(x) {
  fit_gene_modeles_avecGC(x, arbre)
})

modeles30GC


meilleurmodeleGC <- as.data.frame(table(modeles30GC))
colnames(meilleurmodeleGC) <- c("modele", "nombre_de_genes")
meilleurmodeleGC

ggplot(meilleurmodeleGC, aes(x = modele, y = nombre_de_genes)) +
  geom_bar(stat = "identity", fill = "#E5B80B") +
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")





# figure 2a de l'article avec les tests
library(phylolm)


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

modeles30GCecr <- apply(ecrevisses14[1:30, ], 1, function(x) {
  fit_gene_modeles_avecGCbis(x, arbre_ecrevisses)
})


saveRDS(modeles100GCecr, file = "modeles100GCecr.rds")

meilleurmodeleGCbis <- as.data.frame(table(modeles30GCecr))

colnames(meilleurmodeleGCbis) <- c("modele", "nombre_de_genes")

meilleurmodeleGCbis$modele <- factor(
  meilleurmodeleGCbis$modele,
  levels = c("BM", "OU", "EB", "GC")
)

ggplot(meilleurmodeleGCbis,
       aes(x = modele, y = nombre_de_genes)) +
  geom_bar(stat = "identity", fill = "#E5B80B") +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")

plotecrevisse = ggplot(meilleurmodeleGCbis,
                      aes(x = modele, y = nombre_de_genes)) +
  geom_bar(stat = "identity", fill = "#E5B80B") +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  xlab("modèle évolutif") +
  ylab("nombre de gènes")

#ggsave("meilleurmodeleGCbis.jpeg", plot = plotecrevisse)






