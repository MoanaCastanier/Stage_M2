library(dplyr)
library(ggtree)
library(phytools)
library(ape)

# make rooted and binary phylogeny
prepa_tree <- function (tree, outgroup){
  # Raciner 
  tree<- root(tree, outgroup, resolve.root = TRUE)  
  
  # résoudre les polytomies
  tree <- multi2di(tree, random = TRUE)
  
  # Supprimer branche de longueur = 0
  tree$edge.length[which(tree$edge.length == 0)]<-0.00001
  tree <- ladderize (tree, T)
  
  print(is.rooted(tree))
  print(is.binary(tree))
  return(tree) 
}

Echantillonnage_apisum <- read.csv2("~/Documents/Echantillonnage_apisum_clear.csv", head=TRUE, sep=",")
phylo_buchnera <- read.tree("~/Documents/phylo_snippy/phylo_buchnera_without_outgroup.nwk")

phylo_buchnera <- prepa_tree(phylo_buchnera, c("wamizawa", "CS041_S97", "CS042_S98") )

vecteur_loc <- setNames(Echantillonnage_apisum$zone_geographique,
                        Echantillonnage_apisum [["read_file_name"]])
mod_ARD_loc <- fitMk(phylo_buchnera, vecteur_loc, model="ARD", fixedQ=NULL)
ancr (mod_ARD_loc, type = "joint")


vecteur_serratia <- setNames(Echantillonnage_apisum$presence_serratia,
                             Echantillonnage_apisum [["read_file_name"]])
mod_ARD_serratia <- fitMk(phylo_buchnera, vecteur_serratia, model="ARD", fixedQ=NULL)
ancr (mod_ARD_serratia, type ="joint")