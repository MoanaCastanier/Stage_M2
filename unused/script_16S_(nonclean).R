library(readxl)
library(dplyr)
library(tibble)
#library(seqinr)
#library(ggplot2)
#library(tidyr)
#library(ggtree)
#library(ggnewscale)
#library(RColorBrewer)
#library(patchwork)

ApisumDonnees16Straitees <- read_excel("ApisumDonnees16Straitees.xlsx")
ApisumDonnees16Straitees <- ApisumDonnees16Straitees [-44,]
ApisumDonnees16S2021 <- read.csv2("ApisumDonnees16S2021.csv", header = TRUE, sep=".")

####### Homogeneisation des noms : 

ApisumDonnees16Straitees <- ApisumDonnees16Straitees %>%
  mutate(across(where(is.character), ~gsub(";", "_", .)))
ApisumDonnees16S2021  <- ApisumDonnees16S2021 %>%
  mutate(across(where(is.character), ~gsub(";", "_", .)))

colnames(ApisumDonnees16Straitees) <- gsub("\\.", "_", colnames(ApisumDonnees16Straitees))
colnames(ApisumDonnees16S2021) <- gsub("\\.", "_", colnames(ApisumDonnees16S2021))

ApisumDonnees16Straitees <- ApisumDonnees16Straitees %>%
  column_to_rownames(var = "rdp_tax_and_bootstrap")
rownames(ApisumDonnees16Straitees) <- gsub("\\.", "_", rownames(ApisumDonnees16Straitees))

ApisumDonnees16Straitees <- ApisumDonnees16Straitees %>%
  mutate(across(everything(), ~gsub("\\.", ",", as.character(.))))
ApisumDonnees16S2021 <- ApisumDonnees16S2021 %>%
  mutate(across(everything(), ~gsub("\\.", ",", as.character(.))))


# Passer en pourcentage pour 2021 : 

# Fonction pour calculer des pourcentage (colonnes) : 
calculer_pourcentages <- function(your_dataframe, colonnes_individus) {
  
  # Convertir les colonnes en numeric si nécessaire
  for(col in colonnes_individus) {
    if(is.list(your_dataframe[[col]])) {
      your_dataframe[[col]] <- as.numeric(unlist(your_dataframe[[col]]))
    } else {
      your_dataframe[[col]] <- as.numeric(your_dataframe[[col]])
    }
  }
  
  # 1. Calculer read_total pour chaque colonne
  read_totals <- colSums(your_dataframe[, colonnes_individus], na.rm = TRUE)
  
  # 2. CORRECTION : Créer une nouvelle ligne remplie de NA
  ligne_read_total <- your_dataframe[1, ]
  ligne_read_total[1, ] <- NA  # Mettre NA partout
  
  # Remplir uniquement les colonnes individus avec les totaux
  for(col in colonnes_individus) {
    ligne_read_total[[col]] <- read_totals[col]
  }
  
  rownames(ligne_read_total) <- "read_total"
  
  # 3. Ajouter la ligne au dataframe
  your_dataframe_avec_total <- rbind(your_dataframe, ligne_read_total)
  
  # 4. Calculer les pourcentages
  your_dataframe_pourcentages <- your_dataframe_avec_total
  for(col in colonnes_individus) {
    total <- as.numeric(your_dataframe_avec_total["read_total", col])
    lignes_a_convertir <- rownames(your_dataframe_avec_total) != "read_total"
    
    your_dataframe_pourcentages[lignes_a_convertir, col] <- 
      (as.numeric(your_dataframe_avec_total[lignes_a_convertir, col]) / total) * 100
  }
  
  return(your_dataframe_pourcentages)
}
#Fonction, rajouter ligne = somme (%read < 0.5%), remplacer par 0 ensuite 
traiter_seq_non_considere <- function(your_dataframe, colonnes_individus, seuil = 0.5) {
  
  # 1. Identifier les lignes où TOUTES les valeurs < seuil
  lignes_normales <- !(rownames(your_dataframe) %in% c("read_total", "seq_non_considere"))
  
  lignes_inferieures_seuil <- apply(your_dataframe[lignes_normales, colonnes_individus], 1, 
                                    function(x) all(x < seuil, na.rm = TRUE))
  
  noms_lignes_a_sommer <- names(lignes_inferieures_seuil[lignes_inferieures_seuil])
  
  # 2. CORRECTION : Créer la ligne seq_non_considere remplie de NA
  ligne_seq_non_considere <- your_dataframe[1, ]
  ligne_seq_non_considere[1, ] <- NA  # Mettre NA partout
  
  # Remplir uniquement les colonnes individus
  ligne_seq_non_considere[, colonnes_individus] <- 
    colSums(your_dataframe[noms_lignes_a_sommer, colonnes_individus, drop = FALSE], na.rm = TRUE)
  
  rownames(ligne_seq_non_considere) <- "seq_non_considere"
  
  your_dataframe_avec_seq <- rbind(your_dataframe, ligne_seq_non_considere)
  
  # 3. Remplacer par 0 les valeurs < seuil (sauf lignes spéciales)
  lignes_a_modifier <- !(rownames(your_dataframe_avec_seq) %in% c("read_total", "seq_non_considere"))
  
  for(col in colonnes_individus) {
    your_dataframe_avec_seq[lignes_a_modifier, col] <- 
      ifelse(your_dataframe_avec_seq[lignes_a_modifier, col] < seuil, 
             0, 
             your_dataframe_avec_seq[lignes_a_modifier, col])
  }
  
  return(your_dataframe_avec_seq)
}
#Fonction, colonne somme ligne avec possibilité de supprimer si =0
ajouter_total_blast <- function(your_dataframe, colonnes_a_sommer, supprimer_zeros = TRUE) {
  
  # 1. Créer la colonne total_blast (somme des colonnes spécifiées)
  your_dataframe$total_blast <- rowSums(your_dataframe[, colonnes_a_sommer], na.rm = TRUE)
  
  # 2. Optionnel : supprimer les lignes où total_blast = 0
  if(supprimer_zeros) {
    # Toujours garder les lignes spéciales
    lignes_speciales <- rownames(your_dataframe) %in% c("read_total", "seq_non_considere")
    lignes_a_garder <- (your_dataframe$total_blast != 0) | lignes_speciales
    
    your_dataframe <- your_dataframe[lignes_a_garder, ]
  }
  
  return(your_dataframe)
}


# Convertir le tibble en data.frame
ApisumDonnees16S2021 <- as.data.frame(ApisumDonnees16S2021)
colonnes_acyp <- grep("ACYP", colnames(ApisumDonnees16S2021), value = TRUE)
# Convertir les colonnes ACYP en numeric
for(col in colonnes_acyp) {
  ApisumDonnees16S2021[[col]] <- as.numeric(unlist(ApisumDonnees16S2021[[col]]))
}

# Appliquer les fonctions corrigées
ApisumDonnees16S2021 <- calculer_pourcentages(ApisumDonnees16S2021, colonnes_acyp)
ApisumDonnees16S2021 <- traiter_seq_non_considere(ApisumDonnees16S2021, colonnes_acyp, seuil = 0.5)
ApisumDonnees16S2021 <- ajouter_total_blast(ApisumDonnees16S2021, colonnes_acyp, supprimer_zeros = TRUE)


####### Identifier les individus avec > 1000 reads + créer tableau associé :

#Apisum16S 2023
colonnes_individus <- grep("^ACYP", colnames(ApisumDonnees16Straitees), value = TRUE)
colonnes_info_environnement <- setdiff(colnames(ApisumDonnees16Straitees), colonnes_individus)

valeur_read_total <- as.numeric(ApisumDonnees16Straitees["read_total", colonnes_individus])
names(valeur_read_total ) <- colonnes_individus

# Séparer avec seuil = 1000 reads
individus_sup_1000 <- names(valeur_read_total [!is.na(valeur_read_total ) & valeur_read_total  >= 1000])
individus_inf_1000 <- names(valeur_read_total [!is.na(valeur_read_total ) & valeur_read_total  < 1000])

Apisum_Donnee16S_sup_1000_reads <- ApisumDonnees16Straitees[, c(colonnes_info_environnement, individus_sup_1000)]
Apisum_Donnee16S_inf_1000_reads <- ApisumDonnees16Straitees[, c(colonnes_info_environnement, individus_inf_1000)]

#Faire les moyenne des individus
moyenner_replicats <- function(df, colonnes_individus) {
  
  # Identifier les paires de colonnes (colonne normale et son _BIS)
  colonnes_sans_bis <- colonnes_individus[!grepl("_BIS", colonnes_individus)]
  colonnes_bis <- colonnes_individus[grepl("_BIS", colonnes_individus)]
  
  # Créer un nouveau dataframe avec les colonnes non-individus
  colonnes_non_individus <- setdiff(colnames(df), colonnes_individus)
  df_moyenne <- df[, colonnes_non_individus, drop = FALSE]
  
  # Pour chaque colonne sans _BIS
  for(col in colonnes_sans_bis) {
    # Chercher si un _BIS existe
    col_bis <- paste0(col, "_BIS")
    
    if(col_bis %in% colonnes_bis) {
      # Si _BIS existe, faire la moyenne
      nom_nouvelle_col <- paste0(col, "_moyenne")
      df_moyenne[[nom_nouvelle_col]] <- (as.numeric(df[[col]]) + as.numeric(df[[col_bis]])) / 2
    } else {
      # Si pas de _BIS, garder la colonne telle quelle avec _moyenne
      nom_nouvelle_col <- paste0(col)
      df_moyenne[[nom_nouvelle_col]] <- as.numeric(df[[col]])
    }
  }
  
  return(df_moyenne)
}

# UTILISATION :
ApisumDonnees16S2021<- moyenner_replicats(ApisumDonnees16S2021, colonnes_acyp)

#Fonction : 
# Fonction pour ajouter une ligne de présence/absence d'une bactérie
ajouter_presence_bacterie <- function(your_dataframe, bacterie, colonnes_individus, colonne_recherche = "blast_taxonomy") {
  
  # Vérifier que la colonne de recherche existe
  if (!colonne_recherche %in% names(your_dataframe)) {
    stop(paste("La colonne", colonne_recherche, "n'existe pas dans les données"))
  }
  
  # Créer le nom de la nouvelle ligne
  nom_ligne <- paste0("presence_", tolower(bacterie))
  
  # Créer une nouvelle ligne avec la même structure que le dataframe
  nouvelle_ligne <- your_dataframe[1, , drop = FALSE]
  nouvelle_ligne[1, ] <- NA
  rownames(nouvelle_ligne) <- nom_ligne
  
  # Pour chaque colonne individu
  for(col in colonnes_individus) {
    # Vérifier si la bactérie est présente dans au moins une ligne pour cet individu
    presence <- any(
      !is.na(your_dataframe[[col]]) & 
        as.numeric(your_dataframe[[col]]) > 0 & 
        !is.na(your_dataframe[[colonne_recherche]]) &
        grepl(bacterie, your_dataframe[[colonne_recherche]], ignore.case = TRUE),
      na.rm = TRUE
    )
    
    nouvelle_ligne[[col]] <- ifelse(presence, "présent", "absent")
  }
  
  # Ajouter la ligne au dataframe
  your_dataframe <- rbind(your_dataframe, nouvelle_ligne)
  
  # Afficher un résumé
  presents <- sum(nouvelle_ligne[, colonnes_individus] == "présent", na.rm = TRUE)
  absents <- sum(nouvelle_ligne[, colonnes_individus] == "absent", na.rm = TRUE)
  
  cat("✓ Ligne", nom_ligne, "créée\n")
  cat("  Recherche dans la colonne:", colonne_recherche, "\n")
  cat("  Présent dans", presents, "individus\n")
  cat("  Absent dans", absents, "individus\n")
  
  return(your_dataframe)
}

#####ajout de regroupe meme individu 

colonnes_individus <- grep("^ACYP", colnames(Apisum_Donnee16S_sup_1000_reads), value = TRUE)
Apisum_Donnee16S_sup_1000_reads <- ajouter_presence_bacterie(
  Apisum_Donnee16S_sup_1000_reads, 
  "Serratia",
  colonnes_individus = colonnes_individus,
  colonne_recherche = "blast_taxonomy"
)

Apisum_Donnee16S_sup_1000_reads <- ajouter_presence_bacterie(
  Apisum_Donnee16S_sup_1000_reads, 
  "Buchnera",
  colonnes_individus = colonnes_individus,
  colonne_recherche = "blast_taxonomy"
)

colonnes_acyp <- grep("ACYP", colnames(ApisumDonnees16S2021), value = TRUE)
ApisumDonnees16S2021 <- ajouter_presence_bacterie(
  ApisumDonnees16S2021, 
  "Serratia",
  colonnes_individus = colonnes_acyp,
  colonne_recherche = "blast_taxonomy"
)

ApisumDonnees16S2021 <- ajouter_presence_bacterie(
  ApisumDonnees16S2021, 
  "Buchnera",
  colonnes_individus = colonnes_acyp,
  colonne_recherche = "blast_taxonomy"
)

write.table(Apisum_Donnee16S_sup_1000_reads, "Apisum_Donnee16S_sup_1000_reads.csv", 
            sep = ".", dec = ",", row.names = TRUE, quote = FALSE)
write.table(ApisumDonnees16S2021, "ApisumDonnees16S2021.csv", 
            sep = ".", dec = ",", row.names = TRUE, quote = FALSE)