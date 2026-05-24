#library
library(readxl)
library(dplyr)
library(tibble)
library(seqinr)
library(ggplot2)
library(tidyr)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(patchwork)
library(ggtreeExtra)
library(phytools)
library(nnet)
library(ape)
library(vegan)
library (nlme)
library(lmtest)
library(geiger)
library(rrnni)

options(scipen = 999)

###### Fonction #########

# Fonction : Calcul proportion de nucléotides non-ACGT
calcul_pourcentage_non_acgt <- function(sequences) {
  # Concaténer toutes les séquences en un seul vecteur
  tous_nucleotides <- unlist(sequences)
  
  # Nombre total de nucléotides
  total <- length(tous_nucleotides)
  
  # Compter les nucléotides qui ne sont pas a, c, g, t
  non_standard <- sum(!(tous_nucleotides %in% c("a", "c", "g", "t")))
  
  # Calculer la proportion
  proportion <- non_standard / total
  
  # Retourner les résultats
  return(list(
    non_standard = non_standard,
    total = total,
    proportion = proportion,
    pourcentage = proportion * 100
  ))
}

# Fonction : Analyse qualités des séquences pour tous les fichier fasta d'un dossier
analyser_qualite_fasta <- function(dossier, 
                                   pattern = "\\.fasta$") {
  
  # Vérifier que le dossier existe
  if (!dir.exists(dossier)) {
    stop("Le dossier spécifié n'existe pas : ", dossier)
  }
  
  # Lister tous les fichiers .fasta dans le dossier
  fichiers_fasta <- list.files(path = dossier, 
                               pattern = pattern, 
                               full.names = TRUE)
  
  # Vérifier qu'il y a des fichiers
  if (length(fichiers_fasta) == 0) {
    warning("Aucun fichier FASTA trouvé dans le dossier : ", dossier)
    return(NULL)
  }
  
  cat("Nombre de fichiers trouvés :", length(fichiers_fasta), "\n")
  
  # Créer un data frame vide pour stocker les résultats
  resultats_tableau <- data.frame(
    nom_sequence = character(),
    non_acgt = numeric(),
    total_nucleotides = numeric(),
    pourcentage_non_acgt = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Boucle sur tous les fichiers FASTA
  for (i in seq_along(fichiers_fasta)) {
    fichier <- fichiers_fasta[i]
    cat("Traitement du fichier", i, "/", length(fichiers_fasta), ":", basename(fichier), "\n")
    
    # Lire le fichier FASTA
    sequences <- read.fasta(fichier)
    
    # Calculer les proportions
    resultats <- calcul_pourcentage_non_acgt(sequences)
    
    # Extraire le nom du fichier (sans le chemin)
    nom_fichier <- basename(fichier)
    
    # Ajouter une ligne au tableau
    resultats_tableau <- rbind(resultats_tableau, 
                               data.frame(
                                 nom_sequence = nom_fichier,
                                 non_acgt = resultats$non_standard,
                                 total_nucleotides = resultats$total,
                                 pourcentage_non_acgt = resultats$pourcentage
                               ))
  }
  
  return(resultats_tableau)
}

analyser_qualite_fa <- function(dossier, 
                                pattern = "\\.fa$") {
  
  # Vérifier que le dossier existe
  if (!dir.exists(dossier)) {
    stop("Le dossier spécifié n'existe pas : ", dossier)
  }
  
  # Lister tous les fichiers .fasta dans le dossier
  fichiers_fasta <- list.files(path = dossier, 
                               pattern = pattern, 
                               full.names = TRUE)
  
  # Vérifier qu'il y a des fichiers
  if (length(fichiers_fasta) == 0) {
    warning("Aucun fichier FASTA trouvé dans le dossier : ", dossier)
    return(NULL)
  }
  
  cat("Nombre de fichiers trouvés :", length(fichiers_fasta), "\n")
  
  # Créer un data frame vide pour stocker les résultats
  resultats_tableau <- data.frame(
    nom_sequence = character(),
    non_acgt = numeric(),
    total_nucleotides = numeric(),
    pourcentage_non_acgt = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Vecteur pour stocker les fichiers vides ou problématiques
  fichiers_vides <- character()
  
  # Boucle sur tous les fichiers FASTA
  for (i in seq_along(fichiers_fasta)) {
    fichier <- fichiers_fasta[i]
    cat("Traitement du fichier", i, "/", length(fichiers_fasta), ":", basename(fichier), "\n")
    
    # Utiliser tryCatch pour gérer les erreurs
    resultats <- tryCatch({
      # Lire le fichier FASTA
      sequences <- read.fasta(fichier)
      
      # Vérifier si le fichier est vide
      if (is.null(sequences) || length(sequences) == 0) {
        fichiers_vides <- c(fichiers_vides, basename(fichier))
        cat("  -> Fichier vide ou sans séquence, passage au suivant\n")
        NULL
      } else {
        # Calculer les proportions
        calcul_pourcentage_non_acgt(sequences)
      }
    }, error = function(e) {
      # En cas d'erreur, ajouter à la liste des fichiers problématiques
      fichiers_vides <- c(fichiers_vides, basename(fichier))
      cat("  -> Erreur lors de la lecture :", e$message, "\n")
      NULL
    })
    
    # Si des résultats ont été obtenus, les ajouter au tableau
    if (!is.null(resultats)) {
      nom_fichier <- basename(fichier)
      resultats_tableau <- rbind(resultats_tableau, 
                                 data.frame(
                                   nom_sequence = nom_fichier,
                                   non_acgt = resultats$non_standard,
                                   total_nucleotides = resultats$total,
                                   pourcentage_non_acgt = resultats$pourcentage
                                 ))
    }
  }
  
  # Afficher le résumé des fichiers vides
  cat("\n=== RÉSUMÉ ===\n")
  cat("Fichiers traités avec succès :", nrow(resultats_tableau), "\n")
  cat("Fichiers vides ou problématiques :", length(fichiers_vides), "\n")
  
  if (length(fichiers_vides) > 0) {
    cat("\nListe des fichiers vides/problématiques :\n")
    for (f in fichiers_vides) {
      cat("  -", f, "\n")
    }
  }
  
  # Retourner les résultats avec la liste des fichiers vides en attribut
  attr(resultats_tableau, "fichiers_vides") <- fichiers_vides
  
  return(resultats_tableau)
}

# Fonction : fusion des colonnes avec correspondance partielle (trouver col_nom_cible dans col_nom_source malgré ajout)
ajout_colonnes <- function(df_cible, 
                           df_source, 
                           col_nom_cible, 
                           col_nom_source, 
                           col_a_copier,
                           nouveau_nom_col = col_a_copier) {
  
  # Extraire les vecteurs de noms
  noms_cible <- df_cible[[col_nom_cible]]
  noms_source <- df_source[[col_nom_source]]
  valeurs_source <- df_source[[col_a_copier]]
  
  # Suivre les indices source utilisés
  indices_utilises <- c()
  
  # Créer la nouvelle colonne avec correspondance partielle
  df_cible[[nouveau_nom_col]] <- sapply(noms_cible, function(nom) {
    # Chercher ce nom dans les noms sources
    correspondance <- grep(nom, noms_source, fixed = TRUE)
    
    if (length(correspondance) > 0) {
      indices_utilises <<- c(indices_utilises, correspondance[1])
      return(valeurs_source[correspondance[1]])
    } else {
      return(NA)
    }
  })
  
  # Identifier les lignes non appariées
  indices_non_utilises <- setdiff(seq_len(nrow(df_source)), indices_utilises)
  noms_non_apparier <- noms_source[indices_non_utilises]
  
  # Afficher les noms non appariés
  if (length(noms_non_apparier) > 0) {
    cat("\n=== Colonnes non appariées de", col_a_copier, "===\n")
    cat(paste(noms_non_apparier, collapse = "\n"))
    cat("\n\nTotal:", length(noms_non_apparier), "sur", length(noms_source), "\n\n")
  } else {
    cat("\n=== Toutes les colonnes de", col_a_copier, "ont été appariées ===\n\n")
  }
  
  return(df_cible)
}

# Fonction : calcul des fréquences d'une colonne, filtre possible
calculer_frequences <- function(your_dataframe, your_col, filtre_col = NULL, filtre_valeurs = NULL) {
  
  # Vérifier que la colonne existe
  if (!your_col %in% names(your_dataframe)) {
    stop(paste("La colonne", your_col, "n'existe pas dans les données"))
  }
  
  # Vérifier que la colonne type_sequence existe
  if (!"type_sequence" %in% names(your_dataframe)) {
    stop("La colonne 'type_sequence' n'existe pas dans les données")
  }
  
  # Appliquer le filtre si spécifié
  if (!is.null(filtre_col) && !is.null(filtre_valeurs)) {
    if (!filtre_col %in% names(your_dataframe)) {
      stop(paste("La colonne de filtre", filtre_col, "n'existe pas dans les données"))
    }
    your_dataframe <- your_dataframe %>%
      filter(.data[[filtre_col]] %in% filtre_valeurs)
    
    cat("✓ Filtre appliqué:", filtre_col, "=", paste(filtre_valeurs, collapse = ", "), "\n")
    cat("  Nombre de lignes après filtre:", nrow(your_dataframe), "\n")
  }
  
  # Calculer les fréquences
  frequences <- your_dataframe %>%
    group_by(across(all_of(your_col)), type_sequence) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(your_col))) %>%
    mutate(
      total = sum(n),
      frequence = n / total,
      pourcentage = frequence * 100
    ) %>%
    ungroup()
  
  return(frequences)
}

# Fonction : Creation d'histogramme basé sur une colonne, filtre possible
creer_histogrammes <- function(your_dataframe, 
                               your_col, 
                               palette_fixe, 
                               titre_principal = NULL, 
                               sauvegarder = FALSE, filtre_col = NULL, filtre_valeurs = NULL) {
  
  # Appliquer le filtre si spécifié
  if (!is.null(filtre_col) && !is.null(filtre_valeurs)) {
    if (!filtre_col %in% names(your_dataframe)) {
      stop(paste("La colonne de filtre", filtre_col, "n'existe pas dans les données"))
    }
    your_dataframe <- your_dataframe %>%
      filter(.data[[filtre_col]] %in% filtre_valeurs)
    
    cat("✓ Filtre appliqué:", filtre_col, "=", paste(filtre_valeurs, collapse = ", "), "\n")
    cat("  Nombre de lignes après filtre:", nrow(your_dataframe), "\n")
  }
  
  # Calculer les fréquences
  frequences <- calculer_frequences(your_dataframe, your_col)
  
  # Titre par défaut avec info sur le filtre
  if (is.null(titre_principal)) {
    titre_principal <- paste("Fréquence des séquences par", your_col)
    if (!is.null(filtre_col)) {
      titre_principal <- paste0(titre_principal, "\n(", filtre_col, " = ", 
                                paste(filtre_valeurs, collapse = ", "), ")")
    }
  }
  
  # Créer le graphique avec facettes
  p <- ggplot(frequences, aes(x = type_sequence, y = frequence, fill = type_sequence)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(pourcentage, 1), "%\n(n=", n, ")")), 
              vjust = -0.3, size = 3) +
    facet_wrap(as.formula(paste("~", your_col)), scales = "free_y") +
    labs(
      title = titre_principal,
      x = "Type de séquence",
      y = "Pourcentage d'individu",
      fill = "Type de séquence"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 2.5, face = "bold")
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = palette_fixe)
  
  print(p)
  
  # Sauvegarder si demandé
  if (sauvegarder) {
    suffixe <- if (!is.null(filtre_col)) paste0("_", filtre_col, "_", paste(filtre_valeurs, collapse = "_")) else ""
    nom_fichier <- paste0("histogrammes_", your_col, suffixe, ".png")
    ggsave(nom_fichier, plot = p, width = 12, height = 8, dpi = 300)
    cat("✓ Graphique sauvegardé:", nom_fichier, "\n")
  }
  
  return(frequences)
}



# Fonction : Creation tableau récapitulatif de colonne, filtre possible 
creer_tableau_recap <- function(your_dataframe, your_col, sauvegarder = TRUE, 
                                filtre_col = NULL, filtre_valeurs = NULL) {
  
  # Appliquer le filtre si spécifié
  if (!is.null(filtre_col) && !is.null(filtre_valeurs)) {
    if (!filtre_col %in% names(your_dataframe)) {
      stop(paste("La colonne de filtre", filtre_col, "n'existe pas dans les données"))
    }
    your_dataframe <- your_dataframe %>%
      filter(.data[[filtre_col]] %in% filtre_valeurs)
    
    cat("✓ Filtre appliqué:", filtre_col, "=", paste(filtre_valeurs, collapse = ", "), "\n")
    cat("  Nombre de lignes après filtre:", nrow(your_dataframe), "\n")
  }
  
  recap <- your_dataframe %>%
    group_by(across(all_of(your_col))) %>%
    summarise(
      total_individus = n(),
      buchnera_seul = sum(type_sequence == "Buchnera uniquement"),
      les_deux = sum(type_sequence == "Buchnera + Serratia"),
      .groups = "drop"
    ) %>%
    mutate(
      pct_buchnera_seul = round(buchnera_seul / total_individus * 100, 1),
      pct_les_deux = round(les_deux / total_individus * 100, 1)
    )
  
  print(recap)
  
  if (sauvegarder) {
    suffixe <- if (!is.null(filtre_col)) paste0("_", filtre_col, "_", paste(filtre_valeurs, collapse = "_")) else ""
    nom_fichier <- paste0("recap_sequences_", your_col, suffixe, ".csv")
    write.csv(recap, nom_fichier, row.names = FALSE)
    cat("✓ Tableau sauvegardé:", nom_fichier, "\n")
  }
  
  return(recap)
}

# Fonction : creation d'une carte selon colonne 
carte_points <- function(df, country_col = "Country", titre = NULL, sauvegarder = FALSE) {
  
  if (!requireNamespace("maps", quietly = TRUE)) {
    install.packages("maps")
    library(maps)
  }
  
  # Calculer statistiques par pays
  stats_pays <- df %>%
    group_by(.data[[country_col]]) %>%
    summarise(
      total = n(),
      buchnera_seul = sum(type_sequence == "Buchnera uniquement"),
      serratia_seul = sum(type_sequence == "Serratia uniquement"),
      les_deux = sum(type_sequence == "Buchnera + Serratia"),
      .groups = "drop"
    ) %>%
    mutate(
      pct_buchnera = buchnera_seul / total * 100,
      pct_serratia = serratia_seul / total * 100,
      pct_deux = les_deux / total * 100
    )
  
  # Carte du monde
  world_map <- map_data("world")
  
  # Correspondances
  correspondances_maps <- c(
    "United States of America" = "USA",
    "United Kingdom" = "UK",
    "Czechia" = "Czech Republic"
  )
  
  stats_pays <- stats_pays %>%
    mutate(country_map = ifelse(.data[[country_col]] %in% names(correspondances_maps),
                                correspondances_maps[.data[[country_col]]],
                                .data[[country_col]]))
  
  # Centroïdes
  centroides <- world_map %>%
    group_by(region) %>%
    summarise(lon = mean(long), lat = mean(lat), .groups = "drop")
  
  # Joindre
  data_carte <- stats_pays %>%
    left_join(centroides, by = c("country_map" = "region")) %>%
    filter(!is.na(lon))
  
  # Titre
  if (is.null(titre)) {
    titre <- "Échantillonnage par pays"
  }
  
  # Carte
  p <- ggplot() +
    geom_polygon(data = world_map, 
                 aes(x = long, y = lat, group = group),
                 fill = "lightgray", color = "white", size = 0.2) +
    geom_point(data = data_carte,
               aes(x = lon, y = lat, size = total, 
                   color = pct_deux,
                   text = paste0(.data[[country_col]], "\n",
                                 "Total: ", total, " individus\n",
                                 "Buchnera seul: ", round(pct_buchnera, 1), "%\n",
                                 "Serratia seul: ", round(pct_serratia, 1), "%\n",
                                 "Les deux: ", round(pct_deux, 1), "%")),
               alpha = 0.8) +
    scale_size_continuous(range = c(4, 10), name = "Nombre d'individus") +
    scale_color_gradient(low = "#B2A5F3", high = "#FF9A47", 
                         name = "% avec les deux séquences") +
    labs(title = titre) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed(1.3)
  
  print(p)
  
  if (sauvegarder) {
    ggsave("carte_points.png", plot = p, width = 14, height = 8, dpi = 300)
    cat("✓ Carte sauvegardée: carte_points.png\n")
  }
  
  return(p)
}

# Fonction : extraire les longueurs de branches des feuilles 
extract_length <- function (tree) {
  tip_edges <- which(tree$edge[,2] <= Ntip(tree))
  
  # Récupérer les longueurs et les noms associés
  tip_labels <- tree$tip.label[tree$edge[tip_edges, 2]]
  tip_lengths <- tree$edge.length[tip_edges]
  
  # Créer le dataframe
  df_tips <- data.frame(
    feuille       = tip_labels,
    longueur_branche = tip_lengths
  )
}

# Creer un df des moyennes des colonnes numérique d'un df donnée
tab_mean_gene <- function (df) {
  rownames(df) <- df$INDIVIDU
  df <- df [,-c(1,3,4)]
  df[] <- lapply(df, as.numeric)
  
  # Calcul des moyennes de toutes les colonnes numériques
  moyennes <- colMeans(df[, sapply(df, is.numeric)], na.rm = TRUE)
  
  # Création du dataframe résumé
  df_resume <- as.data.frame(t(moyennes))
}

# Fonction : produire des df de moyenne, rajout de ratio et possible groupe 
make_summary_tables <- function(tab, 
                                total_nucleotide, 
                                groupe = NULL) {
  
  # Filtrage par groupe si un vecteur d'individus est fourni
  if (!is.null(groupe)) {
    tab    <- tab    %>% filter(INDIVIDU %in% groupe)
  }
  
  # Calcul des moyennes
  mean_mut    <- tab_mean_gene(tab)
  
  # Ajout des ratios SNP et INDEL
  mean_mut <- mean_mut %>%
    mutate(ratio_snp   = SNP   / total_nucleotide,
           ratio_indel = INDEL / total_nucleotide)
  
  # Retourne une liste avec les deux tableaux
  return(mean_mut)
}

# Fonction : nettoyage des arbres et du df lié (= verification des associations de noms)
nettoyage_arbre <- function (tree, 
                             metadata,
                             col_taxon, 
                             col_localisation, 
                             col_plante) {
  
  # Individus dans le tableau mais absents de l'arbre → signalement + suppression
  absents <- metadata[[col_taxon]][!metadata[[col_taxon]] %in% tree$tip.label]
  if (length(absents) > 0) {
    message("Individus dans le tableau mais absents de l'arbre : ",
            paste(absents, collapse = ", "))
  }
  metadata_clean <- metadata %>%
    filter(.data[[col_taxon]] %in% tree$tip.label)
  
  # Cellules vides → NA propre
  metadata_clean <- metadata_clean %>%
    mutate(
      across(c(all_of(col_localisation), all_of(col_plante)),
             ~ na_if(trimws(as.character(.x)), ""))
    )
  
  # Feuilles de l'arbre sans données dans le tableau → ajout avec NA
  tips_manquants <- setdiff(tree$tip.label, metadata_clean[[col_taxon]])
  if (length(tips_manquants) > 0) {
    message("Feuilles sans données dans le tableau : ",
            paste(tips_manquants, collapse = ", "))
    df_manquants <- data.frame(
      taxon = tips_manquants,
      stringsAsFactors = FALSE
    )
    names(df_manquants)[1] <- col_taxon
    metadata_clean <- bind_rows(metadata_clean, df_manquants)
  }
  
  # Rownames pour %<+%
  rownames(metadata_clean) <- metadata_clean[[col_taxon]]
  
  return(metadata_clean)
}

# Fonction : figure phylogénie avec plante hote 
phylo_plante <- function(tree, 
                         metadata_clean, 
                         col_plante,
                         col_taxon) {
  
  p_test <- ggtree(tree) %<+% metadata_clean
  range(p_test$data$x, na.rm = TRUE)  # limites x
  range(p_test$data$y, na.rm = TRUE)  
  
  p <- p <- ggtree(tree,
                   layout = "rectangular",
                   color  = "grey40",
                   size   = 0.4) +
    geom_treescale(
      x        = min(p_test$data$x, na.rm = TRUE),  # droite du plot
      y        = -2,   # légèrement sous le dernier tip (y=1 est le bas)
      width    = 0.00001,
      offset   = 0.5,
      color    = "black",
      linesize = 0.5,
      fontsize = 3
    ) +
    geom_tiplab(
      size     = 2.5,
      offset   = 0.000000001,
      align    = TRUE,
      linetype = NA
    ) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      legend.box       = "vertical",
      legend.spacing.y = unit(0.3, "cm")
    ) +
    ggtitle("Phylogénie – Plante hôte")
  
  p <- p +
    geom_fruit(
      data    = metadata_clean,
      geom    = geom_point,
      mapping = aes(y     = !!sym(col_taxon),
                    fill  = !!sym(col_plante),
                    shape = !!sym(col_plante)),
      size    = 2,
      stroke  = 0.3,
      color   = "black",
      offset  = 0.09,
      axis.params = list(axis = "none")  # "none" au lieu de NULL
    ) +
    scale_fill_manual(
      name     = "Plante hôte",
      values   = palette_plante,
      na.value = "grey60"
    ) +
    scale_shape_manual(
      name     = "Plante hôte",
      values   = formes_plante,
      na.value = 21
    )
  
  print(p)
}

# Fonction : figure phylogénie avec géographie
phylo_loc <- function(tree,
                      metadata_clean,
                      col_taxon,
                      col_localisation,
                      titre = "Phylogénie – Localisation") {
  
  couleurs_texte_loc <- setNames(
    rep("black", length(palette_localisation)),
    names(palette_localisation)
  )
  
  p_test <- ggtree(tree) %<+% metadata_clean
  range(p_test$data$x, na.rm = TRUE)  # limites x
  range(p_test$data$y, na.rm = TRUE)  
  
  p <- ggtree(tree,
              layout = "rectangular",
              color  = "grey40",
              size   = 0.4) %<+% metadata_clean +
    geom_treescale(
      x        = min(p_test$data$x, na.rm = TRUE),  # droite du plot
      y        = -2,   # légèrement sous le dernier tip (y=1 est le bas)
      width    = 0.00001,
      offset   = 0.5,
      color    = "black",
      linesize = 0.5,
      fontsize = 3
    ) +
    
    geom_tiplab(
      aes(label = label,
          fill  = .data[[col_localisation]],
          color = .data[[col_localisation]]),
      geom          = "label",
      size          = 2,
      offset        = 0.000000001,
      align         = TRUE,
      linetype      = NA,
      label.padding = unit(0.05, "lines"),
      label.r       = unit(0, "lines"),
      label.size    = 0
    ) +
    scale_fill_manual(
      name     = "Localisation",
      values   = palette_localisation,
      na.value = "white"
    ) +
    scale_color_manual(
      name   = "Localisation",
      values = couleurs_texte_loc,
      guide  = "none"
    ) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      legend.box       = "vertical",
      legend.spacing.y = unit(0.35, "cm")
    ) +
    ggtitle(titre)
  
  print(p)
}

# Fonction : figure phylogénie avec combinaison géographie et plante
phylo_loc_plante <- function(tree,
                             metadata_clean,
                             col_taxon,
                             col_localisation,
                             col_plante,
                             titre = "Phylogénie – Localisation & Plante hôte",
                             width = 2800, height = 3200, res = 300) {
  
  png(output_file, width = width, height = height, res = res)

  
  couleurs_texte_loc <- setNames(
    rep("black", length(palette_localisation)),
    names(palette_localisation)
  )
  
  p_test <- ggtree(tree) %<+% metadata_clean
  range(p_test$data$x, na.rm = TRUE)  # limites x
  range(p_test$data$y, na.rm = TRUE)  
  
  p <- ggtree(tree,
              layout = "rectangular",
              color  = "grey40",
              size   = 0.4) %<+% metadata_clean +
    
    geom_treescale(
      x        = min(p_test$data$x, na.rm = TRUE),  # droite du plot
      y        = -2,   # légèrement sous le dernier tip (y=1 est le bas)
      width    = 0.0001,
      offset   = 0.5,
      color    = "black",
      linesize = 0.5,
      fontsize = 3
    ) +
    
    geom_tiplab(
      aes(label = label,
          fill  = .data[[col_localisation]],
          color = .data[[col_localisation]]),
      geom          = "label",
      size          = 2,
      offset        = 0.000000001,
      align         = TRUE,
      linetype      = NA,
      label.padding = unit(0.05, "lines"),
      label.r       = unit(0, "lines")
    ) +
    scale_fill_manual(
      name     = "Localisation",
      values   = palette_localisation,
      na.value = "white"
    ) +
    scale_color_manual(
      name   = "Localisation",
      values = couleurs_texte_loc,
      guide  = "none"
    ) +
    
    ggnewscale::new_scale_fill() +
    
    # ── data = NULL : utilise les données déjà attachées via %<+% ────────────
    geom_fruit(
      data    = NULL,
      geom    = geom_point,
      mapping = aes(fill  = .data[[col_plante]],
                    shape = .data[[col_plante]]),
      size    = 2.5,
      stroke  = 0.3,
      color   = "black",
      offset  = 0.09,
      axis.params = list(axis = "none")
    ) +
    scale_fill_manual(
      name     = "Plante hôte",
      values   = palette_plante,
      na.value = "grey60"
    ) +
    scale_shape_manual(
      name     = "Plante hôte",
      values   = formes_plante,
      na.value = 21
    ) +
    
    theme(
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      legend.box       = "vertical",
      legend.spacing.y = unit(0.35, "cm")
    ) +
    ggtitle(titre)
  
  print(p)
  dev.off()
  message("Figure sauvegardée : ", output_file)
}

#Fonction : figure de la co phylogénie 
cophylo_loc <- function(tree_left, tree_right,
                        metadata_left, metadata_right,
                        col_localisation,
                        titre = "Cophylogénie – Localisation") {
  
  # ── 1. Arbres COMPLETS (pas de pruning) ───────────────────────────────────
  # On garde les arbres entiers, les liens ne concernent que les tips communs
  tl <- tree_left
  tr <- tree_right
  if (!is.rooted(tl)) tl <- unroot(tl)
  if (!is.rooted(tr)) tr <- unroot(tr)
  
  # Tips communs uniquement pour les liens
  tips_communs <- intersect(tl$tip.label, tr$tip.label)
  
  # Association : tip gauche -> tip droit pour les liens
  assoc <- cbind(tips_communs, tips_communs)
  
  obj <- cophylo(tl, tr, assoc = assoc, rotate = FALSE)
  
  # ── 2. Couleurs ───────────────────────────────────────────────────────────
  get_cols <- function(labels, metadata) {
    cols <- palette_localisation[metadata[[col_localisation]][
      match(labels, metadata$read_file_name)]]
    cols[is.na(cols)] <- "grey85"
    cols
  }
  
  cols_left  <- get_cols(obj$trees[[1]]$tip.label, metadata_left)
  cols_right <- get_cols(obj$trees[[2]]$tip.label, metadata_right)
  
  # ── 3. Plot de base sans labels ni points ─────────────────────────────────
  plot(obj,
       fsize    = 0.0001,
       link.col = "grey70",
       link.lwd = 0.6,
       ftype    = "off",
       points   = FALSE)        # ← supprime les ronds noirs
  
  title(titre, line = 1)
  
  # ── 4. Coordonnées ────────────────────────────────────────────────────────
  pp <- get("last_plot.cophylo", envir = .PlotPhyloEnv)
  
  n_left  <- Ntip(obj$trees[[1]])
  n_right <- Ntip(obj$trees[[2]])
  
  x_left  <- pp[[1]]$xx[seq_len(n_left)]
  y_left  <- pp[[1]]$yy[seq_len(n_left)]
  x_right <- pp[[2]]$xx[seq_len(n_right)]
  y_right <- pp[[2]]$yy[seq_len(n_right)]
  
  labels_left  <- obj$trees[[1]]$tip.label
  labels_right <- obj$trees[[2]]$tip.label
  
  # ── 5. Calcul dynamique de cex pour éviter les chevauchements ─────────────
  # On ajuste cex selon le nombre de tips et la hauteur disponible
  usr    <- par("usr")          # c(x1, x2, y1, y2) en coordonnées user
  height <- usr[4] - usr[3]
  
  n_max  <- max(n_left, n_right)
  # Hauteur disponible par tip (en coordonnées user)
  h_per_tip <- height / n_max
  
  # On teste différents cex jusqu'à ce que strheight < h_per_tip
  cex_label <- 0.55
  while (strheight("A", cex = cex_label) * 1.3 > h_per_tip && cex_label > 0.15) {
    cex_label <- cex_label - 0.02
  }
  
  pad <- strwidth("i", cex = cex_label) * 0.4   # padding proportionnel à la taille du texte
  
  # ── 6. Fonction dessin rectangle + texte ──────────────────────────────────
  draw_labels <- function(labels, x_tips, y_tips, cols, side = "left") {
    for (i in seq_along(labels)) {
      w  <- strwidth(labels[i],  cex = cex_label)
      h  <- strheight(labels[i], cex = cex_label) * 1.3
      
      if (side == "left") {
        xl <- x_tips[i] + pad
        xr <- xl + w + pad
        rect(xl, y_tips[i] - h/2, xr, y_tips[i] + h/2,
             col = cols[i], border = NA)
        text(xl + pad/2, y_tips[i], labels[i],
             cex = cex_label, adj = c(0, 0.5), col = "black")
      } else {
        xr <- x_tips[i] - pad
        xl <- xr - w - pad
        rect(xl, y_tips[i] - h/2, xr, y_tips[i] + h/2,
             col = cols[i], border = NA)
        text(xr - pad/2, y_tips[i], labels[i],
             cex = cex_label, adj = c(1, 0.5), col = "black")
      }
    }
  }
  
  draw_labels(labels_left,  x_left,  y_left,  cols_left,  side = "left")
  draw_labels(labels_right, x_right, y_right, cols_right, side = "right")
  
  # ── 7. Légende ────────────────────────────────────────────────────────────
  loc_presentes <- sort(unique(na.omit(c(
    metadata_left[[col_localisation]][match(tips_communs, metadata_left$read_file_name)],
    metadata_right[[col_localisation]][match(tips_communs, metadata_right$read_file_name)]
  ))))
  
  legend("topleft",
         legend     = loc_presentes,
         fill       = palette_localisation[loc_presentes],
         border     = NA,
         ncol       = ceiling(length(loc_presentes) / 3),
         cex        = 0.65,
         title      = "Localisation",
         title.font = 2,
         bty        = "n",
         xpd        = TRUE)
}

# Fonction : preparation des arbres pour la reconstruction = binaire, raciné et pas de lb=0 
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

# Fonction : établir les modele pour faire de la reconstruction le long d'un arbre 
creation_modele <- function (df, 
                             col_name, 
                             character, 
                             tree){
  
  tree_tips <- tree$tip.label
  # Aligner le dataframe sur l'ordre des tips de l'arbre
  sortedData <- df[match(tree_tips, df[[col_name]]), ]
  
  chara_reconstruit <- sortedData[[character]]
  names(chara_reconstruit) <- tree_tips
  
  # Reconstruction ancestrale avec les 3 modèles
  ERreconstruction  <- ace(chara_reconstruit, tree, type = "discrete", model = "ER")
  SYMreconstruction <- ace(chara_reconstruit, tree, type = "discrete", model = "SYM")
  ARDreconstruction <- ace(chara_reconstruit, tree, type = "discrete", model = "ARD")
  
  n_etats <- length(unique(chara_reconstruit))
  
  df_ER_ARD  <- n_etats * (n_etats - 1) - 1
  df_ER_SYM  <- n_etats * (n_etats - 1) / 2 - 1
  df_SYM_ARD <- n_etats * (n_etats - 1) / 2
  
  # Calcul des p-values
  pval_ER_ARD  <- 1 - pchisq(2 * abs(ERreconstruction$loglik  - ARDreconstruction$loglik),  df_ER_ARD)
  pval_ER_SYM  <- 1 - pchisq(2 * abs(ERreconstruction$loglik  - SYMreconstruction$loglik),  df_ER_SYM)
  pval_SYM_ARD <- 1 - pchisq(2 * abs(SYMreconstruction$loglik - ARDreconstruction$loglik),  df_SYM_ARD)
  
  cat("\nTests LRT (p-value < 0.05 = le modèle complexe est significativement meilleur) :\n")
  cat("  ER vs ARD  :", pval_ER_ARD,  "\n")
  cat("  ER vs SYM  :", pval_ER_SYM,  "\n")
  cat("  SYM vs ARD :", pval_SYM_ARD, "\n")
  
  # Sélection du meilleur modèle
  if (pval_ER_ARD < 0.05 && pval_SYM_ARD < 0.05) {
    best_model <- ARDreconstruction
    best_name  <- "ARD"
  } else if (pval_ER_SYM < 0.05 && pval_SYM_ARD >= 0.05) {
    best_model <- SYMreconstruction
    best_name  <- "SYM"
  } else {
    best_model <- ERreconstruction
    best_name  <- "ER"
  }
  
  message("\nMeilleur modèle sélectionné : ", best_name)
  
  return(best_model)
}

# Fonction : Figure avec application d'un modele pour reconstruire l'état au noeud

applis_modele <- function (best_reconstruction, 
                           tree, df, 
                           col_name, character,
                           fig_title,
                           palette_fixe,
                           save = FALSE, filename) {
  
  tree_tips <- tree$tip.label
  # Aligner le dataframe sur l'ordre des tips de l'arbre
  sortedData <- df[match(tree_tips, df[[col_name]]), ]
  
  chara_reconstruit <- sortedData[[character]]
  names(chara_reconstruit) <- tree_tips
  
  # Palette de couleurs (1 couleur par état géographique)
  etats     <- sort(unique(chara_reconstruit))
  n_etats   <- length(etats)
  couleurs <- palette_fixe[etats]
  
  # Visualisation
  
  if (save) {
    pdf(filename, width = 20, height = 30)  # ← utilise le nom choisi
  }
  
  plotTree(tree, lwd = 1, setEnv = TRUE, offset = 0.5, fsize = 0.4)
  
  # Récupère les limites réelles
  pp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xlim <- pp$x.lim
  ylim <- pp$y.lim
  
  # Juste à droite de la légende (bas gauche) 
  add.scale.bar(
    x      = xlim[2] * 0.85,   # bas droite
    y      = ylim[1],
    length = 0.0001,
    lwd    = 2,
    cex    = 0.7
  )
  
  # Noeuds internes : camemberts des probabilités ancestrales
  n_tips  <- length(tree_tips)
  n_nodes <- tree$Nnode
  nodelabels(
    node    = (n_tips + 1):(n_tips + n_nodes),
    pie     = best_reconstruction$lik.anc,
    piecol  = couleurs,
    cex     = 0.3
  )
  
  # Tips : état observé
  tiplabels(
    pie    = to.matrix(chara_reconstruit, etats),
    piecol = couleurs,
    cex    = 0.15
  )
  
  # Légende
  legend("bottomleft", legend = etats, fill = couleurs, cex = 0.7, fig_title)
  
  if (save) {
    dev.off()
    message("Fichier sauvegardé : ", filename)  # ← affiche le nom utilisé
  }
  
}


###### Charger Jeu de données propre ####

Echantillonnage_apisum <- read.csv2("~/Documents/Moana/Echantillonnage_apisum.csv", head=TRUE, sep=",")
Echantillonnage_apisum <- Echantillonnage_apisum [Echantillonnage_apisum$read_file_name != "CS189_S180",]
Echantillonnage_apisum <- Echantillonnage_apisum[Echantillonnage_apisum$qualite_sequence_buchnera != "mauvaise_qualite", ]

#Individu avec des mauvaise qualite supprimer 
Echantillonnage_apisum_filtre <- read.csv2("~/Documents/Moana/Echantillonnage_apisum_filtre.csv", head=TRUE, sep=",")
Echantillonnage_apisum_filtre <- Echantillonnage_apisum_filtre [Echantillonnage_apisum_filtre$read_file_name != "CS189_S180",]

tab_codant_buchnera <- read.csv2("~/Documents/Moana/taux_evolution/tab_resum_codant_buchnera.csv", head=TRUE, sep=",")
tab_noncodant_buchnera <- read.csv2("~/Documents/Moana/taux_evolution/tab_resum_noncodant_buchnera.csv", head=TRUE, sep=",")

tab_noncodant_serratia <- read.csv2("~/Documents/Moana/taux_evolution/resum_noncodant_serratia_ET.csv", head=TRUE, sep=",")
tab_codant_serratia <- read.csv2("~/Documents/Moana/taux_evolution/resum_codant_serratia_ET.csv", head=TRUE, sep=",") 

longueur_gene_buchnera <- read.csv2 ("~/Documents/Moana/taux_evolution/longueur_gene_buchnera.csv", head = TRUE, sep=",")
longueur_gene_serratia <- read.csv2 ("~/Documents/Moana/taux_evolution/longueur_gene_serratia.csv", head = TRUE, sep=",")

phylo_buchnera <- read.tree("~/Documents/Moana/phylo_snippy/phylo_buchnera_without_outgroup.nwk")
phylo_buchnera_incomplet<- read.tree("~/Documents/Moana/phylo_snippy/phylo_buchnera_incomplet.nwk")
phylo_serratia <- read.tree("~/Documents/Moana/phylo_snippy/phylo_serratia_without_outgroup.nwk")

data_clim <- read.csv2("~/Documents/Moana/presabs_serratia.csv", head=TRUE, sep=",")

palette_presence_serratia <- c(
  "present"     = "#FF9A47",  
  "absent" = "#B2A5F3"
)

palette_presence_serratia_histo <- c(
  "Buchnera + Serratia" = "#FF9A47",  
  "Buchnera uniquement" = "#B2A5F3" 
)

palette_resum_localisation <- c(
  "Amerique"     = "#FFCC66", 
  "Asie_central" = "#778F00",  
  "Asie_est"     = "#4EACF9",  
  "Asie_ouest"   = "#B2A5F3", 
  "Europe"       = "#FF667A" 
)

palette_localisation <- c(
  "Japan"  = "#00CCCC",
  "China" = "#3366FF" ,
  "Kazakhstan"   = "#9AB87A",
  "Chili"    = "#EBE06C",
  "Canada" = "#FFCC66",
  "USA" = "#FFFF33",
  "France"   =  "#E31C00",
  "Spain" = "#FF4100",
  "Germany"  = "#960808",
  "UK" = "#FF7700",
  "Israel" = "#CC3399",
  "Cyprus" = "#FFCCFF",
  "Iran" = "#9933CC",
  "Tunisia"  = "#CC6699" 
)

palette_plante <- c(
  "Melilotus_suaveolens" = "#FF6600" ,
  "Melilotus_officinalis"  = "#FF6600" ,
  "Melilotus_sp"  = "#FF6600" ,
  "Melilotus_alba"  =  "#FF6600",
  "Lathyrus_pratensis"  = "#9900FF",
  "Lathyrus_tuberosus"  = "#9900FF",
  "Lathyrus"    =  "#9900FF",
  "Ononis_spinosa"  =  "#CC0099",
  "Pisum_sativum"   =  "#66CC66",
  "Medicago_sativa" =   "#660000",
  "Medicago_lupulina"   = "#660000" ,
  "Vicia_cracca"  = "#FFCA3A",
  "Vicia_sp"  =  "#FFCA3A",
  "Trifolium_pratense"  =  "#FF99FF" ,
  "Trifolium sp"   = "#FF99FF" ,
  "Genista_sagittalis"  =  "#33CCCC",
  "Genista_tinctoria"  = "#33CCCC",
  "Lotus_corniculatus"  =  "#FF9999",
  "Lotus_pedunculatus"  = "#FF9999" ,
  "Onobrychis"    = "#CCCC99" ,
  "Cytisus" =  "orange",
  "Cytisus_scoparius"   = "orange",
  "Securigera"  = "#3A42D0" ,
  "Securigera_varia"  =  "#3A42D0",
  "Sophora_flavescens"  =  "#996600",
  "Macrosiphum silenum"   =  "#99CC00",
  "Macrosiphum funestum"  =  "#99CC00"
)

formes_plante <- c(
  "Melilotus_suaveolens" = 21,
  "Melilotus_officinalis" = 22,  # carré
  "Melilotus_sp" = 23,            # diamant
  "Melilotus_alba" = 24,          # triangle
  "Lathyrus_pratensis" = 21,
  "Lathyrus_tuberosus" = 22,
  "Lathyrus" = 24,
  "Genista_sagittalis" = 22,
  "Genista_tinctoria" = 21,
  "Macrosiphum funestum" = 22,
  "Macrosiphum silenum" = 21,
  "Ononis_spinosa" = 21,
  "Pisum_sativum" = 21,
  "Lotus_corniculatus" = 21,
  "Lotus_pedunculatus" = 22,
  "Medicago_sativa" = 21,
  "Medicago_lupulina" = 22,
  "Vicia_cracca" = 22,
  "Vicia_sp" = 21,
  "Trifolium_pratense" = 21,
  "Trifolium sp" = 22,
  "Sophora_flavescens" = 21,
  "Cytisus_scoparius" = 21,
  "Cytisus" = 22,
  "Securigera_varia" = 22,
  "Securigera" = 21,
  "Onobrychis" = 21
)

palette_resum_plante <- c(
  "Melilotus_sp"  = "#FF6600" ,
  "Lathyrus_sp"    =  "#9900FF",
  "Ononis_spinosa"  =  "#CC0099",
  "Pisum_sativum"   =  "#66CC66",
  "Medicago_sp" =   "#660000",
  "Vicia_cracca"  = "#FFC",
  "Vicia_sp"  =  "#FFCA3A",
  "Trifolium_pratense"  =  "#FF99FF" ,
  "Genista_sp"  =  "#33CCCC",
  "Lotus_sp"  =  "#FF9999",
  "Onobrychis_sp"    = "#CCCC99" ,
  "Cytisus_sp" =  "orange",
  "Securigera_sp"  = "#3A42D0" 
)

###### Visualiser les données #####

table (Echantillonnage_apisum$qualite_sequence_buchnera)
table (Echantillonnage_apisum$qualite_sequence_serratia)

#verifier ou sont les Serratia nulle 
serratia_nul <- Echantillonnage_apisum [Echantillonnage_apisum$qualite_sequence_serratia == "mauvaise_qualite",]
table (serratia_nul$Country)

Echantillonnage_apisum <- Echantillonnage_apisum %>%
  mutate(type_sequence = case_when(
    qualite_sequence_buchnera %in% c("bonne_qualite", "mauvaise_qualite") & qualite_sequence_serratia %in% c("bonne_qualite", "mauvaise_qualite") ~ "Buchnera + Serratia",
    qualite_sequence_buchnera%in% c("bonne_qualite", "mauvaise_qualite") & qualite_sequence_serratia == "absent"   ~ "Buchnera uniquement",
  ))

histogramme_loc <- creer_histogrammes(Echantillonnage_apisum , "Country", palette_presence_serratia_histo, sauvegarder = TRUE)
histogramme_plante <- creer_histogrammes(Echantillonnage_apisum , "resum_biotype", palette_presence_serratia_histo, sauvegarder = TRUE)

creer_tableau_recap(Echantillonnage_apisum , "Country", sauvegarder = FALSE)
creer_tableau_recap(Echantillonnage_apisum , "biotype", sauvegarder = FALSE)

carte_points(Echantillonnage_apisum, country_col = "Country", sauvegarder = TRUE)

###### Analyse stat ######

Echantillonnage_apisum$profondeur_serratia <- as.numeric(Echantillonnage_apisum$profondeur_serratia)
str (Echantillonnage_apisum)
summary(Echantillonnage_apisum)

Echantillonnage_apisum_prof <- Echantillonnage_apisum %>% filter (!is.na (profondeur_serratia))
class(Echantillonnage_apisum_prof$qualite_sequence_serratia)
Echantillonnage_apisum_prof$qualite_sequence_serratia <- as.factor(Echantillonnage_apisum_prof$qualite_sequence_serratia)

#profondeur 
tapply(Echantillonnage_apisum_prof$profondeur_serratia, 
       Echantillonnage_apisum_prof$qualite_sequence_serratia, 
       summary)
# Voir si profondeur varie entre qualité (attention fait avec tout le monde meme ceux ou bubu de mauvaise qualité)
model_prof <- glm (profondeur_serratia ~ qualite_sequence_serratia, data = Echantillonnage_apisum_prof)
summary (model_prof)

#presence Serratia 

#besoin de la version numerique
Echantillonnage_apisum <- Echantillonnage_apisum %>%
  mutate(presence_serratia_num = case_when(
    qualite_sequence_serratia %in% c("bonne_qualite", "mauvaise_qualite")  ~ 1,
    qualite_sequence_serratia == "absent" ~ 0,
  ))

Echantillonnage_apisum$zone_geographique <- relevel(
  factor(Echantillonnage_apisum$zone_geographique), 
  ref = "Europe"
)

Echantillonnage_apisum$resum_biotype <- relevel(
  factor(Echantillonnage_apisum$resum_biotype), 
  ref = "Lathyrus_sp"
)

#climat : 
phylo_buchnera <- prepa_tree(phylo_buchnera, c("wamizawa", "CS041_S97", "CS042_S98") )

Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                               "ID_climat", "ID_climat", 
                                               "bio02")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                               "ID_climat", "ID_climat", 
                                               "bio05")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                               "ID_climat", "ID_climat", 
                                               "bio06")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                               "ID_climat", "ID_climat", 
                                               "bio10")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                               "ID_climat", "ID_climat", 
                                               "bio11")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                               "ID_climat", "ID_climat", 
                                               "bio17")

Echantillonnage_apisum <- Echantillonnage_apisum [Echantillonnage_apisum$read_file_name != "CS124_S169",]
str (Echantillonnage_apisum)
cols_environnemental <- c("bio02", "bio05", "bio06", "bio10", "bio11", "bio17")

#transformation numerique
Echantillonnage_apisum <- Echantillonnage_apisum %>%
  mutate(across(all_of(cols_environnemental), ~ as.numeric(.)))

Echantillonnage_apisum$presence_serratia_num <- as.numeric(Echantillonnage_apisum$presence_serratia_num)


#test signal phylo : 

lambda_pagel <- function (variable, tree) {
  var <-setNames(variable,
                   Echantillonnage_apisum [["read_file_name"]])
  
  lambda_var <- phylosig(tree, var,
                           method="lambda",test=TRUE)
  return(lambda_var)
}

lambda_bio02 <- lambda_pagel (Echantillonnage_apisum$bio02, phylo_buchnera)
lambda_bio05 <- lambda_pagel (Echantillonnage_apisum$bio05, phylo_buchnera)
lambda_bio06 <- lambda_pagel (Echantillonnage_apisum$bio06, phylo_buchnera)
lambda_bio10 <- lambda_pagel (Echantillonnage_apisum$bio10, phylo_buchnera)
lambda_bio11 <- lambda_pagel (Echantillonnage_apisum$bio11, phylo_buchnera)
lambda_bio17 <- lambda_pagel (Echantillonnage_apisum$bio17, phylo_buchnera)
print (lambda_bio02)
print (lambda_bio05)
print (lambda_bio06)
print (lambda_bio10)
print (lambda_bio11)
print (lambda_bio17)

library(phylosignalDB)

zone_geo_df <- data.frame(zone_geo = Echantillonnage_apisum$zone_geographique, row.names = Echantillonnage_apisum$read_file_name)
zone_geo_df[,1] <- as.factor(zone_geo_df[,1])
zone_geo_dist <- gower_dist(x = zone_geo_df, type = list(factor = 1))
phylosignal_M(zone_geo_dist, phylo_buchnera, reps = 999) # reps=999 better

biotype_df <- data.frame(biotype = Echantillonnage_apisum$resum_biotype, row.names = Echantillonnage_apisum$read_file_name)
biotype_df[,1] <- as.factor(biotype_df[,1])
biotype_dist <- gower_dist(x = biotype_df, type = list(factor = 1))
phylosignal_M(biotype_dist, phylo_buchnera, reps = 999) # reps=999 better

# Phylo Buchnera car on utilise aussi les individus pas ouf de serratia 

mod_brownian_complet <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio06 + bio10 + bio11 + bio17 + resum_biotype, correlation = corBrownian(phy = phylo_buchnera, form = ~read_file_name), data = Echantillonnage_apisum, method = "ML")
mod_OU_complet <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio06 + bio10 + bio11 + bio17 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
#mod_nul_complet <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio06 + bio10 + bio11 + bio17 + resum_biotype,data = Echantillonnage_apisum, method = "ML")

anova (mod_brownian_complet, mod_OU_complet)

mod_OU_biotype <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio06 + bio10 + bio11 + bio17, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
mod_OU_bio17 <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio06 + bio10 + bio11 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
mod_OU_bio11 <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio06 + bio10 + bio17 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
mod_OU_bio10 <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio06 + bio11 + bio17 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
mod_OU_bio06 <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio05 + bio10 + bio11 + bio17 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
mod_OU_bio05 <- gls (presence_serratia_num ~ zone_geographique  + bio02 + bio06 + bio10 + bio11 + bio17 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
mod_OU_bio02 <- gls (presence_serratia_num ~ zone_geographique  + bio05 + bio06 + bio10 + bio11 + bio17 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
mod_OU_zonegéographique <- gls (presence_serratia_num ~ bio02 + bio05 + bio06 + bio10 + bio11 + bio17 + resum_biotype, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")

anova (mod_OU_complet, mod_OU_biotype)
anova (mod_OU_complet, mod_OU_bio17)
anova (mod_OU_complet, mod_OU_bio11)
anova (mod_OU_complet, mod_OU_bio10)
anova (mod_OU_complet, mod_OU_bio06)
anova (mod_OU_complet, mod_OU_bio05)
anova (mod_OU_complet, mod_OU_bio02)
anova (mod_OU_complet, mod_OU_zonegéographique)


mod_OU_final <- gls (presence_serratia_num ~ zone_geographique + bio02 + bio05 + bio06 + bio10 + bio11 + bio17, correlation = corMartins(value = 1, phy = phylo_buchnera, form = ~read_file_name, fixed = FALSE), data = Echantillonnage_apisum, method = "ML")
anova (mod_OU_complet, mod_OU_final)
summary (mod_OU_final)

###boxplot climat 

cols_environnemental <- c("bio02", "bio05", "bio06", "bio10", "bio11", "bio17")

png("box_plot_climat_bio17.png", width = 1200, height = 800, res = 150)

ggplot(Echantillonnage_apisum, aes(x= presence_serratia, y= bio02, fill = presence_serratia))+ 
  geom_boxplot(alpha=0.8)+ 
  scale_fill_manual(values = c("present" = "#FF9A47", "absent" = "#B2A5F3")) +
  xlab("presence Serratia")
dev.off()

ggplot(Echantillonnage_apisum, aes(x= presence_serratia, y= bio05, fill = presence_serratia))+ 
  geom_boxplot(alpha=0.8)+ 
  scale_fill_manual(values = c("present" = "#FF9A47", "absent" = "#B2A5F3")) +
  xlab("presence Serratia")

ggplot(Echantillonnage_apisum, aes(x= presence_serratia, y= bio06, fill = presence_serratia))+ 
  geom_boxplot(alpha=0.8)+ 
  scale_fill_manual(values = c("present" = "#FF9A47", "absent" = "#B2A5F3")) +
  xlab("presence Serratia")

ggplot(Echantillonnage_apisum, aes(x= presence_serratia, y= bio10, fill = presence_serratia))+ 
  geom_boxplot(alpha=0.8)+ 
  scale_fill_manual(values = c("present" = "#FF9A47", "absent" = "#B2A5F3")) +
  xlab("presence Serratia")

ggplot(Echantillonnage_apisum, aes(x= presence_serratia, y= bio11, fill = presence_serratia))+ 
  geom_boxplot(alpha=0.8)+ 
  scale_fill_manual(values = c("present" = "#FF9A47", "absent" = "#B2A5F3")) +
  xlab("presence Serratia")

ggplot(Echantillonnage_apisum, aes(x= presence_serratia, y= bio17, fill = presence_serratia))+ 
  geom_boxplot(alpha=0.8)+ 
  scale_fill_manual(values = c("present" = "#FF9A47", "absent" = "#B2A5F3")) +
  xlab("presence Serratia")

dev.off()

make_boxplot <- function(var) {
  ggplot(Echantillonnage_apisum, aes(x = presence_serratia, y = .data[[var]], fill = presence_serratia)) +
    geom_boxplot(alpha = 1) +
    scale_fill_manual(values = c("present" = "#FF9A47", "absent" = "#B2A5F3")) +
    xlab("Présence Serratia") +
    ylab(var) +
    theme(legend.position = "none")  # légende unique à la fin
}

# Générer tous les plots
plots <- lapply(cols_environnemental, make_boxplot)

# Assembler en grille 2 colonnes × 3 lignes + légende commune
figure_finale <- wrap_plots(plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Sauvegarder
png("figure_bio_boxplots.png", width = 3600, height = 4200, res = 600)
print(figure_finale)
dev.off()

###### ACP #####
Echantillonnage_apisum <- Echantillonnage_apisum %>%
  mutate(presence_serratia_num = case_when(
    qualite_sequence_serratia %in% c("bonne_qualite", "mauvaise_qualite")  ~ 1,
    qualite_sequence_serratia == "absent" ~ 0,
  ))

phylo_buchnera <- prepa_tree(phylo_buchnera, c("wamizawa", "CS041_S97", "CS042_S98") )

Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                          "ID_climat", "ID_climat", 
                                          "bio02")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                          "ID_climat", "ID_climat", 
                                          "bio05")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                          "ID_climat", "ID_climat", 
                                          "bio06")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                          "ID_climat", "ID_climat", 
                                          "bio10")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                          "ID_climat", "ID_climat", 
                                          "bio11")
Echantillonnage_apisum <- ajout_colonnes (Echantillonnage_apisum, data_clim, 
                                          "ID_climat", "ID_climat", 
                                          "bio17")

Echantillonnage_apisum <- Echantillonnage_apisum [Echantillonnage_apisum$read_file_name != "CS124_S169",]
str (Echantillonnage_apisum)
cols_environnemental <- c("bio02", "bio05", "bio06", "bio10", "bio11", "bio17")

#transformation numerique
Echantillonnage_apisum <- Echantillonnage_apisum %>%
  mutate(across(all_of(cols_environnemental), ~ as.numeric(.)))
Echantillonnage_apisum$presence_serratia_num <- as.numeric(Echantillonnage_apisum$presence_serratia_num)
Echantillonnage_apisum["log_bio17"] <-log1p(Echantillonnage_apisum$bio17)

#matrice : 

matrice_climat <- as.matrix( Echantillonnage_apisum[,c("read_file_name","bio02", "bio05", "bio06", "bio10", "bio11", "bio17")])
rownames(matrice_climat) <- matrice_climat[,"read_file_name"]
matrice_climat <- matrice_climat[,-1]
str(matrice_climat)
matrice_climat <- matrix(as.numeric(matrice_climat), 
                     nrow = nrow(matrice_climat), 
                     ncol = ncol(matrice_climat),
                     dimnames = dimnames(matrice_climat))

pca_climat <- phyl.pca(phylo_buchnera, matrice_climat, method="BM", mode="cov")
biplot(pca_climat)
as.princomp(pca_climat)
as.prcomp(pca_climat)
eig <- pca_climat$Eval
inertie <- round((diag(eig) / sum(diag(eig))) * 100, 2)
print(inertie)

# Scores et loadings
scores   <- pca_climat$S
loadings <- pca_climat$L

# Pourcentages
eig     <- diag(pca_climat$Eval)
pc1_pct <- round((eig[1] / sum(eig)) * 100, 1)
pc2_pct <- round((eig[2] / sum(eig)) * 100, 1)

# Mise à l'échelle des loadings pour qu'ils s'affichent bien
scale_factor <- max(abs(scores)) / max(abs(loadings)) * 0.7
x_all <- c(scores[, 1], loadings[, 1] * scale_factor)
y_all <- c(scores[, 2], loadings[, 2] * scale_factor)

pdf("ACP_climat.pdf", width = 20, height = 30)
plot(scores[, 1], scores[, 2],
     col  = palette_presence_serratia,
     pch  = 16,
     xlab = paste0("PC1 (", pc1_pct, "%)"),
     ylab = paste0("PC2 (", pc2_pct, "%)"),
     main = "Phylogenetic PCA - Climat",
     xlim = range(x_all) * 1.2,
     ylim = range(y_all) * 1)

text(scores[, 1], scores[, 2],
     labels = rownames(scores),
     col    = palette_presence_serratia,
     cex    = 0.5,
     pos    = 3)

# Flèches variables
arrows(0, 0,
       loadings[, 1] * scale_factor,
       loadings[, 2] * scale_factor,
       col    = "darkgreen",
       length = 0.1)

text(loadings[, 1] * scale_factor,
     loadings[, 2] * scale_factor,
     labels = rownames(loadings),
     col    = "darkgreen",
     cex    = 0.8)

legend("bottomright",
       legend = c("Presente", "Absente"),
       col    = c("#FF9A47", "#B2A5F3"),
       pch    = 16)
dev.off()


##### Regression ####

sub_serratia<- extract_length(phylo_serratia)
sub_buchnera<- extract_length(phylo_buchnera_incomplet)

sub_total <- ajout_colonnes(
  df_cible = sub_serratia,
  df_source = sub_buchnera,
  col_nom_cible = "feuille",
  col_nom_source = "feuille",
  col_a_copier = "longueur_branche",
  nouveau_nom_col = "longueur_branche_buchnera"
)

sub_total <- ajout_colonnes(
  df_cible = sub_total,
  df_source = Echantillonnage_apisum,
  col_nom_cible = "feuille",
  col_nom_source = "read_file_name",
  col_a_copier = "Country",
  nouveau_nom_col = "country"
)

pdf("regression_substitution.pdf", width = 15, height = 15)
ggplot (sub_total) +
  aes (x = longueur_branche_buchnera, y = longueur_branche) +
  geom_point(aes(colour = country), alpha = 0.8) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  #geom_smooth (method = "lm")+
  labs (x = "Taux de substituion par site Buchnera", y = "Taux de substitution par site Serratia", colour = "Country") +
  scale_colour_manual(values = palette_localisation) +
  theme_light ()
dev.off ()

##### Tableau moran ####

# Etablir la longueur des nucléotide codant/non codant 

longueur_gene_buchnera$Length <- as.numeric(longueur_gene_buchnera$Length)
longueur_gene_serratia$Length <- as.numeric(longueur_gene_serratia$Length)


longueur_gene_serratia_filtre <- longueur_gene_serratia[!(longueur_gene_serratia[,1] == "gene" & longueur_gene_serratia[,2] == "gene"), ]

long_gene_b <-longueur_gene_buchnera [,5]
total_gene_buchnera <- sum (long_gene_b)
total_noncodant_buchnera <- 642122 - total_gene_buchnera

long_gene_s <-longueur_gene_serratia_filtre [,5]
total_gene_serratia <- sum (long_gene_s)

long_gene_s_non_et <-longueur_gene_serratia [,5]
total_gene_serratia_non_et <- sum (long_gene_s_non_et)
total_noncodant_serratia <- 2736352 - total_gene_serratia_non_et

total_noncodant_buchnera <- total_noncodant_buchnera/1000
total_noncodant_serratia <-total_noncodant_serratia/1000
total_gene_serratia <- total_gene_serratia/1000
total_gene_buchnera <- total_gene_buchnera/1000

# Tableau brut, ensemble des individus : 

mean_noncodant_buchnera <- make_summary_tables(tab_noncodant_buchnera, total_noncodant_buchnera)
mean_codant_buchnera <- make_summary_tables(tab_codant_buchnera, total_gene_buchnera)
mean_noncodant_serratia <- make_summary_tables(tab_noncodant_serratia, total_noncodant_serratia)
mean_codant_serratia <- make_summary_tables(tab_codant_serratia, total_gene_serratia)

# Tableau par groupe : 

groupe_2_name <- c("CS036_S36", "CS206", "CS011_S11", "CS027_S27", "CS021_S21", "CS228", "CS231", "CS030_S30", "CS238", "CS019_S19", 
                   "CS233", "CS020_S20", "CS016_S16", "CS035_S35", "CS028_S28", "CS234", "CS203", "CS224", "CS227", "CS204", "CS200", 
                   "CS007_S7", "CS008_S8", "CS223", "CS029_S29", "CS010_S10", "CS218", "CS001_S1", "CS009_S9", "CS022_S22", "CS229", 
                   "CS237", "CS114_S125", "CS048_S44", "CS113_S124", "CS046_S43", "CS0054_S48", "CS052_S102", "CS004_S4", "CS236", "CS032_32", "CS230", "CS235", "CS037_S37", "CS212", "CS214", 
                   "CS222", "CS207", "CS215", "CS003_S3", "CS195", "CS213")
groupe_1_name <- c( "CS006_S6", "CS002_S2", "CS209", "CS210", "CS193", "CS208", "CS159_S148", "CS162_S150", "CS163_S151", "CS163_S151", "CS225", "CS202",
                    "CS221", "CS013_S13", "CS219", "CS198", "CS211", "CS205", "CS199", "CS018_S18", "CS197", "CS161_S172", "CS164_S152", "CS094_S177", "CS087_S67", 
                    "CS093_S192", "CS101_S118", "CS099_S71", "CS102_S119", "CS104_S72", "CS176_S159", "CS190_S166", "CS005_S5", "CS014_S14")
groupe_3_name <- c("CS144_S82", "CS185_S165", "CS118_S126", "CS090_S114", "CS078_S110", "CS098_S117", "CS129_S128", "CS179_S91", "CS076_S108", "CS119_S186", "CS089_S69", "CS170_S87", "CS023_S23")

#les buchnera et non codant
mean_noncodant_buchnera_groupe1 <- make_summary_tables(tab_noncodant_buchnera, total_noncodant_buchnera, groupe = groupe_1_name)
mean_noncodant_buchnera_groupe2 <- make_summary_tables(tab_noncodant_buchnera, total_noncodant_buchnera, groupe = groupe_2_name)
mean_noncodant_buchnera_groupe3 <- make_summary_tables(tab_noncodant_buchnera, total_noncodant_buchnera, groupe = groupe_3_name)
# les buchnera et codant 
mean_codant_buchnera_groupe1 <- make_summary_tables(tab_codant_buchnera, total_gene_buchnera, groupe = groupe_1_name)
mean_codant_buchnera_groupe2 <- make_summary_tables(tab_codant_buchnera, total_gene_buchnera, groupe = groupe_2_name)
mean_codant_buchnera_groupe3 <- make_summary_tables(tab_codant_buchnera, total_gene_buchnera, groupe = groupe_3_name)
# les serratia codant
mean_codant_serratia_groupe1 <- make_summary_tables(tab_codant_serratia, total_gene_serratia, groupe = groupe_1_name)
mean_codant_serratia_groupe2 <- make_summary_tables(tab_codant_serratia, total_gene_serratia, groupe = groupe_2_name)
mean_codant_serratia_groupe3 <- make_summary_tables(tab_codant_serratia, total_gene_serratia, groupe = groupe_3_name)
#les serratia et non codant 
mean_noncodant_serratia_groupe1 <- make_summary_tables(tab_noncodant_serratia, total_noncodant_serratia, groupe = groupe_1_name)
mean_noncodant_serratia_groupe2 <- make_summary_tables(tab_noncodant_serratia, total_noncodant_serratia, groupe = groupe_2_name)
mean_noncodant_serratia_groupe3 <- make_summary_tables(tab_noncodant_serratia, total_noncodant_serratia, groupe = groupe_3_name)

library(knitr)
library(tinytex)
# Récupère automatiquement tous les objets dont le nom commence par "mean_"
noms <- ls(pattern = "^mean_")
summary_table <- bind_rows(mget(noms), .id = "Jeu_de_données")
kable(summary_table, digits = 4)

summary_table <- summary_table %>%
  mutate(
    organisme = case_when(
      grepl("buchnera", `Jeu_de_données`) ~ "Buchnera",
      grepl("serratia", `Jeu_de_données`) ~ "Serratia",
      TRUE ~ "Inconnu"
    ),
    region = case_when(
      grepl("noncodant", `Jeu_de_données`)                                 ~ "Non codant",
      grepl("codant", `Jeu_de_données`) & !grepl("non", `Jeu_de_données`) ~ "Codant",
      TRUE ~ "Inconnu"
    ),
    groupe = case_when(
      grepl("groupe1$",    `Jeu_de_données`) ~ "Groupe 1",
      grepl("groupe2$",    `Jeu_de_données`) ~ "Groupe 2",
      grepl("groupe3$",    `Jeu_de_données`) ~ "Groupe 3",
      TRUE ~ "Global"
    )
  ) %>%
  mutate(groupe = factor(groupe, levels = c(
    "Global", "Groupe 1", "Groupe 2",
    "Groupe 3"
  ))) %>%
  arrange(organisme, region, groupe)

# Formatage dynamique des petits nombres 
# Affiche jusqu'au premier chiffre significatif non nul + 1 décimale supplémentaire
format_small <- function(x) {
  sapply(x, function(v) {
    if (is.na(v) || v == 0) return("0")
    if (abs(v) < 0.0001) {
      # Notation scientifique : ex. 2.71e-06
      formatC(v, format = "e", digits = 2)
    } else {
      n_dec <- ceiling(-log10(abs(v))) + 1
      n_dec <- max(2, min(n_dec, 8))
      formatC(v, format = "f", digits = n_dec)
    }
  })
}

# Vérifie les valeurs brutes (affiche dans la console)
message("Valeurs brutes ratio_indel :")
print(summary_table$ratio_indel)

# Applique le formatage
# IMPORTANT : pas de round() sur les ratios, format_small() gère seul les décimales
summary_table <- summary_table %>%
  mutate(
    SNP           = round(SNP,   2),
    INDEL         = round(INDEL, 2),
    ratio_snp_f   = format_small(ratio_snp),
    ratio_indel_f = format_small(ratio_indel)
  )

message("Valeurs formatées ratio_indel :")
print(summary_table$ratio_indel_f)

# ── Fonctions utilitaires LaTeX ────────────────────────────────────────────────
latex_escape <- function(x) {
  x <- gsub("&",  "\\\\&",  x)
  x <- gsub("%",  "\\\\%",  x)
  x <- gsub("_",  "\\\\_",  x)
  x <- gsub("#",  "\\\\#",  x)
  x
}

# Formate un entier avec séparateur de milliers (espace fine LaTeX)
format_nucleotides <- function(n) {
  # formatC avec big.mark donne des virgules ; on remplace par \, (espace fine LaTeX)
  s <- formatC(as.integer(n), format = "d", big.mark = ",")
  gsub(",", "\\\\,", s)
}

make_row <- function(row, is_global) {
  # Colonnes affichées : Groupe | SNP | INDEL | SNP/nt | INDEL/nt
  # MNP et COMPLEX supprimés
  vals <- c(
    as.character(row$groupe),
    as.character(row$SNP),
    as.character(row$INDEL),
    as.character(row$ratio_snp_f),
    as.character(row$ratio_indel_f)
  )
  line <- paste(latex_escape(vals), collapse = " & ")
  if (is_global) {
    line <- paste0("\\color{mutedgray}\\textit{",
                   paste(latex_escape(vals), collapse = "} & \\color{mutedgray}\\textit{"),
                   "}")
  }
  paste0(line, " \\\\")
}

# ── Construit le corps LaTeX ───────────────────────────────────────────────────
build_latex_body <- function(df) {
  lines <- c()
  current_section <- ""
  
  # En-tête du tableau — 5 colonnes (Groupe, SNP, INDEL, SNP/nt, INDEL/nt)
  header <- paste0(
    "\\begin{tabular}{lrrrr}\n",
    "\\toprule\n",
    "\\textbf{Groupe} & \\textbf{SNP} & \\textbf{INDEL} & ",
    "\\textbf{SNP / nucl\\'{e}otide} & \\textbf{INDEL / nucl\\'{e}otide} \\\\\n",
    "\\midrule"
  )
  lines <- c(lines, header)
  
  for (i in seq_len(nrow(df))) {
    row     <- df[i, ]
    section <- paste(row$organisme, "--", row$region)
    org     <- row$organisme
    reg     <- row$region
    
    if (section != current_section) {
      if (current_section != "") lines <- c(lines, "\\midrule")
      
      # ── Choix du total de nucléotides selon organisme × région ──────────────
      total_nt <- switch(
        paste(org, reg),
        "Buchnera Codant"       = total_gene_buchnera,
        "Buchnera Non codant"   = total_noncodant_buchnera,
        "Serratia Codant"       = total_gene_serratia,
        "Serratia Non codant"   = total_noncodant_serratia,
        NA
      )
      
      # Libellé du total : "Nucl. codants : X" ou "Nucl. non codants : X"
      if (!is.na(total_nt)) {
        nt_label <- sprintf(
          "\\quad {\\small\\color{mutedgray}(nucl.~%s~:~%s~kb)}",
          ifelse(reg == "Codant", "codants", "non codants"),
          format_nucleotides(total_nt)
        )
      } else {
        nt_label <- ""
      }
      
      color <- ifelse(org == "Buchnera", "buchnera", "serratia")
      
      # En-tête de section sur 5 colonnes
      lines <- c(lines, sprintf(
        "\\multicolumn{5}{l}{\\cellcolor{%sbg}\\textbf{\\textcolor{%s}{%s --- %s}}%s} \\\\",
        color, color,
        latex_escape(org), latex_escape(reg),
        nt_label
      ))
      lines <- c(lines, "\\addlinespace[2pt]")
      current_section <- section
    }
    
    lines <- c(lines, make_row(row, row$groupe == "Global"))
  }
  
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  return(paste(lines, collapse = "\n"))
}

# ── Construit le document LaTeX complet ───────────────────────────────────────
latex_body <- build_latex_body(summary_table)

latex_doc <- sprintf('\\documentclass[11pt,a4paper]{article}
\\usepackage[margin=2cm]{geometry}
\\usepackage{booktabs}
\\usepackage{xcolor}
\\usepackage{colortbl}
\\usepackage{array}
\\usepackage[T1]{fontenc}
\\usepackage[utf8]{inputenc}
 
\\definecolor{buchnera}{HTML}{185FA5}
\\definecolor{serratia}{HTML}{0F6E56}
\\definecolor{buchnerabg}{HTML}{E6F1FB}
\\definecolor{serratiabg}{HTML}{E1F5EE}
\\definecolor{mutedgray}{HTML}{888888}
 
\\begin{document}
 
\\begin{center}
  {\\Large\\textbf{Mutations par organisme et r\\\'egion}}\\\\[4pt]
  {\\small\\color{mutedgray} Moyennes par groupe phylog\\\'en\\\'etique --- SNP, INDEL et ratios normalis\\\'es}\\\\[2pt]
  {\\small\\color{mutedgray} %s}
\\end{center}
 
\\vspace{10pt}
 
%s
 
\\end{document}',
                     format(Sys.Date(), "%d %B %Y"),
                     latex_body
)

# ── Compile en PDF ─────────────────────────────────────────────────────────────
tex_file <- "_tmp_tableau_mutations.tex"
writeLines(latex_doc, tex_file)

if (requireNamespace("tinytex", quietly = TRUE)) {
  tinytex::pdflatex(tex_file)
} else {
  system(paste("pdflatex", tex_file))
}

pdf_out <- sub(".tex", ".pdf", tex_file)
if (file.exists(pdf_out)) {
  file.rename(pdf_out, "tableau_mutations.pdf")
  message("PDF généré : tableau_mutations.pdf")
} else {
  message("Erreur : PDF non généré.")
}

file.remove(Sys.glob("_tmp_tableau_mutations.*"))



###### Traiter arbre #####

#produire le df propre associer a l'arbre 
echantillon_apisum_propre_serratia <- nettoyage_arbre (phylo_serratia, Echantillonnage_apisum, "read_file_name", "Country", "biotype")
echantillon_apisum_propre_buchnera <- nettoyage_arbre(phylo_buchnera, Echantillonnage_apisum, "read_file_name", "Country", "biotype")

phylo_serratia_localisation <- phylo_loc (phylo_serratia,echantillon_apisum_propre_serratia, "read_file_name", "Country", "Phylogénie Serratia")
phylo_serratia_plante <- phylo_plante (phylo_serratia, echantillon_apisum_propre_serratia, "biotype", "read_file_name")
phylo_serratia_loc_plante <- phylo_loc_plante (phylo_serratia, echantillon_apisum_propre_serratia, "read_file_name", "Country", "biotype", "Phylogénie Serratia")

phylo_buchnera_localisation <- phylo_loc (phylo_buchnera, echantillon_apisum_propre_buchnera, "read_file_name", "Country", "Phylogénie Buchnera")
phylo_buchnera_plante <- phylo_plante (phylo_buchnera, echantillon_apisum_propre_buchnera, "biotype", "read_file_name")
phylo_buchnera_loc_plante <- phylo_loc_plante (phylo_buchnera, echantillon_apisum_propre_buchnera, "read_file_name", "Country", "biotype", "Phylogénie Serratia")

# Export PDF
ggsave("phylo_buchnera_localisation.pdf",
       plot   = phylo_buchnera_localisation,
       width  = 14,
       height = 25,
       units  = "in")

ggsave("phylo_buchnera_plante.pdf",
       plot   = phylo_buchnera_plante,
       width  = 14,
       height = 25,
       units  = "in")

ggsave("phylo_buchnera_localisation_plante.pdf",
       plot   = phylo_buchnera_loc_plante,
       width  = 14,
       height = 25,
       units  = "in")

ggsave("phylo_serratia_localisation.pdf",
       plot   = phylo_serratia_localisation,
       width  = 14,
       height = 14,
       units  = "in")

ggsave("phylo_serratia_plante.pdf",
       plot   = phylo_serratia_plante,
       width  = 14,
       height = 14,
       units  = "in")

ggsave("phylo_serratia_localisation_plante.pdf",
       plot   = phylo_serratia_loc_plante,
       width  = 14,
       height = 14,
       units  = "in")

############### Coplot ###############

#Création de l'objet  
communs <- intersect(phylo_buchnera$tip.label, phylo_serratia$tip.label)
assoc   <- cbind(communs, communs)
p <- cophylo(phylo_buchnera, phylo_serratia, assoc = assoc, rotate = FALSE)

#Création de la figure
#Arbre Serratia reroot pour limiter croisement sur la figure (fonction rotate de cophylo pas active du a diff de feuille)
phylo_serratia_reroot <- read.tree("~/Documents/Moana/phylo_snippy/phylo_serratia_reroot.nwk")
pdf("cophylogenie_localisation.pdf", width = 20, height = 30)
plot_phylo_test <- cophylo_loc(phylo_buchnera,  phylo_serratia_reroot,
                               echantillon_apisum_propre_buchnera,
                               echantillon_apisum_propre_serratia,
                               "Country",
                               "Phylo commune")
dev.off()

#que individu avec serratia 
individus_serratia <- Echantillonnage_apisum_filtre$read_file_name[Echantillonnage_apisum_filtre$presence_serratia == "present"]
phylo_buchnera_only_with_serratia <- keep_tips(phylo_buchnera_incomplet, individus_serratia)
echantillonnage_apisum_propre_serratia <- nettoyage_arbre (phylo_serratia_reroot, Echantillonnage_apisum_filtre, "read_file_name", "Country", "biotype")
echantillonnage_apisum_propre_buchnera <- nettoyage_arbre(phylo_buchnera_only_with_serratia, Echantillonnage_apisum_filtre, "read_file_name", "Country", "biotype")
pdf("cophylogenie_localisation.pdf", width = 20, height = 20)
plot_phylo_test <- cophylo_loc(phylo_buchnera_only_with_serratia,  phylo_serratia_reroot,
                               echantillonnage_apisum_propre_buchnera,
                               echantillonnage_apisum_propre_serratia,
                               "Country",
                               "Phylo commune")
dev.off()

########### Reconstruction historique ##########

phylo_buchnera <- prepa_tree(phylo_buchnera, c("wamizawa", "CS041_S97", "CS042_S98") )
modele_reconstruction_zone_geo <- creation_modele (Echantillonnage_apisum, "read_file_name", "zone_geographique", phylo_buchnera)
applis_modele(modele_reconstruction_zone_geo, phylo_buchnera, Echantillonnage_apisum, "read_file_name", "zone_geographique", "reconstruction geographie", palette_resum_localisation, save = TRUE, "reconstruction_localisation.pdf")

modele_reconstruction_presence_serratia <- creation_modele (Echantillonnage_apisum, "read_file_name", "presence_serratia", phylo_buchnera)
applis_modele(modele_reconstruction_presence_serratia, phylo_buchnera, Echantillonnage_apisum, "read_file_name", "presence_serratia", "reconstruction presence serratia", palette_presence_serratia, save = TRUE, "reconstruction_presence.pdf")

modele_reconstruction_resum_biotype <- creation_modele (Echantillonnage_apisum, "read_file_name", "resum_biotype", phylo_buchnera)
applis_modele(modele_reconstruction_resum_biotype, phylo_buchnera, Echantillonnage_apisum, "read_file_name", "resum_biotype", "reconstruction biotype", palette_plante, save = FALSE)




vecteur_geo <- setNames(Echantillonnage_apisum$zone_geographique,
         Echantillonnage_apisum [["read_file_name"]])
mod_ARD <- fitMk(phylo_buchnera, vecteur_geo, model="ARD", fixedQ=NULL)
mod_ARD
ancr (mod_ARD, tree = phylo_buchnera, type = "joint")
vecteur_presence <- setNames(Echantillonnage_apisum$presence_serratia,
                        Echantillonnage_apisum [["read_file_name"]])
mod_ARD <- fitMk(phylo_buchnera, vecteur_presence, model="ARD", fixedQ=NULL)
mod_ARD
ancr (mod_ARD, type ="joint")


#### PARAFit #####

#  faire la matrice d'association : 

individus_serratia <- Echantillonnage_apisum_filtre$read_file_name[Echantillonnage_apisum_filtre$presence_serratia == "present"]
tous_individus <- Echantillonnage_apisum_filtre$read_file_name

matrice_association <- outer(tous_individus, individus_serratia , FUN = function(x, y) as.integer(x == y))
rownames(matrice_association) <- tous_individus
colnames(matrice_association) <-individus_serratia

# Matrices de distances patristiques
buchnera.dis <- as.dist(cophenetic(phylo_buchnera_incomplet))
serratia.dis <- as.dist(cophenetic(phylo_serratia))

# Mettre data dans le bon ordre (lignes = Buchnera, colonnes = Serratia)
matrice_parafit <- matrice_association[match(phylo_buchnera_incomplet$tip.label, 
                                             rownames(matrice_association)), ]
matrice_parafit <- t(matrice_parafit)
matrice_parafit <- matrice_parafit[match(phylo_serratia$tip.label, 
                                         rownames(matrice_parafit)), ]
matrice_parafit <- t(matrice_parafit)

# Lancer ParaFit
res_parafit <- parafit(buchnera.dis, serratia.dis, matrice_parafit, 
                       nperm = 999, test.links = TRUE,
                       correction = "cailliez")

# Résultats
print(res_parafit)
association <- res_parafit$link.table
association


### Comparaison génomique  ####

Echantillonnage_apisum_filtre <- Echantillonnage_apisum_filtre %>%
  mutate(type_sequence = case_when(
    qualite_sequence_buchnera == "bonne_qualite" & qualite_sequence_serratia == "bonne_qualite" ~ "Buchnera_Serratia",
    qualite_sequence_buchnera == "bonne_qualite" & qualite_sequence_serratia == "absent"   ~ "Buchnera",
  ))

Couple_serratia_buchnera <- Echantillonnage_apisum_filtre [Echantillonnage_apisum_filtre$type_sequence == "Buchnera_Serratia",]
Uniq_buchnera <- Echantillonnage_apisum_filtre [Echantillonnage_apisum_filtre$type_sequence == "Buchnera",]

table (Couple_serratia_buchnera$zone_geographique)

Couple_serratia_buchnera <- Couple_serratia_buchnera [Couple_serratia_buchnera$pourcentage_non_acgt_serratia <1,]
Couple_serratia_buchnera <- Couple_serratia_buchnera [Couple_serratia_buchnera$pourcentage_non_acgt_buchnera <1,]

couple_zone_chaude <- Couple_serratia_buchnera [Couple_serratia_buchnera$zone_geographique %in% c("Asie_central", "Asie_ouest"),]
couple_zone_stable <- Couple_serratia_buchnera [!Couple_serratia_buchnera$zone_geographique %in% c("Asie_central", "Asie_ouest"),]

sample(couple_zone_chaude$read_file_name, 4)
sample(couple_zone_stable$read_file_name, 3)

Uniq_buchnera <- Uniq_buchnera [Uniq_buchnera$pourcentage_non_acgt_buchnera <1,]

uniq_zone_chaude <- Uniq_buchnera [Uniq_buchnera$zone_geographique %in% c("Asie_central", "Asie_ouest"),]
uniq_zone_stable <- Uniq_buchnera [!Uniq_buchnera$zone_geographique %in% c("Asie_central", "Asie_ouest"),]

sample(uniq_zone_chaude$read_file_name, 4)
sample(uniq_zone_stable$read_file_name, 4)

### Allele privée 

test <- read.table("snps_buchnera_CS032_S32.tab", sep = "\t", header = TRUE)

dossier <- "~/Documents/Moana/taux_evolution/comparaison_kazak/buchnera_tab"
dossier_serratia <- "~/Documents/Moana/taux_evolution/comparaison_kazak/serratia_tab"

groupes <- list (
  BS_stable = c("snps_buchnera_CS052_S102.tab", "snps_buchnera_CS113_S124.tab", "snps_buchnera_CS170_S87.tab", "snps_buchnera_L15Sv14.tab"),
  B_stable = c("snps_buchnera_CS070_S58.tab", "snps_buchnera_CS137_S188.tab", "snps_buchnera_SRR11096966.tab", "snps_buchnera_L14Os02.tab"),
  BS_instable = c("snps_buchnera_CS032_S32.tab", "snps_buchnera_CS237.tab", "snps_buchnera_CS159_S148.tab", "snps_buchnera_CS099_S71.tab"),
  B_instable = c("snps_buchnera_CS226.tab", "snps_buchnera_CS166_S154.tab", "snps_buchnera_CS084_S65.tab", "snps_buchnera_CS180_S92.tab")
)

groupes <- list (
  BS_stable = c("snps_serratia_CS052_S102.tab", "snps_serratia_CS113_S124.tab", "snps_serratia_CS170_S87.tab", "snps_serratia_L15Sv14.tab"),
  BS_instable = c("snps_serratia_CS032_S32.tab", "snps_serratia_CS237.tab", "snps_serratia_CS159_S148.tab", "snps_serratia_CS099_S71.tab")
)

n_ind_par_groupe <- 4


# 2. LECTURE ET FILTRAGE DES FICHIERS

tous_snp <- data.frame()

for (nom_groupe in names(groupes)) {
  fichiers <- groupes[[nom_groupe]]
  
  for (i in seq_along(fichiers)) {
    fichier <- fichiers[i]
    chemin  <- file.path(dossier, fichier)
    
    df <- read.delim(chemin, header = TRUE, stringsAsFactors = FALSE)
    
    # Filtrer : TYPE == "snp" et GENE non vide
    df <- df[df$TYPE == "snp" & !is.na(df$GENE) & df$GENE != "", ]
    
    if (nrow(df) == 0) next
    
    # Garder uniquement les colonnes utiles
    df <- df[, c("POS", "REF", "ALT", "GENE")]
    
    # Ajouter les colonnes individu et groupe
    df$individu <- fichier          # nom du fichier = identifiant individu
    df$groupe   <- nom_groupe       # nom du groupe depuis la liste
    
    tous_snp <- rbind(tous_snp, df)
  }
}

print(colnames(tous_snp))
print(head(tous_snp))

# 3. CALCUL DES FRÉQUENCES PAR GROUPE

# Dédoublonner : un individu ne compte qu'une fois par mutation
mutations_uniques <- unique(tous_snp[, c("groupe", "individu", "POS", "REF", "ALT", "GENE")])

# Cle regroupe en 1 var groupe__POS__ALT (sert a compter mut commune a 1 groupe avec compte qui denombre le nombre de cles identique)
cles <- paste(mutations_uniques$groupe, mutations_uniques$POS, mutations_uniques$ALT, sep = "__")
compte <- tapply(mutations_uniques$individu, cles, function(x) length(unique(x)))

# Reconstruire un dataframe de fréquences
mutations_uniques$cle <- paste(mutations_uniques$groupe, mutations_uniques$POS, mutations_uniques$ALT, sep = "__")
mutations_uniques$n_individus <- as.integer(compte[mutations_uniques$cle])
mutations_uniques$frequence   <- mutations_uniques$n_individus / n_ind_par_groupe
mutations_uniques$cle         <- NULL

# Exclure fréquence == 0.25 (1 seul individu)
freq_par_groupe <- mutations_uniques[mutations_uniques$frequence > 0.25, ]


# 4. IDENTIFIER LES MUTATIONS EXCLUSIVES OU PARTAGÉES


# Compter combien de groupes ont chaque mutation (POS + ALT)
cle_mutation <- paste(freq_par_groupe$POS, freq_par_groupe$ALT, sep = "__")
n_groupes_par_mutation <- tapply(freq_par_groupe$groupe, cle_mutation, function(x) length(unique(x)))

freq_par_groupe$cle_mutation  <- cle_mutation
freq_par_groupe$n_groupes     <- as.integer(n_groupes_par_mutation[cle_mutation])


# 5. DATAFRAME 1 — Mutations exclusives à un groupe


df1_exclusif <- freq_par_groupe[freq_par_groupe$n_groupes == 1, 
                                c("groupe", "POS", "REF", "ALT", "GENE", "frequence")]
df1_exclusif <- df1_exclusif[order(df1_exclusif$groupe, df1_exclusif$GENE, df1_exclusif$POS), ]
rownames(df1_exclusif) <- NULL

cat("=== DATAFRAME 1 : Mutations exclusives à un groupe ===\n")
print(df1_exclusif)


# 6. DATAFRAME 2 — Mutations partagées entre plusieurs groupes


df2_long <- freq_par_groupe[freq_par_groupe$n_groupes > 1, 
                            c("groupe", "POS", "REF", "ALT", "GENE", "frequence")]

# Pivoter en format large : une colonne de fréquence par groupe
noms_groupes  <- names(groupes)
mutations_shared <- unique(df2_long[, c("POS", "REF", "ALT", "GENE")])

for (g in noms_groupes) {
  col_name <- paste0("freq_", g)
  freq_g   <- df2_long[df2_long$groupe == g, c("POS", "ALT", "frequence")]
  
  mutations_shared[[col_name]] <- 0  # valeur par défaut
  
  for (j in seq_len(nrow(mutations_shared))) {
    match_idx <- which(freq_g$POS == mutations_shared$POS[j] & 
                         freq_g$ALT == mutations_shared$ALT[j])
    if (length(match_idx) > 0) {
      mutations_shared[[col_name]][j] <- freq_g$frequence[match_idx[1]]
    }
  }
}

# Calculer la fréquence totale
freq_cols <- paste0("freq_", noms_groupes)
mutations_shared$freq_totale <- rowSums(mutations_shared[, freq_cols])

# Exclure si freq totale == 0.25
df2_partage <- mutations_shared[mutations_shared$freq_totale > 0.25, ]
df2_partage <- df2_partage[order(df2_partage$GENE, df2_partage$POS), ]
rownames(df2_partage) <- NULL

cat("\n=== DATAFRAME 2 : Mutations partagées entre groupes ===\n")
print(df2_partage)


# 7. EXPORT

write.csv(df1_exclusif, "mutations_exclusives.csv", row.names = FALSE)
write.csv(df2_partage,  "mutations_partagees.csv",  row.names = FALSE)






### RDA ######

data_clim <- read.csv2("~/Documents/Moana/presabs_serratia.csv", head=TRUE, sep=",")
str (data_clim)

cols_latitude_longitude <- c("x", "y")

cols_environnemental <- data_clim %>%
  select(where(is.character) & !all_of(cols_latitude_longitude)) %>%
  colnames()

#transformation numerique
data_clim <- data_clim %>%
  mutate(across(all_of(cols_environnemental), ~ as.numeric(.)))
data_clim$presence_serratia <- as.numeric(data_clim$presence_serratia)

#transformation a log
cols_log <- c("bio01", "bio06", "bio09", "bio11")  
cols_a_logger <- setdiff(cols_environnemental, cols_log)
for (col in cols_a_logger) {
  data_clim[[paste0("log_", col)]] <- log1p(data_clim[[col]])
}

# vecteur de var clim
vars_clim_log <- grep("^log_", colnames(data_clim), value = TRUE)
vars_clim_log <- c (vars_clim_log, "bio01", "bio06", "bio09", "bio11")

# RDA
Y <- data_clim[, "presence_serratia", drop = FALSE]
X <- data_clim[, vars_clim_log]
rda_clim <- rda(Y ~ ., data = X, scale = TRUE)

summary(rda_clim)

# plot
scores_sites <- as.data.frame(scores(rda_clim, display = "sites"))
scores_sites$presence_serratia <- data_clim$presence_serratia  

plot(rda_clim, type = "n", scaling = 2) 
points(scores_sites$RDA1, scores_sites$PC1,
       col = ifelse(scores_sites$presence_serratia == 1, "red", "blue"),
       pch = 16, cex = 1.2)
text(rda_clim, display = "bp", col = "darkgreen", arrow.mul = 1)
# Légende
legend("topright",
       legend = c("Présence Serratia", "Absence Serratia"),
       col = c("red", "blue"),
       pch = 16)

#proportion variance expliquée
prop_explained <- rda_clim$CCA$tot.chi / (rda_clim$CCA$tot.chi + rda_clim$CA$tot.chi)
paste("Proportion de variance expliquée par les variables :", 
      round(prop_explained * 100, 2), "%")

cor_matrix <- cor(data_clim[, vars_clim_log], use = "complete.obs")
heatmap(cor_matrix, symm = TRUE)


#pas correlé 


str (data_clim)

cols_environnemental <- c("bio01", "bio02", "bio03", "bio05", "bio06", "bio8", "npp")

#transformation numerique
data_clim <- data_clim %>%
  mutate(across(all_of(cols_environnemental), ~ as.numeric(.)))
data_clim$presence_serratia <- as.numeric(data_clim$presence_serratia)

#transformation a log
data_clim[["log_npp"]] <- log1p(data_clim[["npp"]])

# vecteur de var clim
vars_clim_log <- c ("log_npp", "bio01", "bio06", "bio09", "bio11")

# RDA
Y <- data_clim[, "presence_serratia", drop = FALSE]
X <- data_clim[, vars_clim_log]
rda_clim <- rda(Y ~ ., data = X, scale = TRUE)

summary(rda_clim)

# plot
scores_sites <- as.data.frame(scores(rda_clim, display = "sites"))
scores_sites$presence_serratia <- data_clim$presence_serratia  

plot(rda_clim, type = "n", scaling = 2) 
points(scores_sites$RDA1, scores_sites$PC1,
       col = ifelse(scores_sites$presence_serratia == 1, "red", "blue"),
       pch = 16, cex = 1.2)
text(rda_clim, display = "bp", col = "darkgreen", arrow.mul = 1)
# Légende
legend("topright",
       legend = c("Présence Serratia", "Absence Serratia"),
       col = c("red", "blue"),
       pch = 16)

#proportion variance expliquée
prop_explained <- rda_clim$CCA$tot.chi / (rda_clim$CCA$tot.chi + rda_clim$CA$tot.chi)
paste("Proportion de variance expliquée par les variables :", 
      round(prop_explained * 100, 2), "%")

cor_matrix <- cor(data_clim[, vars_clim_log], use = "complete.obs")
heatmap(cor_matrix, symm = TRUE)
















####ŋßðdr#####

Echantillonnage_apisum <- read.csv2("~/Documents/Moana/Echantillonnage_apisum.csv", head=TRUE, sep=",")
Echantillonnage_apisum <- Echantillonnage_apisum [Echantillonnage_apisum$read_file_name != "CS189_S180",]
Echantillonnage_apisum <- Echantillonnage_apisum[Echantillonnage_apisum$qualite_sequence_buchnera != "mauvaise_qualite", ]

phylo_buchnera <- read.tree("~/Documents/Moana/phylo_snippy/phylo_buchnera_without_outgroup.nwk")
individus_a_supprimer <- c("wamizawa", "CS041_S97", "CS042_S98")
phylo_buchnera_filtre <- drop.tip(phylo_buchnera, individus_a_supprimer)

echantillon_apisum_propre_buchnera_filtre <- nettoyage_arbre(phylo_buchnera_filtre, Echantillonnage_apisum, "read_file_name", "Country", "biotype")
phylo_serratia_reroot <- read.tree("~/Documents/Moana/phylo_snippy/phylo_serratia_reroot.nwk")

cophylo_loc <- function(tree_left, tree_right,
                        metadata_left, metadata_right,
                        col_localisation,
                        titre = "Cophylogénie – Localisation",
                        output_file = "cophylogenie.png",
                        width = 2800, height = 3200, res = 300) {
  
  # ── 1. Arbres COMPLETS ────────────────────────────────────────────────────
  tl <- tree_left
  tr <- tree_right
  if (!is.rooted(tl)) tl <- unroot(tl)
  if (!is.rooted(tr)) tr <- unroot(tr)
  
  tips_communs <- intersect(tl$tip.label, tr$tip.label)
  assoc <- cbind(tips_communs, tips_communs)
  obj <- cophylo(tl, tr, assoc = assoc, rotate = FALSE)
  
  # ── 2. Couleurs ───────────────────────────────────────────────────────────
  get_cols <- function(labels, metadata) {
    cols <- palette_localisation[metadata[[col_localisation]][
      match(labels, metadata$read_file_name)]]
    cols[is.na(cols)] <- "grey85"
    cols
  }
  
  cols_left  <- get_cols(obj$trees[[1]]$tip.label, metadata_left)
  cols_right <- get_cols(obj$trees[[2]]$tip.label, metadata_right)
  
  # ── 3. Sauvegarde PNG ─────────────────────────────────────────────────────
  png(output_file, width = width, height = height, res = res)
  
  # ── 4. Plot ───────────────────────────────────────────────────────────────
  plot(obj,
       fsize    = 0.0001,
       link.col = "black",   # ← noir
       link.lwd = 1,       # ← un peu plus épais pour mieux voir
       ftype    = "off",
       points   = FALSE,
       tip.rad  = 0)
  
  title(titre, line = 1)
  
  # ── 5. Coordonnées ────────────────────────────────────────────────────────
  pp <- get("last_plot.cophylo", envir = .PlotPhyloEnv)
  
  n_left  <- Ntip(obj$trees[[1]])
  n_right <- Ntip(obj$trees[[2]])
  
  x_left  <- pp[[1]]$xx[seq_len(n_left)]
  y_left  <- pp[[1]]$yy[seq_len(n_left)]
  x_right <- pp[[2]]$xx[seq_len(n_right)]
  y_right <- pp[[2]]$yy[seq_len(n_right)]
  
  labels_left  <- obj$trees[[1]]$tip.label
  labels_right <- obj$trees[[2]]$tip.label
  
  # ── 6. Calcul dynamique de cex ────────────────────────────────────────────
  usr       <- par("usr")
  height_u  <- usr[4] - usr[3]
  n_max     <- max(n_left, n_right)
  h_per_tip <- height_u / n_max
  
  # ── 6. Calcul dynamique de cex ────────────────────────────────────────────
  cex_label <- 0.5
  pad <- strwidth("i", cex = cex_label) * 0.4   # ← pad actuel (petit)
  tip_offset <- strwidth("i", cex = cex_label) * 2  # ← nouveau : espace branche → label
  
  # ── 7. Labels ─────────────────────────────────────────────────────────────
  draw_labels <- function(labels, x_tips, y_tips, cols, side = "left") {
    for (i in seq_along(labels)) {
      w <- strwidth(labels[i],  cex = cex_label)
      h <- strheight(labels[i], cex = cex_label) * 1.3
      
      if (side == "left") {
        xl <- x_tips[i] + tip_offset*2          # ← décalage depuis la branche
        xr <- xl + w + pad
        rect(xl, y_tips[i] - h/2, xr, y_tips[i] + h/2,
             col = cols[i], border = NA)
        text(xl + pad/2, y_tips[i], labels[i],
             cex = cex_label, adj = c(0, 0.5), col = "black")
      } else {
        xr <- x_tips[i] - tip_offset          # ← décalage depuis la branche
        xl <- xr - w - pad
        rect(xl, y_tips[i] - h/2, xr, y_tips[i] + h/2,
             col = cols[i], border = NA)
        text(xr - pad/2, y_tips[i], labels[i],
             cex = cex_label, adj = c(1, 0.5), col = "black")
      }
    }
  }
  
  draw_labels(labels_left,  x_left,  y_left,  cols_left,  side = "left")
  draw_labels(labels_right, x_right, y_right, cols_right, side = "right")
  
  # ── 8. Légende ────────────────────────────────────────────────────────────
  loc_presentes <- sort(unique(na.omit(c(
    metadata_left[[col_localisation]][match(tips_communs, metadata_left$read_file_name)],
    metadata_right[[col_localisation]][match(tips_communs, metadata_right$read_file_name)]
  ))))
  
  legend("topleft",
         legend     = loc_presentes,
         fill       = palette_localisation[loc_presentes],
         border     = NA,
         ncol       = ceiling(length(loc_presentes) / 5),
         cex        = 1,
         title      = "Localisation",
         title.font = 2,
         bty        = "n",
         xpd        = TRUE)
  
  dev.off()
  message("Figure sauvegardée : ", output_file)
}

cophylo_loc(
  tree_left        = phylo_buchnera,
  tree_right       = phylo_serratia_reroot,
  metadata_left    = echantillon_apisum_propre_buchnera,
  metadata_right   = echantillon_apisum_propre_serratia,
  col_localisation = "Country",
  output_file      = "cophylogenie_avec_racine.png",
  width  = 15000,
  height = 18500,
  res    = 1000
)

phylo_loc_plante <- function(tree,
                             metadata_clean,
                             col_taxon,
                             col_localisation,
                             col_plante,
                             titre = "Phylogénie – Localisation & Plante hôte",
                             save_plot = FALSE,
                             filename  = "phylogenie.png",
                             width     = 3000,
                             height    = 4000,
                             res       = 300
) {
  couleurs_texte_loc <- setNames(
    rep("black", length(palette_localisation)),
    names(palette_localisation)
  )
  
  p_test <- ggtree(tree) %<+% metadata_clean
  
  p <- ggtree(tree,
              layout = "rectangular",
              color  = "grey40",
              size   = 0.4) %<+% metadata_clean +
    
    geom_treescale(
      x        = min(p_test$data$x, na.rm = TRUE),
      y        = -2,
      width    = 0.0001,
      offset   = 0.5,
      color    = "black",
      linesize = 0.5,
      fontsize = 3
    ) +
    
    geom_tiplab(
      aes(label = label,
          fill  = .data[[col_localisation]],
          color = .data[[col_localisation]]),
      geom          = "label",
      size          = 2,
      offset        = 0.000000001,
      align         = TRUE,
      linetype      = NA,
      label.padding = unit(0.05, "lines"),
      label.r       = unit(0, "lines")
    ) +
    scale_fill_manual(
      name     = "Localisation",
      values   = palette_localisation,
      na.value = "white"
    ) +
    scale_color_manual(
      name   = "Localisation",
      values = couleurs_texte_loc,
      guide  = "none"
    ) +
    
    ggnewscale::new_scale_fill() +
    
    geom_fruit(
      data    = NULL,
      geom    = geom_point,
      mapping = aes(fill  = .data[[col_plante]],
                    shape = .data[[col_plante]]),
      size    = 2.5,
      stroke  = 0.3,
      color   = "black",
      offset  = 0.09,
      axis.params = list(axis = "none")
    ) +
    scale_fill_manual(
      name     = "Plante hôte",
      values   = palette_plante,
      na.value = "grey60"
    ) +
    scale_shape_manual(
      name     = "Plante hôte",
      values   = formes_plante,
      na.value = 21
    ) +
    
    theme(
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 11),
      legend.box       = "vertical",
      legend.spacing.y = unit(0.35, "cm")
    ) +
    ggtitle(titre)
  
  # ── Affichage interactif ───────────────────────────────────────────────────
  print(p)
  
  # ── Export PNG via png() ───────────────────────────────────────────────────
  if (save_plot) {
    png(filename = filename,
        width    = width,
        height   = height,
        res      = res)
    print(p)
    dev.off()
    message("Plot sauvegardé : ", filename,
            " (", width, "×", height, " px, ", res, " dpi)")
  }
  
  invisible(p)
}


phylo_serratia_loc_plante <- phylo_loc_plante (phylo_serratia, echantillon_apisum_propre_serratia, "read_file_name", "Country", "biotype", "Phylogénie Serratia")
phylo_buchnera_loc_plante <- phylo_loc_plante (phylo_buchnera, echantillon_apisum_propre_buchnera, "read_file_name", "Country", "biotype", "Phylogénie Serratia")

phylo_buchnera_loc_plante <- phylo_loc_plante(
  phylo_buchnera, echantillon_apisum_propre_buchnera,
  "read_file_name", "Country", "biotype", "Phylogénie Buchnera",
  save_plot = TRUE,
  filename  = "phylo_buchnera_avec_racine.png",
  width     = 9000,
  height    = 9999,
  res       = 750
)

phylo_serratia_loc_plante <- phylo_loc_plante(
  phylo_serratia, echantillon_apisum_propre_serratia,
  "read_file_name", "Country", "biotype", "Phylogénie Serratia",
  save_plot = TRUE,
  filename  = "phylo_serratia.png",
  width     = 9000,
  height    = 9999,
  res       = 750
)


cophylo_loc_plante <- function(tree_left, tree_right,
                               metadata_left, metadata_right,
                               col_localisation,
                               col_plante,
                               titre = "Cophylogénie – Localisation & Plante hôte",
                               output_file = "cophylogenie.png",
                               width = 2800, height = 3200, res = 300) {
  
  # ── 1. Arbres COMPLETS ────────────────────────────────────────────────────
  tl <- tree_left
  tr <- tree_right
  if (!is.rooted(tl)) tl <- unroot(tl)
  if (!is.rooted(tr)) tr <- unroot(tr)
  
  tips_communs <- intersect(tl$tip.label, tr$tip.label)
  assoc <- cbind(tips_communs, tips_communs)
  obj <- cophylo(tl, tr, assoc = assoc, rotate = FALSE)
  
  # ── 2. Couleurs localisation ──────────────────────────────────────────────
  get_cols <- function(labels, metadata) {
    cols <- palette_localisation[metadata[[col_localisation]][
      match(labels, metadata$read_file_name)]]
    cols[is.na(cols)] <- "grey85"
    cols
  }
  
  cols_left  <- get_cols(obj$trees[[1]]$tip.label, metadata_left)
  cols_right <- get_cols(obj$trees[[2]]$tip.label, metadata_right)
  
  # ── 3. Couleurs & formes plante hôte ─────────────────────────────────────
  get_plante <- function(labels, metadata) {
    plantes <- metadata[[col_plante]][match(labels, metadata$read_file_name)]
    list(
      cols   = ifelse(is.na(plantes), "grey60", palette_plante[plantes]),
      formes = ifelse(is.na(plantes), 21,       formes_plante[plantes])
    )
  }
  
  plante_left  <- get_plante(obj$trees[[1]]$tip.label, metadata_left)
  plante_right <- get_plante(obj$trees[[2]]$tip.label, metadata_right)
  
  # ── 4. Sauvegarde PNG ─────────────────────────────────────────────────────
  png(output_file, width = width, height = height, res = res)
  
  # ── 5. Plot cophylo de base ───────────────────────────────────────────────
  plot(obj,
       fsize    = 0.0001,
       link.col = "black",
       link.lwd = 1.4,
       lwd = 1.4,
       ftype    = "off",
       points   = FALSE,
       tip.rad  = 0)
  
  title(titre, line = 1)
  
  # ── 6. Coordonnées des feuilles ───────────────────────────────────────────
  pp <- get("last_plot.cophylo", envir = .PlotPhyloEnv)
  
  n_left  <- Ntip(obj$trees[[1]])
  n_right <- Ntip(obj$trees[[2]])
  
  x_left  <- pp[[1]]$xx[seq_len(n_left)]
  y_left  <- pp[[1]]$yy[seq_len(n_left)]
  x_right <- pp[[2]]$xx[seq_len(n_right)]
  y_right <- pp[[2]]$yy[seq_len(n_right)]
  
  labels_left  <- obj$trees[[1]]$tip.label
  labels_right <- obj$trees[[2]]$tip.label
  
  # ── 7. Paramètres texte ───────────────────────────────────────────────────
  cex_label  <- 0.5
  pad        <- strwidth("i", cex = cex_label) * 0.4
  tip_offset <- strwidth("i", cex = cex_label) * 2
  symbol_size <- strheight("i", cex = cex_label) * 1.5
  
  # ── Calcul des x d'alignement ─────────────────────────────────────────────
  x_align_left  <- max(x_left)  + tip_offset * 2          # bord gauche aligné
  x_align_right <- min(x_right) - tip_offset               # bord droit aligné
  
  # ── 8. Labels + symboles plante ───────────────────────────────────────────
  draw_labels <- function(labels, x_tips, y_tips, cols_loc,
                          cols_plante, formes_plante_vec, side = "left") {
    for (i in seq_along(labels)) {
      w <- strwidth(labels[i],  cex = cex_label)
      h <- strheight(labels[i], cex = cex_label) * 1.3
      
      if (side == "left") {
        xl <- x_align_left            # ← tous alignés sur le même x
        xr <- xl + w + pad
        rect(xl, y_tips[i] - h/2, xr, y_tips[i] + h/2,
             col = cols_loc[i], border = NA)
        text(xl + pad/2, y_tips[i], labels[i],
             cex = cex_label, adj = c(0, 0.5), col = "black")
        
        x_sym <- xr + symbol_size
        points(x_sym, y_tips[i],
               pch = formes_plante_vec[i],
               bg  = cols_plante[i],
               col = "black",
               cex = 0.8)
        
      } else {
        xr <- x_align_right           # ← tous alignés sur le même x
        xl <- xr - w - pad
        rect(xl, y_tips[i] - h/2, xr, y_tips[i] + h/2,
             col = cols_loc[i], border = NA)
        text(xr - pad/2, y_tips[i], labels[i],
             cex = cex_label, adj = c(1, 0.5), col = "black")
        
        x_sym <- xl - symbol_size
        points(x_sym, y_tips[i],
               pch = formes_plante_vec[i],
               bg  = cols_plante[i],
               col = "black",
               cex = 0.8)
      }
    }
  }
  
  # ← CES DEUX LIGNES MANQUENT :
  draw_labels(labels_left,  x_left,  y_left,
              cols_left,  plante_left$cols,  plante_left$formes,  side = "left")
  draw_labels(labels_right, x_right, y_right,
              cols_right, plante_right$cols, plante_right$formes, side = "right")
  
  # ── 9. Légende localisation ───────────────────────────────────────────────
  loc_presentes <- sort(unique(na.omit(c(
    metadata_left[[col_localisation]][match(tips_communs, metadata_left$read_file_name)],
    metadata_right[[col_localisation]][match(tips_communs, metadata_right$read_file_name)]
  ))))
  
  leg_loc <- legend("topleft",
         legend     = loc_presentes,
         fill       = palette_localisation[loc_presentes],
         border     = NA,
         ncol       = ceiling(length(loc_presentes) / 5),
         cex        = 0.8,
         title      = "Localisation",
         title.font = 2,
         bty        = "n",
         xpd        = TRUE)
  
  # ── 10. Légende plante hôte ───────────────────────────────────────────────
  plantes_presentes <- sort(unique(na.omit(c(
    metadata_left[[col_plante]][match(tips_communs, metadata_left$read_file_name)],
    metadata_right[[col_plante]][match(tips_communs, metadata_right$read_file_name)]
  ))))
  
  legend(x      = leg_loc$rect$left,              # ← même x que localisation
         y      = leg_loc$rect$top - leg_loc$rect$h - strheight("i") * 2,  # ← juste en dessous
         legend = plantes_presentes,
         pch    = formes_plante[plantes_presentes],
         pt.bg  = palette_plante[plantes_presentes],
         col    = "black",
         pt.cex = 1.2,
         cex    = 0.8,
         title      = "Plante hôte",
         title.font = 2,
         bty        = "n",
         xpd        = TRUE)
  
  dev.off()
  message("Figure sauvegardée : ", output_file)
}

cophylo_loc_plante(
  tree_left        = phylo_buchnera,
  tree_right       = phylo_serratia_reroot,
  metadata_left    = echantillon_apisum_propre_buchnera,
  metadata_right   = echantillon_apisum_propre_serratia,
  col_localisation = "Country",
  col_plante = "biotype",
  titre = "Cophylogénie – Localisation & Plante hôte",
  output_file      = "cophylogenie_testbranche_moche.png",
  width  = 9000,
  height = 9999,
  res    = 750
)



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

# Fonction : établir les modele pour faire de la reconstruction le long d'un arbre 
creation_modele <- function (df, 
                             col_name, 
                             character, 
                             tree){
  
  tree_tips <- tree$tip.label
  # Aligner le dataframe sur l'ordre des tips de l'arbre
  sortedData <- df[match(tree_tips, df[[col_name]]), ]
  
  chara_reconstruit <- sortedData[[character]]
  names(chara_reconstruit) <- tree_tips
  
  # Reconstruction ancestrale avec les 3 modèles
  ERreconstruction  <- ace(chara_reconstruit, tree, type = "discrete", model = "ER")
  SYMreconstruction <- ace(chara_reconstruit, tree, type = "discrete", model = "SYM")
  ARDreconstruction <- ace(chara_reconstruit, tree, type = "discrete", model = "ARD")
  
  n_etats <- length(unique(chara_reconstruit))
  
  df_ER_ARD  <- n_etats * (n_etats - 1) - 1
  df_ER_SYM  <- n_etats * (n_etats - 1) / 2 - 1
  df_SYM_ARD <- n_etats * (n_etats - 1) / 2
  
  # Calcul des p-values
  pval_ER_ARD  <- 1 - pchisq(2 * abs(ERreconstruction$loglik  - ARDreconstruction$loglik),  df_ER_ARD)
  pval_ER_SYM  <- 1 - pchisq(2 * abs(ERreconstruction$loglik  - SYMreconstruction$loglik),  df_ER_SYM)
  pval_SYM_ARD <- 1 - pchisq(2 * abs(SYMreconstruction$loglik - ARDreconstruction$loglik),  df_SYM_ARD)
  
  cat("\nTests LRT (p-value < 0.05 = le modèle complexe est significativement meilleur) :\n")
  cat("  ER vs ARD  :", pval_ER_ARD,  "\n")
  cat("  ER vs SYM  :", pval_ER_SYM,  "\n")
  cat("  SYM vs ARD :", pval_SYM_ARD, "\n")
  
  # Sélection du meilleur modèle
  if (pval_ER_ARD < 0.05 && pval_SYM_ARD < 0.05) {
    best_model <- ARDreconstruction
    best_name  <- "ARD"
  } else if (pval_ER_SYM < 0.05 && pval_SYM_ARD >= 0.05) {
    best_model <- SYMreconstruction
    best_name  <- "SYM"
  } else {
    best_model <- ERreconstruction
    best_name  <- "ER"
  }
  
  message("\nMeilleur modèle sélectionné : ", best_name)
  
  return(best_model)
}

# Fonction : Figure avec application d'un modele pour reconstruire l'état au noeud

applis_modele <- function (best_reconstruction, 
                           tree, df, 
                           col_name, character,
                           fig_title,
                           palette_fixe,
                           save     = FALSE, 
                           filename = "reconstruction.png",
                           width    = 2800,
                           height   = 3200, 
                           res      = 300
){
  
  tree_tips  <- tree$tip.label
  sortedData <- df[match(tree_tips, df[[col_name]]), ]
  
  chara_reconstruit <- sortedData[[character]]
  names(chara_reconstruit) <- tree_tips
  
  etats    <- sort(unique(chara_reconstruit))
  couleurs <- palette_fixe[etats]
  
  if (save) png(filename, width = width, height = height, res = res)
  
  # Labels désactivés dans plotTree, redessinés manuellement ensuite
  plotTree(tree, lwd = 2, setEnv = TRUE, fsize = 0.0001, ftype = "off")
  
  pp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xlim <- pp$x.lim
  ylim <- pp$y.lim
  
  # ── Labels alignés ────────────────────────────────────────────────────────
  n_tips    <- length(tree_tips)
  x_tips    <- pp$xx[seq_len(n_tips)]
  y_tips    <- pp$yy[seq_len(n_tips)]
  
  cex_label <- 0.4
  x_align   <- max(x_tips) + strwidth("i", cex = cex_label) * 2
  
  segments(x_tips, y_tips, x_align, y_tips,
           col = "grey70", lwd = 0.3, lty = "dotted")
  
  text(x_align, y_tips,
       labels = tree_tips,
       cex    = cex_label,
       adj    = c(0, 0.5),
       xpd    = TRUE)
  # ─────────────────────────────────────────────────────────────────────────
  
  add.scale.bar(
    x      = xlim[2] * 0.85,
    y      = ylim[1],
    length = 0.0001,
    lwd    = 2,
    cex    = 0.7
  )
  
  # Noeuds internes : camemberts des probabilités ancestrales (non modifiés)
  n_nodes <- tree$Nnode
  nodelabels(
    node    = (n_tips + 1):(n_tips + n_nodes),
    pie     = best_reconstruction$lik.anc,
    piecol  = couleurs,
    cex     = 0.3
  )
  
  # Tips : état observé
  tiplabels(
    pie    = to.matrix(chara_reconstruit, etats),
    piecol = couleurs,
    cex    = 0.2
  )
  
  legend("bottomleft", legend = etats, fill = couleurs, cex = 2, fig_title)
  
  if (save) {
    dev.off()
    message("Figure sauvegardée : ", filename)
  }
}

applis_modele_sans_outgroup <- function (best_reconstruction, 
                                         tree,                        # arbre pour la reconstruction (complet)
                                         tree_display = NULL,         # arbre pour l'affichage (filtré), optionnel
                                         df, 
                                         col_name, character,
                                         fig_title,
                                         palette_fixe,
                                         save     = FALSE, 
                                         filename = "reconstruction.png",
                                         width    = 2800,
                                         height   = 3200, 
                                         res      = 300) {
  
  # Arbre d'affichage = arbre complet si non spécifié
  tree_plot <- if (!is.null(tree_display)) tree_display else tree
  
  tree_tips      <- tree$tip.label
  tree_plot_tips <- tree_plot$tip.label
  
  sortedData <- df[match(tree_tips, df[[col_name]]), ]
  chara_reconstruit <- sortedData[[character]]
  names(chara_reconstruit) <- tree_tips
  
  etats    <- sort(unique(chara_reconstruit))
  couleurs <- palette_fixe[etats]
  
  # ── Mapping des noeuds : arbre complet → arbre filtré ────────────────────
  if (!is.null(tree_display)) {
    
    # Pour chaque noeud interne de l'arbre filtré,
    # on identifie les tips descendants puis on trouve
    # le noeud correspondant dans l'arbre complet
    n_nodes_plot <- tree_plot$Nnode
    n_tips_plot  <- length(tree_plot_tips)
    lik_display  <- matrix(NA, 
                           nrow = n_nodes_plot, 
                           ncol = length(etats),
                           dimnames = list(NULL, etats))
    
    for (i in seq_len(n_nodes_plot)) {
      node_id <- n_tips_plot + i
      
      # Tips descendants dans l'arbre filtré
      desc_tips <- tree_plot_tips[
        phangorn::Descendants(tree_plot, node_id, type = "tips")[[1]]
      ]
      
      # Noeud correspondant dans l'arbre complet
      # = MRCA des mêmes tips dans l'arbre complet
      if (length(desc_tips) >= 2) {
        node_complet <- getMRCA(tree, desc_tips)
      } else {
        node_complet <- which(tree$tip.label == desc_tips[1])
      }
      
      # Index dans lik.anc (noeuds internes seulement, donc - n_tips_complet)
      n_tips_complet <- length(tree_tips)
      idx <- node_complet - n_tips_complet
      
      if (!is.na(idx) && idx > 0 && idx <= nrow(best_reconstruction$lik.anc)) {
        lik_display[i, ] <- best_reconstruction$lik.anc[idx, etats]
      }
    }
    
  } else {
    lik_display <- best_reconstruction$lik.anc
  }
  # ─────────────────────────────────────────────────────────────────────────
  
  if (save) png(filename, width = width, height = height, res = res)
  
  plotTree(tree_plot, lwd = 2, setEnv = TRUE, fsize = 0.0001, ftype = "off")
  
  pp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xlim <- pp$x.lim
  ylim <- pp$y.lim
  
  # Labels alignés
  n_tips_plot <- length(tree_plot_tips)
  x_tips      <- pp$xx[seq_len(n_tips_plot)]
  y_tips      <- pp$yy[seq_len(n_tips_plot)]
  
  cex_label <- 0.4
  x_align   <- max(x_tips) + strwidth("i", cex = cex_label) * 2
  
  segments(x_tips, y_tips, x_align, y_tips,
           col = "grey70", lwd = 0.3, lty = "dotted")
  
  text(x_align, y_tips,
       labels = tree_plot_tips,
       cex    = cex_label,
       adj    = c(0, 0.5),
       xpd    = TRUE)
  
  add.scale.bar(
    x      = xlim[2] * 0.85,
    y      = ylim[1],
    length = 0.0001,
    lwd    = 2,
    cex    = 0.7
  )
  
  n_nodes_plot <- tree_plot$Nnode
  nodelabels(
    node    = (n_tips_plot + 1):(n_tips_plot + n_nodes_plot),
    pie     = lik_display,
    piecol  = couleurs,
    cex     = 0.3
  )
  
  # Tips : état observé (uniquement ceux présents dans l'arbre filtré)
  chara_plot <- chara_reconstruit[tree_plot_tips]
  tiplabels(
    pie    = to.matrix(chara_plot, etats),
    piecol = couleurs,
    cex    = 0.2
  )
  
  legend("bottomleft", legend = etats, fill = couleurs, cex = 2, fig_title)
  
  if (save) {
    dev.off()
    message("Figure sauvegardée : ", filename)
  }
}

applis_modele_modif_pourcentage <- function (best_reconstruction, 
                                             tree, df, 
                                             col_name, character,
                                             fig_title,
                                             palette_fixe,
                                             save     = FALSE, 
                                             filename = "reconstruction.png",
                                             width    = 2800,
                                             height   = 3200, 
                                             res      = 300) {
  
  tree_tips  <- tree$tip.label
  sortedData <- df[match(tree_tips, df[[col_name]]), ]
  
  chara_reconstruit <- sortedData[[character]]
  names(chara_reconstruit) <- tree_tips
  
  etats    <- sort(unique(chara_reconstruit))
  couleurs <- palette_fixe[etats]
  
  # Simplification des camemberts ancestraux
  lik_matrix <- best_reconstruction$lik.anc
  lik_simple <- lik_matrix
  
  for (i in seq_len(nrow(lik_matrix))) {
    max_prob  <- max(lik_matrix[i, ])
    max_state <- which.max(lik_matrix[i, ])
    if (max_prob > 0.60) {
      lik_simple[i, ]          <- 0
      lik_simple[i, max_state] <- 1
    }
  }
  
  if (save) png(filename, width = width, height = height, res = res)
  
  # fsize = 0 pour ne pas afficher les labels natifs de plotTree
  plotTree(tree, lwd = 2, setEnv = TRUE, fsize = 0.0001, ftype = "off")
  
  pp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xlim <- pp$x.lim
  ylim <- pp$y.lim
  
  # ── Coordonnées des tips ──────────────────────────────────────────────────
  n_tips <- length(tree_tips)
  
  # xx et yy contiennent d'abord les tips (1:n_tips) puis les noeuds
  x_tips <- pp$xx[seq_len(n_tips)]
  y_tips <- pp$yy[seq_len(n_tips)]
  
  # x d'alignement : le plus à droite de tous les tips + un décalage
  cex_label <- 0.4
  x_align   <- max(x_tips) + strwidth("i", cex = cex_label) * 2
  
  # ── Labels alignés ────────────────────────────────────────────────────────
  # Petite ligne de connexion entre le bout de branche et le label
  segments(x_tips, y_tips, x_align, y_tips,
           col = "grey70", lwd = 0.3, lty = "dotted")
  
  text(x_align, y_tips,
       labels = tree_tips,
       cex    = cex_label,
       adj    = c(0, 0.5),   # aligné à gauche sur x_align
       xpd    = TRUE)
  # ─────────────────────────────────────────────────────────────────────────
  
  add.scale.bar(
    x      = xlim[2] * 0.85,
    y      = ylim[1],
    length = 0.0001,
    lwd    = 2,
    cex    = 0.7
  )
  
  n_nodes <- tree$Nnode
  nodelabels(
    node    = (n_tips + 1):(n_tips + n_nodes),
    pie     = lik_simple,
    piecol  = couleurs,
    cex     = 0.3
  )
  
  tiplabels(
    pie    = to.matrix(chara_reconstruit, etats),
    piecol = couleurs,
    cex    = 0.2
  )
  
  legend("bottomleft", legend = etats, fill = couleurs, cex =2, fig_title)
  
  if (save) {
    dev.off()
    message("Figure sauvegardée : ", filename)
  }
}

applis_modele_modif_pourcentage_sansoutgroup <- function (best_reconstruction, 
                                                          tree,
                                                          tree_display = NULL,
                                                          df, 
                                                          col_name, character,
                                                          fig_title,
                                                          palette_fixe,
                                                          save     = FALSE, 
                                                          filename = "reconstruction.png",
                                                          width    = 2800,
                                                          height   = 3200, 
                                                          res      = 300) {
  
  tree_plot      <- if (!is.null(tree_display)) tree_display else tree
  tree_tips      <- tree$tip.label
  tree_plot_tips <- tree_plot$tip.label
  
  sortedData <- df[match(tree_tips, df[[col_name]]), ]
  chara_reconstruit <- sortedData[[character]]
  names(chara_reconstruit) <- tree_tips
  
  etats    <- sort(unique(chara_reconstruit))
  couleurs <- palette_fixe[etats]
  
  # ── Mapping des noeuds : arbre complet → arbre filtré ────────────────────
  if (!is.null(tree_display)) {
    n_nodes_plot <- tree_plot$Nnode
    n_tips_plot  <- length(tree_plot_tips)
    lik_display  <- matrix(NA, 
                           nrow = n_nodes_plot, 
                           ncol = length(etats),
                           dimnames = list(NULL, etats))
    
    for (i in seq_len(n_nodes_plot)) {
      node_id   <- n_tips_plot + i
      desc_tips <- tree_plot_tips[
        phangorn::Descendants(tree_plot, node_id, type = "tips")[[1]]
      ]
      
      if (length(desc_tips) >= 2) {
        node_complet <- getMRCA(tree, desc_tips)
      } else {
        node_complet <- which(tree$tip.label == desc_tips[1])
      }
      
      n_tips_complet <- length(tree_tips)
      idx <- node_complet - n_tips_complet
      
      if (!is.na(idx) && idx > 0 && idx <= nrow(best_reconstruction$lik.anc)) {
        lik_display[i, ] <- best_reconstruction$lik.anc[idx, etats]
      }
    }
    
  } else {
    lik_display <- best_reconstruction$lik.anc
  }
  # ─────────────────────────────────────────────────────────────────────────
  
  # ── Simplification des camemberts à 60% ──────────────────────────────────
  lik_simple <- lik_display
  
  for (i in seq_len(nrow(lik_display))) {
    max_prob  <- max(lik_display[i, ], na.rm = TRUE)
    max_state <- which.max(lik_display[i, ])
    
    if (!is.na(max_prob) && max_prob > 0.60) {
      lik_simple[i, ]          <- 0
      lik_simple[i, max_state] <- 1
    }
  }
  # ─────────────────────────────────────────────────────────────────────────
  
  if (save) png(filename, width = width, height = height, res = res)
  
  plotTree(tree_plot, lwd = 2, setEnv = TRUE, fsize = 0.0001, ftype = "off")
  
  pp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xlim <- pp$x.lim
  ylim <- pp$y.lim
  
  n_tips_plot <- length(tree_plot_tips)
  x_tips      <- pp$xx[seq_len(n_tips_plot)]
  y_tips      <- pp$yy[seq_len(n_tips_plot)]
  
  cex_label <- 0.4
  x_align   <- max(x_tips) + strwidth("i", cex = cex_label) * 2
  
  segments(x_tips, y_tips, x_align, y_tips,
           col = "grey70", lwd = 0.3, lty = "dotted")
  
  text(x_align, y_tips,
       labels = tree_plot_tips,
       cex    = cex_label,
       adj    = c(0, 0.5),
       xpd    = TRUE)
  
  add.scale.bar(
    x      = xlim[2] * 0.85,
    y      = ylim[1],
    length = 0.0001,
    lwd    = 2,
    cex    = 0.7
  )
  
  n_nodes_plot <- tree_plot$Nnode
  nodelabels(
    node    = (n_tips_plot + 1):(n_tips_plot + n_nodes_plot),
    pie     = lik_simple,                            # ← camemberts simplifiés
    piecol  = couleurs,
    cex     = 0.3
  )
  
  chara_plot <- chara_reconstruit[tree_plot_tips]
  tiplabels(
    pie    = to.matrix(chara_plot, etats),
    piecol = couleurs,
    cex    = 0.2
  )
  
  legend("bottomleft", legend = etats, fill = couleurs, cex = 2, fig_title)
  
  if (save) {
    dev.off()
    message("Figure sauvegardée : ", filename)
  }
}


# modele sans modif de camembert 
phylo_buchnera <- prepa_tree(phylo_buchnera, c("wamizawa", "CS041_S97", "CS042_S98") )
individus_a_supprimer <- c("wamizawa", "CS041_S97", "CS042_S98")
phylo_buchnera_filtre <- drop.tip(phylo_buchnera, individus_a_supprimer)
phylo_buchnera_filtre <- prepa_tree(phylo_buchnera_filtre, c("L6Lap03", "CS148_S140", "L14Lap03", "L10Lap02") )


# geographie : 
modele_reconstruction_zone_geo <- creation_modele (Echantillonnage_apisum, "read_file_name", "zone_geographique", phylo_buchnera)
applis_modele(modele_reconstruction_zone_geo, phylo_buchnera, 
              Echantillonnage_apisum, "read_file_name", 
              "zone_geographique", "reconstruction geographie", 
              palette_resum_localisation, save = TRUE, "reconstruction_localisation_complete.png",
              width  = 9000,
              height = 9999,
              res    = 750)

applis_modele_sans_outgroup(modele_reconstruction_zone_geo, phylo_buchnera, 
                            phylo_buchnera_filtre,
                            Echantillonnage_apisum, "read_file_name", 
                            "zone_geographique", "reconstruction geographie", 
                            palette_resum_localisation, save = TRUE, 
                            "reconstruction_localisation_complete_sans_outgroup.png",
                            width  = 9000,
                            height = 9999,
                            res    = 750)



#modele avec modif de camembert : 

modele_reconstruction_presence_serratia <- creation_modele (Echantillonnage_apisum, "read_file_name", "presence_serratia", phylo_buchnera)
applis_modele_modif_pourcentage (modele_reconstruction_presence_serratia, phylo_buchnera, 
                                 Echantillonnage_apisum,  "read_file_name", "presence_serratia", 
                                 "reconstruction presence serratia",  palette_presence_serratia,  
                                 save = TRUE,   "reconstruction_presence_serratia_modif.png",
                                 width  = 9000,
                                 height = 9999,
                                 res    = 750)

applis_modele_modif_pourcentage_sansoutgroup (modele_reconstruction_presence_serratia, phylo_buchnera, phylo_buchnera_filtre, 
                                              Echantillonnage_apisum,  "read_file_name", "presence_serratia", 
                                              "reconstruction presence serratia",  palette_presence_serratia,  
                                              save = TRUE,   "reconstruction_presence_serratia_modif_sans_outgroup.png",
                                              width  = 9000,
                                              height = 9999,
                                              res    = 750)

# plante hôte : 

#modele avec modif de camembert : 

modele_reconstruction_biotype <- creation_modele (Echantillonnage_apisum, "read_file_name", "resum_biotype", phylo_buchnera)
applis_modele_modif_pourcentage (modele_reconstruction_biotype, phylo_buchnera, 
                                 Echantillonnage_apisum,  "read_file_name", "resum_biotype", 
                                 "reconstruction biotype",  palette_resum_plante,  
                                 save = TRUE,   "reconstruction_biotype_modif.png",
                                 width  = 9000,
                                 height = 9999,
                                 res    = 700)

applis_modele_modif_pourcentage_sansoutgroup (modele_reconstruction_biotype, phylo_buchnera, phylo_buchnera_filtre, 
                                              Echantillonnage_apisum,  "read_file_name", "resum_biotype", 
                                              "reconstruction biotype",  palette_resum_plante,  
                                              save = TRUE,   "reconstruction_resum_biotype_modif_sans_outgroup.png",
                                              width  = 9000,
                                              height = 9999,
                                              res    = 700)



applis_modele(modele_reconstruction_biotype, phylo_buchnera, 
              Echantillonnage_apisum, "read_file_name", 
              "resum_biotype", "reconstruction biotype", 
              palette_resum_plante, save = TRUE, "reconstruction_resum_biotype_complete.png",
              width  = 9000,
              height = 9999,
              res    = 700)

applis_modele_sans_outgroup(modele_reconstruction_biotype, phylo_buchnera, 
                            phylo_buchnera_filtre,
                            Echantillonnage_apisum, "read_file_name", 
                            "resum_biotype", "reconstruction biotype", 
                            palette_resum_plante, save = TRUE, 
                            "reconstruction_resum_biotype_complete_sans_outgroup.png",
                            width  = 9000,
                            height = 9999,
                            res    = 700)



carte_points <- function(df, 
                         lat_col     = "Lat", 
                         long_col    = "Long",
                         country_col = "Country",
                         titre       = NULL, 
                         sauvegarder = FALSE,
                         filename    = "carte_points.png",
                         width       = 2800,
                         height      = 2000,
                         res         = 300) {
  
  if (!requireNamespace("maps", quietly = TRUE)) {
    install.packages("maps")
    library(maps)
  }
  
  world_map <- map_data("world")
  
  if (is.null(titre)) titre <- "Échantillonnage par pays"
  
  data_carte <- df %>%
    filter(!is.na(.data[[long_col]]), !is.na(.data[[lat_col]])) %>%
    group_by(.data[[lat_col]], .data[[long_col]], .data[[country_col]]) %>%
    summarise(
      total      = n(),
      n_serratia = sum(type_sequence == "Buchnera + Serratia"),
      .groups    = "drop"
    ) %>%
    mutate(pct_serratia = n_serratia / total * 100)
  
  p <- ggplot() +
    geom_polygon(data = world_map, 
                 aes(x = long, y = lat, group = group),
                 fill = "lightgray", color = "white", linewidth = 0.2) +
    geom_point(data = data_carte,
               aes(x     = .data[[long_col]], 
                   y     = .data[[lat_col]],
                   size  = total,
                   color = pct_serratia,
                   text  = paste0(.data[[country_col]], "\n",
                                  "Total: ", total, " individus\n",
                                  "Buchnera + Serratia: ", round(pct_serratia, 1), "%")),
               alpha = 0.85) +
    scale_size_continuous(
      range  = c(2, 8), 
      name   = "Nombre d'individus",
      breaks = c(5, 20, 40, 60, 80)
    ) +
    scale_color_gradient(
      low    = "#B2A5F3",
      high   = "#FF9A47",
      name   = "% Buchnera + Serratia",
      limits = c(0, 100)
    ) +
    labs(title = titre) +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      panel.grid      = element_blank(),
      axis.text       = element_blank(),
      axis.title      = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed(1.3)
  
  print(p)
  
  # ── Enregistrement PNG ────────────────────────────────────────────────────
  if (sauvegarder) {
    png(filename, width = width, height = height, res = res)
    print(p)
    dev.off()
    cat("✓ Carte sauvegardée :", filename, "\n")
  }
  # ─────────────────────────────────────────────────────────────────────────
  
  return(invisible(p))
}


Echantillonnage_apisum <- Echantillonnage_apisum %>%
  mutate(type_sequence = case_when(
    qualite_sequence_buchnera %in% c("bonne_qualite", "mauvaise_qualite") & qualite_sequence_serratia %in% c("bonne_qualite", "mauvaise_qualite") ~ "Buchnera + Serratia",
    qualite_sequence_buchnera%in% c("bonne_qualite", "mauvaise_qualite") & qualite_sequence_serratia == "absent"   ~ "Buchnera uniquement",
  ))
Echantillonnage_apisum$Long_regroup<- as.numeric(Echantillonnage_apisum$Long_regroup)
Echantillonnage_apisum$Lat_regroup <- as.numeric(Echantillonnage_apisum$Lat_regroup)

carte_points(Echantillonnage_apisum,
             lat_col     = "Lat_regroup",
             long_col    = "Long_regroup",
             country_col = "Country",
             titre       = "Fréquence de Serratia par site",
             sauvegarder = TRUE,
             filename    = "carte_points.png",
             width       = 6800,
             height      = 6000,
             res         = 800)
