###############################################################################
################ Verification de la qualité des séquences #####################
###############################################################################

library(seqinr)

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

###### Traitement des données séquences #########

Verification_sequence_serratia <- analyser_qualite_fa("~/Documents/Moana/Seq_symbiote/ConsensusMars2026_/Serratia_chromosome")
Verification_sequence_serratia$pourcentage_non_acgt<- as.numeric(Verification_sequence_serratia$pourcentage_non_acgt)
Verification_sequence_serratia[,'qualite_sequence_serratia'] <-NA
# seuil de non acgt = 10%
Verification_sequence_serratia$qualite_sequence_serratia <- case_when(
  Verification_sequence_serratia$pourcentage_non_acgt > 50 ~ "absent",
  Verification_sequence_serratia$pourcentage_non_acgt < 10 ~ "bonne_qualité",
  TRUE ~ "mauvaise_qualité"  # pour tous les autres cas (entre 10 et 50)
)

Verification_sequence_buchnera <- analyser_qualite_fa("~/Documents/Moana/Seq_symbiote/ConsensusMars2026_/Buchnera_chromosome")
Verification_sequence_buchnera$pourcentage_non_acgt<- as.numeric(Verification_sequence_buchnera$pourcentage_non_acgt)
Verification_sequence_buchnera[,'qualite_sequence_buchnera'] <-NA
# seuil de non acgt = 10%, seuil ou conscidéré comme absent = 50%
Verification_sequence_buchnera$qualite_sequence_buchnera <- case_when(
  Verification_sequence_buchnera$pourcentage_non_acgt > 50 ~ "absent",
  Verification_sequence_buchnera$pourcentage_non_acgt < 10 ~ "bonne_qualité",
  TRUE ~ "mauvaise_qualité"  # pour tous les autres cas (entre 10 et 50)
)

table(Verification_sequence_buchnera$qualite_sequence_buchnera)
table(Verification_sequence_serratia$qualite_sequence_serratia)

write.csv(Verification_sequence_buchnera, "Verification_sequence_buchnera.csv", row.names = FALSE)
write.csv(Verification_sequence_serratia, "Verification_sequence_serratia.csv", row.names = FALSE)


## Comparaison méthode assemblage : 

# Verification des 4 méthodes : 
Verification_sequence_buchnera <- analyser_qualite_fasta("~/Documents/Moana/Seq_symbiote/diff_assemblage/Buchnera_bowtie")
Verification_sequence_buchnera_Denovo <- analyser_qualite_fasta("~/Documents/Moana/Seq_symbiote/diff_assemblage/Buchnera_Denovo")
Verification_sequence_buchnera_Marjo <- analyser_qualite_fasta("~/Documents/Moana/Seq_symbiote/diff_assemblage/Buchnera_Marjo")
Verification_sequence_buchnera_reassembly <- analyser_qualite_fasta("~/Documents/Moana/Seq_symbiote/diff_assemblage/Buchnera_reassembly")
write.csv(Verification_sequence_buchnera, "Verification_sequence_buchnera_bowtie_partielle.csv", row.names = FALSE)
write.csv(Verification_sequence_buchnera_Denovo , "Verification_sequence_buchnera_Denovo.csv", row.names = FALSE)
write.csv(Verification_sequence_buchnera_Marjo , "Verification_sequence_buchnera_Marjo.csv", row.names = FALSE)
write.csv(Verification_sequence_buchnera_reassembly, "Verification_sequence_buchnera_reassembly.csv", row.names = FALSE)


# Ajout des qualité des séquences a Dataset global : 
Verification_sequence_serratia <- read.csv2 ("Verification_sequence_serratia.csv", sep=",")
Verification_sequence_buchnera <- read.csv2 ("Verification_sequence_buchnera.csv", sep=",")
Echantillonnage_apisum<- read.csv2("Echantillonnage_apisum.csv", sep=",")

Echantillonnage_apisum <- ajout_colonnes(
  df_cible = Echantillonnage_apisum,
  df_source = Verification_sequence_serratia,
  col_nom_cible = "read_file_name",
  col_nom_source = "nom_sequence",
  col_a_copier = "pourcentage_non_acgt",
  nouveau_nom_col = "pourcentage_non_acgt_serratia"
)

Echantillonnage_apisum <- ajout_colonnes(
  df_cible = Echantillonnage_apisum,
  df_source = Verification_sequence_serratia,
  col_nom_cible = "read_file_name",
  col_nom_source = "nom_sequence",
  col_a_copier = "qualite_sequence_serratia",
  nouveau_nom_col = "qualite_sequence_serratia"
)

Echantillonnage_apisum <- ajout_colonnes(
  df_cible = Echantillonnage_apisum,
  df_source = Verification_sequence_buchnera,
  col_nom_cible = "read_file_name",
  col_nom_source = "nom_sequence",
  col_a_copier = "pourcentage_non_acgt",
  nouveau_nom_col = "pourcentage_non_acgt_buchnera"
)

Echantillonnage_apisum <- ajout_colonnes(
  df_cible = Echantillonnage_apisum,
  df_source = Verification_sequence_buchnera,
  col_nom_cible = "read_file_name",
  col_nom_source = "nom_sequence",
  col_a_copier = "qualite_sequence_buchnera",
  nouveau_nom_col = "qualite_sequence_buchnera"
)

table(Echantillonnage_apisum$qualite_sequence_buchnera)
table(Echantillonnage_apisum$qualite_sequence_serratia)

# Sauvegarder
write.csv(Echantillonnage_apisum , "Echantillonnage_apisum.csv", row.names = FALSE)


####### Isoler les séquences de bonne qualité #####

trier_fichiers_par_qualite <- function(df, 
                                       dossier_source, 
                                       dossier_cible, 
                                       col_id, 
                                       col_qualite, 
                                       valeur_qualite = "bonne_qualité",
                                       verbose = TRUE) {
  
  # 1. Vérifications de sécurité
  if (!dir.exists(dossier_source)) {
    stop(paste0("Erreur : Le dossier source '", dossier_source, "' n'existe pas."))
  }
  
  if (!all(c(col_id, col_qualite) %in% names(df))) {
    cols_manquantes <- c(col_id, col_qualite)[!c(col_id, col_qualite) %in% names(df)]
    stop(paste("Erreur : Les colonnes suivantes sont introuvables dans le dataframe :", 
               paste(cols_manquantes, collapse = ", ")))
  }
  
  # 2. Création du dossier de destination
  if (!dir.exists(dossier_cible)) {
    dir.create(dossier_cible, recursive = TRUE)
    if (verbose) message(paste("Dossier cible créé :", dossier_cible))
  }
  
  # 3. Filtrage des individus
  # Remplacement de str_trim() + tolower() + filter() + pull() du tidyverse
  qualite_col <- tolower(trimws(df[[col_qualite]]))
  valeur_ref  <- tolower(trimws(valeur_qualite))
  
  individus_cibles <- unique(df[[col_id]][qualite_col == valeur_ref])
  # Suppression des éventuels NA
  individus_cibles <- individus_cibles[!is.na(individus_cibles)]
  
  if (length(individus_cibles) == 0) {
    warning(paste0("Aucun individu trouvé avec la qualité '", valeur_qualite,
                   "' dans la colonne '", col_qualite, "'."))
    return(invisible(NULL))
  }
  
  if (verbose) {
    message(paste(length(individus_cibles), "individus sélectionnés pour la qualité '",
                  valeur_qualite, "'."))
  }
  
  # 4. Liste des fichiers disponibles dans le dossier source
  fichiers_disponibles <- list.files(dossier_source, pattern = "\\.(fa|fasta|fna)$",
                                     full.names = TRUE, ignore.case = TRUE)
  
  if (length(fichiers_disponibles) == 0) {
    warning(paste0("Aucun fichier .fa/.fasta/.fna trouvé dans le dossier '", dossier_source, "'."))
    return(invisible(NULL))
  }
  
  compteur_copies    <- 0
  compteur_manquants <- 0
  
  # 5. Boucle de copie
  for (id in individus_cibles) {
    id_str <- as.character(id)
    
    # Remplacement de str_detect(..., fixed()) par grepl(..., fixed = TRUE)
    fichiers_a_copier <- fichiers_disponibles[
      grepl(id_str, basename(fichiers_disponibles), fixed = TRUE)
    ]
    
    if (length(fichiers_a_copier) > 0) {
      for (fichier in fichiers_a_copier) {
        fichier_dest <- file.path(dossier_cible, basename(fichier))
        
        if (!file.exists(fichier_dest)) {
          file.copy(fichier, fichier_dest)
          if (verbose) message(paste("  [OK] Copié :", basename(fichier)))
          compteur_copies <- compteur_copies + 1
        } else {
          if (verbose) message(paste("  [SKIP] Existe déjà :", basename(fichier)))
        }
      }
    } else {
      if (verbose) message(paste("  [WARN] Aucun fichier trouvé pour l'individu :", id_str))
      compteur_manquants <- compteur_manquants + 1
    }
  }
  
  # 6. Résumé final
  message("\n=== Terminé ===")
  message(paste("Fichiers copiés     :", compteur_copies))
  if (compteur_manquants > 0) {
    message(paste("Individus sans fichier :", compteur_manquants))
  }
  message(paste("Résultat dans :", dossier_cible))
  
  return(invisible(list(copies = compteur_copies, manquants = compteur_manquants)))
}

trier_fichiers_par_qualite(
  df = Echantillonnage_apisum,
  dossier_source = "~/Documents/Moana/Seq_symbiote/ConsensusMars2026_/Buchnera_chromosome",
  dossier_cible = "~/Documents/Moana/Seq_symbiote/ConsensusMars2026_/Buchnera_chromosome_good",
  col_id = "read_file_name",       # Remplacez par le vrai nom de la colonne ID
  col_qualite = "qualite_sequence_buchnera",
  valeur_qualite = "bonne_qualité"
)

trier_fichiers_par_qualite(
  df = Echantillonnage_apisum,
  dossier_source = "~/Documents/Moana/Seq_symbiote/ConsensusMars2026_/Serratia_chromosome",
  dossier_cible = "~/Documents/Moana/Seq_symbiote/ConsensusMars2026_/Serratia_chromosome_good",
  col_id = "read_file_name",       # Remplacez par le vrai nom de la colonne ID
  col_qualite = "qualite_sequence_serratia",
  valeur_qualite = "bonne_qualité"
)