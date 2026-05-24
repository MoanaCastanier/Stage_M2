###### Traitement des annotation ####

# FONCTION 1: LISTE DES GÈNES ET LEUR FICHIER

liste_genes_fichiers <- function(fichier_csv) {
  # Charger les données
  data <- read.csv(fichier_csv, stringsAsFactors = FALSE)
  
  # Créer une liste unique de chaque gène avec le(s) fichier(s) où on le trouve
  genes_fichiers <- data %>%
    select(gene, fichier) %>%
    distinct() %>%
    arrange(gene, fichier)
  
  cat(paste("\nNombre total de combinaisons gène-fichier:", nrow(genes_fichiers), "\n"))
  cat(paste("Nombre de gènes uniques:", n_distinct(genes_fichiers$gene), "\n"))
  cat(paste("Nombre de fichiers uniques:", n_distinct(genes_fichiers$fichier), "\n\n"))
  
  return(genes_fichiers)
}

# FONCTION 2: GÈNES PRÉSENTS DANS TOUS LES FICHIERS

genes_dans_tous_fichiers <- function(fichier_csv) {
  # Charger les données
  data <- read.csv(fichier_csv, stringsAsFactors = FALSE)
  
  # Compter le nombre total de fichiers
  nb_fichiers_total <- n_distinct(data$fichier)
  
  # Trouver les gènes présents dans tous les fichiers
  genes_tous <- data %>%
    group_by(gene) %>%
    summarise(nb_fichiers = n_distinct(fichier), .groups = "drop") %>%
    filter(nb_fichiers == nb_fichiers_total) %>%
    arrange(gene) %>%
    select(gene)
  
  cat(paste("\nNombre total de fichiers:", nb_fichiers_total, "\n"))
  cat(paste("Gènes présents dans TOUS les fichiers:", nrow(genes_tous), "\n"))
  
  if (nrow(genes_tous) > 0) {
    cat("\nListe des gènes:\n")
    print(genes_tous, n = Inf)
  } else {
    cat("Aucun gène présent dans tous les fichiers.\n")
  }
  cat("\n")
  
  return(genes_tous)
}

# FONCTION 3: GÈNES PRÉSENTS CHEZ TOUS LES INDIVIDUS 

genes_tous_individus_par_fichier <- function(fichier_csv) {
  # Charger les données
  data <- read.csv(fichier_csv, stringsAsFactors = FALSE)
  
  # Définir le nombre d'individus attendus par fichier
  # Par défaut 9, sauf pour le fichier spécifique qui en a 7
  fichier_special <- "Alignment_buchnera_09003A_Ar_Po_58_ArPo28_ArPo31_CS001_S1_a_CS003_S3.fasta"
  
  # Compter le nombre d'individus distincts par fichier
  nb_individus_par_fichier <- data %>%
    group_by(fichier) %>%
    summarise(nb_individus_reel = n_distinct(individu), .groups = "drop") %>%
    mutate(nb_individus_attendu = ifelse(fichier == fichier_special, 7, 9))
  
  cat("\nNombre d'individus par fichier:\n")
  print(nb_individus_par_fichier, n = Inf)
  cat("\n")
  
  # Pour chaque combinaison fichier-gène, compter le nombre d'individus
  genes_count <- data %>%
    group_by(fichier, gene) %>%
    summarise(nb_individus_gene = n_distinct(individu), .groups = "drop")
  
  # Joindre avec le nombre attendu et filtrer
  genes_tous_individus <- genes_count %>%
    left_join(nb_individus_par_fichier, by = "fichier") %>%
    filter(nb_individus_gene == nb_individus_attendu) %>%
    select(gene, fichier, nb_individus = nb_individus_attendu) %>%
    arrange(fichier, gene)
  
  cat(paste("Nombre total de combinaisons gène-fichier où le gène est présent chez TOUS les individus:", 
            nrow(genes_tous_individus), "\n\n"))
  
  # Résumé par fichier
  resume <- genes_tous_individus %>%
    group_by(fichier) %>%
    summarise(nb_genes = n(), .groups = "drop") %>%
    arrange(desc(nb_genes))
  
  cat("Résumé par fichier:\n")
  print(resume, n = Inf)
  cat("\n")
  
  return(genes_tous_individus)
}

# Substitution 

consecutif_genes_fichiers_substitution <- liste_genes_fichiers("genes_differences_consecutives.csv")

consecutif_tous_fichiers_substitution <- genes_dans_tous_fichiers("genes_differences_consecutives.csv")

consecutif_tous_individus_substitution <- genes_tous_individus_par_fichier("genes_differences_consecutives.csv")