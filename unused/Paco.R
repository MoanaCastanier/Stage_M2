###############################################################################
######################      Script PACO       #################################
################################################################################

### Créer la matrice d'association : 
Echantillonnage_apisum <- read.csv2("Echantillonnage_apisum.csv", head=TRUE, sep=",")
Echantillonnage_apisum <- Echantillonnage_apisum [Echantillonnage_apisum$read_file_name != "CS189_S180",]
Echantillonnage_apisum <- Echantillonnage_apisum %>%
  mutate(presence_serratia = case_when(
    qualite_sequence_serratia %in% c("bonne_qualite", "mauvaise_qualite")  ~ 1,
    qualite_sequence_serratia == "absent" ~ 0,
  ))

# Filtrer les individus de mauvaise qualité
Echantillonnage_apisum_filtre <- Echantillonnage_apisum[Echantillonnage_apisum$qualite_sequence_serratia != "mauvaise_qualite" & 
                                                          Echantillonnage_apisum$qualite_sequence_buchnera != "mauvaise_qualite", ]
individus_serratia <- Echantillonnage_apisum_filtre$read_file_name[Echantillonnage_apisum_filtre$presence_serratia == 1]
tous_individus <- Echantillonnage_apisum_filtre$read_file_name

matrice_association <- outer(tous_individus, individus_serratia , FUN = function(x, y) as.integer(x == y))
rownames(matrice_association) <- tous_individus
colnames(matrice_association) <-individus_serratia 

rownames(matrice_association) <- paste0(tous_individus, "_B")
colnames(matrice_association) <- paste0(individus_serratia, "_S")

phylo_buchnera_incomplet<- read.tree("~/Documents/Moana/phylo_snippy/phylo_buchnera_incomplet.nwk")
phylo_serratia <- read.tree("~/Documents/Moana/phylo_snippy/phylo_serratia_without_outgroup.nwk")
matrix_buchnera <- cophenetic(phylo_buchnera_incomplet)
matrix_serratia <- cophenetic (phylo_serratia)
rownames(matrix_buchnera) <- paste0(rownames(matrix_buchnera), "_B")
colnames(matrix_buchnera) <- paste0(colnames(matrix_buchnera), "_B")

rownames(matrix_serratia) <- paste0(rownames(matrix_serratia), "_S")
colnames(matrix_serratia) <- paste0(colnames(matrix_serratia), "_S")


D <- prepare_paco_data(matrix_buchnera, matrix_serratia, matrice_association)
D <- add_pcoord(D, correction = "cailliez")
final_paco <- PACo(D, nperm = 10000, method = "r0")
print (final_paco$gof)


#Function to generate a random binary 0/1 matrix
generate_constrained_matrix <- function(rows, cols) {
  
  # On tire sans remise `cols` lignes parmi `rows`
  # Chaque colonne j reçoit un 1 dans la ligne selected_rows[j]
  selected_rows <- sample(rows, cols, replace = FALSE)
  
  m <- matrix(0L, nrow = rows, ncol = cols)
  for (j in seq_len(cols)) {
    m[selected_rows[j], j] <- 1L
  }
  m
}

p_value <- 0
nperm <- 1000
ss_obs <- final_paco$gof$ss
ss_perm_values <- c()

for (p in 1:nperm) {
  
  #Generate the random matrix
  random_matrix <- generate_constrained_matrix(length(tous_individus), length(individus_serratia))
  rownames(random_matrix) <- tous_individus
  colnames(random_matrix) <- individus_serratia
  
  rownames(random_matrix) <- paste0(tous_individus, "_B")
  colnames(random_matrix) <- paste0(individus_serratia, "_S")
  
  P <- prepare_paco_data(matrix_buchnera, matrix_serratia, random_matrix)
  P <- add_pcoord(P, correction = "cailliez")
  permuted_paco <- PACo(P, nperm = 1, method = "r0")
  ss_perm <- permuted_paco$gof$ss
  # Stocker le ss_perm de cette permutation
  ss_perm_values <- c(ss_perm_values, ss_perm)
  
  if (ss_perm <= ss_obs) {
    p_value <- p_value + 1
  }
}

