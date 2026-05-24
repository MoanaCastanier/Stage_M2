###############################################################################
################# Rendre les noms des arbres propres ##########################
###############################################################################


#Clean les noms : 
input_file  <- "phylo_buchnera_incomplet.nwk"   # chemin vers votre fichier Newick
output_file <- "phylo_buchnera_incomplet.nwk"  # fichier de sortie

# --- LECTURE ---
tree_text <- readLines(input_file, warn = FALSE)
tree_text <- paste(tree_text, collapse = "")

# --- REMPLACEMENT ---
# Supprime le préfixe  : Buchnera_cons_Buchnera_cons_mod_
# Supprime le suffixe  : _Cons_Bowtie_snippy  (avec ou sans _snippy)
tree_simplified <- gsub(
  pattern     = "Buchnera_([^:,()]+?)(?:_Sp)?_BAMmarjo_snippy",
  replacement = "\\1",
  x           = tree_text,
  perl        = TRUE
)

# --- REMPLACEMENT ---
# Supprime le préfixe  : Buchnera_cons_Buchnera_cons_mod_
# Supprime le suffixe  : _Cons_Bowtie_snippy  (avec ou sans _snippy)
tree_simplified_2 <- gsub(
  pattern     = "Buchnera_([^:,()]+?)(?:_Sp)?_BAMbowtie(?:_Moana)?_snippy",
  replacement = "\\1",
  x           = tree_simplified,
  perl        = TRUE
)

# --- ÉCRITURE ---
writeLines(tree_simplified_2, output_file)

#Clean les noms : 
input_file  <- "phylo_serratia_without_outgroup.nwk"   # chemin vers votre fichier Newick
output_file <- "phylo_serratia_without_outgroup.nwk"  # fichier de sortie

# --- LECTURE ---
tree_text <- readLines(input_file, warn = FALSE)
tree_text <- paste(tree_text, collapse = "")

# --- REMPLACEMENT ---
# Supprime le préfixe  : Buchnera_cons_Buchnera_cons_mod_
# Supprime le suffixe  : _Cons_Bowtie_snippy  (avec ou sans _snippy)
tree_simplified <- gsub(
  pattern     = "Serratia_([^:,()]+?)(?:_Sp)?_BAMmarjo_snippy",
  replacement = "\\1",
  x           = tree_text,
  perl        = TRUE
)

# --- REMPLACEMENT ---
# Supprime le préfixe  : Buchnera_cons_Buchnera_cons_mod_
# Supprime le suffixe  : _Cons_Bowtie_snippy  (avec ou sans _snippy)
tree_simplified_2 <- gsub(
  pattern     = "Serratia_([^:,()]+?)(?:_Sp)?_BAMbowtie(?:_Moana)?_snippy",
  replacement = "\\1",
  x           = tree_simplified,
  perl        = TRUE
)

# --- ÉCRITURE ---
writeLines(tree_simplified_2, output_file)