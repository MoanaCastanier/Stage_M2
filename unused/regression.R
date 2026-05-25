##### Regression ####
Echantillonnage_apisum <- read.csv2("~/Documents/Moana/Echantillonnage_apisum.csv", head=TRUE, sep=",")
Echantillonnage_apisum <- Echantillonnage_apisum [Echantillonnage_apisum$read_file_name != "CS189_S180",]
Echantillonnage_apisum <- Echantillonnage_apisum[Echantillonnage_apisum$qualite_sequence_buchnera != "mauvaise_qualite", ]

phylo_buchnera_incomplet<- read.tree("~/Documents/Moana/phylo_snippy/phylo_buchnera_incomplet.nwk")
phylo_serratia <- read.tree("~/Documents/Moana/phylo_snippy/phylo_serratia_without_outgroup.nwk")

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