# Analyses statistiques - portage de bactéries résistantes chez les petits
# mammifères sauvages urbains
#
# Etude transversale conduite dans une métropole française (19 sites, 99 individus)
# Espèces cibles : Rattus norvegicus et Apodemus sylvaticus
# Bactéries ciblées : Entérobactérales BLSE, S. aureus résistant macrolides/pénicilline
#                    et Acinetobacter baumannii
#
# Le code suit la structure suivante :
#   1. Chargement des packages
#   2. Import et nettoyage des données
#   3. Statistiques descriptives
#   4. Tests bivariés (Fisher, Kruskal-Wallis, Wilcoxon)
#   5. Modélisation par GLM Firth (niveau site)
#   6. GLMM (niveau individu, annexe)
#   7. Diagnostics de colinéarité et EPV


# ===========================================================================
# 1. Packages
# ===========================================================================

library(tidyverse)    # manipulation et visualisation des données
library(janitor)      # nettoyage automatique des noms de colonnes
library(lme4)         # modèles mixtes (GLMM)
library(MuMIn)        # sélection de modèles par AICc (dredge)
library(performance)  # ICC et diagnostics de modèles
library(knitr)        # tableaux formatés
library(brglm2)       # GLM avec correction de Firth (biais réduit)


# ===========================================================================
# 2. Import et mise en forme des données
# ===========================================================================

# Le fichier source est un CSV avec séparateur ";" exporté depuis Excel.
# clean_names() transforme tous les noms en snake_case minuscules.
Data_indiv <- read.csv("data_individual_records.csv", sep = ";") |>
  clean_names()

# On renomme certaines variables pour plus de lisibilté dans la suite
Data_indiv <- Data_indiv |>
  rename(
    HFP        = hfp_100,        # indice d'empreinte humaine (Human Footprint)
    site       = loc_number,     # identifiant numérique du site
    Loctype    = loc_type,       # type de site (City, Sewage, Park, Zoo)
    BLSE       = entero_blse,    # portage entérobactérales BLSE (0/1)
    Staph_macro = staph_macro    # portage S. aureus résistant macrolides (0/1)
  ) |>
  mutate(
    site    = as.character(site),
    Loctype = as.factor(Loctype),
    Sex     = as.factor(sex),
    Season  = factor(season, levels = c("Autumn", "Winter", "Spring")),
    Species = as.factor(species),

    # On regroupe les types de sites en deux grandes catégories pour les
    # tests de Fisher : milieu urbain dense vs espaces verts
    Habitat = case_when(
      Loctype %in% c("City", "Sewage") ~ "Urbain_dense",
      Loctype %in% c("Park", "Zoo")    ~ "Parc",
      TRUE                              ~ NA_character_
    ) |> as.factor(),

    # Les variables continues sont centrées-réduites pour la modélisation,
    # ce qui facilite la comparaison des effets et la convergence des GLMM
    Weight_scaled = as.numeric(scale(weight)),
    HFP_scaled    = as.numeric(scale(HFP))
  )

# Sous-ensemble principal : Rattus norvegicus uniquement
# Les analyses bivariées et multivariées principales portent sur cette espèce
# car les effectifs des autres espèces sont trop faibles pour la modélisation.
Datarats <- Data_indiv |>
  filter(species == "Rattus norvegicus")

cat("Effectif total :", nrow(Data_indiv), "\n")
cat("Rattus norvegicus :", nrow(Datarats), "\n")
cat("BLSE positifs (Rn) :", sum(Datarats$BLSE, na.rm = TRUE), "\n")
cat("Staph macro positifs (Rn) :", sum(Datarats$Staph_macro, na.rm = TRUE), "\n")


# ===========================================================================
# 3. Statistiques descriptives
# ===========================================================================

# Effectifs par espèce
Data_indiv |>
  count(Species, name = "n") |>
  kable(caption = "Effectifs par espèce", booktabs = TRUE)

# Effectifs de R. norvegicus par site et saison
Datarats |>
  count(site, Loctype, Season, HFP) |>
  arrange(Loctype, site) |>
  kable(
    caption = "Effectifs de Rattus norvegicus par site et saison",
    col.names = c("Site", "Type", "Saison", "HFP-100", "n"),
    booktabs = TRUE
  )

# Prévalences brutes avec IC 95% (méthode de Wald)
tibble(
  Bactérie  = c("Entérobactérales BLSE", "S. aureus résistants macrolides"),
  n_positifs = c(
    sum(Datarats$BLSE, na.rm = TRUE),
    sum(Datarats$Staph_macro, na.rm = TRUE)
  ),
  n_total   = nrow(Datarats),
  Prévalence = paste0(round(n_positifs / n_total * 100, 1), "%"),
  IC_95 = paste0(
    "[",
    round((n_positifs / n_total) -
      1.96 * sqrt((n_positifs / n_total) * (1 - n_positifs / n_total) / n_total), 3) * 100,
    "–",
    round((n_positifs / n_total) +
      1.96 * sqrt((n_positifs / n_total) * (1 - n_positifs / n_total) / n_total), 3) * 100,
    "%]"
  )
) |>
  kable(caption = "Prévalences brutes chez Rattus norvegicus", booktabs = TRUE)

# Distribution des variables explicatives
cat("\nSexe :\n")
table(Datarats$Sex) |> kable(col.names = c("Sexe", "n"), booktabs = TRUE)

cat("\nSaison :\n")
table(Datarats$Season) |> kable(col.names = c("Saison", "n"), booktabs = TRUE)

cat("\nPoids (g) :\n")
summary(Datarats$weight) |> print()


# ===========================================================================
# 4. Tests bivariés
# ===========================================================================

# ---- 4a. Fisher : Habitat (Urbain dense vs Parc) et portage BLSE / Staph ----

fisher_BLSE <- fisher.test(table(Datarats$Habitat, Datarats$BLSE))
fisher_SM   <- fisher.test(table(Datarats$Habitat, Datarats$Staph_macro))

# Tableau récapitulatif proportion de sites positifs par type d'habitat
datasite <- Datarats |>
  group_by(site, Loctype, Habitat, HFP, HFP_scaled) |>
  summarise(
    n_total    = n(),
    n_pos_BLSE = sum(BLSE, na.rm = TRUE),
    n_pos_SM   = sum(Staph_macro, na.rm = TRUE),
    .groups = "drop"
  )

datasite |>
  group_by(Habitat) |>
  summarise(
    Sites_positifs_BLSE = sum(n_pos_BLSE > 0),
    Sites_total         = n(),
    Prop_BLSE = paste0(round(Sites_positifs_BLSE / Sites_total * 100, 1), "%"),
    Sites_positifs_SM   = sum(n_pos_SM > 0),
    Prop_SM   = paste0(round(Sites_positifs_SM / Sites_total * 100, 1), "%")
  ) |>
  kable(caption = "Proportion de sites positifs par type d'habitat", booktabs = TRUE)

# Résultats des tests de Fisher (OR et IC 95%)
tibble(
  Bactérie = c("Entérobactérales BLSE", "S. aureus macrolides-R"),
  OR = c(
    "Inf (séparation complète)",
    as.character(round(as.numeric(fisher_SM$estimate), 2))
  ),
  IC_95 = c(
    paste0("[", round(fisher_BLSE$conf.int[1], 2), " ; Inf]"),
    paste0("[", round(as.numeric(fisher_SM$conf.int[1]), 2),
           " ; ", round(as.numeric(fisher_SM$conf.int[2]), 2), "]")
  ),
  p_value = c(
    ifelse(fisher_BLSE$p.value < 0.001, "< 0,001",
           as.character(round(fisher_BLSE$p.value, 3))),
    ifelse(fisher_SM$p.value < 0.001, "< 0,001",
           as.character(round(fisher_SM$p.value, 3)))
  )
) |>
  kable(
    caption = "Test exact de Fisher : Ville vs Parc. OR = odds-ratio.
    La séparation complète pour les BLSE (0 positif en Parc) empêche
    l'estimation d'un OR fini.",
    booktabs = TRUE
  )

# ---- 4b. Facteurs individuels : sexe, saison, poids -------------------------

# BLSE
tbl_sexe_BLSE <- Datarats |>
  filter(!is.na(Sex)) |>
  group_by(Sex) |>
  summarise(n = n(), n_pos = sum(BLSE, na.rm = TRUE),
            prev = round(100 * n_pos / n, 1))
print(tbl_sexe_BLSE)
fisher.test(table(Datarats$Sex, Datarats$BLSE))

# Test de Fisher avec p-value simulée quand certaines cellules sont trop petites
fisher.test(table(Datarats$Season, Datarats$BLSE), simulate.p.value = TRUE)

Datarats <- Datarats |>
  mutate(Weight_q = ntile(weight, 4))

kruskal.test(weight ~ BLSE, data = Datarats)

# S. aureus résistant macrolides (même structure)
fisher.test(table(Datarats$Sex, Datarats$Staph_macro))
fisher.test(table(Datarats$Season, Datarats$Staph_macro), simulate.p.value = TRUE)
kruskal.test(weight ~ Staph_macro, data = Datarats)

# ---- 4c. Apodemus sylvaticus : portage A. baumannii -------------------------

Data_apo <- Data_indiv |>
  filter(Species == "Apodemus sylvaticus") |>
  mutate(ACBA_pos = acba == 1)

# Sexe vs portage
tab_sex_ACBA <- table(Data_apo$Sex, Data_apo$ACBA_pos)
fisher.test(tab_sex_ACBA)

# Poids vs portage : test de Wilcoxon (Mann-Whitney) car effectifs petits
wilcox.test(weight ~ ACBA_pos, data = Data_apo, exact = FALSE)

# Résumé descriptif porteurs vs non-porteurs
Data_apo |>
  group_by(ACBA_pos) |>
  summarise(
    n            = n(),
    n_males      = sum(Sex == "M", na.rm = TRUE),
    n_females    = sum(Sex == "F", na.rm = TRUE),
    mean_weight  = mean(weight, na.rm = TRUE),
    median_weight = median(weight, na.rm = TRUE),
    min_weight   = min(weight, na.rm = TRUE),
    max_weight   = max(weight, na.rm = TRUE)
  )


# ===========================================================================
# 5. GLM Firth au niveau du site (variable réponse : site positif oui/non)
# ===========================================================================

# Le GLM classique est inutilisable ici à cause de la séparation quasi-complète
# (tous les sites BLSE positifs sont en milieu dense, aucun en parc).
# La correction de Firth (pénalisation de la vraisemblance) permet d'obtenir
# des estimateurs moins biaisés et des IC plus fiables dans ce contexte.

datasite <- datasite |>
  mutate(
    positif    = as.integer(n_pos_BLSE > 0),  # portage BLSE au niveau site
    positif_SM = as.integer(n_pos_SM > 0)      # portage Staph macro au niveau site
  )

# Modèle pour les BLSE
mod_firth <- glm(
  positif ~ HFP_scaled,
  data   = datasite,
  family = binomial(link = "logit"),
  method = brglmFit
)
summary(mod_firth)

# OR + IC 95% (Wald, acceptable avec Firth sur cet effectif)
exp(cbind(
  OR     = coef(mod_firth),
  IC_low  = coef(mod_firth) - 1.96 * sqrt(diag(vcov(mod_firth))),
  IC_high = coef(mod_firth) + 1.96 * sqrt(diag(vcov(mod_firth)))
))

# Modèle pour S. aureus macrolides
mod_firth_SM <- glm(
  positif_SM ~ HFP_scaled,
  data   = datasite,
  family = binomial(link = "logit"),
  method = brglmFit
)
summary(mod_firth_SM)

exp(cbind(
  OR     = coef(mod_firth_SM),
  IC_low  = coef(mod_firth_SM) - 1.96 * sqrt(diag(vcov(mod_firth_SM))),
  IC_high = coef(mod_firth_SM) + 1.96 * sqrt(diag(vcov(mod_firth_SM)))
))

# Graphique prévalence BLSE ~ HFP (bulle par site, taille proportionnelle à l'effectif)
datasite |>
  mutate(prevalence = n_pos_BLSE / n_total) |>
  ggplot(aes(x = HFP, y = prevalence)) +
  geom_point(aes(size = n_total, color = Loctype), alpha = 0.8) +
  geom_text(
    aes(label = ifelse(n_pos_BLSE > 0, paste0(n_pos_BLSE, "/", n_total), "")),
    vjust = -1, size = 3
  ) +
  scale_y_continuous(
    labels = \(x) paste0(round(x * 100), "%"),
    limits = c(-0.05, 1.05)
  ) +
  scale_size_continuous(name = "Effectif\npar site", range = c(2, 8)) +
  scale_color_brewer(palette = "Set2", name = "Type de site") +
  labs(
    x = "Indice HFP-100",
    y = "Prévalence BLSE par site",
    caption = "Taille des points proportionnelle à l'effectif.
    Etiquettes : n positifs / n total pour les sites avec au moins un positif.
    La modélisation utilise la présence/absence au niveau site (GLM Firth)."
  )


# ===========================================================================
# 6. GLMM au niveau individuel (annexe)
# ===========================================================================

# Cette section tente de modéliser le portage individuel avec le site comme
# effet aléatoire. Elle documente les limites de l'approche dans ce jeu
# de données (colinéarité HFP/site, faible EPV).

Datarats_annexe <- Data_indiv |>
  filter(
    species == "Rattus norvegicus",
    !is.na(Weight_scaled),
    sex %in% c("M", "F"),
    !is.na(season),
    !is.na(BLSE)
  ) |>
  mutate(site = as.factor(site))

# Modèle global : tous les prédicteurs + effet aléatoire site
glmm_complet <- glmer(
  BLSE ~ HFP_scaled + Weight_scaled + sex + season + (1 | site),
  data   = Datarats_annexe,
  family = binomial(link = "logit"),
  na.action = na.fail
)
summary(glmm_complet)
icc(glmm_complet)  # si ICC proche de 1, le site absorbe toute la variance

# Sélection de modèles par AICc sur toutes les combinaisons de prédicteurs
options(na.action = "na.fail")
modeles_complet <- dredge(glmm_complet, rank = "AICc")
modeles_complet |>
  as.data.frame() |>
  select(any_of(c("HFP_scaled", "Weight_scaled", "sex", "season",
                  "df", "AICc", "delta", "weight"))) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  kable(booktabs = TRUE)
options(na.action = "na.omit")

# Meilleur modèle retenu après sélection
glmm_best <- glmer(
  BLSE ~ HFP_scaled + Weight_scaled + (1 | site),
  data   = Datarats_annexe,
  family = binomial(link = "logit"),
  na.action = na.fail
)
summary(glmm_best)
icc(glmm_best)

# Même démarche pour S. aureus macrolides
Datarats_annexe_SM <- Data_indiv |>
  filter(
    species == "Rattus norvegicus",
    !is.na(Weight_scaled),
    sex %in% c("M", "F"),
    !is.na(season),
    !is.na(Staph_macro)
  ) |>
  mutate(site = as.factor(site))

glmm_complet_SM <- glmer(
  Staph_macro ~ HFP_scaled + Weight_scaled + sex + season + (1 | site),
  data   = Datarats_annexe_SM,
  family = binomial(link = "logit"),
  na.action = na.fail
)
summary(glmm_complet_SM)
icc(glmm_complet_SM)

options(na.action = "na.fail")
dredge(glmm_complet_SM, rank = "AICc") |>
  as.data.frame() |>
  select(any_of(c("HFP_scaled", "Weight_scaled", "sex", "season",
                  "df", "AICc", "delta", "weight"))) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  kable(booktabs = TRUE)
options(na.action = "na.omit")


# ===========================================================================
# 7. Diagnostics : colinéarité HFP / site et EPV
# ===========================================================================

# ---- 7a. Colinéarité HFP ~ site ------------------------------------------
# HFP est une valeur fixe par site (une seule valeur raster par coordonnée).
# Inclure HFP comme effet fixe ET site comme effet aléatoire dans le même
# modèle revient à modéliser deux fois la même source de variance.

Datarats |>
  group_by(site) |>
  summarise(
    n        = n(),
    HFP_mean = round(mean(HFP_scaled, na.rm = TRUE), 3),
    HFP_sd   = round(sd(HFP_scaled, na.rm = TRUE), 4),
    .groups  = "drop"
  ) |>
  kable(
    caption = "SD intra-site du HFP centré-réduit : doit être 0 si le HFP
    est strictement constant par site (confirme la colinéarité structurelle).",
    col.names = c("Site", "N individus", "HFP moyen", "HFP SD intra-site"),
    booktabs = TRUE
  )

r2_hfp_site <- summary(lm(HFP_scaled ~ as.factor(site), data = Datarats))$r.squared
cat("R² de HFP expliqué par le site :", round(r2_hfp_site, 4), "\n")
# Un R² = 1 confirme que HFP est entièrement déterminé par le site

# ---- 7b. Events Per Variable (EPV) ----------------------------------------
# Règle heuristique : au moins 10 événements (porteurs) par prédicteur inclus
# dans le modèle. En dessous, le modèle risque d'être sur-ajusté.

n_events    <- sum(Datarats$BLSE, na.rm = TRUE)
n_predictors <- 3  # poids, sexe, saison
epv          <- n_events / n_predictors

cat("Événements (BLSE positifs) :", n_events, "\n")
cat("Nombre de prédicteurs :", n_predictors, "\n")
cat("EPV :", round(epv, 1), "(seuil recommandé : >= 10)\n")
# Si EPV < 10, on ne peut pas inclure plusieurs prédicteurs simultanement

# ---- 7c. Diagnostic GEE (non retenu - trop peu de clusters) ----------------
# Une approche GEE aurait été une alternative aux GLMM, mais le nombre de
# clusters (sites) est insuffisant pour que les corrections de variance soient
# fiables (règle : au moins 40 clusters recommandés pour GEE robuste).

diag_clusters <- Datarats |>
  group_by(site) |>
  summarise(
    n_indiv = n(),
    n_pos   = sum(BLSE, na.rm = TRUE),
    .groups = "drop"
  )

cat("Nombre de clusters :", nrow(diag_clusters), "\n")
cat("Taille moyenne :", round(mean(diag_clusters$n_indiv), 1), "\n")
cat("Range :", min(diag_clusters$n_indiv), "–", max(diag_clusters$n_indiv), "\n")
