#### RIPDIV DIVERSITY ####
library(here)
library(dplyr)

#### USER SETTTINGS ####

# directories
dir_data <- "..."
dir_out <- "..."

# set parameter for community data analyses
i.param <- 1

parameter <- case_when(i.param == 1 ~ "bac",
                       i.param == 2 ~ "fun",
                       i.param == 3 ~ "phago",
                       i.param == 4 ~ "photo")

save_output <- FALSE
overall_stats <- TRUE

# input filenames
filename_meta <- paste0("sample_id_", parameter, ".txt")
filename_counts <- paste0("counts_", parameter, ".txt")
filename_taxo <- paste0("taxa_", parameter, ".txt")
filename_env <- paste0("envfac_", parameter, ".csv")
filename_bacteria <- "bacteria_stats.csv"
filename_metab <- "metab_stats.csv"
filename_eea <- "eea_stats.csv"
filename_pigments <- "pigments_stats.csv"
filename_afdm <- "afdm_stats.csv"

# output table filenames
filename_alphadiversity <- paste0("alphadiv_", parameter, ".txt")
filename_alphadiversity_patchmean <- paste0("alphadiv_mean_", parameter, ".txt")
filename_gammadiversity <- paste0("gammadiv_", parameter, ".txt")
filename_dominant <- paste0("dominant_", parameter, ".txt")
filename_permanova <- paste0("beta_permanova_", parameter, ".txt")
filename_aovp_add <- paste0("aovp_add", ".txt")
filename_emmeans_inter <- paste0("emmeans_hedges_inter", ".txt")
filename_emmeans_add <- paste0("emmeans_hedges_add", ".txt")

# output plot filenames
filename_alphadiversity_plot <- paste0("Richness_", parameter, ".png")
filename_shannon_plot <- paste0("Shannon_", parameter, ".png")
filename_pigments_plot <- paste0("Pigments", ".png")
filename_pigments_legend <- paste0("Pigments_legend", ".png")
filename_betadiversity_plot <- paste0("NMDS_", parameter, ".png")
filename_betadiversity_legend <- paste0("NMDS_legend_", parameter, ".png")
filename_dbRDA_plot <- paste0("dbRDA_", parameter, ".png")
filename_dbRDA_legend <- paste0("dbRDA_legend_", parameter, ".png")
filename_effectsize_plot <- paste0("effectsizes", ".png")
filename_effectsize_facet_plot <- paste0("effectsizes_facets", ".png")
filename_corr_plot <- paste0("correlations", ".png")
filename_corr_plot_p <- paste0("correlations_signif", ".png")
filename_venn_plot <- paste0("Venn_", parameter, ".png")
filename_venn_plot2 <- paste0("Venn2_", parameter, ".png")

# plot labels
if (parameter == "bac") {
  y_label <- "Bacteria"
} else if (parameter == "fun") {
  y_label <- "Fungi"
} else if (parameter == "phago") {
  y_label <- "Phagotrophic protists"
} else if (parameter == "photo") {
  y_label <- "Autotrophic protists"
}

# colours & shapes
mycolor <- c("#C2DF23FF", "#ffffff", "#3F4788FF", "#ffffff") # define color for bg for treatment
myedge <- c("#C2DF23FF", "#1E9B8AFF", "#3F4788FF", "#481568FF") # define edge for col for treatment
myshape <- c(21, 22, 23, 24, 25) # define shapes for pch for ST

### END OF USER SETTINGS ####

#### LIBRARIES ####
library(phyloseq)
library(ggcorrplot)
library(vegan)
library(tidyverse)
library(data.table)
library(scales)
library(reshape2)
library(ggvenn)
library(ggVennDiagram)
library(ggrepel)
library(ggpubr)
library(readr)
library(devtools) # rtools must be installed as well
library(metagenomeSeq)
library(effectsize)
library(lmPerm)
library(emmeans)
library(stringi)
library(AICcPermanova)
library(corrplot)
library(sinkr)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
scores <- vegan::scores
here <- here::here

options(scipen = 999)

#### FUNCTIONS ####

# function veganotu to convert phyloseq object to vegan matrix
veganotu <- function(physeq) {
  require("vegan")
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

annotation_compass <- function(label, position = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"), fontsize = 20, padding = grid::unit(c(0.5, 0.5), "line"), ...) {
  position <- match.arg(position)
  x <- switch(position,
    N = 0.5,
    NE = 1,
    E = 1,
    SE = 1,
    S = 0.5,
    SW = 0,
    W = 0,
    NW = 0
  )
  y <- switch(position,
    N = 1,
    NE = 1,
    E = 0.5,
    SE = 0,
    S = 0,
    SW = 0,
    W = 0.5,
    NW = 1
  )
  hjust <- switch(position,
    N = 0.5,
    NE = 1,
    E = 1,
    SE = 1,
    S = 0.5,
    SW = 0,
    W = 0,
    NW = 0
  )
  vjust <- switch(position,
    N = 1,
    NE = 1,
    E = 0.5,
    SE = 0,
    S = 0,
    SW = 0,
    W = 0.5,
    NW = 1
  )
  f1 <- switch(position,
    N = 0,
    NE = -1,
    E = -1,
    SE = -1,
    S = 0,
    SW = 1,
    W = 1,
    NW = 1
  )
  f2 <- switch(position,
    N = -1,
    NE = -1,
    E = 0,
    SE = 1,
    S = 1,
    SW = 1,
    W = 0,
    NW = -1
  )
  annotation_custom(grid::textGrob(label,
    gp = gpar(col = "black", fontsize = fontsize),
    x = grid::unit(x, "npc") + f1 * padding[1],
    y = grid::unit(y, "npc") + f2 * padding[2],
    hjust = hjust, vjust = vjust, ...
  ))
}

# ggplot theme
mytheme <- theme(
  plot.margin = margin(t = 15, r = 10, b = 10, l = 10, unit = "pt"),
  plot.title = element_text(size = rel(1.5), face = "bold", hjust = 1),
  axis.text.x = element_text(size = 16, angle = 0), # , hjust = 0.8, vjust = 0.8),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16),
  axis.line = element_line(linewidth = 0.5),
  axis.ticks.x.bottom = element_line(linewidth = 0.5),
  axis.ticks.y.left = element_line(linewidth = 0.5),
  axis.ticks.length = unit(0.2, "cm"),
  panel.grid.major = element_blank(), # line(colour="grey", size = 0.3),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
  plot.background = element_rect(fill = "white", colour = "white"),
  strip.background.x = element_rect(fill = "white"),
  strip.background.y = element_rect(fill = "white"),
  strip.text = element_text(size = 18),
  # legend.key.width = unit(0.1, "cm"), # 0.3
  # legend.key.height = unit(0.1, "cm"), # 0.3
  legend.key = element_blank(),
  title = element_text(size = 10, face = "plain")
)


#### IMPORT DATA ####
setwd(dir_data)

# Metadata or Sample_id
meta <- read.delim(filename_meta, header = TRUE, row.names = 1)
meta <- meta %>%
  mutate(treatment = paste(transport, layer, sep = "-")) %>%
  relocate(stream, transport, layer, treatment, replicate)
head(meta)
metas <- sample_data(meta)

# ASV count
counts <- read.table(filename_counts, header = TRUE, row.names = 1, check.names = FALSE)
head(counts)
colSums(counts != 0)
tcounts <- as.data.frame(t(counts))
OTU <- otu_table(tcounts, taxa_are_rows = F)

# ASV taxa classification
tax <- read.delim(filename_taxo, header = TRUE, row.names = 1)
head(tax)
taxa <- tax_table(as.matrix(tax))

# merge to phyloseq object
merged_data <- phyloseq(OTU, taxa, metas) # sample names/OTU id must be similar
merged_data

# environmental variables
env <- read.csv(filename_env, header = TRUE, sep = ";")
env <- env %>%
  mutate(treatment = paste(transport, layer, sep = "-")) %>%
  relocate(stream, transport, layer, treatment, replicate)
head(env)

# Organic matter
afdm <- read.table(here(filename_afdm), sep = ";", header = TRUE, stringsAsFactors = TRUE)
afdm <- afdm %>%
  select(stream:AFDM) %>%
  mutate(treatment = paste(transport, layer, sep = "-")) %>%
  relocate(stream, transport, layer, treatment, replicate) %>%
  group_by(stream, transport, layer) %>%
  mutate(
    AFDM = ifelse(is.na(AFDM), mean(AFDM, na.rm = TRUE), AFDM),
    across(everything(), ~ ifelse(is.nan(.x), NA, .x))
  ) %>%
  ungroup()
head(afdm)

# metabolism
metab <- read.table(here(filename_metab), sep = ";", header = TRUE, stringsAsFactors = TRUE)
metab <- metab %>%
  # filter(replicate < 4) %>%
  select(stream:NEP) %>%
  complete(stream, transport, layer, replicate) %>%
  mutate(treatment = paste(transport, layer, sep = "-")) %>%
  relocate(stream, transport, layer, treatment, replicate) %>%
  mutate(CR = -CR) %>%
  group_by(stream, transport, layer) %>%
  mutate(
    across(.cols = c(CR, NEP), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
    across(everything(), ~ ifelse(is.nan(.x), NA, .x))
  ) %>%
  ungroup()
head(metab)

# enzymes
eea <- read.table(here(filename_eea), sep = ";", header = TRUE, stringsAsFactors = TRUE)
eea <- eea %>%
  select(stream:PP) %>%
  mutate(treatment = paste(transport, layer, sep = "-")) %>%
  relocate(stream, transport, layer, treatment, replicate) %>%
  select(-BX) %>%
  group_by(stream, transport, layer) %>%
  mutate(
    across(.cols = c(BG, AP, NAG, LAP, PO, PP), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
    across(everything(), ~ ifelse(is.nan(.x), NA, .x))
  ) %>%
  ungroup()
head(eea)

# pigments
pigments <- read.table(here(filename_pigments), sep = ";", header = TRUE, stringsAsFactors = TRUE)
pigments <- pigments %>%
  mutate(treatment = paste(transport, layer, sep = "-")) %>%
  relocate(stream, transport, layer, treatment, replicate) %>%
  select(-c(Phaeo.a, Phaeo.b, Chlorophyllide)) %>%
  group_by(stream, transport, layer) %>%
  mutate(
    across(.cols = c(
      Chl.a, Fucoxanthin, Lutein, Chl.b, Diadinoxanthin, Zeaxanthin,
      b.Carotin, Violaxanthin, Alloxanthin
    ), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
    across(everything(), ~ ifelse(is.nan(.x), NA, .x))
  ) %>%
  ungroup()
head(pigments)

# bacterial cell abundance and production
bacteria <- read.table(here(filename_bacteria), sep = ";", header = TRUE, stringsAsFactors = TRUE)
bacteria <- bacteria %>%
  mutate(treatment = paste(transport, layer, sep = "-")) %>%
  relocate(stream, transport, layer, treatment, replicate) %>%
  select(-c(BB, bspBP)) %>%
  group_by(stream, transport, layer) %>%
  mutate(
    across(.cols = c(BP, BA, cspBP, X16S, ITS), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)),
    across(everything(), ~ ifelse(is.nan(.x), NA, .x))
  ) %>%
  ungroup()
head(bacteria)

bacteria %>%
  group_by(transport) %>%
  summarize(m_BA = mean(BA, na.rm = T))

# join all environmental/response variables
env_all <- env %>%
  left_join(afdm, by = c("stream", "transport", "layer", "treatment", "replicate")) %>%
  left_join(metab, by = c("stream", "transport", "layer", "treatment", "replicate")) %>%
  left_join(eea, by = c("stream", "transport", "layer", "treatment", "replicate")) %>%
  left_join(pigments, by = c("stream", "transport", "layer", "treatment", "replicate")) %>%
  left_join(bacteria, by = c("stream", "transport", "layer", "treatment", "replicate")) %>%
  select(-c(X16S.copies, ITS.copies, om, ba, bp, cbp, chl.a, fuco, cr, ncp, ncp.om)) %>%
  rename(MV = mig.vel, DO = do, sorting = sort.c, OM = AFDM, NCP = NEP, "16S" = X16S) %>%
  group_by(stream, transport, layer, treatment) %>%
  mutate(
    across(
      .cols = c(
        MV, DO, sorting, OM, CR, NCP, BG, AP, NAG, LAP, PO, PP,
        Chl.a, Fucoxanthin, Lutein, Chl.b, Diadinoxanthin, Zeaxanthin,
        b.Carotin, Violaxanthin, Alloxanthin, BP, BA, cspBP, "16S", ITS
      ),
      ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)
    ),
    across(everything(), ~ ifelse(is.nan(.x), NA, .x))
  ) %>%
  ungroup()

row.names(env_all) <- env_all$id

#################### CODE BELOW: ONLY DO ONCE ##################################
#### pANOVA AND EFFECT SIZES (HEDGES' G) ####

if (overall_stats == TRUE){
  
  # convert data.frames to long format
  env_all_long <- env_all %>%
    pivot_longer(-c(stream, transport, layer, treatment, replicate, id), names_to = "variable", values_to = "value") %>%
    relocate(stream, transport, layer, treatment, replicate) %>%
    mutate(group = case_when(
      variable %in% c("CR", "NCP", "BP", "cspBP") ~ "metabolism",
      variable %in% c("BG", "AP", "NAG", "LAP", "PO", "PP") ~ "resources",
      variable %in% c(
        "Chl.a", "Fucoxanthin", "Lutein", "Chl.b", "Diadinoxanthin",
        "Zeaxanthin", "b.Carotin", "Violaxanthin", "Alloxanthin"
      ) ~ "pigments",
      variable %in% c("MV", "DO", "sorting", "OM") ~ "environment",
      variable %in% c("BA", "16S", "ITS") ~ "biomass",
      TRUE ~ NA_character_
    )) %>%
    mutate(variable = factor(variable, levels = unique(variable))) %>%
    mutate(group = factor(group, levels = unique(group))) %>%
    relocate(stream:id, group, variable, value) %>%
    droplevels() %>%
    as.data.table()
  
  aovdat <- env_all_long
  head(aovdat)
  
  # set contrasts
  options(contrasts = c("contr.sum", "contr.poly")) # required for Anova Type III SS
  
  # create lists for results
  template_list <- vector("list", length = length(levels(aovdat$variable)))
  aovdatlist <- template_list
  aovlist_inter <- template_list
  aovlist_add <- template_list
  emmlist_inter <- template_list
  emmlist_add <- template_list
  hedgeslist_inter <- template_list
  hedgeslist_add <- template_list
  bestmodellist <- template_list
  
  template_names <- levels(aovdat$variable)
  names(aovdatlist) <- template_names
  names(aovlist_inter) <- template_names
  names(aovlist_add) <- template_names
  names(emmlist_inter) <- template_names
  names(emmlist_add) <- template_names
  names(hedgeslist_inter) <- template_names
  names(hedgeslist_add) <- template_names
  names(bestmodellist) <- template_names
  
  # set individual contrasts (only contrasts of interests)
  MigU <- c(1, 0, 0, 0)
  StaU <- c(0, 1, 0, 0)
  MigS <- c(0, 0, 1, 0)
  StaS <- c(0, 0, 0, 1)
  
  custom_contrasts <- list(
    "migrating underlying - stationary underlying" = MigU - StaU,
    "migrating underlying - migrating superficial" = MigU - MigS,
    "stationary underlying - stationary superficial" = StaU - StaS,
    "migrating superficial - stationary superficial" = MigS - StaS
  )
  
  # calculate permANOVA and post-hoc test (selected pairwise comparisons) # here with aovp (lmPerm), can alternatively be done with aovperm (permuco)
  mylevs <- levels(aovdat$variable)
  
  for (m in mylevs) {
    g <- aovdat %>%
      filter(variable == m) %>%
      distinct(group, variable) %>%
      pull(group)
  
    aovdat_sub <- aovdat %>%
      filter(variable == m) %>%
      droplevels() %>%
      mutate(treatment = paste(transport, layer, sep = "-")) %>%
      mutate(treatment = as.factor(treatment))
    aovdatlist[[m]] <- aovdat_sub
  
    aovdat_sub_summary <- aovdat_sub %>%
      group_by(stream, transport, treatment, layer, group, variable) %>%
      filter(!is.na(value)) %>%
      summarize(replicates = n()) %>%
      relocate(stream, transport, layer, treatment, group, variable, replicates) %>%
      as.data.table()
  
    # permANOVA
    set.seed(111)
    model_inter <- aovp(value ~ transport * layer + Error(stream), aovdatlist[[m]], perm = "Exact", seqs = F, maxIter = 99999) # lmPerm
    set.seed(111)
    model_add <- aovp(value ~ transport + layer + Error(stream), aovdatlist[[m]], perm = "Exact", seqs = F, maxIter = 99999) # lmPerm
  
    AIC_model_inter <- AIC(model_inter$Within)
    AIC_model_add <- AIC(model_add$Within)
  
    if (AIC_model_inter < AIC_model_add | abs(AIC_model_inter - AIC_model_add) > 2) {
      print("Best model: Interaction")
      bestmodellist[[m]] <- "interaction"
    } else if (AIC_model_inter > AIC_model_add) {
      print("Best model: Addition")
      bestmodellist[[m]] <- "addition"
    }
  
    df_inter <- as.data.frame(summary(model_inter$Within)[[1]])
    df_inter$factor <- row.names(df_inter)
    df_inter <- df_inter %>%
      mutate(variable = m) %>%
      relocate(variable, factor)
    row.names(df_inter) <- 1:nrow(df_inter)
    aovlist_inter[[m]] <- df_inter
    
    df_add <- as.data.frame(summary(model_add$Within)[[1]])
    df_add$factor <- row.names(df_add)
    df_add <- df_add %>%
      mutate(variable = m) %>%
      relocate(variable, factor)
    row.names(df_add) <- 1:nrow(df_add)
    aovlist_add[[m]] <- df_add
  
    # emmeans
    emm_inter <- emmeans(model_inter, specs = pairwise ~ transport * layer) # calculate estimated marginal means
    emm_transport <- emmeans(model_add, specs = pairwise ~ transport) # calculate estimated marginal means
    emm_layer <- emmeans(model_add, specs = pairwise ~ layer) # calculate estimated marginal means
  
    # main contrasts (p-values)
    emm.contrasts_inter <- contrast(emm_inter, # custom contrasts (only selected pairwise comparisons)
      method = custom_contrasts,
      # method = "pairwise",
      adjust = "mvt"
    ) %>% as.data.frame()
  
    emm.contrasts_transport <- contrast(emm_transport, # custom contrasts (only selected pairwise comparisons)
      method = "pairwise",
      adjust = "mvt"
    ) %>% as.data.frame()
  
    emm.contrasts_layer <- contrast(emm_layer, # custom contrasts (only selected pairwise comparisons)
      method = "pairwise",
      adjust = "mvt"
    ) %>% as.data.frame()
  
    # calculate Hedge's g with package "effectsize"
    # addition model
    treatment_contrasts <- emm.contrasts_inter$contrast
    aovdat_superficial <- aovdat_sub %>%
      filter(layer == "superficial")
    aovdat_underlying <- aovdat_sub %>%
      filter(layer == "underlying")
    aovdat_migrating <- aovdat_sub %>%
      filter(transport == "migrating")
    aovdat_stationary <- aovdat_sub %>%
      filter(transport == "stationary")
  
    if (all(aovdat_superficial$value == 0)) {
      hedges_superficial <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_superficial <- hedges_g(value ~ transport, data = aovdat_superficial)
    }
  
    if (all(aovdat_underlying$value == 0)) {
      hedges_underlying <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_underlying <- hedges_g(value ~ transport, data = aovdat_underlying)
    }
  
    if (all(aovdat_migrating$value == 0)) {
      hedges_migrating <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_migrating <- hedges_g(value ~ layer, data = aovdat_migrating)
    }
  
    if (all(aovdat_stationary$value == 0)) {
      hedges_stationary <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_stationary <- hedges_g(value ~ layer, data = aovdat_stationary)
    }
  
    df_hedges_inter <- data.frame(
      "contrast" = treatment_contrasts,
      "hedges_g" = c(
        hedges_underlying$Hedges_g,
        hedges_migrating$Hedges_g,
        hedges_stationary$Hedges_g,
        hedges_superficial$Hedges_g
      ),
      "CI_lower" = c(
        hedges_underlying$CI_low,
        hedges_migrating$CI_low,
        hedges_stationary$CI_low,
        hedges_superficial$CI_low
      ),
      "CI_upper" = c(
        hedges_underlying$CI_high,
        hedges_migrating$CI_high,
        hedges_stationary$CI_high,
        hedges_superficial$CI_high
      )
    )
    df_hedges_inter <- df_hedges_inter %>%
      mutate("interpretation" = interpret_hedges_g(df_hedges_inter$hedges_g, rules = "gignac2016")) %>%
      left_join(emm.contrasts_inter, by = "contrast") %>%
      mutate(signif = case_when(p.value < 0.05 ~ TRUE, TRUE ~ FALSE))
  
    # transport model
    treatment_contrasts_transport <- emm.contrasts_transport$contrast
  
    # layer model
    treatment_contrasts_layer <- emm.contrasts_layer$contrast
  
    if (all(aovdat_sub$value == 0)) {
      hedges_transport <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
      hedges_layer <- data.frame(Hedges_g = NA, CI_low = NA, CI_high = NA)
    } else {
      hedges_transport <- hedges_g(value ~ transport, data = aovdat_sub)
      hedges_layer <- hedges_g(value ~ layer, data = aovdat_sub)
    }
  
    df_hedges_transport <- data.frame(
      "contrast" = treatment_contrasts_transport,
      "hedges_g" = hedges_transport$Hedges_g,
      "CI_lower" = hedges_transport$CI_low,
      "CI_upper" = hedges_transport$CI_high
    )
  
    df_hedges_layer <- data.frame(
      "contrast" = treatment_contrasts_layer,
      "hedges_g" = hedges_layer$Hedges_g,
      "CI_lower" = hedges_layer$CI_low,
      "CI_upper" = hedges_layer$CI_high
    )
  
    df_hedges_transport <- df_hedges_transport %>%
      mutate("interpretation" = interpret_hedges_g(df_hedges_transport$hedges_g, rules = "gignac2016")) %>%
      left_join(emm.contrasts_transport, by = "contrast") %>%
      mutate(signif = case_when(p.value < 0.05 ~ TRUE, TRUE ~ FALSE))
  
    df_hedges_layer <- df_hedges_layer %>%
      mutate("interpretation" = interpret_hedges_g(df_hedges_layer$hedges_g, rules = "gignac2016")) %>%
      left_join(emm.contrasts_layer, by = "contrast") %>%
      mutate(signif = case_when(p.value < 0.05 ~ TRUE, TRUE ~ FALSE))
  
    # combine both factor of addition model
    emm.contrasts_add <- emm.contrasts_transport %>%
      rbind(emm.contrasts_layer)
    df_hedges_add <- df_hedges_transport %>%
      rbind(df_hedges_layer)
  
    # write to output list
    emmlist_inter[[m]] <- emm.contrasts_inter
    hedgeslist_inter[[m]] <- df_hedges_inter
    emmlist_add[[m]] <- emm.contrasts_add
    hedgeslist_add[[m]] <- df_hedges_add
  
    # add group
    emmlist_inter[[m]]$variable <- m
    hedgeslist_inter[[m]]$variable <- m
    emmlist_add[[m]]$variable <- m
    hedgeslist_add[[m]]$variable <- m
  
    # add variable
    emmlist_inter[[m]]$group <- g
    hedgeslist_inter[[m]]$group <- g
    emmlist_add[[m]]$group <- g
    hedgeslist_add[[m]]$group <- g
  }
  
  # combine to single dt
  aovtable <- rbindlist(aovlist_add)
  bestmodeltable <- unlist(bestmodellist)
  dt_bestmodel <- data.table(parameter = names(bestmodeltable), best_model = bestmodeltable)
  
  emm_inter <- rbindlist(emmlist_inter)
  emm_add <- rbindlist(emmlist_add)
  
  hedges_inter <- rbindlist(hedgeslist_inter) %>%
    mutate(variable = factor(variable, levels = mylevs)) %>%
    mutate(group = factor(group)) %>%
    mutate(contrast = factor(contrast)) %>%
    arrange(group)
  hedges_add <- rbindlist(hedgeslist_add) %>%
    mutate(variable = factor(variable, levels = mylevs)) %>%
    mutate(group = factor(group)) %>%
    mutate(contrast = factor(contrast)) %>%
    mutate(contrast = recode_factor(contrast,
      "migrating - stationary" = "transport",
      "superficial - underlying" = "layer"
    )) %>%
    arrange(group)
  
  if (save_output == TRUE) {
    write.table(aovtable, here(dir_out, filename_aovp_add), row.names = FALSE)
    write.table(hedges_inter, here(dir_out, filename_emmeans_inter), row.names = FALSE)
    write.table(hedges_add, here(dir_out, filename_emmeans_add), row.names = FALSE)
  }
}

#### PLOT: HEDGES' G AND CI ####

if (overall_stats == TRUE) {

  # only transport
  y_reordered <- hedges_add %>%
    select(variable, group) %>%
    pull(variable) %>%
    unique() %>%
    as.character()
  # y_reordered <- y_reordered[!(y_reordered %in% c("OM", "BA"))]
  # y_reordered <- c(y_reordered[1:11], "OM", y_reordered[12:13], "BA", y_reordered[14:length(y_reordered)])
  
  hedges_add$variable <- factor(hedges_add$variable, levels = y_reordered)
  
  gr <- hedges_add %>%
    group_by(group) %>%
    summarise(n = length(unique(variable))) %>%
    mutate(nn = cumsum(n)) %>%
    mutate(group = str_to_title(group))
  gr$group <- recode_factor(gr$group, "Resources" = "Resource\nacquisition")
  hedges_add$dummy <- as.factor(paste(hedges_add$contrast, hedges_add$signif, sep = "-"))
  hedges_add$dummy <- factor(hedges_add$dummy, levels = rev(levels(hedges_add$dummy)))
  
  eff_plot <- ggplot(hedges_add) +
    geom_errorbarh(height = 1 / 4, aes(y = variable, xmin = CI_lower, xmax = CI_upper, colour = contrast, group = contrast), position = position_dodge(width = 0.8), show.legend = FALSE) +
    geom_point(aes(y = variable, x = hedges_g, fill = dummy, colour = contrast, group = contrast), size = 4, shape = 21, position = position_dodge(width = 0.8)) +
    geom_vline(xintercept = 0, linetype = "solid", size = 0.5) +
    geom_hline(yintercept = gr$nn[gr$group != "Biomass"] + 0.5, linetype = "longdash", size = 0.5) +
    geom_text(data = gr, aes(x = 3, y = nn - 0.7, label = group), size = 5, hjust = 1) +
    # scale_y_discrete(labels = c("Richness", "Shannon", expression(Chl~italic(a)), "Fuxoxanthin", "Lutein", expression(Chl~italic(b)),
    #                             "Diadinoxanthin", "Zeaxanthin", expression(beta*"-"*Carotin),
    #                             "Violaxanthin", "Alloxanthin", "OM", "BP", "cell-sp. BP", "BA",
    #                             "16S copies", "ITS copies",
    #                             "CR", "NCP", "BG", "AP", "NAG", "LAP", "PO", "PP")) +
    scale_y_discrete(labels = c(
      "Migr. Velocity", "Diss. Oxygen", "Sorting", "Org. Matter", "CR", "NCP", "Bact. Production", "cell-spec. BP",
      "BG", "AP", "NAG", "LAP", "PO", "PP", expression(Chl ~ italic(a)),
      "Fucoxanthin", "Lutein", expression(Chl ~ italic(b)), "Diadinoxanthin",
      "Zeaxanthin", expression(beta * "-" * Carotin), "Violaxanthin", "Alloxanthin",
      "Bact. Abundance", "16S copies", "ITS copies"
    )) +
    # scale_fill_manual(values = c("white", "black"), labels = c("No", "Yes"), guide = guide_legend(override.aes = list(shape = 21))) +
    scale_fill_manual(
      values = c(
        "transport-TRUE" = myedge[1], "transport-FALSE" = "white",
        "layer-TRUE" = myedge[3], "layer-FALSE" = "white"
      ),
      breaks = c("transport-TRUE", "transport-FALSE"), labels = c("Yes", "No"), guide = guide_legend(order = 2, override.aes = list(shape = 21, fill = c("black", "white")))
    ) +
    scale_colour_manual(values = c(myedge[1], myedge[3]), 
                        labels = c("Transport", "Layer"),
                        guide = guide_legend(order = 1, override.aes = list(shape = 21, fill = c(myedge[1], myedge[3])))) +
    labs(x = "Hedges' g effect size", y = NULL, fill = "Significant", colour = "Factor") +
    mytheme +
    # theme(legend.position = "bottom")
    theme(
      legend.position = c(0.82, 0.68),
      legend.box.background = element_rect(colour = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.direction = "vertical",
      legend.box = "horizontal",
      legend.spacing.x = unit(0, "pt")
    )
  eff_plot
  
  if (save_output == TRUE) {
    ggsave(here(dir_out, filename_effectsize_plot), eff_plot, width = 20, height = 18, units = "cm")
  }
}

#### PLOT: CORRELATION ####

if (overall_stats == TRUE) {

  aovdat_wide <- aovdat %>%
    select(-group) %>%
    pivot_wider(id_cols = c(stream:id), names_from = "variable", values_from = "value") %>%
    select(-c(stream:id)) %>%
    filter(!is.na(rowSums(.))) %>%
    relocate(all_of(y_reordered))
  
  cor_matrix <- cor(aovdat_wide, method = "spearman")
  set.seed(111)
  cor_pmatrix <- cor_pmat(aovdat_wide, method = "spearman")
  
  colnames(cor_matrix) <- c(
    "Migr. Velocity", "Diss. Oxygen", "Sorting", "Org. Matter", "CR", "NCP",
    "Bact. Production", "cell-sp. BP", "BG", "AP", "NAG", "LAP", "PO", "PP",
    "Chl a", "Fucoxanthin", "Lutein", "Chl b",
    "Diadinoxanthin", "Zeaxanthin", ":beta-Carotin", "Violaxanthin", "Alloxanthin",
    "Bact. Abundance", "16S copies", "ITS copies"
  )
  row.names(cor_matrix) <- c(
    "Migr. Velocity", "Diss. Oxygen", "Sorting", "Org. Matter", "CR", "NCP",
    "Bact. Production", "cell-sp. BP", "BG", "AP", "NAG", "LAP", "PO", "PP",
    "Chl a", "Fucoxanthin", "Lutein", "Chl b",
    "Diadinoxanthin", "Zeaxanthin", ":beta-Carotin", "Violaxanthin", "Alloxanthin",
    "Bact. Abundance", "16S copies", "ITS copies"
  )
  
  colnames(cor_pmatrix) <- colnames(cor_matrix)
  row.names(cor_pmatrix) <- row.names(cor_matrix)
  
  png(
    file = here(dir_out, filename_corr_plot_p),
    width = 25, height = 25, units = "cm", res = 600
  )
  corr_plot_p <- corrplot(cor_matrix,
    method = "color", tl.col = "black", insig = "blank", p.mat = cor_pmatrix,
    type = "lower", outline = "grey", order = "original", col = brewer.pal(n = 8, name = "RdBu"),
    addCoef.col = "black", diag = FALSE, number.cex = 0.8, number.digits = 2, number.font = 1
  )
  dev.off()
  
  png(
    file = here(dir_out, filename_corr_plot),
    width = 25, height = 25, units = "cm", res = 600
  )
  corr_plot <- corrplot(cor_matrix,
    method = "color", tl.col = "black", insig = "blank",
    type = "lower", outline = "grey", order = "original", col = brewer.pal(n = 8, name = "RdBu"),
    addCoef.col = "black", diag = FALSE, number.cex = 0.8, number.digits = 2, number.font = 1
  )
  dev.off()
}

#### PLOT: PIGMENT COMPOSITION ####

pigments_wide <- pigments %>%
  select(-Chl.a) %>%
  pivot_longer(!c("stream", "transport", "layer", "treatment", "replicate"), names_to = "variable", values_to = "value") %>%
  group_by(stream, transport, layer, treatment, replicate) %>%
  mutate(rel_val = value/sum(value)*100) %>%
  ungroup()
pigments_sum <- pigments_wide %>%
  group_by(transport, layer) %>%
  summarise(sum_val = sum(value))

pigments_wide_mean <- pigments_wide %>%
  group_by(transport, layer, variable) %>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            sd_vak = sd(value, na.rm = TRUE),
            mean_rel = mean(rel_val, na.rm = TRUE),
            sd_rel = sd(rel_val, na.rm = TRUE)) %>%
  arrange(transport, -mean_rel) %>%
  mutate(variable = as.factor(variable))
pigments_wide_mean$variable <- factor(pigments_wide_mean$variable, 
                                      levels = c("Fucoxanthin", "Chl.b", "Lutein",
                                                 "Diadinoxanthin", "b.Carotin", "Zeaxanthin",
                                                 "Violaxanthin", "Alloxanthin"))

transport.labs <- c("migrating" = "Migrating", "stationary" = "Stationary")

color_values <- rev(RColorBrewer::brewer.pal(8, "Paired"))

pigment_plot <- ggplot(pigments_wide_mean, aes(x = layer, y = mean_val, fill = forcats::fct_rev(variable))) +
  geom_bar(stat = "identity", colour = "black", size = 0.3, width = 0.7) +
  xlab("") +
  ylab(expression(paste("Pigments (", mu, "g g DM"^"-1" * ")"))) +
  scale_x_discrete(
    limits = c("superficial", "underlying"),
    labels = c("Superficial", "Underlying")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) +
  facet_grid(. ~ transport, labeller = labeller(transport = transport.labs)) +
  scale_fill_manual(
    name = "Pigments", values = color_values,
    breaks = levels(forcats::fct_rev(pigments_wide_mean$variable)),
    labels = rev(c("Fucoxanthin", "Chl b", "Lutein",
               "Diadinoxanthin", ~ beta * "-Carotin", "Zeaxanthin",
               "Violaxanthin", "Alloxanthin"))
    ) +
  # limits = c("h", "b", "c", "d", "e", "f", "g", "a")) +
  guides(fill = guide_legend(reverse = F)) +
  mytheme
pigment_plot

pigment_plot_legend <- get_legend(pigment_plot + theme(
  legend.position = "right",
  legend.box.margin = margin(0, 0, 0, 12), # create some space to the left of the legend
  legend.title = element_blank()
))
as_ggplot(pigment_plot_legend)

pigment_plot <- pigment_plot +
  theme(legend.position = "none")

if (save_output == TRUE) {
  ggsave(here(dir_out, filename_pigments_plot), pigment_plot, width = 20, height = 18, units = "cm")
  ggsave(here(dir_out, filename_pigments_legend), pigment_plot_legend, width = 10, height = 10, units = "cm")
}

pigments_wide_mean %>%
  filter(variable == "Fucoxanthin")
pigments_wide %>%
  

############## CODE BELOW: SPECIFIC FOR EACH PARAMETER #########################
#### RAREFACTION ####

# Rarefy to even sequencing depth and calculate indices:
data_to_rarefy <- merged_data # phyloseq data we need to rarefy
data_stats <- meta

raref_data <- rarefy_even_depth(data_to_rarefy, rngseed = TRUE)

pf1 <- plot_richness(raref_data, measures = "Observed") # ACE, Shannon, Simpson InvSimpson, Fisher, ..

if (parameter %in% c("bac", "fun")) {
  rare_step <- 1000
} else {
  rare_step <- 5
}
rarefaction_plot <- rarecurve(tcounts, step = rare_step, label = FALSE)

#### ALPHADIVERSITY ####

# Calculate Alpha-diversity Indices for each sample
alphadiv <- estimate_richness(raref_data, split = TRUE)
alphadiv <- cbind(data_stats, alphadiv) # combine design data with indices
head(alphadiv)

# export as txt
if (save_output == TRUE) {
  write.table(alphadiv, here(dir_out, filename_alphadiversity), row.names = F)
}

# summarize alphadiv data for table with diversity indices
colnames(alphadiv)
alpha_mean <- alphadiv[, c(6:7, 9, 11:14)] # select indices (columns) to analyze
head(alpha_mean)
sum_alpha <- aggregate(alpha_mean, list(meta$stream, meta$transport, meta$layer, meta$treatment), function(x) c(M = mean(x), SE = sd(x) / sqrt(length(x))))
sum_alpha <- sum_alpha %>%
  rename(stream = Group.1, transport = Group.2, layer = Group.3, treatment = Group.4)
head(sum_alpha)

if (save_output == TRUE) {
  write.table(sum_alpha, here(dir_out, filename_alphadiversity_patchmean), row.names = FALSE)
}

set.seed(111)
alpha_model_add <- aovp(Shannon ~ transport + layer + Error(stream), alphadiv, perm = "Exact", seqs = F, maxIter = 99999) # lmPerm
summary(alpha_model_add)

alpha_emm_transport <- emmeans(alpha_model_add, specs = pairwise ~ transport) # calculate estimated marginal means
alpha_emm.contrasts_transport <- contrast(alpha_emm_transport, # custom contrasts (only selected pairwise comparisons)
                                    method = "pairwise",
                                    adjust = "mvt") %>% as.data.frame()

alpha_effectsize <- hedges_g(Observed ~ transport, data = alphadiv)
alpha_effectsize

#### GAMMA DIVERSITY ####

gammadiv <- estimate_richness(raref_data, split = FALSE)
gammadiv <- cbind(data_stats, gammadiv) # combine design data with indices
head(gammadiv)

if (save_output == TRUE) {
  write.table(gammadiv, here(dir_out, filename_gammadiversity), row.names = FALSE)
}

#### PLOT OBSERVED RICHNESS (RAREFIED DATA) ####

transport_labs <- c("Migrating ripple", "Stationary")
names(transport_labs) <- c("migrating", "stationary")
alphadiv$treatment <- factor(alphadiv$treatment)

if (parameter == "bac") {
  richness_breaks <- seq(0, 3000, 500)
} else if (parameter == "fun") {
  richness_breaks <- seq(0, 800, 200)
} else if (parameter == "phago") {
  richness_breaks <- seq(0, 90, 20)
} else if (parameter == "photo") {
  richness_breaks <- seq(0, 200, 40)
}

richness_plot <- ggplot(alphadiv, aes(x = layer, y = Observed, group = treatment)) +
  geom_boxplot(outlier.colour = NA, colour = c("black"), fill = myedge, size = 0.7) +
  geom_point(aes(shape = stream),
    size = 3, fill = "gray20", colour = "gray60",
    position = position_jitter(width = 0.2, height = 0)
  ) +
  #geom_text(aes(label = row.names(alphadiv)), position = position_jitter(width = 0.2, height = 0)) +
  #labs(x = NULL, y = paste(y_label, "(# ASV)")) + 
  labs(x = NULL, y = paste0("Observed richness ", "(", y_label, ")")) +
  scale_y_continuous(breaks = richness_breaks, limits = c(min(richness_breaks), max(richness_breaks))) +
  scale_x_discrete(
    limits = c("superficial", "underlying"),
    labels = c("Superficial", "Underlying")
  ) +
  theme(legend.position = "none") +
  scale_shape_manual("Stream",
    values = c(21, 22, 23, 24, 25),
    limits = c("Verlorenwasser", "Pulsnitz", "Plane", "Oker", "Spree"),
    labels = c("Verlorenwasser", "Pulsnitz", "Plane", "Oker", "Spree")
  ) +
  facet_grid(. ~ transport, labeller = labeller(transport = transport_labs)) +
  mytheme

if (parameter != "photo") {
  richness_plot <- richness_plot +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_x_discrete(labels = NULL, breaks = NULL)
}

if (parameter != "bac") {
  richness_plot <- richness_plot +
    theme(strip.text.x = element_blank())
}
richness_plot

if (save_output == TRUE) {
  ggsave(here(dir_out, filename_alphadiversity_plot), richness_plot, width = 16, height = 13, units = "cm")
}

#### PLOT SHANNON INDEX ####

shannon_plot <- ggplot(alphadiv, aes(x = layer, y = Shannon)) +
  geom_boxplot(outlier.colour = NA, colour = c("black"), fill = myedge, size = 0.7) +
  geom_point(aes(shape = stream),
    size = 3, fill = "gray20", colour = "gray60",
    position = position_jitter(width = 0.2, height = 0)
  ) +
  labs(x = NULL, y = y_label) +
  scale_y_continuous(breaks = seq(0, 8, 1), limits = c(0, 8)) +
  scale_x_discrete(
    limits = c("superficial", "underlying"),
    labels = c("Superficial", "Underlying")
  ) +
  theme(legend.position = "none") +
  scale_shape_manual("Stream",
    values = c(21, 22, 23, 24, 25),
    limits = c("Verlorenwasser", "Pulsnitz", "Plane", "Oker", "Spree"),
    labels = c("Verlorenwasser", "Pulsnitz", "Plane", "Oker", "Spree")
  ) +
  facet_grid(. ~ transport, labeller = labeller(transport = transport_labs)) +
  mytheme

if (save_output == TRUE) {
  ggsave(here(dir_out, filename_shannon_plot), shannon_plot, width = 16, height = 13, units = "cm")
}

# dominant taxa per sample
d <- dominant(merged_data, level = NULL)
head(d)
length(d)

if (save_output == TRUE) {
  write.table(d, here(dir_out, filename_dominant), row.names = TRUE)
}

#### BETADIVERSITY ####

# # data standardization: cumulative sum scaled (CSS) by library metagenomeSeq
# # convert OTU table into package format
# metaSeqObject <- newMRexperiment(counts) # rows: taxa; columns: 60 samples
# # CSS normalization
# metaSeqObject_css <- cumNorm(metaSeqObject , p = cumNormStatFast(metaSeqObject))
# # convert CSS normalized data into data.frame-formatted OTU table (log transformed data)
# counts_css <- data.frame(MRcounts(metaSeqObject_css, norm = TRUE, log = FALSE))
# tcounts_css <- as.data.frame(t(counts_css))

# treatment factors transport, layer, stream
transp <- as.factor(meta$transport)
names(transp) <- row.names(meta)
layer <- as.factor(meta$layer)
names(layer) <- row.names(meta)
stream <- as.factor(meta$stream)
names(stream) <- row.names(meta)

comm <- tcounts
# comm <- tcounts_css
str(comm)

# remove zero columns (ASVs that do not occur in any sample)
ncol_before <- ncol(comm)
comm <- comm[, colSums(comm, na.rm = TRUE) != 0]
ncol_after <- ncol(comm)
print(paste(ncol_before - ncol_after, "ASVs removed"))

# Hellinger transformation (= asymetrical transformation to reduce problems with double zero data, best transformation for density data)
com_hel <- hellinger(comm) # decostand(com, method = "hellinger")

# create distance matrix
Fdis.hel <- vegdist(com_hel, method = "bray", na.rm = TRUE) # use bray

# betadisper (= multivariate version of Levene's test for homogeneity of variances)
Fbeta.hel <- betadisper(Fdis.hel, transp)
set.seed(111)
permutest(Fbeta.hel) # p > 0.05 means that variance (of the distance to centroid) is not significantly different across groups, so we meet the assumption of multivariate variance dispersion (aka. multivariate homoscedasticity)
# boxplot(Fbeta.hel)
# plot(Fbeta.hel, hull = FALSE, ellipse = TRUE) ## sd ellipse

Fbeta.hel <- betadisper(Fdis.hel, layer)
set.seed(111)
permutest(Fbeta.hel) # p > 0.05 means that variance (of the distance to centroid) is not significantly different across groups, so we meet the assumption of multivariate variance dispersion (aka. multivariate homoscedasticity)
# boxplot(Fbeta.hel)
# plot(Fbeta.hel, hull = FALSE, ellipse = TRUE) ## sd ellipse

# So by transforming after Hellinger, we obtain a smoother variance within groups that reduce the overall distinction of 'distance to centroids'

#### PERMANOVA ####

# CAUTION: order of factors matters!! (Type 2 SS)
# overall test
set.seed(111)
permanova_otu_all <- adonis2(com_hel ~ transport * layer, strata = stream, data = meta, permutations = 999, method = "bray", by = NULL)
# test by terms
set.seed(111)
permanova_otu_interact <- adonis2(com_hel ~ transport * layer, strata = stream, data = meta, permutations = 999, method = "bray")
set.seed(111)
permanova_otu_add <- adonis2(com_hel ~ transport + layer, strata = stream, data = meta, permutations = 999, method = "bray")

# calculate AIC corrected for small sample size (AICc)
permanova_interact_AIC <- AICc_permanova2(permanova_otu_interact)
permanova_interact_AIC
permanova_add_AIC <- AICc_permanova2(permanova_otu_add)
permanova_add_AIC

permanova_otu_interact
permanova_otu_add
if (save_output == TRUE) {
  write.table(permanova_otu_add, here(dir_out, filename_permanova), row.names = FALSE)
}

#### NMDS ####

# set dissimilarity index
disind <- "bray"

com_nmds <- metaMDS(com_hel, k = 2, distance = disind, autotransform = FALSE, trace = TRUE)
com_nmds_ndim <- com_nmds$ndim
com_nmds_stress <- com_nmds$stress

stressplot(com_nmds)

# calling enviromental parameters as vectors, it shows the correlation of vectors to the B-diversity
env_selected <- env_all %>%
  select(MV:ITS)
row.names(env_selected) <- row.names(env_all)

# only rows that match counts file
env_selected_metamatch <- env_selected %>%
  filter(row.names(.) %in% row.names(meta))
set.seed(111)
ef <- envfit(com_nmds, env_selected_metamatch, perm = 999, na.rm = TRUE) # env factors fit
ef

ef_df <- as.data.frame(scores(ef, display = "vectors"))
ef_df <- cbind(ef_df, "R2" = (ef$vectors)$r, "p" = (ef$vectors)$pvals)
names(ef_df) <- c("NMDS1", "NMDS2", "R2", "p")
ef_df_signif <- ef_df %>%
  filter(p < 0.05)
vector_lengths <- sqrt(ef_df_signif$NMDS1^2+ef_df_signif$NMDS2^2)
vector_lengths

#### PLOT NMDS ####

# NMDS scores for samples ("sites")
nmds_site_scores <- as.data.frame(scores(com_nmds, "sites")) # using the scores function from vegan to extract the site scores and convert to a data.frame
nmds_site_scores$site <- rownames(nmds_site_scores) # create a column of site names, from the rownames of data.scores
nmds_site_scores <- merge(nmds_site_scores, cbind(site = row.names(meta), meta), by = "site")
nmds_site_scores$treatment <- factor(nmds_site_scores$treatment, levels = c("migrating-superficial", "migrating-underlying", "stationary-superficial", "stationary-underlying"))
nmds_site_scores$stream <- factor(nmds_site_scores$stream, levels = c("Verlorenwasser", "Pulsnitz", "Plane", "Oker", "Spree"))

head(nmds_site_scores)

# NMDS scores for species ("species")
nmds_species_scores <- as.data.frame(scores(com_nmds, "species")) # using the scores function from vegan to extract the species scores and convert to a data.frame
nmds_species_scores$species <- rownames(nmds_species_scores)
head(nmds_species_scores)

# scores for environmental vectors ("vectors")
ordiArrowMul(ef)
multiplier <- case_when(parameter == "bac" ~ 2,
                        parameter == "fun" ~ 6,
                        parameter == "phago" ~ 1,
                        parameter == "photo" ~ 1.5,
                        TRUE ~ 1)
# if (parameter == "fun") {
#   multiplier <- 6
# } else {
#   multiplier <- 2
# }
# nmds_envir_scores <- as.data.frame(scores(ef, "vectors")) * multiplier
nmds_envir_scores <- ef_df_signif[c("NMDS1", "NMDS2")] * multiplier
head(nmds_envir_scores)
row.names(nmds_envir_scores)

if (parameter == "bac") {
  envir_labels <- c(
    "MV", "DO", "Sorting", "OM", "CR", "NCP",
    "AP", "NAG", "LAP", "PO", "PP",
    "Chl~italic(a)", "Fucoxanthin", "Lutein", "Chl~italic(b)",
    "Diadinoxanthin", "Zeaxanthin", paste("\u03B2","Carotin", sep ="-"), "Violaxanthin", "Alloxanthin",
    "BP", "BA", "cell-sp.~BP", "16*S"
  )
} else if (parameter == "fun") {
  envir_labels <- c(
    "Migration~velocity", "Sorting", "PP", "Bact.~abundance"
  )
} else if (parameter == "phago") {
  envir_labels <- c(
    "DO", "Sorting", "OM", "NCP",
    "AP", "NAG", "LAP", "PO", "PP",
    "Chl~italic(a)", "Fucoxanthin", "Lutein", "Chl~italic(b)",
    "Diadinoxanthin", "Zeaxanthin", paste("\u03B2","Carotin", sep ="-"),
    "BP", "BA"
  )
} else if (parameter == "photo") {
  envir_labels <- c(
    "DO", "Sorting", "CR", "NCP",
    "NAG", "LAP", "PO", "PP",
    "Chl~italic(a)", "Fucoxanthin",
    "Diadinoxanthin", paste("\u03B2","Carotin", sep ="-"),
    "BP", "BA"
  )
}

# initiate new plot
plotdat <- nmds_site_scores
specdat <- nmds_species_scores

plot.new()
ord_nmds <- invisible(ordiellipse(com_nmds, plotdat[["treatment"]], choices = c(1, 2), display = "sites", kind = "ehull", conf = 0.95, label = T))

df_ell_nmds <- data.frame()
for (g in levels(plotdat[, which(names(plotdat) %in% "treatment")])) {
  df_ell_nmds <- rbind(df_ell_nmds, cbind(as.data.frame(with(
    plotdat[plotdat[["treatment"]] == g, ],
    veganCovEllipse(ord_nmds[[g]]$cov, ord_nmds[[g]]$center, ord_nmds[[g]]$scale)
  )), treatment = g))
}

fill_labels <- levels(plotdat$treatment)
fill_labels <- str_replace(fill_labels, "-", " ")

beta_plot <- ggplot() +
  geom_point(data = plotdat, aes(x = NMDS1, y = NMDS2, shape = stream, fill = treatment), size = 3.5) + # add the point markers
  geom_path(data = df_ell_nmds, aes(x = NMDS1, y = NMDS2, colour = treatment), linetype = 1, show.legend = F) +
  scale_y_continuous(breaks = breaks_pretty(n = 10)) +
  scale_x_continuous(breaks = breaks_pretty(n = 10)) +
  scale_colour_manual(values = myedge) +
  # scale_fill_manual(values=c(mycolours, "grey20"))+
  scale_fill_manual(values = myedge, labels = fill_labels) +
  scale_shape_manual(values = myshape, na.translate = FALSE) + # breaks=c("0", "21", "105"), labels=c("0 days", "21 days", "105 days"))+
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  # guides(fill = guide_legend(override.aes=list(shape=22)))+
  labs(title = NULL, fill = NULL, shape = NULL) +
  guides(fill = guide_legend(override.aes = list(fill = myedge, shape = 21))) +
  # coord_equal()+
  annotation_compass(fontsize = 16, label = paste0("Stress = ", round(com_nmds_stress, 3), ", ", "k = ", com_nmds_ndim), position = "SE") +
  # annotation_compass(fontsize = 18, label = "p = 0.15", position = "SE") +
  annotation_compass(label = y_label, position = "NW", fontsize = 20) +
  geom_segment(data = nmds_envir_scores, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = "grey30", show.legend = F, alpha = 0.8) +
  geom_text_repel(data = nmds_envir_scores, max.overlaps = 50, point.padding = NA, nudge_x = 0.02, nudge_y = 0.02, aes(x = NMDS1, y = NMDS2, label = envir_labels), parse = TRUE, colour = "black", size = 4, show.legend = F) +
  mytheme +
  # theme(legend.key = element_rect(colour="black"))+
  # theme(legend.position = "none", panel.background = element_rect(colour = "black")) +
  theme(aspect.ratio = 1)
beta_plot

beta_plot_legend <- get_legend(beta_plot + theme(
  legend.position = "right",
  legend.box.margin = margin(0, 0, 0, 12), # create some space to the left of the legend
  legend.title = element_blank()
))
as_ggplot(beta_plot_legend)

beta_plot <- beta_plot +
  theme(legend.position = "none")

if (save_output == TRUE) {
  ggsave(here(dir_out, filename_betadiversity_plot), beta_plot, width = 20, height = 20, units = "cm")
  ggsave(here(dir_out, filename_betadiversity_legend), beta_plot_legend, width = 10, height = 10, units = "cm")
}

#### BIOENV/BVSTEP ####

# standardize data
env_selected_scaled <- as.data.frame(scale(env_selected))
env_selected_scaled_metamatch <- as.data.frame(scale(env_selected_metamatch))

# select for most dominant environmental variables ## CAUTION: takes too long for >10 variables, better use BVSTEP (below)
# sol <- bioenv(com_hel, env_selected_scaled,
#   index = "bray", upto = ncol(env), trace = FALSE, partial = NULL,
#   method = "spearman",
#   metric = c("euclidean", "mahalanobis", "manhattan", "gower")
# )
# summary(sol)
# 
# sol_summary <- summary(sol)
# dominant_env <- unlist(str_split(sol_summary$variables[which.max(sol_summary$correlation)], pattern = " "))
# dominant_env
# 
# dominant_env_selected_scaled <- env_selected_scaled[colnames(env_selected_scaled) %in% dominant_env]

# alternative to bioenv: bvstep routine by  Clarke and Ainsworth (1993)
all_samples <- row.names(env_selected_scaled)
env_selected_scaled_metamatch_nona <- env_selected_scaled_metamatch %>%
  drop_na()
#dropped_samples <- all_samples[!(all_samples %in% row.names(env_selected_scaled_nona))]
com_hel_nona <- com_hel %>%
  filter(row.names(com_hel) %in% row.names(env_selected_scaled_metamatch_nona))
remaining_samples <- row.names(com_hel_nona)

# code transport and layer as dummy variables
transp_nona <- transp[names(transp) %in% remaining_samples]
layer_nona <- layer[names(layer) %in% remaining_samples]

dummy_transp <- scale(as.numeric(transp_nona == "migrating"))
dummy_layer <- scale(as.numeric(layer_nona == "superficial"))

env_selected_scaled_metamatch_nona_transp_layer <- cbind("Transport" = dummy_transp, "Layer" = dummy_layer, env_selected_scaled_metamatch_nona)

# step 1
set.seed(111)
res_bvstep1 <- bvStep(com_hel_nona, env_selected_scaled_metamatch_nona_transp_layer,
  fix.dist.method = "bray", var.dist.method = "euclidean",
  scale.fix = FALSE, scale.var = FALSE,
  max.rho = 0.95, min.delta.rho = 0.001,
  random.selection = TRUE,
  prop.selected.var = 0.2,
  num.restarts = 100,
  output.best = 10,
  var.always.include = NULL
)
res_bvstep1

if (parameter == "bac") {
  step1_vars <- c(26)
} else if (parameter == "fun") {
  step1_vars <- c(26)
} else if (parameter == "phago") {
  step1_vars <- c(4)
} else if (parameter == "photo") {
  step1_vars <- c(14)
}

# step 2
set.seed(111)
res_bvstep2 <- bvStep(com_hel_nona, env_selected_scaled_metamatch_nona_transp_layer,
  fix.dist.method = "bray", var.dist.method = "euclidean",
  scale.fix = FALSE, scale.var = FALSE,
  max.rho = 0.95, min.delta.rho = 0.001,
  random.selection = TRUE,
  prop.selected.var = 0.3,
  num.restarts = 100,
  output.best = 10,
  var.always.include = step1_vars
)
res_bvstep2

if (parameter == "bac") {
  step2_vars <- c(3, 26)
} else if (parameter == "fun") {
  step2_vars <- c(3, 26)
} else if (parameter == "phago") {
  step2_vars <- c(1, 4, 25)
} else if (parameter == "photo") {
  step2_vars <- c(14, 16)
}

# step 3
set.seed(111)
res_bvstep3 <- bvStep(com_hel_nona, env_selected_scaled_metamatch_nona_transp_layer,
  fix.dist.method = "bray", var.dist.method = "euclidean",
  scale.fix = FALSE, scale.var = FALSE,
  max.rho = 0.95, min.delta.rho = 0.001,
  random.selection = TRUE,
  prop.selected.var = 0.3,
  num.restarts = 100,
  output.best = 10,
  var.always.include = step2_vars
)
res_bvstep3

highest_corr <- res_bvstep3$best.model.rho
highest_corr
dominant_env <- unlist(str_split(res_bvstep3$best.model.vars, pattern = ","))
dominant_env

#dominant_env_selected_scaled <- env_selected_scaled_metamatch[colnames(env_selected_scaled_metamatch) %in% dominant_env]
dominant_env_selected_scaled_nona <- env_selected_scaled_metamatch_nona_transp_layer[colnames(env_selected_scaled_metamatch_nona_transp_layer) %in% dominant_env]

#### MANTEL TEST ####
# (is also included in bioenv result, no need to do here actually)
Fdis.hel_nona <- vegdist(com_hel_nona, method = "bray", na.rm = TRUE) # use bray

set.seed(111)
res_mantel <- mantel(Fdis.hel_nona, vegdist(dominant_env_selected_scaled_nona, method = "euc"), method = "spearman", permutations = 9999)
mantel_corr <- res_mantel$statistic
mantel_corr
mantel_p <- res_mantel$signif
mantel_p

#### dbRDA ####
stream_condition <- meta$stream[row.names(meta) %in% remaining_samples]

if (parameter == "bac") {
  dbRDA <- dbrda(com_hel_nona ~ MV + Fucoxanthin + cspBP + Condition(stream_condition), dominant_env_selected_scaled_nona, dist = "bray")
} else if (parameter == "fun") {
  dbRDA <- dbrda(com_hel_nona ~ MV + PO + cspBP + Condition(stream_condition), dominant_env_selected_scaled_nona, dist = "bray")
} else if (parameter == "phago") {
  dbRDA <- dbrda(com_hel_nona ~ Transport + DO + Fucoxanthin + BA + Condition(stream_condition), dominant_env_selected_scaled_nona, dist = "bray")
} else if (parameter == "photo") {
  dbRDA <- dbrda(com_hel_nona ~ PP + Fucoxanthin + Condition(stream_condition), dominant_env_selected_scaled_nona, dist = "bray")
}

summary(dbRDA)
set.seed(111)
anova.cca(dbRDA) # overall test of the significant of the analysis
set.seed(111)
anova.cca(dbRDA, by = "axis", perm.max = 500) # test axes for significance

how <- how(nperm = 200, plots = Plots(strata = stream_condition))
set.seed(111)
anova.cca(dbRDA, by = "terms", permutations = how)

#### PLOT dbRDA ####

dbRDA_summary <- summary(dbRDA)
dbRDA1_explainedvar_fitted <- round(dbRDA_summary$concont$importance[2, "dbRDA1"] * 100, 1)
dbRDA2_explainedvar_fitted <- round(dbRDA_summary$concont$importance[2, "dbRDA2"] * 100, 1)
dbRDA1_explainedvar_total <- round(dbRDA_summary$cont$importance[2, "dbRDA1"] * 100, 1)
dbRDA2_explainedvar_total <- round(dbRDA_summary$cont$importance[2, "dbRDA2"] * 100, 1)

dbRDA_sites <- as.data.frame(scores(dbRDA)$sites)
dbRDA_env <- as.data.frame(scores(dbRDA)$biplot)

row.names(dbRDA_env)

# if (parameter == "bac") {
#   row.names(dbRDA_env) <- c("MV", "Fucoxanthin", "cspBP")
# } else if (parameter == "fun") {
#   row.names(dbRDA_env) <- c("MV", "PO", "cspBP")
# } else if (parameter == "phago") {
#   row.names(dbRDA_env) <- c("DO", "Fucoxanthin", "BA")
# } else if (parameter == "photo") {
#   row.names(dbRDA_env) <- c("PP", "Fucoxanthin")
# }

# NMDS scores for samples ("sites")
dbRDA_sites$id <- rownames(dbRDA_sites) # create a column of site names, from the rownames of data.scores
dbRDA_sites <- merge(dbRDA_sites, cbind(id = row.names(meta), meta), by = "id")
dbRDA_sites$treatment <- factor(dbRDA_sites$treatment, levels = c("migrating-superficial", "migrating-underlying", "stationary-superficial", "stationary-underlying"))
dbRDA_sites$stream <- factor(dbRDA_sites$stream, levels = c("Verlorenwasser", "Pulsnitz", "Plane", "Oker", "Spree"))

head(dbRDA_sites)

if (parameter == "bac") {
  x_lims <- c(-2, 2.6)
  y_lims <- c(-3.5, 3)
} else if (parameter == "fun") {
  x_lims <- c(-4, 1.5)
  y_lims <- c(-4.3, 2)
} else if (parameter == "phago") {
  x_lims <- c(-1.7, 1.7)
  y_lims <- c(-1.6, 2.1)
} else if (parameter == "photo") {
  x_lims <- c(-2, 2.3)
  y_lims <- c(-3, 2.5)
}

dbrda_plot <- ggplot(dbRDA_sites, aes(x = dbRDA1, y = dbRDA2)) +
  geom_point(data = dbRDA_sites, aes(x = dbRDA1, y = dbRDA2, shape = stream, fill = treatment), size = 3.5) + # add the point markers
  geom_text(data = dbRDA_sites, aes(x = dbRDA1, y = dbRDA2, label = id), size = 5) + # add the point markers
  #scale_y_continuous(breaks = seq(-5, 5, by = 0.5)) +
  #scale_x_continuous(breaks = seq(-5, 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 0.5), limits = y_lims) +
  scale_x_continuous(breaks = seq(-5, 5, by = 0.5), limits = x_lims) +
  scale_colour_manual(values = myedge) +
  # scale_fill_manual(values=c(mycolours, "grey20"))+
  scale_fill_manual(values = myedge, labels = fill_labels) +
  scale_shape_manual(values = myshape, na.translate = FALSE) + # breaks=c("0", "21", "105"), labels=c("0 days", "21 days", "105 days"))+
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  labs(
    title = NULL, fill = NULL, shape = NULL,
    x = paste("dbRDA1", paste0("(", dbRDA1_explainedvar_fitted, "% of fitted, ", dbRDA1_explainedvar_total, "% of total variation)")),
    y = paste("dbRDA2", paste0("(", dbRDA2_explainedvar_fitted, "% of fitted, ", dbRDA2_explainedvar_total, "% of total variation)"))
  ) +
  guides(fill = guide_legend(override.aes = list(fill = myedge, shape = 21))) +
  # annotation_compass(fontsize = 16, label = paste0("Stress = ", round(com_nmds_stress, 3), ", ", "k = ", com_nmds_ndim), position = "SW") +
  # annotation_compass(fontsize = 18, label = "p = 0.15", position = "SE") +
  annotation_compass(label = y_label, position = "NW", fontsize = 20) +
  geom_segment(
    data = dbRDA_env, aes(x = 0, xend = dbRDA1, y = 0, yend = dbRDA2),
    arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = "black", show.legend = F
  ) +
  #geom_text_repel(data = dbRDA_env, seed = 1, point.padding = NA, nudge_x = 0.05, nudge_y = 0.06, aes(x = dbRDA1, y = dbRDA2, label = row.names(dbRDA_env)), colour = "black", size = 6, show.legend = F) +
  coord_fixed() +
  mytheme +
  theme(aspect.ratio = 1) +
  theme(legend.position = "none")
dbrda_plot

if (save_output == TRUE) {
  ggsave(here(dir_out, filename_dbRDA_plot), dbrda_plot, width = 20, height = 20, units = "cm")
}

#### Venn Diagram ####

# split count data by treatment
id_migS <- meta %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "migrating", layer == "superficial") %>%
  pull(id)
id_migU <- meta %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "migrating", layer == "underlying") %>%
  pull(id)
id_staS <- meta %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "stationary", layer == "superficial") %>%
  pull(id)
id_staU <- meta %>%
  mutate(id = row.names(.)) %>%
  filter(transport == "stationary", layer == "underlying") %>%
  pull(id)

# list containing ASV in each treatment
ASV_migS <- counts[, c(id_migS)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_migU <- counts[, c(id_migU)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_staS <- counts[, c(id_staS)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

ASV_staU <- counts[, c(id_staU)] %>%
  filter(rowSums(.) !=0) %>%
  mutate(id = row.names(.)) %>%
  pull(id)

venn_list <- list("Mig-S" = ASV_migS, 
                  "Mig-U" = ASV_migU, 
                  "Sta-S" = ASV_staS, 
                  "Sta-U" = ASV_staU)

# set breaks and colours
if (parameter == "bac") {
  venn_breaks <- seq(0, 5000, 1000)
} else if (parameter == "fun") {
  venn_breaks <- seq(0, 2500, 500)
} else if (parameter == "phago") {
  venn_breaks <- seq(0, 150, 25)
} else if (parameter == "photo") {
  venn_breaks <- seq(0, 200, 50)
}
venn_nbr <- length(venn_breaks) - 1
venn_palette <- "Heat"
venn_colours <- rev(sequential_hcl(venn_palette, n = venn_nbr))

venn_plot <- ggVennDiagram(venn_list, 
                           label_alpha = 0, 
                           label_size = 6,
                           set_size = 6,
                           category.names = c("Mig-S","Mig-U","Sta-S", "Sta-U")) + 
  scale_fill_stepsn(breaks = venn_breaks,
                     limits = c(min(venn_breaks), max(venn_breaks)),
                     colours = venn_colours) +
  #scale_fill_distiller(palette = "YlOrBr", direction = 1) +
  #scale_fill_continuous(palette = sequential_hcl(n = 6, palette = "Grays")) +
  labs(fill = NULL) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  scale_y_continuous(expand = expansion(mult = .1)) +
  #guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme(legend.text = element_text(size = 14))
venn_plot

if (save_output == TRUE) {
  ggsave(here(dir_out, filename_venn_plot), venn_plot, width = 20, height = 20, units = "cm")
}

venn_plot2 <- ggvenn(
  venn_list,
  fill_color = myedge,
  stroke_color = FALSE,
  stroke_size = 0.5, set_name_size = 6,
  text_size = 6, text_color = "black"
)
venn_plot2

if (save_output == TRUE) {
  ggsave(here(dir_out, filename_venn_plot2), venn_plot2, width = 20, height = 20, units = "cm")
}
