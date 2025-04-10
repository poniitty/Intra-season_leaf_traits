library(tidyverse)
library(mgcv)
library(MetBrewer)

spe_cols <- as.character(met.brewer("Redon", 11, type = "discrete", return_hex = TRUE))
gro_cols <- as.character(met.brewer("Hokusai3", 3, type = "discrete", return_hex = TRUE))

p_value_to_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.10 ~ ".",
    TRUE ~ ""
  )
}

Sys.setlocale("LC_TIME", "English")

d <- read_csv("data/Seasonal_traits_microclim.csv") %>% 
  mutate(Elevation = case_match(
    site,
    "RA009" ~ "Low",
    "MI259" ~ "Mid",
    "RA068" ~ "High"
  ))

d <- d %>% 
  select(site, Elevation, date, species, leaf_area, SLA, LDMC, d_weight, BITM, ExG) %>% 
  pivot_longer(cols = leaf_area:ExG, names_to = "trait") %>% 
  mutate(doy = yday(date)) %>% 
  relocate(doy, .after = date) %>% 
  mutate(Elevation = factor(Elevation, levels = c("Low","Mid","High"))) %>% 
  relocate(Elevation, .after = site) %>% 
  mutate(trait = case_match(
    trait,
    "leaf_area" ~ "Leaf area",
    "d_weight" ~ "Leaf dry mass",
    "BITM" ~ "Brigthness index (BITM)",
    "ExG" ~ "Greenness index (ExG)",
    "SLA" ~ "Specific leaf area (SLA)",
    "LDMC" ~ "Leaf dry matter content (LDMC)",
    .default = as.character(trait)
  )) %>% 
  mutate(trait = factor(trait, levels = c("Specific leaf area (SLA)",
                                          "Leaf area",
                                          "Greenness index (ExG)",
                                          "Leaf dry matter content (LDMC)",
                                          "Leaf dry mass",
                                          "Brigthness index (BITM)"))) %>% 
  mutate(species = factor(species, levels = c("Bistorta vivipara",
                                              "Solidago virgaurea",
                                              "Trollius europaeus",
                                              "Viola biflora",
                                              "Betula nana",
                                              "Salix lapponum",
                                              "Salix reticulata",
                                              "Vaccinium myrtillus",
                                              "Vaccinium uliginosum",
                                              "Dryas octopetala",
                                              "Vaccinium vitis-idaea")))

specs <- unique(d$species) %>% sort()
traits <- unique(d$trait) %>% sort()

############################################
# Single model

r2s <- tibble()
all_preds <- tibble()
for(tr in traits){
  
  dt <- d %>% 
    filter(trait == tr) %>% 
    mutate(species = factor(species),
           site = factor(site))
  
  m <- gam(value ~ s(doy, by = species, k = 4) + s(species, bs = 're') + site, data = dt, method = "REML")
  summary(m)
  
  ms <- summary(m)
  
  p_df <- ms$s.table %>% 
    as.data.frame() %>% 
    rownames_to_column("species") %>% 
    mutate(species = gsub("s\\(doy):species","",species)) %>% 
    mutate(trait = tr) %>% 
    relocate(trait) %>% 
    filter(species != "s(species)")
  
  prds <- dt %>% 
    mutate(prds = predict(m, dt, type = "response", 
                          se.fit = FALSE)) %>% 
    group_by(species, trait, doy) %>% 
    summarise(across(c(value,prds), mean), .groups = "drop")
  
  prds_df <- expand_grid(trait = tr,
                         dt %>% select(site, species) %>% distinct(),
                         doy = seq(min(dt$doy), max(dt$doy), length.out = 100))
  
  prds_df <- bind_cols(prds_df,
                       predict(m,
                               prds_df, 
                               type = "response", se.fit = TRUE) %>% 
                         as.data.frame()) %>% 
    group_by(species, trait, doy) %>% 
    summarise(across(c(fit,se.fit), mean), .groups = "drop")
  
  effect <- prds_df %>%
    group_by(species, trait) %>%
    slice_min(fit, n = 1) %>%
    bind_rows(prds_df %>%
                group_by(species, trait) %>%
                slice_max(fit, n = 1)) %>%
    ungroup() %>% 
    arrange(species, doy) %>% 
    mutate(fitmin = fit - se.fit,
           fitmax = fit + se.fit) %>% 
    group_by(species, trait) %>% 
    summarise(maxdiff = max(fitmax) - min(fitmin),
              mindiff = max(fitmin) - min(fitmax),
              .groups = "drop") %>% 
    ungroup
  
  vls <- prds_df %>% 
    group_by(species, trait) %>% 
    summarise(medianv = median(fit),
              maxv = max(fit),
              minv = min(fit), .groups = "drop")
    
  
  r2si <- prds %>%
    group_by(species, trait) %>% 
    summarise(r2 = cor(value, prds)^2,
              .groups = "drop")
  
  
  r2s <- bind_rows(r2s,
                   full_join(p_df,
                             r2si, by = join_by(trait, species)) %>% 
                     full_join(.,
                               vls, by = join_by(trait, species)) %>% 
                     full_join(.,
                               effect, by = join_by(trait, species)) %>% 
                     rename(pvalue = `p-value`))
  
  all_preds <- bind_rows(all_preds,
                         prds_df)
}

r2s <- r2s %>% 
  mutate(significance = p_value_to_stars(pvalue))

r2s %>% 
  pivot_wider(id_cols = species, names_from = trait, values_from = pvalue)

r2s %>% 
  pivot_wider(id_cols = species, names_from = trait, values_from = significance)

r2s %>% 
  mutate(pvalue = as.character(round(pvalue, 4))) %>% 
  mutate(pvalue = ifelse(pvalue == "0", "< 0.001", pvalue)) %>% 
  mutate(pvalue = paste(pvalue, " ", significance)) %>% 
  pivot_wider(id_cols = species, names_from = trait, values_from = pvalue) %>% 
  write_csv2("output/species_pvalues.csv")


all_preds <- all_preds %>% 
  mutate(
    `Growth form` = case_when(
      species %in% c("Bistorta vivipara","Solidago virgaurea","Trollius europaeus","Viola biflora") ~ "Forb",
      species %in% c("Betula nana","Salix lapponum","Salix reticulata","Vaccinium myrtillus","Vaccinium uliginosum") ~ "Deciduous shrub",
      species %in% c("Dryas octopetala","Vaccinium vitis-idaea") ~ "Evergreen shrub"
    )
  ) %>% 
  mutate(`Growth form` = factor(`Growth form`, levels = c("Forb",
                                                          "Deciduous shrub",
                                                          "Evergreen shrub"))) %>% 
  mutate(trait = factor(trait, levels = c("Specific leaf area (SLA)",
                                          "Leaf area",
                                          "Greenness index (ExG)",
                                          "Leaf dry matter content (LDMC)",
                                          "Leaf dry mass",
                                          "Brigthness index (BITM)")))


r2s <- r2s %>% 
  mutate(trait = factor(trait, levels = c("Specific leaf area (SLA)",
                                          "Leaf area",
                                          "Greenness index (ExG)",
                                          "Leaf dry matter content (LDMC)",
                                          "Leaf dry mass",
                                          "Brigthness index (BITM)"))) %>% 
  mutate(species = factor(species, levels = c("Bistorta vivipara",
                                              "Solidago virgaurea",
                                              "Trollius europaeus",
                                              "Viola biflora",
                                              "Betula nana",
                                              "Salix lapponum",
                                              "Salix reticulata",
                                              "Vaccinium myrtillus",
                                              "Vaccinium uliginosum",
                                              "Dryas octopetala",
                                              "Vaccinium vitis-idaea")))


#############################################
# Plots

r2s %>% 
  # filter(!trait %in% c("BITM","ExG")) %>%
  mutate(cvar = (maxv-minv)/medianv,
         cvar_max = (maxdiff)/medianv,
         cvar_min = (mindiff)/medianv) %>% 
  mutate(cvar_min = ifelse(cvar_min < 0, 0, cvar_min)) %>% 
  ggplot(aes(x = species, y = cvar, col = species)) +
  geom_errorbar(aes(ymin = cvar_min, ymax = cvar_max), linewidth = 1, show.legend = FALSE) +
  geom_point(aes(col = species), size = 4) +
  scale_color_manual(values = spe_cols) +
  facet_wrap(vars(trait), scales = "free_y") +
  ylab("Seasonal variation compared to median (Range-to-median ratio)") +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(aspect.ratio=1) +
  theme(legend.text = element_text(face = "italic"))
ggsave("figures/seasonal_CV.pdf", width = 9, height = 5.5)
ggsave("figures/seasonal_CV.jpg", width = 9, height = 5.5)

# SLA
all_preds %>% 
  filter(trait == "Specific leaf area (SLA)") %>% 
  ggplot(aes(y = fit, x = doy)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill = `Growth form`), 
              linetype=0, alpha=0.8) +
  geom_point(data = d %>% filter(trait == "Specific leaf area (SLA)"),
             aes(y = value), col = "gray70", size = 1) +
  scale_fill_manual(values = gro_cols) +
  facet_wrap(vars(species), scales = "free_y") +
  ylab("Specific leaf area (SLA)") +
  xlab("Day of year") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  geom_text(data = r2s %>% filter(trait == "Specific leaf area (SLA)") %>% mutate(r2 = paste0(round(r2, 2)," ", significance)),
            aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1.8, label = r2), size = 3) +
  theme(aspect.ratio=1)
ggsave("figures/SLA.pdf", width = 9, height = 7)
ggsave("figures/SLA.jpg", width = 9, height = 7)


# LDMC
all_preds %>% 
  filter(trait == "Leaf dry matter content (LDMC)") %>% 
  ggplot(aes(y = fit, x = doy)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill = `Growth form`), 
              linetype=0, alpha=0.8) +
  geom_point(data = d %>% filter(trait == "Leaf dry matter content (LDMC)"),
             aes(y = value), col = "gray70", size = 1) +
  scale_fill_manual(values = gro_cols) +
  facet_wrap(vars(species), scales = "free_y") +
  ylab("Leaf dry matter content (LDMC)") +
  xlab("Day of year") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  geom_text(data = r2s %>% filter(trait == "Leaf dry matter content (LDMC)") %>% mutate(r2 = paste0(round(r2, 2)," ", significance)),
            aes(x = -Inf, y = -Inf, hjust = -0.4, vjust = -0.8, label = r2), size = 3) +
  theme(aspect.ratio=1)
ggsave("figures/LDMC.pdf", width = 9, height = 7)
ggsave("figures/LDMC.jpg", width = 9, height = 7)

# Leaf area
all_preds %>% 
  filter(trait == "Leaf area") %>% 
  ggplot(aes(y = fit, x = doy)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill = `Growth form`), 
              linetype=0, alpha=0.8) +
  geom_point(data = d %>% filter(trait == "Leaf area"),
             aes(y = value), col = "gray70", size = 1) +
  scale_fill_manual(values = gro_cols) +
  facet_wrap(vars(species), scales = "free_y") +
  ylab("Leaf area") +
  xlab("Day of year") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  geom_text(data = r2s %>% filter(trait == "Leaf area") %>% mutate(r2 = paste0(round(r2, 2)," ", significance)),
            aes(x = -Inf, y = Inf, hjust = -0.2, vjust = 1.8, label = r2), size = 3) +
  theme(aspect.ratio=1)
ggsave("figures/Leaf_area.pdf", width = 9, height = 7)
ggsave("figures/Leaf_area.jpg", width = 9, height = 7)

# Leaf dry mass
all_preds %>% 
  filter(trait == "Leaf dry mass") %>% 
  ggplot(aes(y = fit, x = doy)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill = `Growth form`), 
              linetype=0, alpha=0.8) +
  geom_point(data = d %>% filter(trait == "Leaf dry mass"),
             aes(y = value), col = "gray70", size = 1) +
  scale_fill_manual(values = gro_cols) +
  facet_wrap(vars(species), scales = "free_y") +
  ylab("Leaf dry mass") +
  xlab("Day of year") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  geom_text(data = r2s %>% filter(trait == "Leaf dry mass") %>% mutate(r2 = paste0(round(r2, 2)," ", significance)),
            aes(x = -Inf, y = -Inf, hjust = -0.4, vjust = -0.8, label = r2), size = 3) +
  theme(aspect.ratio=1)
ggsave("figures/Leaf_dry_mass.pdf", width = 9, height = 7)
ggsave("figures/Leaf_dry_mass.jpg", width = 9, height = 7)

# BITM
all_preds %>% 
  filter(trait == "Brigthness index (BITM)") %>% 
  ggplot(aes(y = fit, x = doy)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill = `Growth form`), 
              linetype=0, alpha=0.8) +
  geom_point(data = d %>% filter(trait == "Brigthness index (BITM)"),
             aes(y = value), col = "gray70", size = 1) +
  scale_fill_manual(values = gro_cols) +
  facet_wrap(vars(species), scales = "free_y") +
  ylab("Brigthness index (BITM)") +
  xlab("Day of year") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  geom_text(data = r2s %>% filter(trait == "Brigthness index (BITM)") %>% mutate(r2 = paste0(round(r2, 2)," ", significance)),
            aes(x = -Inf, y = Inf, hjust = -0.4, vjust = 1.8, label = r2), size = 3) +
  theme(aspect.ratio=1)
ggsave("figures/BITM.pdf", width = 9, height = 7)
ggsave("figures/BITM.jpg", width = 9, height = 7)

# ExG
all_preds %>% 
  filter(trait == "Greenness index (ExG)") %>% 
  ggplot(aes(y = fit, x = doy)) +
  geom_line(linewidth = 0.6) +
  geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit, fill = `Growth form`), 
              linetype=0, alpha=0.8) +
  geom_point(data = d %>% filter(trait == "Greenness index (ExG)"),
             aes(y = value), col = "gray70", size = 1) +
  scale_fill_manual(values = gro_cols) +
  facet_wrap(vars(species), scales = "free_y") +
  ylab("Greenness index (ExG)") +
  xlab("Day of year") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text = element_text(face = "italic"),
        legend.title = element_blank()) +
  geom_text(data = r2s %>% filter(trait == "Greenness index (ExG)") %>% mutate(r2 = paste0(round(r2, 2)," ", significance)),
            aes(x = -Inf, y = -Inf, hjust = -0.4, vjust = -0.8, label = r2), size = 3) +
  theme(aspect.ratio=1)
ggsave("figures/ExG.pdf", width = 9, height = 7)
ggsave("figures/ExG.jpg", width = 9, height = 7)

# Plots by species
for(sp in specs){
  
  gg <- all_preds %>% 
    filter(species == sp) %>% 
    ggplot(aes(y = fit, x = doy)) +
    geom_line(linewidth = 0.6) +
    geom_ribbon(aes(ymin = fit - se.fit, ymax = fit + se.fit), 
                linetype=0, alpha=0.7, fill = "lightblue") +
    geom_point(data = d %>% filter(species == sp),
               aes(y = value), col = "gray70", size = 1) +
    facet_wrap(vars(trait), scales = "free_y") +
    ylab("Trait value") +
    xlab("Day of year") +
    theme_bw() +
    geom_text(data = r2s %>% filter(species == sp) %>%  mutate(r2 = paste0(round(r2, 2)," ", significance)),
              aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1.8, label = r2), size = 3) +
    ggtitle(sp) +
    theme(aspect.ratio=1) +
    theme(plot.title = element_text(face = "italic"))
  
  ggsave(paste0("figures/xxx_species_",gsub(" ","_",sp),".pdf"), plot = gg, width = 7, height = 5)
  ggsave(paste0("figures/xxx_species_",gsub(" ","_",sp),".jpg"), plot = gg, width = 7, height = 5)
  
}


r2s %>% 
  mutate(cvar = (maxv-minv)/medianv) %>% 
  group_by(trait) %>% 
  summarise(r2 = mean(r2),
            cvar = mean(cvar)) %>% arrange(cvar)

r2s %>% 
  mutate(cvar = (maxv-minv)/medianv) %>% 
  # filter(!trait %in% c("BITM","ExG")) %>% 
  group_by(species) %>% 
  summarise(r2 = mean(r2),
            cvar = mean(cvar)) %>% 
  arrange(desc(cvar))

all_preds %>% 
  mutate(fit = ifelse(trait %in% c("Leaf area", "Leaf dry mass", "Specific leaf area (SLA)"), log(fit), fit)) %>% 
  ggplot(aes(y = fit, x = doy, group = species, color = species)) +
  scale_color_manual(values = spe_cols) +
  geom_line(size = 1.5, alpha = 0.7) +
  theme_bw() +
  facet_wrap(vars(trait), scales = "free_y") +
  ylab("Trait value") + xlab("Day of year") +
  theme(legend.text = element_text(face = "italic"))
ggsave("figures/seasonal_trends_combined.pdf", width = 10, height = 5)
ggsave("figures/seasonal_trends_combined.jpg", width = 10, height = 5)



######################################################
# Functional groups

d <- d %>% 
  mutate(
    growth_form = case_when(
      species %in% c("Bistorta vivipara","Solidago virgaurea","Trollius europaeus","Viola biflora") ~ "Forb",
      species %in% c("Betula nana","Salix lapponum","Salix reticulata","Vaccinium myrtillus","Vaccinium uliginosum") ~ "Deciduous shrub",
      species %in% c("Dryas octopetala","Vaccinium vitis-idaea") ~ "Evergreen shrub"
    )
  ) %>% 
  mutate(growth_form = factor(growth_form, levels = c("Forb",
                                                          "Deciduous shrub",
                                                          "Evergreen shrub")))

p_df_all <- tibble()
for(tr in traits){
  
  dt <- d %>% 
    filter(trait == tr) %>% 
    mutate(species = factor(species),
           site = factor(site),
           growth_form = factor(growth_form))
  
  m <- gam(value ~ s(doy, by = species, k = 4) + s(doy, by = growth_form, k = 4) + s(species, bs = 're') + site, data = dt, method = "REML")
  
  ms <- summary(m)
  
  p_df <- ms$s.table %>% 
    as.data.frame() %>% 
    rownames_to_column("growth_form") %>% 
    filter(grepl("growth_form", growth_form)) %>% 
    mutate(growth_form = gsub("s\\(doy):growth_form","",growth_form)) %>% 
    mutate(trait = tr) %>% 
    relocate(trait) %>% 
    rename(pvalue = `p-value`)
  
  p_df_all <- bind_rows(p_df_all, p_df)
  
}

p_df_all <- p_df_all %>% 
  mutate(significance = p_value_to_stars(pvalue))

p_df_all %>% 
  pivot_wider(id_cols = growth_form, names_from = trait, values_from = pvalue)


p_df_all %>% 
  pivot_wider(id_cols = growth_form, names_from = trait, values_from = significance)

p_df_all %>% 
  mutate(pvalue = as.character(round(pvalue, 4))) %>% 
  mutate(pvalue = ifelse(pvalue == "0", "< 0.001", pvalue)) %>% 
  mutate(pvalue = paste(pvalue, " ", significance)) %>% 
  pivot_wider(id_cols = growth_form, names_from = trait, values_from = pvalue) %>% 
  write_csv2("output/growthform_pvalues.csv")
