# devtools::install_github("BlakeRMills/MetBrewer")
library(tidyverse)
library(MetBrewer)

Sys.setlocale("LC_TIME", "English")

d <- read_csv("data/Seasonal_traits_microclim.csv") %>% 
  mutate(Elevation = case_match(
    site,
    "RA009" ~ "Low",
    "MI259" ~ "Mid",
    "RA068" ~ "High"
  ))

d <- d %>% 
  mutate(date = as.Date(ifelse(date == "2024-07-24", as.Date("2024-07-23"), date)),
         date = as.Date(ifelse(date == "2024-08-05", as.Date("2024-08-06"), date)),
         date = as.Date(ifelse(date == "2024-09-09", as.Date("2024-09-10"), date)))

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

d <- d %>% 
  group_by(date, doy, species, trait) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

traits <- unique(d$trait) %>% sort()

cols <- as.character(met.brewer("Redon", 11, type = "discrete", return_hex = TRUE))

d %>% 
  filter(date != "2024-09-17") %>% 
  group_by(trait, date) %>% 
  arrange(trait, date, desc(value)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>% 
  ggplot(aes(y = rank, x = doy, group = species, color = species)) +
  # geom_point() +
  geom_line(size = 3, alpha = 0.7) +
  scale_color_manual(values = cols) +
  scale_y_reverse() +
  theme_bw() +
  facet_wrap(vars(trait)) +
  theme(legend.text = element_text(face = "italic"))
ggsave("figures/seasonal_ranking.pdf", width = 9, height = 5)
ggsave("figures/seasonal_ranking.jpg", width = 9, height = 5)

for(tr in traits){
  
  dt <- d %>% 
    filter(trait == tr)
  
  crs <- dt %>% 
    pivot_wider(id_cols = species, names_from = date, values_from = value) %>% 
    select(-species) %>% 
    cor(., method = "spearman", use = "pairwise.complete.obs")
  
  diag(crs) <- NA
  
  print(tr)
  print(mean(as.numeric(crs), na.rm = TRUE))
  
}





d %>% 
  mutate(value = ifelse(trait %in% c("Leaf area", "Leaf dry mass", "SLA"), log(value), value)) %>% 
  ggplot(aes(y = value, x = doy, group = species, color = species)) +
  # geom_point() +
  geom_line(size = 1) +
  theme_bw() +
  facet_wrap(vars(trait), scales = "free_y")

d %>% 
  mutate(value = ifelse(trait %in% c("Leaf area", "Leaf dry mass", "SLA"), log(value), value)) %>% 
  ggplot(aes(y = value, x = doy, group = species, color = species)) +
  # geom_point() +
  geom_smooth(method = "gam", size = 1, se = FALSE) +
  theme_bw() +
  facet_wrap(vars(trait), scales = "free_y")
