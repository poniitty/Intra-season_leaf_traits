library(tidyverse)

Sys.setlocale("LC_TIME", "English")

d <- read_csv("data/Seasonal_traits_microclim.csv")

d <- d %>% 
  select(site, date, species, leaf_area, R, G, B) %>% 
  mutate(date = as.Date(ifelse(date == "2024-07-24", as.Date("2024-07-23"), date)),
         date = as.Date(ifelse(date == "2024-08-05", as.Date("2024-08-06"), date)),
         date = as.Date(ifelse(date == "2024-09-09", as.Date("2024-09-10"), date))) %>% 
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

unique(d$date) %>% sort

d %>% 
  group_by(site, date, species) %>% 
  summarise(across(c(leaf_area, R, G, B), ~mean(.x, na.rm = TRUE)))

d$rgb <- rgb(d$R, d$G, d$B)

d %>% 
  ggplot(aes(y = species, x = as.character(date))) +
  geom_point(size = 12, col = d$rgb) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(face = "italic")) +
  ylab(NULL) + xlab(NULL) +
  scale_y_discrete(limits=rev)
ggsave("figures/leaf_colors.pdf", width = 7, height = 5)
ggsave("figures/leaf_colors.jpg", width = 7, height = 5)
