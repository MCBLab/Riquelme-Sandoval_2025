library(dplyr)
library(ggplot2)
library(tidyverse)
library(vegan)

setwd("/Documents and Settings/diego.coelho/Documents/brain_dev/Riquelme-Sandoval_2025/")

data <- readxl::read_excel("data/Scatter Plot.xlsx", sheet = 2)

write.csv(data, "data/scatter_plot.csv", quote = F, row.names = F)

control_means <- data %>% filter(tipo == "control") %>% summarise(
  mean_prolongamento = mean(prolongamento, na.rm = TRUE),
  mean_index = mean(Index_complex, na.rm = TRUE)
)
data <- data %>% mutate(
  Prolongamento_norm = 50 * prolongamento / control_means$mean_prolongamento,
  Index_complex_norm = 50 * Index_complex / control_means$mean_index
)

data_filt <- data #%>% filter(tipo %in% c("control", "DAG"))

X <- data_filt %>%
  select(Prolongamento_norm, Index_complex_norm)

dist_mat <- vegdist(X, method = "euclidean")

bd <- betadisper(dist_mat, data_filt$tipo, type = "centroid")
dispersion <- anova(bd)

result <- adonis2(dist_mat ~ tipo, data = data_filt, permutations = 10000)

pt <- permutest(bd, pairwise = TRUE, permutations = 10000)

pvals <- as.data.frame(pt$pairwise$permuted)
pvals$Group1 <- rownames(pvals)

pvals_long <- pvals %>%
  pivot_longer(-Group1, names_to = "Group2", values_to = "p") %>%
  filter(!is.na(p)) %>% filter(grepl("control-", Group1))

# Heatmap with p-values

ggplot(pvals_long, aes(Group1, Group2, fill = p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", p)), size = 5) +
  scale_fill_viridis_c(direction = -1) +
  theme_classic() +
  labs(
    title = "Pairwise PERMDISP (betadisper)",
    fill = "p-value",
    x = "",
    y = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot clusters: prolongamento vs. Index_complex, colored by tipo
data %>% filter(tipo %in% c("control", "cortico", "CB1_NEG", "cortico_CB1_NEG")) %>%
ggplot(aes(x = Prolongamento_norm, y = Index_complex_norm)) +
  geom_point(size = 1) +
  stat_density_2d(
    aes(fill = tipo, alpha = after_stat(level)), 
    geom = "polygon", show.legend = F,
    color = NA
  ) +
  theme_classic() +
  labs(y = "Complexity Index", x = "Process Index") +
  facet_wrap(.~ tipo) +
  labs(title = "Clusters by tipo (log y)", x = "Prolongamento", y = "Index_complex (log scale)") +
  scale_y_log10(limits = c(1,1000)) +
  scale_x_continuous(limits = c(0,150)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red")


data %>% filter(tipo %in% c("control", "cortico", "CB1_NEG", "cortico_CB1_NEG")) %>%
  ggplot(aes(x = Prolongamento_norm, y = Index_complex_norm)) +
  geom_point(size = 1) +
  geom_density_2d(
    aes(color = tipo),
    linetype = "solid", show.legend = F
  ) +
  theme_classic() +
  labs(y = "Complexity Index", x = "Process Index") +
  facet_wrap(.~ tipo) +
  labs(title = "Clusters by tipo (log y)", x = "Prolongamento", y = "Index_complex (log scale)") +
  scale_y_log10(limits = c(1,1000)) +
  scale_x_continuous(limits = c(0,150)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red")




nmds <- metaMDS(dist_mat, k = 2, trymax = 100)
scores_df <- as.data.frame(scores(nmds))
scores_df$tipo <- data_filt$tipo

ggplot(scores_df, aes(NMDS1, NMDS2, color = tipo)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_classic() +
  labs(
    title = "NMDS (metaMDS) by tipo",
    subtitle = paste("Stress =", round(nmds$stress, 3))
  )
