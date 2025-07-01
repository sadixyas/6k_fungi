# my starting directory for all this analysis
setwd("~/bigdata/6000fungiproject")

# load libraries
library(DBI)
library(duckdb)
library(dplyr)
library(tidyr)    
library(ggplot2)

# connect to DuckDB
con <- dbConnect(duckdb(), dbdir="functionalDB/function.duckdb", read_only=TRUE)

# fetch total protein counts & avg length per genome
total_stats <- dbGetQuery(con, "
  SELECT locustag, COUNT(*) AS total_n, AVG(length) AS avg_length
  FROM gene_proteins
  GROUP BY locustag
")

# fetch secreted protein counts & avg length per genome
secreted_stats <- dbGetQuery(con, "
  SELECT locustag, COUNT(*) AS secreted_n, AVG(length) AS avg_secreted_length
  FROM secreted_proteins
  GROUP BY locustag
")

# fetch PHYLUM from species table
species_meta <- dbReadTable(con, "species") %>%
  rename(locustag = LOCUSTAG) %>%
  select(locustag, PHYLUM) %>%
  distinct()

# fetch trophicMode from funguild, renaming species_prefix to locustag
funguild_meta <- dbReadTable(con, "funguild") %>%
  rename(locustag = species_prefix) %>%
  select(locustag, trophicMode) %>%
  distinct()

#I need to find a shorter way to connect all these tables. I tried doing it with the combined table but could not load it

# disconnect from duckdb
dbDisconnect(con)

# assemble and clean
df <- total_stats %>%
  left_join(secreted_stats, by="locustag") %>%
  replace_na(list(secreted_n = 0, avg_secreted_length = NA_real_)) %>%
  mutate(
    non_secreted_n          = total_n - secreted_n,
    avg_non_secreted_length = (total_n * avg_length -
                                 secreted_n * avg_secreted_length) /
      non_secreted_n
  ) %>%
  left_join(species_meta,  by="locustag") %>%
  left_join(funguild_meta, by="locustag") %>%
  # drop any missing PHYLUM or trophicMode
  filter(!is.na(PHYLUM) & PHYLUM != "",
         !is.na(trophicMode) & trophicMode != "") %>%
  # group trophic modes (since there is no clear rule on the groups, this is how I am doing it)
  mutate(
    Trophic_mode_grouped = case_when(
      trophicMode %in% c(
        "Pathotroph-Pathotroph-Saprotroph",
        "Saprotroph-Pathotroph-Saprotroph",
        "Pathotroph-Saprotroph"
      ) ~ "Pathotroph-Saprotroph",
      trophicMode %in% c(
        "Saprotroph-Saprotroph-Symbiotroph",
        "Symbiotroph-Saprotroph-Symbiotroph",
        "Saprotroph-Symbiotroph"
      ) ~ "Saprotroph-Symbiotroph",
      trophicMode %in% c(
        "Pathotroph-Pathotroph-Saprotroph-Symbiotroph",
        "Pathotroph-Saprotroph-Symbiotroph"
      ) ~ "Pathotroph-Saprotroph-Symbiotroph",
      TRUE ~ trophicMode
    )
  )

# split into subsets for creating different dataframes
df_all <- df
df_sec <- df        %>% filter(secreted_n    > 0)
df_non <- df        %>% filter(non_secreted_n > 0)
# for phylum‐faceted plots, drop two small phyla because they show only one point in scatterplot
df_all_phy <- df_all %>% filter(!PHYLUM %in% c("Cryptomycota", "Sanchytriomycota"))
df_sec_phy <- df_sec %>% filter(!PHYLUM %in% c("Cryptomycota", "Sanchytriomycota"))
df_non_phy <- df_non %>% filter(!PHYLUM %in% c("Cryptomycota", "Sanchytriomycota"))

# compute global Pearson r
r_all_global <- cor(df_all$avg_length,              df_all$total_n,        use="complete.obs")
r_sec_global <- cor(df_sec$avg_secreted_length,     df_sec$secreted_n,     use="complete.obs")
r_non_global <- cor(df_non$avg_non_secreted_length, df_non$non_secreted_n, use="complete.obs")

# compute groupwise Pearson r
r_all_phylum   <- df_all_phy %>%   group_by(PHYLUM)               %>% summarize(r = cor(avg_length, total_n, use="complete.obs"))
r_all_trophic  <- df_all %>%       group_by(Trophic_mode_grouped) %>% summarize(r = cor(avg_length, total_n, use="complete.obs"))
r_sec_phylum   <- df_sec_phy %>%   group_by(PHYLUM)               %>% summarize(r = cor(avg_secreted_length, secreted_n, use="complete.obs"))
r_sec_trophic  <- df_sec %>%       group_by(Trophic_mode_grouped) %>% summarize(r = cor(avg_secreted_length, secreted_n, use="complete.obs"))
r_non_phylum   <- df_non_phy %>%   group_by(PHYLUM)               %>% summarize(r = cor(avg_non_secreted_length, non_secreted_n, use="complete.obs"))
r_non_trophic  <- df_non %>%       group_by(Trophic_mode_grouped) %>% summarize(r = cor(avg_non_secreted_length, non_secreted_n, use="complete.obs"))

# theme and annotation helpers so I won't have to do it for all
base_theme <- theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color="black", fill=NA))

add_global_r <- function(r_val) {
  annotate("text", x=Inf, y=Inf,
           label=paste0("r = ", round(r_val, 2)),
           hjust=1.1, vjust=1.1, size=5)
}

add_group_r <- function(r_df) {
  geom_text(
    data = r_df,
    aes(x = Inf, y = Inf, label = paste0("r = ", round(r, 2))),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.1, size = 4
  )
}

# plots

# All proteins
p_all_global <- ggplot(df_all, aes(avg_length, total_n, color = PHYLUM)) +
  geom_point(size=3, alpha=0.8) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_global_r(r_all_global) +
  labs(title="All Proteins", x="Average protein length (aa)", y="Total proteins") +
  base_theme

p_all_by_phylum <- ggplot(df_all_phy, aes(avg_length, total_n, color = PHYLUM)) +
  geom_point(size=2, alpha=0.7) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_group_r(r_all_phylum) +
  facet_wrap(~PHYLUM, scales="free") +
  labs(title="All Proteins by Phylum", x="Average protein length (aa)", y="Total proteins") +
  base_theme

p_all_by_trophic <- ggplot(df_all, aes(avg_length, total_n, color = PHYLUM)) +
  geom_point(size=2, alpha=0.7) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_group_r(r_all_trophic) +
  facet_wrap(~Trophic_mode_grouped, scales="free") +
  labs(title="All Proteins by Trophic Mode", x="Average protein length (aa)", y="Total proteins") +
  base_theme

# Secreted proteins
p_sec_global <- ggplot(df_sec, aes(avg_secreted_length, secreted_n, color = PHYLUM)) +
  geom_point(size=3, alpha=0.8) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_global_r(r_sec_global) +
  labs(title="Secreted Proteins", x="Average secreted protein length (aa)", y="Secreted protein count") +
  base_theme

p_sec_by_phylum <- ggplot(df_sec_phy, aes(avg_secreted_length, secreted_n, color = PHYLUM)) +
  geom_point(size=2, alpha=0.7) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_group_r(r_sec_phylum) +
  facet_wrap(~PHYLUM, scales="free") +
  labs(title="Secreted Proteins by Phylum", x="Average secreted protein length (aa)", y="Secreted protein count") +
  base_theme

p_sec_by_trophic <- ggplot(df_sec, aes(avg_secreted_length, secreted_n, color = PHYLUM)) +
  geom_point(size=2, alpha=0.7) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_group_r(r_sec_trophic) +
  facet_wrap(~Trophic_mode_grouped, scales="free") +
  labs(title="Secreted Proteins by Trophic Mode", x="Average secreted protein length (aa)", y="Secreted protein count") +
  base_theme

# Non-secreted proteins
p_non_global <- ggplot(df_non, aes(avg_non_secreted_length, non_secreted_n, color = PHYLUM)) +
  geom_point(size=3, alpha=0.8) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_global_r(r_non_global) +
  labs(title="Non-Secreted Proteins", x="Average non-secreted protein length (aa)", y="Non-secreted protein count") +
  base_theme

p_non_by_phylum <- ggplot(df_non_phy, aes(avg_non_secreted_length, non_secreted_n, color = PHYLUM)) +
  geom_point(size=2, alpha=0.7) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_group_r(r_non_phylum) +
  facet_wrap(~PHYLUM, scales="free") +
  labs(title="Non-Secreted Proteins by Phylum", x="Average non-secreted protein length (aa)", y="Non-secreted protein count") +
  base_theme

p_non_by_trophic <- ggplot(df_non, aes(avg_non_secreted_length, non_secreted_n, color = PHYLUM)) +
  geom_point(size=2, alpha=0.7) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") +
  add_group_r(r_non_trophic) +
  facet_wrap(~Trophic_mode_grouped, scales="free") +
  labs(title="Non-Secreted Proteins by Trophic Mode", x="Average non-secreted protein length (aa)", y="Non-secreted protein count") +
  base_theme

# render all plots
print(p_all_global)
print(p_all_by_phylum)
print(p_all_by_trophic)

print(p_sec_global)
print(p_sec_by_phylum)
print(p_sec_by_trophic)

print(p_non_global)
print(p_non_by_phylum)
print(p_non_by_trophic)
