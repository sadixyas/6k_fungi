setwd("~/bigdata/6000fungiproject") #since I have most of my files here it's easier to run from this directory

library("DBI")
library("duckdb")
library("ggplot2")
library("dplyr")

# connect to the duckdb database, if I want to make any changes there then I can write read_only=FALSE
con <- dbConnect(duckdb(), dbdir="functionalDB/function.duckdb", read_only=TRUE)


secreted_counts <- dbGetQuery(con, "
  SELECT locustag, COUNT(*) AS secreted_n
    FROM secreted_proteins
   GROUP BY locustag
")
total_counts <- dbGetQuery(con, "
  SELECT locustag, COUNT(*) AS total_n
    FROM gene_proteins
   GROUP BY locustag
")

secreted_meta <- dbReadTable(con, "secreted_proteins") %>%
                    select(locustag, PHYLUM, trophicMode) %>% distinct()

dbDisconnect(con)

# assemble and compute percentage
df <- merge(secreted_counts, total_counts, by="locustag")
df$percent_secreted <- df$secreted_n / df$total_n * 100

# sort descending
df <- df[order(df$percent_secreted, decreasing=TRUE), ]

df2 <- df %>%
  left_join(secreted_meta, by="locustag") %>%
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

# drop any rows with no PHYLUM
df_phylum <- df2 %>%
  filter(!is.na(PHYLUM) & PHYLUM != "")

df_phylum2 <- df2 %>%
  # drop three extra phyla because they don't have enough data for visualization
  filter(
    !is.na(PHYLUM),
    PHYLUM != "",
    !PHYLUM %in% c("Blastocladiomycota",
                   "Cryptomycota",
                   "Sanchytriomycota")
  )

#  plot
ggplot(df, aes(x = reorder(locustag, percent_secreted),
               y = percent_secreted)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Genome",
    y = "Secreted proteins (%)",
    title = "Proportion of Secreted Proteins per Genome"
  ) +
  theme_minimal(base_size = 14)


ggplot(df2, aes(x = reorder(locustag, percent_secreted), y = percent_secreted)) +
  geom_col(color = "darkblue") +
  coord_flip() +
  labs(
    x = "Genome",
    y = "Secreted proteins (%)",
    title = "Proportion of Secreted Proteins per Genome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.y = element_text(size = 6)
  )

# Histogram by Phylum
# Horizontal, faceted histogram by Phylum
ggplot(df_phylum, aes(x = percent_secreted, fill = PHYLUM)) +
  geom_histogram(color = "black", bins = 30) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~PHYLUM, scales = "fixed", ncol = 2) +
  coord_flip() +
  labs(
    x = "Secreted proteins (%)",
    y = "Number of Genomes",
    title = "Secretome Proportion by Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )

group1 <- c("Ascomycota", "Basidiomycota")
group2 <- c("Chytridiomycota", "Microsporidia")
group3 <- c("Mucoromycota", "Zoopagomycota")

df_AB <- df_phylum2 %>% filter(PHYLUM %in% group1)
df_CM <- df_phylum2 %>% filter(PHYLUM %in% group2)
df_MZ <- df_phylum2 %>% filter(PHYLUM %in% group3)

ggplot(df_AB, aes(x = secreted_n, fill = PHYLUM)) +
  geom_histogram(color = "black", bins = 30) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~PHYLUM, ncol = 1, scales = "fixed") +
  labs(
    x = "Number of Secreted Proteins",
    y = "Number of Genomes",
    title = "Secretome Size in Ascomycota vs Basidiomycota"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")

ggplot(df_CM, aes(x = secreted_n, fill = PHYLUM)) +
  geom_histogram(color = "black", bins = 30) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~PHYLUM, ncol = 1, scales = "fixed") +
  coord_flip() +
  labs(
    x = "Number of Secreted Proteins",
    y = "Number of Genomes",
    title = "Secretome Size in Chytridiomycota vs Microsporidia"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")

ggplot(df_MZ, aes(x = secreted_n, fill = PHYLUM)) +
  geom_histogram(color = "black", bins = 30) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~PHYLUM, ncol = 1, scales = "fixed") +
  coord_flip() +
  labs(
    x = "Number of Secreted Proteins",
    y = "Number of Genomes",
    title = "Secretome Size in Mucoromycota vs Zoopagomycota"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")

# Histogram by Trophic-Mode Group
ggplot(df2, aes(x = percent_secreted, fill = Trophic_mode_grouped)) +
  geom_histogram(color = "black", bins = 30) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~Trophic_mode_grouped, scales = "fixed", ncol = 2) +
  coord_flip() +
  labs(
    x = "Secreted proteins (%)",
    y = "Number of Genomes",
    title = "Secretome Proportion by Trophic Mode"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  )

ggplot(df_phylum2, aes(x = reorder(locustag, secreted_n), y = secreted_n, fill = PHYLUM)) +
  geom_col(color = "black") +
  coord_flip() +
  labs(
    x = "Genome",
    y = "Number of Secreted Proteins",
    title = "Secretome Size per Genome by Phylum"
  ) +
  facet_wrap(~PHYLUM, ncol=1, scales="free_y") +
  theme_minimal()
