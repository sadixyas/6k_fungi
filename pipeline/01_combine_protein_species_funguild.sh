#!/usr/bin/bash -l

#SBATCH -p short
#SBATCH -C ryzen
#SBATCH --mem=128gb
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --out=logs/combine_tables_%A_%a.log


IDB="/bigdata/stajichlab/shared/projects/MosesLab_projects/Fungi_5k/functionalDB/function.duckdb"
ODB="/bigdata/stajichlab/sadikshs/6000fungiproject/functionalDB/function.duckdb"

duckdb "${ODB}" <<SQL
-- Attach the source DuckDB file 
ATTACH '${IDB}' AS idb;

-- Create (or replace) the combined table in the output database
CREATE OR REPLACE TABLE combined_gene_species_funguild AS
SELECT
  gp.*,
  s.* EXCLUDE (LOCUSTAG),       -- drop duplicate column
  f.* EXCLUDE (species_prefix)
FROM idb.gene_proteins AS gp
JOIN idb.species        AS s  ON gp.locustag        = s.LOCUSTAG
JOIN idb.funguild       AS f  ON s.LOCUSTAG   = f.species_prefix
;

-- Clean up
DETACH idb;
SQL

echo "😊 combined_gene_species_funguild is now in '${ODB}'."

