#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH -C ryzen
#SBATCH --mem=128gb
#SBATCH -c 96
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --out=logs/create_secreted_proteins.log

#I don't need to do module load every time because I already have that in my profile

# Define DuckDB paths
IDB="/bigdata/stajichlab/shared/projects/MosesLab_projects/Fungi_5k/functionalDB/function.duckdb"
ODB="/bigdata/stajichlab/sadikshs/6000fungiproject/functionalDB/function.duckdb"

duckdb "$ODB" <<EOF
-- Attach the source database under alias "idb" because my source and destination files are different 
ATTACH '${IDB}' AS idb;

-- Create (or replace) the secreted_proteins table, excluding duplicate protein_id columns
CREATE OR REPLACE TABLE secreted_proteins AS
SELECT
  c.*,
  sp.*  EXCLUDE(species_prefix, protein_id),
  tp.*  EXCLUDE(species_prefix, protein_id),
  tm.*  EXCLUDE(species_prefix, protein_id),
  wp.*  EXCLUDE(species_prefix, protein_id)
FROM combined_gene_species_funguild AS c
  JOIN idb.signalp   AS sp ON c.transcript_id = sp.protein_id
    AND sp.probability > 0.8
  JOIN idb.targetp   AS tp ON c.transcript_id = tp.protein_id
    AND tp.prediction = 'SP'
    AND tp.cleavage_probability > 0.5
  JOIN idb.tmhmm     AS tm ON c.transcript_id = tm.protein_id
    AND tm.PredHel <= 1
    AND (tm.PredHel = 0 OR tm.First60 > 0)
  JOIN wolfpsort_extr AS wp ON c.transcript_id = wp.protein_id
    AND wp.localization = 'extr'
    -- antijoin to drop any transcript_ids found in prosite
  LEFT JOIN idb.prosite AS pr ON c.transcript_id = pr.protein_id
WHERE pr.protein_id IS NULL;

-- Detach the source DB
DETACH idb;
EOF

echo " 😊 secreted_proteins table created in ${ODB}"

