# Epilepsy-Channel Pocket Atlas (pilot)

Do pathogenic missense variants in epilepsy ion channels fall in computationally druggable pockets more often than expected?
Pilot on **Nav1.1 (SCN1A), Nav1.2 (SCN2A), KCNQ2**, with pockets detected by **fpocket on the Vecura platform**.

**Status: exploratory pilot, not a validated result.**

## Pipeline
1. **Pockets**: fpocket (Vecura) on PDB 7DTD (Nav1.1), 6J8E (Nav1.2), 7CR0 (KCNQ2). Druggable = druggability score >= 0.5. Outputs are in `data/pockets/`.
2. **Lining residues**: protein residues with any heavy atom within 4 A of a pocket's alpha-spheres (`scripts/01`, `02`).
3. **Variants**: ClinVar missense (pathogenic/likely pathogenic = case positions; benign/likely benign = sensitivity only) and gnomAD v4 missense seen in >= 2 alleles (controls). Pulled in the browser with `browser/*.js`, stored as hex bitsets of unique residue positions in `data/variants/`.
4. **Tests** (`scripts/03`-`06`): residue-level (unique positions), Fisher exact, logistic regression adjusted for relative solvent accessibility, and permutation within burial bins (and sequence windows).
5. **Positive control**: fpocket's best KCNQ2 pocket lines 7/13 retigabine-contact residues in 7CR2 (`browser/retigabine_contacts.js`).

## Run
```
pip install -r requirements.txt
cd scripts
python3 03_build_table.py && python3 04_clinvar_background_test.py && python3 05_region_matched_test.py && python3 06_gnomad_control_test.py
```

## Pilot results (pathogenic positions, pocket vs elsewhere, gnomAD controls)
| Channel | In pocket | Elsewhere | p (burial-adj.) | p (burial + region) |
|---|---|---|---|---|
| Nav1.1 | 96% | 72% | 0.018 | 0.007 |
| Nav1.2 | 96% | 59% | 0.009 | 0.0004 |
| KCNQ2 | 89% | 63% | 0.0009 | 0.31 |

## Limitations
- Pocket residues are buried and conserved; burial adjustment only partly separates "pocket" from "constrained".
- ClinVar is biased towards tested hotspots. Few pocket residues carry any gnomAD variant.
- Nav structures contain the alpha subunit only (no membrane or beta subunits). KCNQ2 pocket list was assembled from druggable pocket classes, not a full table.
- Numbering between ClinVar/gnomAD and structures is assumed to match; spot-checked only.
- Nav1.7/lamotrigine (8THH) control not yet run. Novelty search was brief.
