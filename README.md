# RNA Ligand-Based Similarity Analysis


> [!NOTE] 
> This repository contains data and scripts for analyzing ligand-based similarity metrics and molecular fingerprints in the context of RNA-targeted small molecule discovery. The content is related to scientific publication
> 
> Dávid Bajusz, Anita Rácz, Filip Stefaniak, *Evaluation of single-template ligand-based methods for  the discovery of small-molecule nucleic acid binders.*
> 

## Repository Structure

- **1-fingerprints/**  
  Contains molecular fingerprint datasets in compressed CSV format. Each file represents a different fingerprinting method (e.g., CDK ECFP, FCFP, MACCS, PubChem, Indigo, MAP4, OpenBabel, RDKit, random bits).

- **3-metrics/**  
  Contains statistical analyses and metrics, including t-tests comparing the performance of different methods and random baselines.  
  - `metrics__ttest.csv`: Summary of statistical tests for various methods.
  - `combined_best_metrics.csv.gz`: combined values for all methods and all similarity metrics
  - `individual_metrics/`: Individual metrics for fingerprints (all subdirectories) and individual ligand-based methods (`tanimoto` directory). The best-performing variants of all methods are in `_best` directory.

- **4-ANOVA/**  
  Contains results of ANOVA analyses and interaction means for different factors:
  - `F1_weighted_means.csv`, `F2_weighted_means.csv`, `F3_weighted_means.csv`: Weighted means for different factors.
  - `F1_F2_interaction_means.csv`, `F2_F3_interaction_means.csv`, `F3_F1_interaction_means.csv`: Interaction means between factors.

- **5-case-study/**  
  Contains case studies data for Glutamine_RS.

## Citation

If you use this repository or its data in your research, please cite the relevant publication or acknowledge the authors.
