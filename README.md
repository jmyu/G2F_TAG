# Reproducible Data&Code for G2F JGRA
Unzip this file first: McNellie_Final_Genotype_g2f_geno_hybrid_imp_v3_numeric_reduced_snpEff.zip
```
Reproducible_Data_Code
¦	README.md
¦	EnvInd_table.txt    
¦	JGRA.function.r
¦	All_Traits_NBT.r
¦	All_Traits_LH82.r
¦	All_Traits_LH185.r
¦	All_Traits_LH198.r
¦	All_Traits_PHB47.r
¦	All_Traits_PHZ51.r
¦	McNellie_Final_Genotype_g2f_geno_hybrid_imp_v3_numeric_reduced_snpEff.csv
+---Pheno_LbE
¦   ¦	LH185_DTA_LbE_table
¦   ¦	LH185_DTS_LbE_table
¦   ¦	LH185_EH_LbE_table
¦   ¦	... 
+---Env_MeanPara
¦       ¦	LH185_diff_EH_envMeanPara_13_36
¦       ¦	LH185_diff_GY_envMeanPara_24_35
¦       ¦	LH185_GDD_DTA_envMeanPara_6_44
¦       ¦	... 
+---result
¦   ¦	LH185_DTA_1to2_out4figure
¦   ¦	LH185_DTA_1to3_out4figure
¦   ¦	LH185_DTA_1to4_out4figure
¦   ¦	...
¦   +---Table
¦       ¦	PredictionAbility_PHB47
¦       ¦	PredictionAbility_PHZ51
¦       ¦	...
+---Pre CERIS Phenotypes
¦   ¦	G2F_Hybrid_Phenotypes_2014_2015.csv
```

##### Testers(1 + 5): NBT, LH82, LH185, LH198, PHB47, PHZ51
##### Traits (6): DTA, DTS, EH, PH, GM, GY

### Datasets
1. **Phenotype** : 
	- Folder `./Pheno_LbE/` includes 36 phenotype files (formated as two-way table)
	- Data source: `../Phenotypes_Post_CERIS/` (NBT_DTA and NBT_DTS reformated)

2. **Environment**: 
	- Folder `./Env_MeanPara/` includes 36 files with the environmental index values for each Tester-Trait
	- Data source: `../env mean para/` (removed all .csv files)
	- **Meta file:** `EnvInd_table.txt` records the environmental indices for 36 Tester-Trait combinations summarized from the 36 file names under `./Env_MeanPara/`

3. **Genotype**:
	- `McNellie_Final_Genotype_g2f_geno_hybrid_imp_v3_numeric_reduced_snpEff.csv`
	- Data source: `../Genotype files/`
 

### Codes

1. JGRA function library: `./JGRA.function.r` 
	- Includes JGRA functions for 1to2, 1to3, and 1to4 predictions
	- Functions are modified based on `../CERIS JGRA/CERIS JGRA functions.R/` to match the data format of this project

2. Run JGRA for each tester:
	- 6 Rscript files named as `ALLTraits_Tester.r`
	- For each Tester-Trait-PredictionScenario, the last reshuffle result will be saved under `./result/` for plotting
	- For each tester, the summarized prediction ability (across and within) will be saved under `./result/Table/`
3. Wrap up results of all testers
	- `Results2Table.r`
	- Output: `PredictiveAbility_across.csv` & `PredictiveAbility_within.csv`


