# E Nielsen Dandoroff PhD code
**This GitHub contains all essential code used in E Nielsen Dandoroff's PhD thesis.**

## Chapter 3 code
* in_silico_tools_threshold_code (_in silico_ predictor cutpoint generation)

## Chapter 4 code
### Remote
* VCF_formatting (bulk convert text files to VCF format)
* BayesDel_CADD_scores (annotate variants with BayesDel addAF and CADD Phred scores)
* vcf_to_csv (bulk convert VCF files to CSV format)
* Script_processing_var_scores (filter variants by BayesDel addAF and CADD Phred scores)

### UKB RAP
#### MGSvar cohort creation
* get-snp-data-20230217 (extract UKB patient IDs for MGSvar cohort pt.1)
* get-snp-data-20230405 (extract UKB patient IDs for MGSvar cohort pt.2)
* MGSvar_script3_END (annotation of MGSvar cohort)

#### pathvar cohort creation
* path_script (extract UKB participant IDs for pathvar cohort)
* pathvar_script_END (annotation of pathvar cohort)

#### UKB cohort creation
* EntireCohort_script (UKB cohort data cleaning and creation of derived variables)

#### Statistical analyses
* EntireCohort_MGSvar_analysis
* EntireCohort_pathvar_Biosynthgenes
* EntireCohort_pathvar_BMgenes
* EntireCohort_pathvar_Bonegenes
* EntireCohort_pathvar_Braingenes
* EntireCohort_pathvar_CAgenes
* EntireCohort_pathvar_CCgenes
* EntireCohort_pathvar_CentCytgenes
* EntireCohort_pathvar_DDRgenes
* EntireCohort_pathvar_GHgenes
* EntireCohort_pathvar_LTUgenes
* EntireCohort_pathvar_MGSgenes
* EntireCohort_pathvar_Mitogenes
* EntireCohort_pathvar_RASgenes
* EntireCohort_pathvar_RNAPgenes
* EntireCohort_pathvar_TFgenes
* EntireCohort_pathvar_TGFbBMPgenes
* EntireCohort_pathvar_THgenes
* EntireCohort_pathvar_WNTSHHgenes

#### Investigation of ICD10 diagnosis codes
* EntireCohort_Diagnoses
