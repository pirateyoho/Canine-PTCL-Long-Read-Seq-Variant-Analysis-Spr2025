# Bioinformatics pipeline for calling structural variants from canine PTCL PacBio long-read sequencing data
## Background
This repository contains scripts that were used for calling large structural variants in PacBio HiFi long-read sequencing of DNA from canine peripheral T-cell lymphoma (PTCL) samples and control samples of whole lymph node from dogs without lymphoma. This repository is intended for internal use by members of the Clinical Hematopathology Laboratory at Colorado State University and their collaborators. 
## References
Adapted from the pbmm2 (https://github.com/PacificBiosciences/pbmm2) and pbsv (https://github.com/PacificBiosciences/pbsv) documentation.
## Raw data
This pipeline utilized bam files of PacBio Hifi long-read sequencing data from 5 canine PTCL samples and 4 control lymph node samples. This data is available from the Avery lab Nas shared drive:
"M:\CHLab data\Sequencing Data\250414_CD4PTCL_PacBioLongReadSeq_Novogene_Owens"
## Software
A conda (version 23.7.4) environment containing the following packages:
* pbmm2 version 1.14.99
* pbsv version 2.9.0
* bcftools version 1.21
* pysam version 0.22.1

R (version 4.4.0) and R Studio (version 2024.09.0) with the following packages loaded:
* ggplot2 (version 3.5.1)
* dplyr (version 1.1.4)
* tidyr (version 1.3.1)
* knitr (version 1.49)
* stringr (version 1.5.1)
  
## Pipeline overview
1. Download reference genome FASTA files for CanFam3.1 and CanFam4 from Ensembl.
2. Build an index of the reference genome with *pbmm2 index*.
3. Align HiFi reads to reference genome with *pbmm2 align*.
4. Identify signatures of structural variation with *pbsv discover*.
5. Call structural variants from structural variant signatures and assign genotypes with *pbsv call*.
6. Filter and annotate variant calls with *svpack* and *bcftools*.
7. Final output: Filtered and annotated structural variant VCF files.
8. Analysis of VCF files in R.
### Sample information
| Sample#  | Group    | Diagnosis                            | Sex | Breed   | Age (yrs) | Tissue     | HiFi reads                  | HiFi Read Length (mean, bp)                           |
| -------- | -------- | ------------------------------------ | --- | ------- | --------- | ---------- | --------------------------- | ----------------------------------------------------- |
| CF215206 | CD4 PTCL | CD4 PTCL                             | FS  | GLDR    | 6         | Lymph node |                   8,063,210 |                                               17,793  |
| CF215393 | CD4 PTCL | CD4 PTCL                             | FS  | MIX     | 6         | Lymph node |                   7,328,965 |                                               19,176  |
| CF215422 | CD4 PTCL | CD4 PTCL                             | FS  | MIX     | 5         | Lymph node |                   7,233,403 |                                               19,354  |
| CF215442 | CD4 PTCL | CD4 PTCL                             | FS  | SEATERR | 3         | Lymph node |                   7,997,232 |                                               16,706  |
| CF217246 | CD4 PTCL | CD4 PTCL                             | F   | MIX     | 6         | Lymph node |                   8,284,937 |                                               16,300  |
| CF54434  | CTRL     | Necrohemorrhagic gastroenterocolitis | FS  | BRDC    | 4.3       | Lymph node |                 11,926,147  |                                                 8,579 |
| CF53469  | CTRL     | Boxer cardiomyopathy                 | FS  | BOX     | 6         | Lymph node |                   6,968,338 |                                               19,083  |
| CF51406  | CTRL     | Hepatocellular carcinoma             | FS  | LAB     | 14.4      | Lymph node |                   7,190,034 |                                               17,273  |
| CF51135  | CTRL     | Chemodectoma and hemangiosarcoma     | MC  | BOX     | 8.4       | Lymph node |                   8,130,492 |                                               16,208  |
