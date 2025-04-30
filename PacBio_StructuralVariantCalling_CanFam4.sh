#!/bin/bash

#SBATCH --job-name=PacBioSV
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=23:59:00
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=edlarsen@colostate.edu
#SBATCH --output=PacBioSV_log_%j.txt

# Conda environment with pbmm2, pbsv, bcftools, and pysam (a depdendency of svpack) should be activated prior to running

# Print versions
pbmm2 --version
pbsv --version
bcftools --version

##################################################
##### DOWNLOAD REFERENCE GENOME FASTA FILES #####
##################################################
# CanFam4
rsync -azvP rsync://ftp.ensembl.org/ensembl/pub/release-109/fasta/canis_lupus_familiarisgsd/dna/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa.gz .

# Check md5sum of FASTA files
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-109/fasta/canis_lupus_familiarisgsd/dna/CHECKSUMS .
grep ".dna.toplevel" CHECKSUMS
sum Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa.gz

# Extract compressed FASTA files
gunzip Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa.gz

############################################################
##### ALIGN HIFI READS TO REFERENCE GENOME WITH PBMM2 #####
############################################################
# Create reference genome index
pbmm2 index Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa CanFam4.mmi --preset HIFI

# Align HiFi reads
infiles="../Data-X202SC25030830-Z01-F001/CF215206/Sequel-Revio/FPAC250144567-1A/CF215206.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF215393/Sequel-Revio/FPAC250144568-1A/CF215393.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF215422/Sequel-Revio/FPAC250144569-1A/CF215422.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF215442/Sequel-Revio/FPAC250144570-1A/CF215442.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF217246/Sequel-Revio/FPAC250144571-1A/CF217246.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF51135/Sequel-Revio/FPAC250144575-1A/CF51135.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF51406/Sequel-Revio/FPAC250144574-1A/CF51406.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF53469/Sequel-Revio/FPAC250144573-1A/CF53469.hifi_reads.bam ../Data-X202SC25030830-Z01-F001/CF54434/Sequel-Revio/FPAC250144572-1A/CF54434.hifi_reads.bam"

for file in $infiles
do
name=$(basename "${file}" | awk -F'.' '{print $(NF-2)}'); echo ${name}
pbmm2 align CanFam4.mmi ${file} ../output/${name}.CanFam4.aligned.bam --preset HIFI --sort --j 16 -J 8 --log-level INFO --sample ${name}
done

echo "pbmm2 steps completed"

###############################################
##### CALL STRUCTURAL VARIANTS WITH PBSV #####
###############################################
# Acquire bed file of tandem repeat locations for reference genome
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam4/bigZips/canFam4.trf.bed.gz'
gunzip canFam4.trf.bed.gz

# Discover signatures of structural variation
for file in $(ls ../output/*.CanFam4.aligned.bam)
do
name=$(basename ${file} .aligned.bam); echo ${name}
pbsv discover --tandem-repeats canFam4.trf.bed ${file} ../output/${name}.svsig.gz
done

# Call structural variants and assign genotypes
# get list of svsig.gz files
Cfam4_infiles=$(ls ../output/*.CanFam4.svsig.gz)

# Run pbsv call jointly for all samples of interest
pbsv call --ccs -j 32 Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa ${Cfam4_infiles} ../output/PTCLandCTRL_StructuralVariants_CanFam4.vcf

echo "pbsv steps completed"

############################################
##### FILTER & ANNOTATE VARIANT CALLS #####
###########################################
# install svpack
git clone https://github.com/PacificBiosciences/svpack.git
cd svpack

# Return PASS SVs at least 50 bp long
./svpack filter --pass-only --min-svlen 50 ../../output/PTCLandCTRL_StructuralVariants_CanFam4.vcf > ../../output/PTCLandCTRL_StructuralVariants_CanFam4.filtered.vcf

echo "svpack filter completed"

### Annotate SVs that impact genes
## Retrieve gff3 files for reference genomes
# CanFam 4
rsync -azvP rsync://ftp.ensembl.org/ensembl/pub/release-109/gff3/canis_lupus_familiarisgsd/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.109.gff3.gz .

# Check md5sum
rsync -azvP rsync://ftp.ensembl.org/ensembl/pub/release-109/gff3/canis_lupus_familiarisgsd/CHECKSUMS .
grep "1.0.109.gff3.gz" CHECKSUMS
sum Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.109.gff3.gz

# Extract compressed gff3 files
gunzip Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.109.gff3.gz

## Add BCSQ tag to variants that impact genes

./svpack consequence ../../output/PTCLandCTRL_StructuralVariants_CanFam4.filtered.vcf Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.109.gff3 > ../../output/PTCLandCTRL_StructuralVariants_CanFam4.filtered.ann.vcf

echo "svpack consequence completed"

## Subset only SVs that impact genes
bcftools view -i 'INFO/BCSQ!=""' ../../output/PTCLandCTRL_StructuralVariants_CanFam4.filtered.ann.vcf -o ../../output/PTCLandCTRL_StructuralVariants_CanFam4.filtered.ann.bcsq.vcf

## Subset only breakend variant calls
# Of all BND that passed filtering
bcftools view -i 'INFO/SVTYPE=="BND"' ../../output/PTCLandCTRL_StructuralVariants_CanFam4.filtered.ann.vcf -o ../../output/PTCLandCTRL_BND_Variants.CanFam4.filtered.ann.vcf
# Of only BND that impact genes
bcftools view -i 'INFO/SVTYPE=="BND"' ../../output/PTCLandCTRL_StructuralVariants_CanFam4.filtered.ann.bcsq.vcf -o ../../output/PTCLandCTRL_BND_Variants.CanFam4.filtered.ann.bcsq.vcf