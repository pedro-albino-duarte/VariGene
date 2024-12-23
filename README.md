# VariGene pipeline

VariGene is a robust and reproducible pipeline to identify genetic variants associated with hereditary diseases.

## Pipeline inputs

This pipeline accepts any sample from a genomic origin that could be single-end or paired-end. For the pipeline to read the files, they must be in either `.fastq` or a compressed format (`.gz`).

## Workflow

1. **Pre-processing:**
    - Data is pre-processed with the `fastp` tool to ensure the quality of the `.fastq` files before proceeding to alignment against the reference genome (`hg19`).
    - This process generates an HTML and JSON files with quality metrics for user evaluation.

2. **Alignment:**
    - The sample is aligned against the reference genome using the `bwa` tool, which utilizes the Burrow-Wheeler method.
    - The reference genome is accessed through the UCSC genome browser, and indexed files are prepared using the `bwa index` command. The alignment generates a SAM file.
    - To streamline the workflow and save computation time, pre-indexed files are hosted on my personal Google Cloud account and provided with each iteration of the workflow, eliminating the need for repeated indexation.

3. **Conversion to BAM:**
    - SAM files are converted to BAM files (binary format) using `samtools view`, reducing storage space and simplifying manipulation.
    - Sorting and indexing are performed using `samtools sort` and `samtools index`.

4. **Variant Calling:**
    - Variants are identified using the `bcftools mpileup` function, which reads the reference genome and the BAM file.
    - The `bcftools call` function is used to assign multiallelic and variant sites (`-mv`), and the outputs are stored in an uncompressed VCF format (`-Ov`).

5. **Filtering Variants:**
    - VCF files are filtered to retain variants with quality scores ≥ 20 using `bcftools view` with the `-i 'QUAL>20'` command.

6. **MultiQC Report:**
    - A function generates a MultiQC report, providing consolidated metrics from `samtools flagstat` and `samtools idxstats`, as well as coverage and variant-calling metrics.
    - Metrics are stored in a directory named `metrics_dir`.

## Pipeline Setup

1. Manually download or clone the Github repository:
   ```bash
   git clone https://github.com/pedro-albino-duarte/VariGene.git
   ```

2. Navigate to the repository folder:
   ```bash
   cd /VariGene
   ```

3. Build the Docker image:
   ```bash
   docker build -t varigene-pedrod-dockerim .
   ```

4. Create and run the Docker container:
   ```bash
   docker run -d --name varigene-container -v $(pwd):/shared varigene-pedrod-dockerim
   ```

5. Access the container:
   ```bash
   docker exec -u root -it varigene-container /bin/bash
   ```

6. *(Optional)* Download a sample file for testing, and move it to the /VariGene folder, in your local computer:
   ```bash
   wget https://storage.googleapis.com/bucket-challenge-pedrod/SRR13065073_1.fastq -P /VariGene
   ```

7. Inside the Docker container, run the pipeline:
   ```bash
   python3 workflow.py <SAMPLE_FASTQ> /data/hg19.fa <OUTPUT_DIRECTORY>
   ```
   Example:
   ```bash
   python3 workflow.py /shared/SRR13065073_1.fastq /data/hg19.fa .
   ```

8. View the generated MultiQC report:
   Open the file `/multiqc_outputs/multiqc_report.html` in your VariGene local directory.

---

## Response to Question 3

### 3.1 Metrics analysis

#### Mapped and Unmapped Reads

To analyze mapped and unmapped reads, the `samtools flagstat` and `samtools idxstats` functions were used on the BAM files generated after aligning the sample `AMOSTRA_A` to the reference genome.

##### Flagstat Metrics
```bash
samtools flagstat AMOSTRA_A.bam > flagstat.txt
```
##### Idxstats Metrics
```bash
samtools idxstats AMOSTRA_A.bam > idxstats.txt
```

The `flagstat.txt` file summarizes the mapped and unmapped reads, while the `idxstats.txt` file shows the number of mapped and unmapped reads per contig. Results indicate 802,979 mapped reads and 0 unmapped reads. Beyond these metrics, these files also provide another metrics that can be useful for the user.

#### Total Number of Identified Variants

The total number of variants was assessed using the `bcftools stats` function:
```bash
bcftools stats AMOSTRA_A.anotada.vcf > vcf_stats.txt
```
Results revealed a total of 83 variants: 75 SNPs and 8 indels.

#### Total Number of Identified Variants per Gene



##### Filtering Variants
1. Filter low-quality variants:

Initially, the bcftools was utilized to remove the low-quality variants, in this case, those who had a quality > 20:

   ```bash
   bcftools filter -i 'QUAL>20' AMOSTRA_A.anotada.vcf -o filtered.vcf
   ```

2. Prepare the BED file:

Afterwards there was the need to use the SEQUENCING_KIT.bed file to extract variants within specific target regions. However, to do so there is the need to remove the header, in order to acomodate the bcftools requirements to read this file. In addition, the bcftools also requires for the the .bed file to have only the first three columns.
   ```bash
   sed '1d' SEQUENCING_KIT.bed > cleaned_SEQUENCING_KIT.bed
   awk '{print $1, $2, $3}' OFS="\t" cleaned_SEQUENCING_KIT.bed > simplified_SEQUENCING_KIT.bed
   ```

3. Compress and index the VCF file:

To accomodate the necessity of the filtered.vcf file to be compressed with bgzip and indexed with tabix, the following code was applied.
   ```bash
   bgzip -c filtered.vcf > filtered.vcf.gz
   tabix -p vcf filtered.vcf.gz
   ```

4. Extract target region variants:

In consequence, the bcftools view function can be used to extract variants within specific target regions by using the following code:

   ```bash
   bcftools view -R simplified_SEQUENCING_KIT.bed filtered.vcf.gz -o filtered_in_regions.vcf
   ```

5. Count variants per gene:

Considering now the file that was generated (filtered.vcf), a new txt file was created to summarize the variants that were found per gene. The results of this analysis are in the file variants_per_gene.txt

   ```bash
   grep -v "^#" filtered_in_regions.vcf | awk -F '\t' '{split($8,info,";"); for (i in info) if (info[i] ~ /^GENE=/) {split(info[i],gene,"="); print gene[2]}}' | sort | uniq -c > variants_per_gene.txt
   ```




### 3.2 Identifying Variants of Clinical Interest



#### Steps
1. Extract variants within clinical interest regions (from `SEQUENCING_KIT.bed`):
   ```bash
   bcftools view -R simplified_sequencing_KIT.bed filtered.vcf.gz -o filtered_in_regions.vcf
   ```

2. Annotate variants using:
   - **SnpEff**: Provides variant effects and potential impacts.
   - **ClinVar Database**: Identifies variants with known clinical significance.

3. Prioritize variants based on:
   - Pathogenicity (ClinVar).
   - Functional impact (e.g., frameshifts or stop mutations).
   - Disease associations (e.g., OMIM).

4. Assess zygosity to determine homozygous or heterozygous states.

5. Validate findings through lab assays.

### 3.3 Analysis of Pathogenic Variant c.244C>T (p.(Gln82*)) in NTHL1

#### Verifying the Variant's Presence

By searching in ClinVar about the pathogenic variant c.244C>T p.(Gln82*) (https://www.ncbi.nlm.nih.gov/clinvar/variation/192319/?oq=%22NTHL1%22[GENE]+%22c.244C%3ET%22[VARNAME]&m=NM_002528.7(NTHL1):c.244C%3ET%20(p.Gln82Ter)%3Fterm=NTHL1%20c.244C%3ET) relative to the gene NTHL1, it was verified that the ClinVar ID relative to this variation was `192319`. Considering this, to verify the presence of this variant in the sample the following command was applied:

```bash
grep "192319" filtered_in_regions.vcf
```
This revealed the variant on chromosome 16, position `2096239`.

#### Determining Zygosity

To verify its zygosity, the column after the GT field (Genotype), which was the 10th column, was assessed:
```bash
grep "192319" filtered_in_regions.vcf | awk -F '\t' '{print $10}'
```
The result was a genotype of `0/1`, indicating heterozygosity.

#### Clinical Relevance
```bash
grep "192319" filtered_in_regions.vcf | awk -F '\t' '{print $8}'
```
Results confirmed:
- **Pathogenicity:** Pathogenic
- **Diseases:** Hereditary cancer-predisposing syndrome, Familial adenomatous polyposis 3, NTHL1-related disorder.

---
This README file details the setup and application of the VariGene pipeline, addressing specific questions and ensuring clarity for reproducibility.

---

Author: Pedro Duarte

Project Copyright: Germano de Sousa
