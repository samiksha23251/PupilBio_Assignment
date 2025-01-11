# ### Task 1 - Data Handling and Statistical Analysis  ###

## Objective
This task involves analyzing CpG coverage data for different tissues, identifying potential PMPs (Pattern Methylation Patterns) with high specificity for tissue differentiation, and estimating thresholds for sequencing reads. The script leverages statistical methods and machine learning to achieve these goals.

---

## 1. Introduction
This script analyzes methylation patterns in CpG sites to identify biomarkers for tissue differentiation. It uses statistical tests, data visualization, and machine learning models to:
- Calculate coverage statistics.
- Visualize tissue-specific patterns.
- Identify PMPs with high specificity.
- Predict tissue types using Random Forest.
- Study how sequencing depth impacts specificity confidence.

---

## 2. Setup and Requirements
To run the script, ensure the following dependencies are installed:
- **Python**: Version 3.7+
- **Required Libraries**: `pandas`, `matplotlib`, `seaborn`, `scikit-learn`, `scipy`

---

## 3. Dataset Preparation
The input dataset (`PupilBioTest_PMP_revA.csv`) should contain:
- Tissue-specific CpG methylation patterns.
- Coverage columns named `000`, `001`, ..., `111`.

### Dataset Cleaning
- Unusual characters in column names are removed for consistency.
- The cleaned dataset is stored in a copy for further analysis.

---

## 4. Analysis Steps

### 4.1 Coverage Statistics Calculation
- **Goal**: Compute the median, mean, standard deviation (std), and coefficient of variation (CV) of single CpG coverage for each tissue.
- **Method**: Aggregate coverage values from `000` to `111` columns and group data by tissue.
- **Output**: A summary table of coverage statistics per tissue.

---

### 4.2 Visualization of Coverage Statistics
**Plots**:
- **Boxplot**: Visualize the distribution of single CpG coverage across tissues.
- **Violin Plot**: Show detailed coverage patterns for each tissue.
- **Bar Plot**: Summarize the median and CV for each tissue.

---

### 4.3 Biomarker (PMP) Identification
- **Goal**: Identify PMPs with high specificity for tissue differentiation.
- **Steps**:
  1. Perform statistical t-tests for each CpG coordinate between two tissues.
  2. Sort results by p-value and filter for significant PMPs (`p-value < 0.05`).
  3. Rank PMPs by feature importance using a Random Forest model.

---

### 4.4 Machine Learning for Tissue Differentiation
- **Model Used**: Random Forest Classifier
- **Steps**:
  1. Prepare data: Features are coverage columns (`000` to `111`), and labels are tissue types.
  2. Split dataset into training (70%) and testing (30%) sets.
  3. Train the model and evaluate using classification metrics (e.g., accuracy, confusion matrix).

---

### 4.5 VRF Analysis
- **VRF Calculation**: Variant Read Fraction (VRF) is computed as:
VRF = Methylated Reads (111) / Single CpG Coverage


- **Goal**: Compare mean VRF values for each PMP across tissues.

---

### 4.6 Sequencing Depth Effects on Specificity Confidence
- **Goal**: Study the impact of varying sequencing depths (100k, 500k, 1M reads) on specificity confidence.
- **Method**:
1. Scale VRF values to simulate different depths.
2. Plot specificity confidence for each depth.

---

### 4.7 Threshold Reads Estimation for Top PMPs
- **Goal**: Estimate the number of reads required to confidently call a tissue type (Tissue #2) at a sequencing depth of 1M reads.
- **Steps**:
1. Select the top 10 PMPs based on feature importance.
2. For each PMP, calculate the read threshold using VRF values.

---

### 4.8 Validation of Hypothesis
- **Goal**: Validate the hypothesis by comparing the specificity of the top 10 PMPs with individual CpG sites.
- **Method**:
1. Filter data for top PMPs.
2. Compare VRF distributions.

---

## 5. Output

### **Tables**:
- Coverage statistics per tissue.
- Significant PMPs with p-values.
- Mean VRF for PMPs across tissues.
- Threshold reads for top PMPs.

### **Plots**:
- Boxplot, violin plot, and bar plot for coverage statistics.
- Confusion matrix for Random Forest predictions.
- Specificity confidence vs sequencing depth.
- Threshold reads for top PMPs.

All plots are saved as PNG files (`conf_matrix.png`, `depth_variations_specificity.png`, etc.).

---

# ### Task 2: NGS Data Analysis  ###

---

## Objective
This task aims to process and analyze raw next-generation sequencing (NGS) data, evaluate the ability to perform quality control, align sequences to a reference genome, and identify somatic mutations.

---

## Dataset
The dataset consists of paired-end FASTQ files, including one sample from normal tissue and one from cancer tissue.

### Normal Tissue:
- **Read 1**: `PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq`
- **Read 2**: `PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq`

### Cancer Tissue:
- **Read 1**: `PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq`
- **Read 2**: `PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq`

---

## Sub-Tasks and Workflow


## Sub-tasks and Workflow

### (1). Quality Control

**Objective:** Perform quality checks on the raw sequencing data to evaluate sequence quality metrics, such as sequence counts, per-base quality, and read duplication levels.

**Tool Used:** FastQC

#### Steps:
1. **Command:** Quality control was performed on all FASTQ files using the following command on a Linux server:
   ```bash
   fastqc ./*.fastq
   ```

2. **Results:** The following output files were generated:
   - PA220KH-lib09-P19-Tumor_S2_L001_R1_001_fastqc.html
   - PA220KH-lib09-P19-Tumor_S2_L001_R2_001_fastqc.html
   - PA221MH-lib09-P19-Norm_S1_L001_R1_001_fastqc.html
   - PA221MH-lib09-P19-Norm_S1_L001_R2_001_fastqc.html

3. **Analysis:** Open the `.html` files in a web browser to view the detailed quality metrics.

### (2). Alignment and Mutation Calling

**Objective:** Align the samples to the human reference genome and identify somatic mutations present in the cancer sample but absent in the normal tissue.

**Tools and Reference Genome:**

- **Alignment Tool:** BWA (Burrows-Wheeler Aligner)
- **Reference Genome:** hg19.fa

#### Sub-tasks and Steps:

##### (a). Alignment:

1. **Command to Align Reads:** The paired-end reads were aligned to the human reference genome using the following command:
   ```bash
   bwa mem $REFERENCE_GENOME $READ1 $READ2 > $SAM_FILE

### Output:
- **Cancer Sample:** PA220KH-lib09-P19-Tumor_S2_L001_aligned.sam
- **Normal Sample:** PA221MH-lib09-P19-Norm_S1_L001_aligned.sam

#### (b).i. Somatic Mutation Identification:

# Step 1: Convert SAM to BAM  
Tool Used: Samtools  
Command: 
```bash
samtools view -Sb $SAM_FILE > $BAM_FILE

# Step 2: Sort BAM Files
Command:
```bash
samtools sort $BAM_FILE -o $SORTED_BAM_FILE

# Step 3: Index BAM Files
Command:
```bash
samtools index $SORTED_BAM_FILE

# Step 4: Call Somatic Mutations  
Tool Used: GATK Mutect2  
Command:
```bash
gatk Mutect2 \
  -R $REFERENCE_GENOME \
  -I $TUMOR_BAM \
  -I $NORMAL_BAM \
  -tumor TUMOR_SAMPLE_ID \
  -normal NORMAL_SAMPLE_ID \
  -O somatic_variants.vcf

# Step 5: Filter Somatic Mutations  
Tool Used: GATK FilterMutectCalls  
Command: 
```bash
gatk FilterMutectCalls -V somatic_variants.vcf -O somatic_variants_filtered.vcf

# Step 6: Extract VCF Data  
Tool Used: bcftools  
Command: 
```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\n' somatic_variants_filtered.vcf > somatic_variants_extracted.tsv

## (b).ii. Custom Code Development:

### Objective:
The goal of this task is to develop custom scripts to perform mutation detection and calculate metrics such as Variant Allele Frequency (VAF). This process refines somatic mutation analysis by filtering out low-quality variants.

### Steps:

1. **Data Loading**:  
   The `somatic_variants_extracted.tsv` file, generated after variant calling, is loaded into a pandas DataFrame.

2. **Renaming Columns**:  
   Columns in the DataFrame are renamed for easier interpretation:
   - **CHROM**: Chromosome
   - **POS**: Position
   - **REF**: Reference Allele
   - **ALT**: Alternate Allele
   - **QUAL**: Quality Score
   - **DP**: Read Depth
   - **AD**: Allele Depth

3. **Calculate Variant Allele Frequency (VAF)**:  
   The function calculates the Variant Allele Frequency (VAF) by dividing the alternate allele depth by the total read depth.

4. **Filter Variants by VAF Threshold**:  
   Variants with a VAF greater than 0.01 are retained for further analysis.

---

## (c). Background Mutation Level Calculation and Confidence Threshold:

### Objective:
This step calculates the median background mutation level using normal tissue data, accounting for sequencing errors or biases. It also determines how many reads per million are required to confidently call a mutation.

### Steps:

1. **Calculate Median Background Mutation Level**:  
   The function computes the background mutation level based on the alternate allele depths (AD) and read depths (DP) for the normal tissue sample. The background mutation level is calculated as the median value of variant allele frequencies (VAF) across the normal tissue data.

2. **Calculate Reads Per Million Required for Confident Mutation Calling**:  
   The function calculates the number of reads per million required to confidently call a mutation based on the median background mutation level. The threshold is inversely proportional to the background mutation level, assuming a simple model with a confidence factor of 99%.
