# FASTQ-to-VCF Pipeline

Germline variant calling pipeline following GATK Best Practices.
Converts paired-end FASTQ files to a high-quality, filtered VCF.

## Quick Start

### 1. Install Dependencies

```bash
conda create -n fq2vcf -c bioconda -c conda-forge python=3.10
conda activate fq2vcf

conda install -c bioconda -c conda-forge \
    bwa-mem2=2.2.1 samtools=1.19 bcftools=1.19 \
    fastp=0.23.4 gatk4=4.5.0.0 picard=3.1.1 \
    multiqc=1.21
```

### 2. Download Reference Genome

```bash
REF_DIR="/data/references"
mkdir -p "$REF_DIR"

# GRCh38 from GATK bundle
wget -P "$REF_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget -P "$REF_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget -P "$REF_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

# Known sites for BQSR
wget -P "$REF_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P "$REF_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Build BWA-MEM2 index (~1 hour)
bwa-mem2 index "$REF_DIR/Homo_sapiens_assembly38.fasta"
```

### 3. Run Pipeline

**Bash script:**
```bash
bash fq2vcf_pipeline.sh \
    --sample NA12878 \
    --fq1 /data/reads/NA12878_R1.fastq.gz \
    --fq2 /data/reads/NA12878_R2.fastq.gz \
    --ref /data/references/Homo_sapiens_assembly38.fasta \
    --known-dbsnp /data/references/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-mills /data/references/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --outdir /data/output/NA12878 \
    --threads 16
```

**Python wrapper (with validation and resumability):**
```bash
python fq2vcf.py \
    --sample NA12878 \
    --fq1 /data/reads/NA12878_R1.fastq.gz \
    --fq2 /data/reads/NA12878_R2.fastq.gz \
    --ref /data/references/Homo_sapiens_assembly38.fasta \
    --known-dbsnp /data/references/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-mills /data/references/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --outdir /data/output/NA12878 \
    --threads 16
```

## Pipeline Steps

| Step | Tool | Purpose |
|------|------|---------|
| 1 | `fastp` | Adapter trimming, quality filtering, dedup |
| 2 | `bwa-mem2` | Alignment to reference genome |
| 3 | `Picard MarkDuplicates` | PCR/optical duplicate flagging |
| 4 | `GATK BQSR` | Base quality score recalibration |
| 5 | `GATK HaplotypeCaller` | Variant calling |
| 6 | `GATK VariantFiltration` | Hard filtering (SNPs + INDELs) |
| 7 | `MultiQC` + `bcftools stats` | QC aggregation |

## Output Files

```
{outdir}/
├── {sample}.final.vcf.gz          # Final filtered VCF
├── {sample}.final.vcf.gz.tbi      # Tabix index
├── 01_fastp/                      # Trimmed FASTQs + QC
├── 02_align/                      # Sorted BAM
├── 03_dedup/                      # Dedup BAM + metrics
├── 04_bqsr/                       # Recalibrated BAM
├── 05_vcf/                        # Raw VCF
├── 06_filter/                     # Filtered VCF
├── 07_qc/                         # QC reports
│   ├── multiqc_report.html        # Aggregated QC
│   └── *.bcftools_stats.txt       # Variant stats
└── {sample}.pipeline.log          # Full log
```

## Filter Thresholds

**SNPs:**
- QD < 2.0 (variant quality by depth)
- FS > 60.0 (strand bias)
- MQ < 40.0 (mapping quality)
- MQRankSum < -12.5 (mapping quality difference)
- ReadPosRankSum < -8.0 (read position bias)
- SOR > 3.0 (strand odds ratio)

**INDELs:**
- QD < 2.0
- FS > 200.0
- ReadPosRankSum < -20.0

## Advanced Options

```bash
# Output gVCF for joint calling (multi-sample)
python fq2vcf.py ... --gvcf

# Skip BQSR (no known sites)
python fq2vcf.py ... --skip-bqsr

# Skip duplicate marking
python fq2vcf.py ... --skip-dedup

# Dry run (validate inputs only)
python fq2vcf.py ... --dry-run
```

## Multi-Sample Joint Calling

After generating per-sample gVCFs:

```bash
# Combine gVCFs
gatk CombineGVCFs \
    -R ref.fasta \
    -V sample1.g.vcf.gz -V sample2.g.vcf.gz \
    -O combined.g.vcf.gz

# Joint genotyping
gatk GenotypeGVCFs \
    -R ref.fasta \
    -V combined.g.vcf.gz \
    -O joint.vcf.gz
```

## System Requirements

- **CPU**: 8+ cores (16+ for WGS)
- **RAM**: 16 GB min, 32 GB recommended (64 GB+ for WGS)
- **Disk**: ~100 GB per WGS sample
- **Java**: JDK 17+ (GATK requirement)
- **OS**: Linux (Ubuntu 20.04+, CentOS 7+)
