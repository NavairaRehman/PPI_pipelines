#!/usr/bin/env bash
# ============================================================================
# fq2vcf_pipeline.sh — FASTQ (paired-end) to high-quality VCF
# ============================================================================
# Germline variant calling pipeline following GATK Best Practices.
#
# Usage:
#   bash fq2vcf_pipeline.sh \
#       --sample SAMPLE_ID \
#       --fq1  /path/to/read1.fastq.gz \
#       --fq2  /path/to/read2.fastq.gz \
#       --ref  /path/to/reference.fasta \
#       --known-dbsnp   /path/to/dbsnp.vcf.gz \
#       --known-mills   /path/to/Mills_and_1000G.indels.vcf.gz \
#       --outdir /path/to/output \
#       [--threads 16] \
#       [--platform ILLUMINA] \
#       [--read-group-name SAMPLE_ID] \
#       [--gvcf]               # output gVCF for joint calling
#       [--skip-bqsr]          # skip BQSR (no known-sites files)
#       [--skip-dedup]         # skip duplicate marking
#       [--help]
#
# Outputs:
#   {outdir}/{sample}.final.vcf.gz      — filtered VCF
#   {outdir}/{sample}.final.vcf.gz.tbi  — tabix index
#   {outdir}/multiqc_report.html        — aggregated QC
#   {outdir}/pipeline.log               — full log
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================
# DEFAULTS
# ============================
THREADS=$(nproc)
PLATFORM="ILLUMINA"
RG_NAME=""
GVCF_MODE=false
SKIP_BQSR=false
SKIP_DEDUP=false
OUTDIR=""
SAMPLE=""
FQ1=""
FQ2=""
REF=""
KNOWN_DBSNP=""
KNOWN_MILLS=""

# ============================
# COLOR OUTPUT
# ============================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log()  { echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $*" | tee -a "${OUTDIR}/pipeline.log" 2>/dev/null || echo "$*"; }
ok()   { echo -e "${GREEN}[OK]${NC} $*" | tee -a "${OUTDIR}/pipeline.log" 2>/dev/null || echo "$*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*" | tee -a "${OUTDIR}/pipeline.log" 2>/dev/null || echo "$*"; }
err()  { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

# ============================
# PARSE ARGUMENTS
# ============================
usage() {
    head -27 "$0" | tail -23
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)     SAMPLE="$2"; shift 2 ;;
        --fq1)        FQ1="$2"; shift 2 ;;
        --fq2)        FQ2="$2"; shift 2 ;;
        --ref)        REF="$2"; shift 2 ;;
        --known-dbsnp) KNOWN_DBSNP="$2"; shift 2 ;;
        --known-mills) KNOWN_MILLS="$2"; shift 2 ;;
        --outdir)     OUTDIR="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --platform)   PLATFORM="$2"; shift 2 ;;
        --read-group-name) RG_NAME="$2"; shift 2 ;;
        --gvcf)       GVCF_MODE=true; shift ;;
        --skip-bqsr)  SKIP_BQSR=true; shift ;;
        --skip-dedup) SKIP_DEDUP=true; shift ;;
        --help|-h)    usage ;;
        *)            err "Unknown option: $1" ;;
    esac
done

# ============================
# VALIDATE INPUTS
# ============================
[[ -z "$SAMPLE" ]]  && err "Missing --sample"
[[ -z "$FQ1" ]]     && err "Missing --fq1"
[[ -z "$FQ2" ]]     && err "Missing --fq2"
[[ -z "$REF" ]]     && err "Missing --ref"
[[ -z "$OUTDIR" ]]  && err "Missing --outdir"

[[ ! -f "$FQ1" ]]  && err "FASTQ R1 not found: $FQ1"
[[ ! -f "$FQ2" ]]  && err "FASTQ R2 not found: $FQ2"
[[ ! -f "$REF" ]]  && err "Reference not found: $REF"

[[ -z "$RG_NAME" ]] && RG_NAME="$SAMPLE"

# Check required tools
for cmd in fastp bwa-mem2 samtools gatk java; do
    command -v "$cmd" &>/dev/null || err "Required tool not found: $cmd"
done

if [[ "$SKIP_BQSR" == false ]] && [[ -z "$KNOWN_DBSNP" || -z "$KNOWN_MILLS" ]]; then
    warn "Known-sites files not fully specified. BQSR will be skipped."
    SKIP_BQSR=true
fi

# ============================
# SETUP
# ============================
mkdir -p "$OUTDIR"/{fastp,align,dedup,bqsr,vcf,qc}

# Initialize log
LOG="${OUTDIR}/pipeline.log"
echo "=== fq2vcf Pipeline ===" > "$LOG"
log "Sample:    $SAMPLE"
log "FASTQ R1:  $FQ1"
log "FASTQ R2:  $FQ2"
log "Reference: $REF"
log "Threads:   $THREADS"
log "Output:    $OUTDIR"
log "gVCF mode: $GVCF_MODE"
log "Skip BQSR: $SKIP_BQSR"
log ""

# Record start time
PIPELINE_START=$(date +%s)

# ============================
# STEP 1: QC & TRIM with fastp
# ============================
log ">>> STEP 1: QC & Trim with fastp"
STEP_START=$(date +%s)

FASTP_R1="${OUTDIR}/fastp/${SAMPLE}.R1.trimmed.fastq.gz"
FASTP_R2="${OUTDIR}/fastp/${SAMPLE}.R2.trimmed.fastq.gz"

fastp \
    --in1 "$FQ1" \
    --in2 "$FQ2" \
    --out1 "$FASTP_R1" \
    --out2 "$FASTP_R2" \
    --json "${OUTDIR}/fastp/${SAMPLE}.fastp.json" \
    --html "${OUTDIR}/fastp/${SAMPLE}.fastp.html" \
    --thread "$THREADS" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --correction \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --overrepresentation_analysis \
    2>&1 | tee -a "$LOG"

STEP_END=$(date +%s)
ok "Step 1 complete ($(( STEP_END - STEP_START ))s)"

# ============================
# STEP 2: ALIGN with BWA-MEM2
# ============================
log ">>> STEP 2: Align with BWA-MEM2"
STEP_START=$(date +%s)

# Check if BWA-MEM2 index exists (look for .01.2bit.amb or .amb)
if [[ ! -f "${REF}.01.2bit.amb" ]] && [[ ! -f "${REF}.amb" ]]; then
    log "Building BWA-MEM2 index..."
    bwa-mem2 index "$REF" 2>&1 | tee -a "$LOG"
fi

# Build RG string
RG="@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:${PLATFORM}\\tLB:${SAMPLE}_lib\\tPU:${SAMPLE}_unit"

BAM_RAW="${OUTDIR}/align/${SAMPLE}.raw.bam"

bwa-mem2 mem \
    -t "$THREADS" \
    -R "$RG" \
    "$REF" \
    "$FASTP_R1" \
    "$FASTP_R2" \
    2>> "$LOG" \
    | samtools view -@ "$((THREADS / 2))" -bS - \
    | samtools sort -@ "$((THREADS / 2))" -m 4G -o "$BAM_RAW" -

samtools index -@ "$((THREADS / 2))" "$BAM_RAW"

STEP_END=$(date +%s)
ok "Step 2 complete ($(( STEP_END - STEP_START ))s)"

# ============================
# STEP 3: MARK DUPLICATES
# ============================
if [[ "$SKIP_DEDUP" == false ]]; then
    log ">>> STEP 3: Mark Duplicates with Picard"
    STEP_START=$(date +%s)

    BAM_DEDUP="${OUTDIR}/dedup/${SAMPLE}.dedup.bam"
    METRICS="${OUTDIR}/dedup/${SAMPLE}.dedup.metrics"

    gatk MarkDuplicates \
        -I "$BAM_RAW" \
        -O "$BAM_DEDUP" \
        -M "$METRICS" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        2>&1 | tee -a "$LOG"

    STEP_END=$(date +%s)
    ok "Step 3 complete ($(( STEP_END - STEP_START ))s)"
else
    log ">>> STEP 3: Skipped (duplicate marking disabled)"
    BAM_DEDUP="$BAM_RAW"
fi

# ============================
# STEP 4: BQSR (Base Quality Score Recalibration)
# ============================
if [[ "$SKIP_BQSR" == false ]]; then
    log ">>> STEP 4: BQSR"
    STEP_START=$(date +%s)

    RECAL_TABLE="${OUTDIR}/bqsr/${SAMPLE}.recal.table"
    BAM_BQSR="${OUTDIR}/bqsr/${SAMPLE}.recal.bam"

    # Build known-sites arguments
    KNOWN_ARGS=""
    [[ -n "$KNOWN_DBSNP" ]] && KNOWN_ARGS="$KNOWN_ARGS --known-sites $KNOWN_DBSNP"
    [[ -n "$KNOWN_MILLS" ]] && KNOWN_ARGS="$KNOWN_ARGS --known-sites $KNOWN_MILLS"

    # Step 4a: Build recalibration model
    gatk BaseRecalibrator \
        -R "$REF" \
        -I "$BAM_DEDUP" \
        $KNOWN_ARGS \
        -O "$RECAL_TABLE" \
        2>&1 | tee -a "$LOG"

    # Step 4b: Apply recalibration
    gatk ApplyBQSR \
        -R "$REF" \
        -I "$BAM_DEDUP" \
        -bqsr "$RECAL_TABLE" \
        -O "$BAM_BQSR" \
        --create-output-bam-index true \
        2>&1 | tee -a "$LOG"

    STEP_END=$(date +%s)
    ok "Step 4 complete ($(( STEP_END - STEP_START ))s)"
else
    log ">>> STEP 4: Skipped (BQSR disabled)"
    BAM_BQSR="$BAM_DEDUP"
fi

# ============================
# STEP 5: VARIANT CALLING
# ============================
log ">>> STEP 5: Variant Calling with GATK HaplotypeCaller"
STEP_START=$(date +%s)

if [[ "$GVCF_MODE" == true ]]; then
    VCF_OUT="${OUTDIR}/vcf/${SAMPLE}.g.vcf.gz"
    gatk HaplotypeCaller \
        -R "$REF" \
        -I "$BAM_BQSR" \
        -O "$VCF_OUT" \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --native-pair-hmm-threads "$THREADS" \
        --min-base-quality-score 20 \
        --standard-min-confidence-threshold-for-calling 30 \
        2>&1 | tee -a "$LOG"
else
    VCF_RAW="${OUTDIR}/vcf/${SAMPLE}.raw.vcf.gz"
    gatk HaplotypeCaller \
        -R "$REF" \
        -I "$BAM_BQSR" \
        -O "$VCF_RAW" \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --native-pair-hmm-threads "$THREADS" \
        --min-base-quality-score 20 \
        --standard-min-confidence-threshold-for-calling 30 \
        2>&1 | tee -a "$LOG"
fi

STEP_END=$(date +%s)
ok "Step 5 complete ($(( STEP_END - STEP_START ))s)"

# ============================
# STEP 6: HARD FILTERING (for non-gVCF output)
# ============================
if [[ "$GVCF_MODE" == false ]]; then
    log ">>> STEP 6: Variant Filtering"
    STEP_START=$(date +%s)

    # Separate SNPs and INDELs for different filter strategies
    VCF_SNP="${OUTDIR}/vcf/${SAMPLE}.raw.snp.vcf.gz"
    VCF_INDEL="${OUTDIR}/vcf/${SAMPLE}.raw.indel.vcf.gz"
    VCF_SNP_FILT="${OUTDIR}/vcf/${SAMPLE}.filtered.snp.vcf.gz"
    VCF_INDEL_FILT="${OUTDIR}/vcf/${SAMPLE}.filtered.indel.vcf.gz"
    VCF_FINAL="${OUTDIR}/${SAMPLE}.final.vcf.gz"

    # 6a: Select SNPs
    gatk SelectVariants \
        -R "$REF" \
        -V "$VCF_RAW" \
        --select-type-to-include SNP \
        -O "$VCF_SNP" \
        2>&1 | tee -a "$LOG"

    # 6b: Select INDELs
    gatk SelectVariants \
        -R "$REF" \
        -V "$VCF_RAW" \
        --select-type-to-include INDEL \
        -O "$VCF_INDEL" \
        2>&1 | tee -a "$LOG"

    # 6c: Apply SNP filters (GATK Best Practices hard filters)
    gatk VariantFiltration \
        -R "$REF" \
        -V "$VCF_SNP" \
        -O "$VCF_SNP_FILT" \
        --filter-expression "QD < 2.0"                              --filter-name "LowQD" \
        --filter-expression "FS > 60.0"                             --filter-name "HighFS" \
        --filter-expression "MQ < 40.0"                             --filter-name "LowMQ" \
        --filter-expression "MQRankSum < -12.5"                     --filter-name "LowMQRankSum" \
        --filter-expression "ReadPosRankSum < -8.0"                 --filter-name "LowReadPosRankSum" \
        --filter-expression "SOR > 3.0"                             --filter-name "HighSOR" \
        2>&1 | tee -a "$LOG"

    # 6d: Apply INDEL filters
    gatk VariantFiltration \
        -R "$REF" \
        -V "$VCF_INDEL" \
        -O "$VCF_INDEL_FILT" \
        --filter-expression "QD < 2.0"                              --filter-name "LowQD" \
        --filter-expression "FS > 200.0"                            --filter-name "HighFS" \
        --filter-expression "ReadPosRankSum < -20.0"                --filter-name "LowReadPosRankSum" \
        2>&1 | tee -a "$LOG"

    # 6e: Merge filtered SNPs and INDELs
    gatk MergeVcfs \
        -I "$VCF_SNP_FILT" \
        -I "$VCF_INDEL_FILT" \
        -O "$VCF_FINAL" \
        2>&1 | tee -a "$LOG"

    # Index final VCF
    bcftools index "$VCF_FINAL"

    STEP_END=$(date +%s)
    ok "Step 6 complete ($(( STEP_END - STEP_START ))s)"
fi

# ============================
# STEP 7: QC STATS
# ============================
log ">>> STEP 7: Generate QC Metrics"
STEP_START=$(date +%s)

# samtools flagstat
samtools flagstat -@ "$THREADS" "$BAM_BQSR" > "${OUTDIR}/qc/${SAMPLE}.flagstat.txt"

# samtools idxstats
samtools idxstats -@ "$THREADS" "$BAM_BQSR" > "${OUTDIR}/qc/${SAMPLE}.idxstats.txt"

# samtools stats
samtools stats -@ "$THREADS" "$BAM_BQSR" > "${OUTDIR}/qc/${SAMPLE}.samtools_stats.txt"

# CollectInsertSizeMetrics (Picard)
gatk CollectInsertSizeMetrics \
    -I "$BAM_BQSR" \
    -O "${OUTDIR}/qc/${SAMPLE}.insert_size_metrics.txt" \
    -H "${OUTDIR}/qc/${SAMPLE}.insert_size_histogram.pdf" \
    2>&1 | tee -a "$LOG" || warn "InsertSizeMetrics failed (non-critical)"

# bcftools stats on final VCF
if [[ "$GVCF_MODE" == false ]]; then
    bcftools stats "$VCF_FINAL" > "${OUTDIR}/qc/${SAMPLE}.bcftools_stats.txt"

    # Quick variant count summary
    log "--- Variant Summary ---"
    bcftools query -f '%TYPE\n' "$VCF_FINAL" | sort | uniq -c | tee -a "$LOG"
    PASS_COUNT=$(bcftools view -f PASS "$VCF_FINAL" | grep -vc '^#' || echo 0)
    TOTAL_COUNT=$(bcftools view "$VCF_FINAL" | grep -vc '^#' || echo 0)
    log "PASS variants: $PASS_COUNT / $TOTAL_COUNT"
fi

# Aggregate with MultiQC
multiqc "$OUTDIR" -o "${OUTDIR}/qc" --force 2>&1 | tee -a "$LOG" || warn "MultiQC failed (non-critical)"

STEP_END=$(date +%s)
ok "Step 7 complete ($(( STEP_END - STEP_START ))s)"

# ============================
# CLEANUP INTERMEDIATE FILES (optional)
# ============================
log ""
log ">>> Pipeline complete!"
PIPELINE_END=$(date +%s)
TOTAL_TIME=$(( PIPELINE_END - PIPELINE_START ))
log "Total time: $(( TOTAL_TIME / 60 ))m $(( TOTAL_TIME % 60 ))s"
log ""
log "=== Output Files ==="

if [[ "$GVCF_MODE" == true ]]; then
    log "  gVCF:       ${OUTDIR}/vcf/${SAMPLE}.g.vcf.gz"
else
    log "  VCF:        ${OUTDIR}/${SAMPLE}.final.vcf.gz"
    log "  VCF index:  ${OUTDIR}/${SAMPLE}.final.vcf.gz.tbi"
fi
log "  QC report:  ${OUTDIR}/qc/multiqc_report.html"
log "  BAM:        ${BAM_BQSR}"
log "  Log:        ${OUTDIR}/pipeline.log"
log ""

# Print flagstat summary
log "=== Alignment Summary ==="
head -5 "${OUTDIR}/qc/${SAMPLE}.flagstat.txt" | tee -a "$LOG"

log ""
log "=== Done ==="

# Exit with specific code
exit 0
