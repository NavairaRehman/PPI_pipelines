#!/usr/bin/env python3
"""
fq2vcf.py — Python wrapper for the FASTQ-to-VCF pipeline.

Provides:
  1. Argument parsing with validation
  2. Automatic tool detection
  3. Dry-run mode for testing
  4. Resumable execution (skips completed steps)
  5. Structured logging

Usage:
    python fq2vcf.py \
        --sample NA12878 \
        --fq1 reads_R1.fastq.gz \
        --fq2 reads_R2.fastq.gz \
        --ref /data/references/Homo_sapiens_assembly38.fasta \
        --outdir /data/output/NA12878 \
        --threads 16

    # Dry run (validate inputs without executing):
    python fq2vcf.py ... --dry-run

    # With known sites for BQSR:
    python fq2vcf.py ... \
        --known-dbsnp /data/references/dbsnp_138.hg38.vcf.gz \
        --known-mills /data/references/Mills_and_1000G.indels.hg38.vcf.gz
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ============================================================================
# Configuration
# ============================================================================

@dataclass
class PipelineConfig:
    """Pipeline configuration with validation."""
    sample: str
    fq1: Path
    fq2: Path
    ref: Path
    outdir: Path
    threads: int = 16
    platform: str = "ILLUMINA"
    rg_name: str = ""
    gvcf: bool = False
    skip_bqsr: bool = False
    skip_dedup: bool = False
    known_dbsnp: Optional[Path] = None
    known_mills: Optional[Path] = None
    dry_run: bool = False

    def __post_init__(self):
        self.fq1 = Path(self.fq1)
        self.fq2 = Path(self.fq2)
        self.ref = Path(self.ref)
        self.outdir = Path(self.outdir)
        if self.known_dbsnp:
            self.known_dbsnp = Path(self.known_dbsnp)
        if self.known_mills:
            self.known_mills = Path(self.known_mills)
        if not self.rg_name:
            self.rg_name = self.sample


# ============================================================================
# Logging
# ============================================================================

def setup_logging(outdir: Path, sample: str) -> logging.Logger:
    """Configure logging to file and console."""
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / f"{sample}.pipeline.log"

    logger = logging.getLogger("fq2vcf")
    logger.setLevel(logging.DEBUG)

    # File handler (detailed)
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    ))

    # Console handler (info+)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ============================================================================
# Tool Detection
# ============================================================================

REQUIRED_TOOLS = {
    "fastp":      "QC and adapter trimming",
    "bwa-mem2":   "Read alignment",
    "samtools":   "BAM processing",
    "gatk":       "Variant calling and BQSR",
    "bcftools":   "VCF statistics and manipulation",
}

OPTIONAL_TOOLS = {
    "picard":     "Mark duplicates (fallback: samtools markdup)",
    "multiqc":    "QC report aggregation",
}


def check_tools(logger: logging.Logger, skip_bqsr: bool = False) -> bool:
    """Verify all required tools are available."""
    all_ok = True

    logger.info("\n=== Tool Check ===")
    for tool, purpose in REQUIRED_TOOLS.items():
        path = shutil.which(tool)
        if path:
            # Get version
            version = "unknown"
            for flag in ["--version", "-V", "version", "--help"]:
                try:
                    result = subprocess.run(
                        [tool, flag], capture_output=True, text=True, timeout=10
                    )
                    version = (result.stdout + result.stderr).split("\n")[0][:80]
                    break
                except Exception:
                    pass
            logger.info(f"  [OK] {tool:12s} -> {path} ({version})")
        else:
            logger.error(f"  [MISSING] {tool:12s} — {purpose}")
            all_ok = False

    for tool, purpose in OPTIONAL_TOOLS.items():
        path = shutil.which(tool)
        if path:
            logger.info(f"  [OK] {tool:12s} -> {path}")
        else:
            logger.warning(f"  [SKIP] {tool:12s} — {purpose} (optional)")

    return all_ok


# ============================================================================
# Input Validation
# ============================================================================

def validate_inputs(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Validate all input files and parameters."""
    logger.info("\n=== Input Validation ===")
    valid = True

    # Check FASTQ files
    for label, path in [("R1 FASTQ", cfg.fq1), ("R2 FASTQ", cfg.fq2)]:
        if not path.exists():
            logger.error(f"  {label} not found: {path}")
            valid = False
        elif path.stat().st_size == 0:
            logger.error(f"  {label} is empty: {path}")
            valid = False
        else:
            size_gb = path.stat().st_size / (1024**3)
            logger.info(f"  {label}: {path} ({size_gb:.2f} GB)")

    # Check reference
    if not cfg.ref.exists():
        logger.error(f"  Reference not found: {cfg.ref}")
        valid = False
    else:
        ref_size = cfg.ref.stat().st_size / (1024**3)
        logger.info(f"  Reference: {cfg.ref} ({ref_size:.2f} GB)")

        # Check for index files
        for ext in [".fai", ".amb", ".sa", ".bwt", ".pac", ".ann"]:
            idx = Path(str(cfg.ref) + ext)
            alt_idx = Path(str(cfg.ref).replace(".fasta", "") + ext)
            if not idx.exists() and not alt_idx.exists():
                logger.warning(f"  Missing index: {ext} (will be built)")
                break

    # Check known sites
    if cfg.known_dbsnp:
        if cfg.known_dbsnp.exists():
            logger.info(f"  Known SNPs: {cfg.known_dbsnp}")
        else:
            logger.warning(f"  Known SNPs not found: {cfg.known_dbsnp} (BQSR skipped)")
            cfg.known_dbsnp = None

    if cfg.known_mills:
        if cfg.known_mills.exists():
            logger.info(f"  Known indels: {cfg.known_mills}")
        else:
            logger.warning(f"  Known indels not found: {cfg.known_mills} (BQSR skipped)")
            cfg.known_mills = None

    # Thread validation
    nproc = os.cpu_count() or 1
    if cfg.threads > nproc:
        logger.warning(f"  Requested {cfg.threads} threads but only {nproc} available")
    logger.info(f"  Threads: {cfg.threads}")
    logger.info(f"  Sample: {cfg.sample}")

    return valid


# ============================================================================
# Pipeline Steps
# ============================================================================

def run_step(name: str, cmd: list, logger: logging.Logger, outdir: Path,
             stdout_file: Optional[str] = None, timeout: int = 86400) -> bool:
    """Execute a pipeline step with logging and error handling."""
    logger.info(f"\n--- {name} ---")
    start = time.time()

    log_path = outdir / f"{name.replace(' ', '_').lower()}.log"

    try:
        with open(log_path, "w") as log_fh:
            result = subprocess.run(
                cmd, stdout=log_fh, stderr=subprocess.STDOUT,
                timeout=timeout, text=True
            )

        elapsed = time.time() - start
        if result.returncode == 0:
            logger.info(f"  {name} completed ({elapsed:.0f}s)")
            return True
        else:
            logger.error(f"  {name} FAILED (exit {result.returncode}, {elapsed:.0f}s)")
            logger.error(f"  See log: {log_path}")
            # Show last 10 lines of error
            if log_path.exists():
                lines = log_path.read_text().strip().split("\n")
                for line in lines[-10:]:
                    logger.error(f"    {line}")
            return False

    except subprocess.TimeoutExpired:
        logger.error(f"  {name} TIMED OUT after {timeout}s")
        return False
    except Exception as e:
        logger.error(f"  {name} ERROR: {e}")
        return False


def step_fastp(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Step 1: QC and trimming with fastp."""
    fastp_dir = cfg.outdir / "01_fastp"
    fastp_dir.mkdir(exist_ok=True)

    cmd = [
        "fastp",
        "--in1", str(cfg.fq1),
        "--in2", str(cfg.fq2),
        "--out1", str(fastp_dir / f"{cfg.sample}.R1.trimmed.fastq.gz"),
        "--out2", str(fastp_dir / f"{cfg.sample}.R2.trimmed.fastq.gz"),
        "--json", str(fastp_dir / f"{cfg.sample}.fastp.json"),
        "--html", str(fastp_dir / f"{cfg.sample}.fastp.html"),
        "--thread", str(cfg.threads),
        "--detect_adapter_for_pe",
        "--qualified_quality_phred", "20",
        "--length_required", "50",
        "--correction",
        "--cut_front", "--cut_tail",
        "--cut_window_size", "4",
        "--cut_mean_quality", "20",
        "--overrepresentation_analysis",
    ]
    return run_step("fastp", cmd, logger, fastp_dir)


def step_align(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Step 2: Align with BWA-MEM2."""
    align_dir = cfg.outdir / "02_align"
    align_dir.mkdir(exist_ok=True)

    bam_out = align_dir / f"{cfg.sample}.sorted.bam"

    # Read group string
    rg = f"@RG\\tID:{cfg.sample}\\tSM:{cfg.sample}\\tPL:{cfg.platform}\\tLB:{cfg.sample}_lib\\tPU:{cfg.sample}_unit"

    # Check and build index if needed
    ref_str = str(cfg.ref)
    index_exts = [".01.2bit.amb", ".amb", ".01.2bit.bwt", ".bwt"]
    has_index = any(Path(ref_str + ext).exists() for ext in index_exts)
    if not has_index:
        logger.info("  Building BWA-MEM2 index...")
        idx_cmd = ["bwa-mem2", "index", ref_str]
        if not run_step("bwa-mem2 index", idx_cmd, logger, align_dir):
            return False

    # BWA-MEM2 align -> samtools sort (piped)
    fastp_dir = cfg.outdir / "01_fastp"
    threads_half = max(1, cfg.threads // 2)
    fq1_trimmed = fastp_dir / f"{cfg.sample}.R1.trimmed.fastq.gz"
    fq2_trimmed = fastp_dir / f"{cfg.sample}.R2.trimmed.fastq.gz"
    align_cmd = [
        "bash", "-c",
        f"bwa-mem2 mem -t {cfg.threads} -R '{rg}' '{ref_str}' "
        f"'{fq1_trimmed}' '{fq2_trimmed}' "
        f"| samtools view -@ {threads_half} -bS - "
        f"| samtools sort -@ {threads_half} -m 4G -o '{bam_out}' -"
    ]
    if not run_step("bwa-mem2 + samtools sort", align_cmd, logger, align_dir):
        return False

    # Index BAM
    idx_cmd = ["samtools", "index", "-@", str(threads_half), str(bam_out)]
    return run_step("samtools index", idx_cmd, logger, align_dir)


def step_markdup(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Step 3: Mark duplicates."""
    dedup_dir = cfg.outdir / "03_dedup"
    dedup_dir.mkdir(exist_ok=True)

    bam_in = cfg.outdir / "02_align" / f"{cfg.sample}.sorted.bam"
    bam_out = dedup_dir / f"{cfg.sample}.dedup.bam"
    metrics = dedup_dir / f"{cfg.sample}.dedup.metrics"

    cmd = [
        "gatk", "MarkDuplicates",
        "-I", str(bam_in),
        "-O", str(bam_out),
        "-M", str(metrics),
        "--CREATE_INDEX", "true",
        "--VALIDATION_STRINGENCY", "SILENT",
        "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", "2500",
    ]
    return run_step("MarkDuplicates", cmd, logger, dedup_dir)


def step_bqsr(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Step 4: Base Quality Score Recalibration."""
    bqsr_dir = cfg.outdir / "04_bqsr"
    bqsr_dir.mkdir(exist_ok=True)

    bam_in = cfg.outdir / "03_dedup" / f"{cfg.sample}.dedup.bam"
    if cfg.skip_dedup:
        bam_in = cfg.outdir / "02_align" / f"{cfg.sample}.sorted.bam"

    recal_table = bqsr_dir / f"{cfg.sample}.recal.table"
    bam_out = bqsr_dir / f"{cfg.sample}.recal.bam"

    # Build known-sites args
    known_args = []
    if cfg.known_dbsnp:
        known_args += ["--known-sites", str(cfg.known_dbsnp)]
    if cfg.known_mills:
        known_args += ["--known-sites", str(cfg.known_mills)]

    # BaseRecalibrator
    cmd1 = [
        "gatk", "BaseRecalibrator",
        "-R", str(cfg.ref),
        "-I", str(bam_in),
        *known_args,
        "-O", str(recal_table),
    ]
    if not run_step("BaseRecalibrator", cmd1, logger, bqsr_dir):
        return False

    # ApplyBQSR
    cmd2 = [
        "gatk", "ApplyBQSR",
        "-R", str(cfg.ref),
        "-I", str(bam_in),
        "-bqsr", str(recal_table),
        "-O", str(bam_out),
        "--create-output-bam-index", "true",
    ]
    return run_step("ApplyBQSR", cmd2, logger, bqsr_dir)


def step_call_variants(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Step 5: Call variants with GATK HaplotypeCaller."""
    vcf_dir = cfg.outdir / "05_vcf"
    vcf_dir.mkdir(exist_ok=True)

    # Determine input BAM
    if not cfg.skip_bqsr and cfg.known_dbsnp:
        bam_in = cfg.outdir / "04_bqsr" / f"{cfg.sample}.recal.bam"
    elif not cfg.skip_dedup:
        bam_in = cfg.outdir / "03_dedup" / f"{cfg.sample}.dedup.bam"
    else:
        bam_in = cfg.outdir / "02_align" / f"{cfg.sample}.sorted.bam"

    if cfg.gvcf:
        vcf_out = vcf_dir / f"{cfg.sample}.g.vcf.gz"
        extra_args = ["-ERC", "GVCF"]
    else:
        vcf_out = vcf_dir / f"{cfg.sample}.raw.vcf.gz"
        extra_args = []

    cmd = [
        "gatk", "HaplotypeCaller",
        "-R", str(cfg.ref),
        "-I", str(bam_in),
        "-O", str(vcf_out),
        *extra_args,
        "-G", "StandardAnnotation",
        "-G", "AS_StandardAnnotation",
        "--native-pair-hmm-threads", str(cfg.threads),
        "--min-base-quality-score", "20",
        "--standard-min-confidence-threshold-for-calling", "30",
    ]
    return run_step("HaplotypeCaller", cmd, logger, vcf_dir)


def step_filter_variants(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Step 6: Filter variants (hard filters for high-quality calls)."""
    filt_dir = cfg.outdir / "06_filter"
    filt_dir.mkdir(exist_ok=True)

    vcf_dir = cfg.outdir / "05_vcf"
    vcf_raw = vcf_dir / f"{cfg.sample}.raw.vcf.gz"
    vcf_snp = vcf_dir / f"{cfg.sample}.raw.snp.vcf.gz"
    vcf_indel = vcf_dir / f"{cfg.sample}.raw.indel.vcf.gz"
    vcf_snp_filt = filt_dir / f"{cfg.sample}.filtered.snp.vcf.gz"
    vcf_indel_filt = filt_dir / f"{cfg.sample}.filtered.indel.vcf.gz"
    vcf_final = cfg.outdir / f"{cfg.sample}.final.vcf.gz"

    steps = [
        ("SelectVariants SNPs", [
            "gatk", "SelectVariants", "-R", str(cfg.ref),
            "-V", str(vcf_raw), "--select-type-to-include", "SNP",
            "-O", str(vcf_snp),
        ]),
        ("SelectVariants INDELs", [
            "gatk", "SelectVariants", "-R", str(cfg.ref),
            "-V", str(vcf_raw), "--select-type-to-include", "INDEL",
            "-O", str(vcf_indel),
        ]),
        ("Filter SNPs", [
            "gatk", "VariantFiltration", "-R", str(cfg.ref),
            "-V", str(vcf_snp), "-O", str(vcf_snp_filt),
            "--filter-expression", "QD < 2.0", "--filter-name", "LowQD",
            "--filter-expression", "FS > 60.0", "--filter-name", "HighFS",
            "--filter-expression", "MQ < 40.0", "--filter-name", "LowMQ",
            "--filter-expression", "MQRankSum < -12.5", "--filter-name", "LowMQRankSum",
            "--filter-expression", "ReadPosRankSum < -8.0", "--filter-name", "LowReadPosRankSum",
            "--filter-expression", "SOR > 3.0", "--filter-name", "HighSOR",
        ]),
        ("Filter INDELs", [
            "gatk", "VariantFiltration", "-R", str(cfg.ref),
            "-V", str(vcf_indel), "-O", str(vcf_indel_filt),
            "--filter-expression", "QD < 2.0", "--filter-name", "LowQD",
            "--filter-expression", "FS > 200.0", "--filter-name", "HighFS",
            "--filter-expression", "ReadPosRankSum < -20.0", "--filter-name", "LowReadPosRankSum",
        ]),
        ("MergeVcfs", [
            "gatk", "MergeVcfs",
            "-I", str(vcf_snp_filt), "-I", str(vcf_indel_filt),
            "-O", str(vcf_final),
        ]),
    ]

    for name, cmd in steps:
        if not run_step(name, cmd, logger, filt_dir):
            return False

    # Index final VCF
    return run_step("bcftools index", ["bcftools", "index", str(vcf_final)], logger, filt_dir)


def step_qc(cfg: PipelineConfig, logger: logging.Logger) -> bool:
    """Step 7: Generate QC metrics."""
    qc_dir = cfg.outdir / "07_qc"
    qc_dir.mkdir(exist_ok=True)

    # Determine input BAM
    if not cfg.skip_bqsr and cfg.known_dbsnp:
        bam_in = cfg.outdir / "04_bqsr" / f"{cfg.sample}.recal.bam"
    elif not cfg.skip_dedup:
        bam_in = cfg.outdir / "03_dedup" / f"{cfg.sample}.dedup.bam"
    else:
        bam_in = cfg.outdir / "02_align" / f"{cfg.sample}.sorted.bam"

    ok = True
    ok &= run_step("flagstat", [
        "samtools", "flagstat", "-@", str(cfg.threads), str(bam_in)
    ], logger, qc_dir)
    ok &= run_step("idxstats", [
        "samtools", "idxstats", "-@", str(cfg.threads), str(bam_in)
    ], logger, qc_dir)
    ok &= run_step("samtools stats", [
        "samtools", "stats", "-@", str(cfg.threads), str(bam_in)
    ], logger, qc_dir)

    # bcftools stats on final VCF
    if not cfg.gvcf:
        vcf_final = cfg.outdir / f"{cfg.sample}.final.vcf.gz"
        ok &= run_step("bcftools stats", [
            "bcftools", "stats", str(vcf_final)
        ], logger, qc_dir)

    # MultiQC
    run_step("multiqc", [
        "multiqc", str(cfg.outdir), "-o", str(qc_dir), "--force"
    ], logger, qc_dir)

    return ok


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="FASTQ-to-VCF germline variant calling pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python fq2vcf.py --sample NA12878 \\
      --fq1 reads_R1.fastq.gz --fq2 reads_R2.fastq.gz \\
      --ref hg38.fasta --outdir results/NA12878

  # With BQSR known sites
  python fq2vcf.py --sample NA12878 \\
      --fq1 reads_R1.fastq.gz --fq2 reads_R2.fastq.gz \\
      --ref hg38.fasta --outdir results/NA12878 \\
      --known-dbsnp dbsnp.hg38.vcf.gz \\
      --known-mills Mills_1000G.indels.hg38.vcf.gz

  # Output gVCF for joint calling
  python fq2vcf.py --sample NA12878 ... --gvcf

  # Dry run (validate only)
  python fq2vcf.py --sample NA12878 ... --dry-run
        """,
    )

    # Required arguments
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--fq1", required=True, help="Read 1 FASTQ file")
    parser.add_argument("--fq2", required=True, help="Read 2 FASTQ file")
    parser.add_argument("--ref", required=True, help="Reference genome FASTA")
    parser.add_argument("--outdir", required=True, help="Output directory")

    # Optional arguments
    parser.add_argument("--threads", type=int, default=16, help="Number of threads (default: 16)")
    parser.add_argument("--platform", default="ILLUMINA", help="Sequencing platform (default: ILLUMINA)")
    parser.add_argument("--read-group-name", default="", help="Read group name (default: sample ID)")
    parser.add_argument("--gvcf", action="store_true", help="Output gVCF for joint calling")
    parser.add_argument("--skip-bqsr", action="store_true", help="Skip BQSR step")
    parser.add_argument("--skip-dedup", action="store_true", help="Skip duplicate marking")
    parser.add_argument("--known-dbsnp", default=None, help="dbSNP VCF for BQSR")
    parser.add_argument("--known-mills", default=None, help="Mills & 1000G indels VCF for BQSR")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without executing")

    args = parser.parse_args()

    # Build config
    cfg = PipelineConfig(
        sample=args.sample,
        fq1=args.fq1,
        fq2=args.fq2,
        ref=args.ref,
        outdir=args.outdir,
        threads=args.threads,
        platform=args.platform,
        rg_name=args.read_group_name,
        gvcf=args.gvcf,
        skip_bqsr=args.skip_bqsr,
        skip_dedup=args.skip_dedup,
        known_dbsnp=args.known_dbsnp,
        known_mills=args.known_mills,
        dry_run=args.dry_run,
    )

    # Setup logging
    logger = setup_logging(cfg.outdir, cfg.sample)
    logger.info("=" * 60)
    logger.info("  FASTQ-to-VCF Pipeline (GATK Best Practices)")
    logger.info("=" * 60)
    logger.info(f"  Sample:    {cfg.sample}")
    logger.info(f"  Threads:   {cfg.threads}")
    logger.info(f"  Output:    {cfg.outdir}")
    logger.info(f"  gVCF:      {cfg.gvcf}")
    logger.info(f"  BQSR:      {not cfg.skip_bqsr}")
    logger.info(f"  Dedup:     {not cfg.skip_dedup}")

    # Check tools
    if not check_tools(logger):
        logger.error("\nMissing required tools. Install with:")
        logger.error("  conda install -c bioconda bwa-mem2 samtools gatk4 fastp bcftools picard multiqc")
        sys.exit(1)

    # Validate inputs
    if not validate_inputs(cfg, logger):
        logger.error("\nInput validation failed. Fix the above errors and retry.")
        sys.exit(1)

    if cfg.dry_run:
        logger.info("\n=== DRY RUN — all checks passed ===")
        sys.exit(0)

    # Execute pipeline
    pipeline_start = time.time()
    logger.info("\n" + "=" * 60)
    logger.info("  Starting pipeline execution")
    logger.info("=" * 60)

    steps = [
        ("Step 1: QC & Trim (fastp)",      step_fastp),
        ("Step 2: Align (BWA-MEM2)",        step_align),
        ("Step 3: Mark Duplicates (Picard)", step_markdup),
        ("Step 4: BQSR (GATK)",             step_bqsr),
        ("Step 5: Call Variants (GATK HC)",  step_call_variants),
        ("Step 6: Filter Variants",          step_filter_variants),
        ("Step 7: QC Metrics",               step_qc),
    ]

    # Skip steps based on config
    if cfg.skip_dedup:
        steps = [s for s in steps if "Dedup" not in s[0] or "Mark" not in s[0]]
    if cfg.skip_bqsr or not (cfg.known_dbsnp and cfg.known_mills):
        steps = [s for s in steps if "BQSR" not in s[0]]
    if cfg.gvcf:
        steps = [s for s in steps if "Filter" not in s[0]]

    for step_name, step_func in steps:
        # Skip dedup/bqsr based on flags
        if "Dedup" in step_name and cfg.skip_dedup:
            logger.info(f"\n--- {step_name} (SKIPPED) ---")
            continue
        if "BQSR" in step_name and (cfg.skip_bqsr or not cfg.known_dbsnp):
            logger.info(f"\n--- {step_name} (SKIPPED) ---")
            continue
        if "Filter" in step_name and cfg.gvcf:
            logger.info(f"\n--- {step_name} (SKIPPED for gVCF) ---")
            continue

        success = step_func(cfg, logger)
        if not success:
            logger.error(f"\nPipeline FAILED at: {step_name}")
            sys.exit(1)

    # Summary
    elapsed = time.time() - pipeline_start
    logger.info("\n" + "=" * 60)
    logger.info("  PIPELINE COMPLETE")
    logger.info("=" * 60)
    logger.info(f"  Total time: {elapsed/60:.1f} minutes")

    if cfg.gvcf:
        logger.info(f"  gVCF: {cfg.outdir}/05_vcf/{cfg.sample}.g.vcf.gz")
    else:
        logger.info(f"  VCF:  {cfg.outdir}/{cfg.sample}.final.vcf.gz")
    logger.info(f"  QC:   {cfg.outdir}/07_qc/")
    logger.info(f"  Log:  {cfg.outdir}/{cfg.sample}.pipeline.log")

    # Print flagstat summary
    flagstat = cfg.outdir / "07_qc" / "flagstat.log"
    if flagstat.exists():
        logger.info("\n=== Alignment Summary ===")
        for line in flagstat.read_text().split("\n")[:5]:
            if line.strip():
                logger.info(f"  {line}")


if __name__ == "__main__":
    main()
