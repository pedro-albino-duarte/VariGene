import os
import subprocess
import argparse
import gzip

def validate_inputs(input_files):
    """
    Validate user input files and determine if they are single-end or paired-end.
    Returns the paths for R1 and R2 files (R2 is None for single-end).
    """
    if len(input_files) == 1:
        # Single-end
        if not input_files[0].endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            raise ValueError("Invalid input format. Provide a valid FASTQ file for single-end input.")
        return input_files[0], None
    elif len(input_files) == 2:
        # Paired-end
        r1, r2 = input_files
        if not (r1.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")) and r2.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))):
            raise ValueError("Invalid input format. Provide valid FASTQ files for paired-end input.")
        return r1, r2
    else:
        raise ValueError("Invalid number of input files. Provide 1 file for single-end or 2 files for paired-end.")

def detect_compression(file_path):
    """
    Check if a file is compressed (GZip format).
    """
    with open(file_path, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'


def index_reference_genome(reference_genome):
    """
    Index the reference genome using bwa, if not already indexed.
    """
    required_files = [f"{reference_genome}.{ext}" for ext in ["bwt", "pac", "ann", "amb", "sa"]]
    if all(os.path.exists(f) for f in required_files):
        print("Reference genome is already indexed. Skipping indexing step.")
        return
    cmd = ["bwa", "index", reference_genome]
    print(f"Indexing reference genome with command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_fastp(input_fastq1, input_fastq2, output_dir):
    """
    Run fastp for quality control of FASTQ files.
    Returns paths to cleaned FASTQ files.
    """
    cleaned_fastq1 = os.path.join(output_dir, "cleaned_1.fastq")
    cleaned_fastq2 = os.path.join(output_dir, "cleaned_2.fastq") if input_fastq2 else None
    fastp_json = os.path.join(output_dir, "fastp_report.json")
    fastp_html = os.path.join(output_dir, "fastp_report.html")

    if input_fastq2:
        # Paired-end
        cmd = [
            "fastp",
            "-i", input_fastq1,
            "-I", input_fastq2,
            "-o", cleaned_fastq1,
            "-O", cleaned_fastq2,
            "--json", fastp_json,
            "--html", fastp_html
        ]
    else:
        # Single-end
        cmd = [
            "fastp",
            "-i", input_fastq1,
            "-o", cleaned_fastq1,
            "--json", fastp_json,
            "--html", fastp_html
        ]

    print(f"Running fastp: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return cleaned_fastq1, cleaned_fastq2

def run_bwa(reference_genome, input_fastq1, input_fastq2, output_dir):
    """
    Run bwa to align FASTQ files against the reference genome.
    Returns the path to the SAM file.
    """
    aligned_sam = os.path.join(output_dir, "aligned.sam")

    if input_fastq2:
        # Paired-end
        cmd = [
            "bwa", "mem", reference_genome, input_fastq1, input_fastq2
        ]
    else:
        # Single-end
        cmd = [
            "bwa", "mem", reference_genome, input_fastq1
        ]

    print(f"Running bwa: {' '.join(cmd)} > {aligned_sam}")

    # Run the command and redirect output to SAM file
    with open(aligned_sam, "w") as sam_file:
        subprocess.run(cmd, stdout=sam_file, check=True)

    return aligned_sam

def convert_sam_to_bam(sam_file, output_dir):
    """
    Convert SAM file to BAM format.
    Returns the path to the BAM file.
    """
    bam_file = os.path.join(output_dir, "aligned.bam")
    cmd = ["samtools", "view", "-Sb", sam_file, "-o", bam_file]
    print(f"Converting SAM to BAM: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return bam_file

def sort_bam(bam_file, output_dir):
    """
    Sort BAM file by reference coordenates.
    Returns the path to the sorted BAM file.
    """
    sorted_bam = os.path.join(output_dir, "sorted.bam")
    cmd = ["samtools", "sort", bam_file, "-o", sorted_bam]
    print(f"Sorting BAM: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return sorted_bam

def index_bam(bam_file):
    """
    Sorted bam files are Indexed.
    """
    cmd = ["samtools", "index", bam_file]
    print(f"Indexing BAM: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_mpileup_and_call_variants(reference_genome, bam_file, output_dir):
    """
    Run samtools mpileup and bcftools call to generate a VCF file.
    Returns the path to the VCF file.
    """
    mpileup_vcf = os.path.join(output_dir, "variants.vcf")
    cmd = [
        "bcftools", "mpileup", "-f", reference_genome, bam_file, "|",
        "bcftools", "call", "-mv", "-Ov", "-o", mpileup_vcf
    ]
    print(f"Running mpileup and variant calling: {' '.join(cmd)}")
    subprocess.run(cmd, shell=True, check=True)
    return mpileup_vcf

def filter_variants(input_vcf, output_dir):
    """
    Filter variants with quality > 20 using bcftools.
    Returns the path to the filtered VCF file.
    """
    filtered_vcf = os.path.join(output_dir, "filtered_variants.vcf")
    cmd = [
        "bcftools", "view", "-i", "'QUAL>20'", input_vcf, "-o", filtered_vcf
    ]
    print(f"Filtering variants: {' '.join(cmd)}")
    subprocess.run(cmd, shell=True, check=True)
    return filtered_vcf

def generate_metrics(output_dir, sorted_bam=None, vcf_file=None):
    """
    Generate various metrics for the files created in the pipeline.
    """
    metrics_dir = os.path.join(output_dir, "metrics")
    os.makedirs(metrics_dir, exist_ok=True)

    # Generate BAM Metrics
    if aligned_bam:
        flagstat_file = os.path.join(metrics_dir, "bam_flagstat.txt")
        idxstats_file = os.path.join(metrics_dir, "bam_idxstats.txt")
        coverage_file = os.path.join(metrics_dir, "bam_coverage.txt")

        # Flagstat
        cmd_flagstat = f"samtools flagstat {aligned_bam} > {flagstat_file}"
        print(f"Generating flagstat metrics for {aligned_bam}")
        subprocess.run(cmd_flagstat, shell=True, check=True)

        # Idxstats
        cmd_idxstats = f"samtools idxstats {aligned_bam} > {idxstats_file}"
        print(f"Generating idxstats metrics for {aligned_bam}")
        subprocess.run(cmd_idxstats, shell=True, check=True)

        # Coverage
        cmd_coverage = f"samtools depth {aligned_bam} > {coverage_file}"
        print(f"Generating coverage metrics for {aligned_bam}")
        subprocess.run(cmd_coverage, shell=True, check=True)

    # Generate VCF Metrics
    if vcf_file:
        vcf_stats = os.path.join(metrics_dir, "vcf_stats.txt")
        cmd_vcf = f"bcftools stats {vcf_file} > {vcf_stats}"
        print(f"Generating VCF stats for {vcf_file}")
        subprocess.run(cmd_vcf, shell=True, check=True)

    cmd_multiqc = [
        "multiqc", metrics_dir,"-o","/shared/multiqc_outputs"
    ]
    subprocess.run(" ".join(cmd_multiqc), shell=True, check=True)

    print(f"Metrics generated and saved in {metrics_dir}")


def main():
    """
    Main function to parse arguments and run the bioinformatics pipeline.
    """
    parser = argparse.ArgumentParser(description="Bioinformatics pipeline for quality control, alignment, and variant calling.")
    parser.add_argument(
        "input_files",
        nargs="+",
        help="Input FASTQ files (1 for single-end, 2 for paired-end). Can be GZipped."
    )
    parser.add_argument(
        "reference_genome",
        help="Path to the reference genome (FASTA format)."
    )
    parser.add_argument(
        "output_dir",
        help="Path to the output directory."
    )
    args = parser.parse_args()

    # Validate inputs
    input_files = args.input_files
    reference_genome = args.reference_genome
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    #Step 1: Run fastp
    print("Step 1: Running fastp")
    cleaned_fastq1, cleaned_fastq2 = run_fastp(input_files[0], input_files[1] if len(input_files) > 1 else None, output_dir)

    # Step 2: Run bwa-mem2
    print("Step 2: Running bwa")
    aligned_sam = run_bwa(reference_genome, cleaned_fastq1, cleaned_fastq2, output_dir)

    # Step 3: Convert SAM to BAM
    print("Step 3: Converting SAM to BAM")
    aligned_bam = convert_sam_to_bam(aligned_sam, output_dir)

    # Step 4: Sort BAM file
    print("Step 4: Sorting BAM file")
    sorted_bam = sort_bam(aligned_bam, output_dir)

    # Step 5: Index BAM file
    print("Step 5: Indexing BAM file")
    index_bam(sorted_bam)

    # Step 6: Generate mpileup and call variants
    print("Step 6: Running mpileup and calling variants")
    mpileup_vcf = run_mpileup_and_call_variants(reference_genome, sorted_bam, output_dir)

    # Step 7: Filter variants by quality
    print("Step 7: Filtering variants with quality > 20")
    filtered_vcf = filter_variants(mpileup_vcf, output_dir)

    # Step 8: Generate metrics
    print("Step 8: Generating metrics for pipeline outputs")
    generate_metrics(
        output_dir=output_dir,
        sorted_bam=sorted_bam,
        vcf_file=filtered_vcf
    )

    print(f"Pipeline completed. Results are in {output_dir}")

if __name__ == "__main__":
    main()
