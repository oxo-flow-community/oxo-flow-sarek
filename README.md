# oxo-flow-sarek — WGS/WES germline and somatic variant calling

[![CI](https://github.com/oxo-flow-community/oxo-flow-sarek/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-sarek/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> ★ Verified · ⇄ Official port of [`nf-core/sarek`](https://github.com/nf-core/sarek) @ `3.10.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

Run the GATK best-practice variant-calling path for whole-genome and
whole-exome sequencing (WGS/WES) data, germline by default: FastQC quality
control, fastp trimming and splitting, BWA-MEM alignment with read-group
metadata, GATK MarkDuplicates with CRAM conversion, base quality score
recalibration (BQSR), single-sample HaplotypeCaller variant calling, CNN 1D
scoring with tranche filtering, VEP annotation, per-sample QC (mosdepth,
samtools stats, bcftools stats, vcftools TsTv/FILTER summaries) and a final
MultiQC report. You get CRAM alignments, recalibration tables, filtered and
annotated VCFs, per-sample QC reports and a MultiQC report.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.11.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/download/v0.11.0/oxo-flow-v0.11.0-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternatively via conda: `conda install -c bioconda oxo-flow-cli` (note: the
conda package may lag behind releases; other platform binaries are available
on the [releases page](https://github.com/Traitome/oxo-flow/releases)).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-sarek.git
cd oxo-flow-sarek
```

### 3. Requirements

- **Reference data (GRCh38, user-provided)** — point the `[config]` block of
  `main.oxoflow` at your bundle:
  - genome FASTA plus `.fai` and `.dict` (`fasta` / `fasta_fai` / `dict`);
  - a BWA index directory (`bwa_index_dir`);
  - GATK bundle known-sites VCFs with `.tbi` (`dbsnp` + `known_indels`):
    `dbsnp_146.hg38.vcf.gz`, `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`,
    `Homo_sapiens_assembly38.known_indels.vcf.gz`;
  - a VEP cache (GRCh38, homo_sapiens, cache version 116) mounted at `/.vep`
    in the container for the annotation step.
- **Input reads** — paired-end FASTQ at `raw/{sample}_R1.fastq.gz` and
  `raw/{sample}_R2.fastq.gz`; supported input size is capped at ~50M read
  pairs per sample (fastp split cap, see Fidelity).
- **Compute** — up to 24 CPUs / 36 GB per rule (BWA-MEM: 24 threads / 30G;
  VEP: 6 threads / 36G).
- **Tools** — Docker containers with pinned images, declared per rule in
  `main.oxoflow` (`[rules.environment]`); Docker is required at runtime (no
  conda environments are used).
- **Disk** — `results/` grows with per-sample CRAMs, VCFs and reports; the
  GRCh38 reference bundle and VEP cache require substantial disk space of
  their own.

## Usage

```bash
# 1. install oxo-flow (see Requirements)
# 2. prepare data: raw/<sample>_R1.fastq.gz / _R2.fastq.gz (see fixtures/)
# 3. point the reference config at your GRCh38 bundle
#    (main.oxoflow [config] fasta/fasta_fai/dict/bwa_index_dir/dbsnp/known_indels)
# 4. preview the plan
oxo-flow dry-run main.oxoflow
# 5. run
oxo-flow run main.oxoflow -j 8
# 6. run a subset
oxo-flow run main.oxoflow -t multiqc --samples first:2
```

The bundled fixtures (`test/fixtures/raw/test_R1.fastq.gz` /
`test_R2.fastq.gz`, from nf-core/test-datasets) let you exercise the DAG from
FastQC through BWA-MEM; the BWA-MEM step itself needs a GRCh38 index
(`config.bwa_index_dir`), and the BQSR/Haplotypecaller/VEP steps need the
corresponding GATK bundle files and a VEP cache.

## Source

Upstream: **[nf-core/sarek](https://github.com/nf-core/sarek)** @ `3.10.0`
(commit `8ccac7ad37b05dd792447763bf9671b719824587`), MIT license. Created
2026-08-15; this workflow may lag behind upstream releases. See
[NOTICE.md](NOTICE.md) for attribution.

## Fidelity

Rows cover every upstream process/rule on the default main execution path.
"not ported" rows carry a reason.

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| FASTQC | `fastqc` | fastqc 0.12.1 | identical command |
| FASTP | `fastp` | fastp 1.1.0 | upstream default trimming/fastp module; TrimGalore is the `--trim_fastq_trimgalore` alternative and is not ported |
| BWA_MEM | `bwa_mem` | bwa 0.7.19 | upstream default aligner is `bwa-mem`, not `bwa-mem2`; read-group flags from sarek.nf; prefix `{meta.id}.{reads[0] token}` = `test.0001` under split_fastq; the same image also carries samtools 1.22.1 (used for `samtools sort`) |
| GATK4_MARKDUPLICATES | `gatk_markduplicates` | gatk4 4.6.2.0 | CRAM-only port (`save_output_as_bam=false`); upstream's `if [[ ${prefix} == *.cram ]]` BAM branch not ported — conversion is unconditional |
| MOSDEPTH (post-MD) | `mosdepth_md` | mosdepth 0.3.14 | ext.prefix `{meta.id}.md`; WGS `--by 500` mode |
| SAMTOOLS_STATS (post-MD) | `samtools_stats_md` | samtools 1.24 | ext.prefix `{meta.id}.md.cram`; 1.24 is the version in the `htslib_samtools` stats/index images (the BWA image carries 1.22.1) |
| GATK4_BASERECALIBRATOR | `gatk_baserecalibrator` | gatk4 4.6.2.0 | known-sites = dbsnp + Mills gold standard + known indels (GRCh38); single whole-genome job, no per-interval scatter (see deviations) |
| GATK4_APPLYBQSR | `gatk_applybqsr` | gatk4 4.6.2.0 | output CRAM (not BAM) per default; single whole-genome job, no per-interval scatter (see deviations) |
| SAMTOOLS_INDEX (recal) | `samtools_index_recal` | samtools 1.24 | indexes the recalibrated CRAM |
| MOSDEPTH (recal) | `mosdepth_recal` | mosdepth 0.3.14 | ext.prefix `{meta.id}.recal` |
| SAMTOOLS_STATS (recal) | `samtools_stats_recal` | samtools 1.24 | ext.prefix `{meta.id}.recal.cram` |
| GATK4_HAPLOTYPECALLER | `gatk_haplotypecaller` | gatk4 4.6.2.0 | default `tools=haplotypecaller,vep` → `call_haplotypecaller=true`; single-sample mode (no `-ERC GVCF`), `--pcr-indel-model CONSERVATIVE`; single whole-genome job, no per-interval scatter (see deviations) |
| GATK4_CNNSCOREVARIANTS | `gatk_cnnscorevariants` | gatk4 4.6.2.0 | VCF_VARIANT_FILTERING_GATK part 1; CNN 1D scoring (module default `--tensor-type 1D`); upstream keeps the `{sample}.cnn.vcf.gz` intermediate unpublished — the port stores it under `results/variant_calling/cnnscorevariants/` for DAG handoff |
| GATK4_FILTERVARIANTTRANCHES | `gatk_filtervarianttranches` | gatk4 4.6.2.0 | VCF_VARIANT_FILTERING_GATK part 2; ext.args `--info-key CNN_1D`, ext.prefix `{meta.id}.haplotypecaller`, known sites (dbsnp + 2 GRCh38 indel sets) passed as `--resource`; produces `{sample}.haplotypecaller.filtered.vcf.gz` |
| BCFTOOLS_STATS | `bcftools_stats` | bcftools 1.23.1 | VCF_QC_BCFTOOLS_VCFTOOLS part 1; runs on the filtered VCF (prefix `{meta.id}.haplotypecaller.filtered`) |
| VCFTOOLS_TSTV_COUNT | `vcftools_tstv_count` | vcftools 0.1.17 | VCF_QC_BCFTOOLS_VCFTOOLS part 2; runs on the filtered VCF |
| VCFTOOLS_TSTV_QUAL | `vcftools_tstv_qual` | vcftools 0.1.17 | VCF_QC_BCFTOOLS_VCFTOOLS part 3; runs on the filtered VCF |
| VCFTOOLS_SUMMARY | `vcftools_filter_summary` | vcftools 0.1.17 | VCF_QC_BCFTOOLS_VCFTOOLS part 4; runs on the filtered VCF |
| ENSEMBLVEP_VEP | `ensemblvep_vep` | ensembl-vep 116.0 | annotates the filtered VCF (`{sample}.haplotypecaller.filtered_VEP.ann.vcf.gz`); requires a VEP cache mounted at `/.vep` in the container (upstream bundles it via `--vep_cache`); `--cache_version 116`, GRCh38 |
| MULTIQC | `multiqc` | multiqc 1.35 | fan-in over all report producers; scans the results dir with the upstream `assets/multiqc_config.yml` |
| NGSCheckMate (BAM_NGSCHECKMATE + BCFTOOLS_MPILEUP) | — not ported | — | sample-identity QC on the CRAM; outside this port's scope (QC fan-in, needs a cohort reference) |
| PREPARE_GENOME (BWA_INDEX, GATK4_CREATESEQUENCEDICTIONARY, SAMTOOLS_FAIDX) | — not ported | — | reference preparation; the port requires a pre-built reference bundle (upstream also accepts pre-built index/dict/fai via `params.bwa`/`dict`/`fasta_fai`, so this is an upstream-compatible shortcut, not a behavior change) |
| BED_PREPARE_INTERVALS (BUILD_INTERVALS, CREATE_INTERVALS_BED, TABIX bgzip/tabix interval split) + per-interval scatter/gather (GATK4_GATHERBQSRREPORTS, CRAM_MERGE_INDEX_SAMTOOLS, GATK4_MERGEVCFS) | — not ported | — | interval preparation and the per-interval scatter/gather of BQSR/ApplyBQSR/HaplotypeCaller; the port runs single whole-genome GATK jobs without `--intervals` (gathers are exact, so results are equivalent, but the per-contig parallelism of the upstream default path is absent) |
| SAREK_UMI_* / Strelka2 / Mutect2 / Manta / DeepVariant / CNVkit / ASCAT / MSIsensor / ... | — not ported | — | UMI workflows, `--tools` alternatives and optional callers; out of scope (default `haplotypecaller,vep` only) |

Deviations (all documented, nothing silently dropped):

- **fastp split parts — supported input cap**: upstream `--split_by_lines`
  produces `0001.`-prefixed outputs; only the first split part (`0001.`) is
  wired to BWA-MEM (the nf-core bwa/mem prefix logic `tokenize('.')[0]` gives
  `{sample}.0001.bam`). The supported input size is therefore capped at one
  split part — ~50M read pairs / 200M lines per sample (the upstream default
  `split_fastq=50,000,000` threshold). Datasets above the threshold are a
  **known unsupported upstream behavior**: upstream aligns every part and
  merges the part BAMs before MarkDuplicates (`BAM_MERGE_INDEX_SAMTOOLS`),
  which this port does not reproduce — **do not run WGS data above ~50M read
  pairs per sample** with this port until the multi-part path is wired.
- **no per-interval scatter/gather**: BQSR, ApplyBQSR and HaplotypeCaller run
  as single whole-genome jobs (no `--intervals`, no
  GATK4_GATHERBQSRREPORTS / CRAM_MERGE_INDEX_SAMTOOLS / GATK4_MERGEVCFS).
  Gathers are exact, so the results are mathematically equivalent to
  upstream's, but wall-time and per-contig parallelism differ substantially.
- **CNN-scored intermediate location**: upstream disables the
  CNNSCOREVARIANTS publishDir (the `{sample}.cnn.vcf.gz` stays in the task
  workdir); oxo-flow hands files between rules through `results/`, so the
  port keeps it under `results/variant_calling/cnnscorevariants/`.
- **markduplicates CRAM branch**: unconditional `samtools view -Ch` +
  `samtools index` replaces upstream's bash conditional; BAM output mode
  (`save_output_as_bam=true`) is not ported.
- **Docker staging**: only the rule workdir is mounted, so reference files are
  copied into the workdir under fixed local names (`reference.fasta`,
  `reference.fasta.fai`, `reference.dict`) — same effect as Nextflow's
  staging of the reference into the task directory.
- **`known_indels`** is a TOML array (2 GRCh38 files); both it and
  `known_indels_tbi` must be updated together when changing references.

## Test

```bash
bash test/run.sh
```

Runs `oxo-flow validate`, `oxo-flow lint` and a `dry-run` smoke check; CI runs
the same script on every push.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md); the upstream MIT license is included verbatim at
[LICENSE.upstream](LICENSE.upstream).
