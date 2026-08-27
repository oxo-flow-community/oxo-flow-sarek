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

Every optional upstream branch is ported alongside the default path and gated
by config flags: reference preparation (`prepare_reference`), the `bwa-mem2`
aligner (`aligner = "bwa-mem2"`), BAM output mode
(`save_output_as_bam = true`), UMI consensus preprocessing
(`umi_read_structure = "3M2S+T"`), optional callers FreeBayes, Strelka2,
Manta, bcftools mpileup, TIDDIT, goleft indexcov (WGS only), DeepVariant and
NGSCheckMate (`call_*` / `tools_ngscheckmate`), and joint germline genotyping
with VQSR (`joint_germline = true`), per-chromosome scatter/gather for the
GATK steps (`scatter_gatk = true`), and per-caller VCF QC + VEP annotation
(auto-enabled with each `call_*` caller). The default configuration
reproduces the upstream default plan exactly — every added branch appears
only as a skip until its flag is set.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
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
  - a BWA index directory (`bwa_index_dir`) — and `bwa_mem2_index_dir` when
    using the `bwa-mem2` aligner (or build both with `prepare_reference = true`);
  - GATK bundle known-sites VCFs with `.tbi` (`dbsnp` + `known_indels`):
    `dbsnp_146.hg38.vcf.gz`, `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`,
    `Homo_sapiens_assembly38.known_indels.vcf.gz`;
  - joint germline only: the 1000G omni2.5 SNP VCF (`known_snps` +
    `known_snps_tbi`) for VQSR;
  - NGSCheckMate only: the SNP position BED (`ngscheckmate_bed`);
  - a VEP cache (GRCh38, homo_sapiens, cache version 112 — must match the
    ensembl-vep 112 binary in envs/vep.yaml; a gtf2vep-built subset cache
    works too) at `vep_dir_cache` (default `/.vep`, as upstream) + `vep_cache_ready = true`
    for the annotation step (upstream fails hard without the cache; the
    port gates the VEP rule on the flag instead).
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

Rows cover every upstream process/rule on the default main execution path plus
every ported optional branch (gated by config flags — all default-off).
"not ported" rows carry a reason + evidence.

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| FASTQC | `fastqc` | fastqc 0.12.1 | identical command |
| FASTP | `fastp` | fastp 1.1.0 | upstream default trimming/splitting module (TrimGalore does not exist in 3.10.0 — `grep -ri trimgalore` over the upstream tree is empty; the `--trim_fastq_trimgalore` param was dropped in 3.10.0) |
| BWA_MEM | `bwa_mem` | bwa 0.7.19 | upstream default aligner is `bwa-mem`, not `bwa-mem2`; read-group flags from sarek.nf; prefix `{meta.id}.{reads[0] token}` = `test.0001` under split_fastq; the same image also carries samtools 1.22.1 (used for `samtools sort`) |
| BWA_MEM2 | `bwa_mem2` | bwa-mem2 2.2.1 | `aligner = "bwa-mem2"` (upstream `params.aligner`); same `-K 100000000 -Y -R` args as BWA_MEM; shares bwa_mem's output path — the two rules are mutually exclusive via `when`, all downstream rules are unchanged; index from `bwa_mem2_index_dir` |
| GATK4_MARKDUPLICATES | `gatk_markduplicates` | gatk4 4.6.2.0 | default CRAM branch (`save_output_as_bam=false`) |
| GATK4_MARKDUPLICATES (BAM branch) | `gatk_markduplicates_bam` | gatk4 4.6.2.0 | `save_output_as_bam = true`; `--CREATE_INDEX true`, no CRAM conversion; `.md.bai` renamed to `.md.bam.bai` (upstream BAM_MERGE_INDEX_SAMTOOLS); downstream rules read `{config.alignment_ext}` / `{config.recal_index_ext}` (set `"bam"` / `"bam.bai"` together) |
| MOSDEPTH (post-MD) | `mosdepth_md` | mosdepth 0.3.14 | ext.prefix `{meta.id}.md`; WGS `--by 500` mode |
| SAMTOOLS_STATS (post-MD) | `samtools_stats_md` | samtools 1.24 | ext.prefix `{meta.id}.md.{alignment_ext}`; 1.24 is the version in the `htslib_samtools` stats/index images (the BWA image carries 1.22.1) |
| GATK4_BASERECALIBRATOR | `gatk_baserecalibrator` | gatk4 4.6.2.0 | known-sites = dbsnp + Mills gold standard + known indels (GRCh38); single whole-genome job by default; per-chromosome jobs under `scatter_gatk = true` (see the scatter/gather row) |
| GATK4_APPLYBQSR | `gatk_applybqsr` | gatk4 4.6.2.0 | output CRAM per default (BAM in the `save_output_as_bam` branch); single whole-genome job by default; per-chromosome jobs under `scatter_gatk = true` (see the scatter/gather row) |
| SAMTOOLS_INDEX (recal) | `samtools_index_recal` | samtools 1.24 | indexes the recalibrated alignment (`.crai` or `.bam.bai`) |
| MOSDEPTH (recal) | `mosdepth_recal` | mosdepth 0.3.14 | ext.prefix `{meta.id}.recal` |
| SAMTOOLS_STATS (recal) | `samtools_stats_recal` | samtools 1.24 | ext.prefix `{meta.id}.recal.{alignment_ext}` |
| GATK4_HAPLOTYPECALLER | `gatk_haplotypecaller` | gatk4 4.6.2.0 | default `tools=haplotypecaller,vep` → `call_haplotypecaller=true`; single-sample mode (no `-ERC GVCF`), `--pcr-indel-model CONSERVATIVE`; gated off when `joint_germline = true` (upstream picks the GVCF branch); single whole-genome job by default; per-chromosome jobs under `scatter_gatk = true` (see the scatter/gather row) |
| GATK4_CNNSCOREVARIANTS | `gatk_cnnscorevariants` | gatk4 4.6.2.0 | VCF_VARIANT_FILTERING_GATK part 1; CNN 1D scoring (module default `--tensor-type 1D`); upstream keeps the `{sample}.cnn.vcf.gz` intermediate unpublished — the port stores it under `results/variant_calling/cnnscorevariants/` for DAG handoff; skipped in joint mode (upstream: no filtering of the joint VCF) |
| GATK4_FILTERVARIANTTRANCHES | `gatk_filtervarianttranches` | gatk4 4.6.2.0 | VCF_VARIANT_FILTERING_GATK part 2; ext.args `--info-key CNN_1D`, ext.prefix `{meta.id}.haplotypecaller`, known sites (dbsnp + 2 GRCh38 indel sets) passed as `--resource`; produces `{sample}.haplotypecaller.filtered.vcf.gz` |
| BCFTOOLS_STATS | `bcftools_stats` | bcftools 1.23.1 | VCF_QC_BCFTOOLS_VCFTOOLS part 1; runs on the filtered VCF (prefix `{meta.id}.haplotypecaller.filtered`) |
| VCFTOOLS_TSTV_COUNT | `vcftools_tstv_count` | vcftools 0.1.17 | VCF_QC_BCFTOOLS_VCFTOOLS part 2; runs on the filtered VCF |
| VCFTOOLS_TSTV_QUAL | `vcftools_tstv_qual` | vcftools 0.1.17 | VCF_QC_BCFTOOLS_VCFTOOLS part 3; runs on the filtered VCF |
| VCFTOOLS_SUMMARY | `vcftools_filter_summary` | vcftools 0.1.17 | VCF_QC_BCFTOOLS_VCFTOOLS part 4; runs on the filtered VCF |
| ENSEMBLVEP_VEP | `ensemblvep_vep` | ensembl-vep 112.0 | annotates the filtered VCF (`{sample}.haplotypecaller.filtered_VEP.ann.vcf.gz`); **gated on `vep_cache_ready`** — upstream fails hard without the cache (bundled at `/.vep` via `--vep_cache`); set the flag after placing a cache at `vep_dir_cache`; `--cache_version 112` (matches the env binary — VEP caches are version-locked), GRCh38 |
| PREPARE_GENOME (BWA_INDEX) | `bwa_index` | bwa 0.7.19 | `prepare_reference = true`; builds `results/reference/bwa/index.{amb,ann,bwt,pac,sa}` — fixed `index` prefix (upstream: fasta basename; irrelevant downstream, the BWA rules find the index by extension). Deviation: upstream does not publish the index unless `save_reference`; the port publishes it because oxo-flow outputs must be tracked |
| PREPARE_GENOME (BWAMEM2_INDEX) | `bwamem2_index` | bwa-mem2 2.2.1 | same gating/prefix note; `results/reference/bwamem2/index.{0123,amb,ann,bwt.2bit.64,pac}` |
| PREPARE_GENOME (GATK4_CREATESEQUENCEDICTIONARY) | `gatk_createsequencedictionary` | gatk4 4.6.2.0 | `--URI` is the fasta basename as upstream; GATK writes `<fasta-basename>.dict`, the port renames it to `reference.dict` for a fixed output path |
| PREPARE_GENOME (SAMTOOLS_FAIDX) | `samtools_faidx` | samtools 1.24 | indexes a workdir copy of the fasta (the `{config.fasta}` path is treated as read-only); output `results/reference/fai/reference.fasta.fai` |
| SAREK_UMI (FGBIO_FASTQTOBAM → SAMTOOLS_BAM2FQ → BWA_MEM → FGBIO_GROUPREADSBYUMI → FGBIO_CALLMOLECULARCONSENSUSREADS → BAM_CONVERT_SAMTOOLS → FASTP) | `fgbio_fastqtobam`, `samtools_bam2fq_umi`, `bwa_mem_umi`, `fgbio_groupreadsbyumi`, `fgbio_callmolecularconsensusreads`, `samtools_bam2fq_consensus`, `fastp_umi` | fgbio 3.1.2, bwa 0.7.19, samtools 1.24, fastp 1.1.0 | `umi_read_structure` set (e.g. `"3M2S+T"`); `fastp_umi` shares the fastp output paths so BWA-MEM downstream is untouched; GroupReadsByUmi histogram/metrics go to `reports/umi/`; see deviations for the BAM_CONVERT_SAMTOOLS collapse |
| BAM_VARIANT_CALLING_FREEBAYES (FREEBAYES_GERMLINE, BCFTOOLS_SORT, TABIX_VC, VCFLIB_VCF_FILTER, TABIX_FILT) | `freebayes`, `bcftools_sort_freebayes`, `tabix_freebayes`, `vcffilter_freebayes`, `tabix_freebayes_filt` | freebayes 1.3.10, bcftools 1.23.1, vcflib 1.0.14 | `call_freebayes = true`; `--min-alternate-fraction 0.1 --min-mapping-quality 1`; QUAL filter threshold from `freebayes_filter` (30, upstream `params.freebayes_filter`); all four VCFs/TBI published to `variant_calling/freebayes/{sample}/` as upstream |
| STRELKA_GERMLINE | `strelka_germline` | strelka 2.9.10 | `call_strelka = true`; ext.prefix `{meta.id}.strelka`; email-check disabled via the upstream `sed`; all six outputs (SNP/INDEL/variants × vcf+tbi) published as upstream |
| MANTA_GERMLINE | `manta_germline` | manta 1.6.0 | `call_manta = true`; ext.prefix `{meta.id}.manta`; only the `diploid_sv` pair is published — candidateSmallIndels/candidateSV stay in the workdir exactly as upstream |
| BCFTOOLS_MPILEUP (germline) | `bcftools_mpileup_call` | bcftools 1.23.1 | `call_mpileup = true`; args `--output-type v --multiallelic-caller`, filter `count(GT=="RR")==0`; the module's own bcftools_stats file stays unpublished (as upstream) |
| TIDDIT_SV + TABIX_BGZIP_TIDDIT_SV | `tiddit_sv`, `tabix_tiddit` | tiddit 3.9.5 | `call_tiddit = true`; `--skip_assembly` (upstream passes an empty bwa index channel for germline); `.ploidies.tab` published as upstream |
| BAM_VARIANT_CALLING_INDEXCOV (SAMTOOLS_REINDEX_BAM + GOLEFT_INDEXCOV) | `samtools_reindex_bam`, `goleft_indexcov` | samtools 1.24, goleft 0.2.4 | WGS only (`!wes && call_indexcov`); per-sample header-only reindex with `-F 3844 -q 30` + `--write-index` over `/dev/null##idx##`; cohort run `--fai --directory indexcov` (no `--extranormalize` — inputs are BAMs, matching upstream's reindex path); bed.gz+tbi published to `variant_calling/indexcov/` |
| RUNDEEPVARIANT | `deepvariant` | deepvariant 1.10.0 | `call_deepvariant = true`; `--model_type=WGS --sample_name {sample}`; vcf + g.vcf pairs published to `variant_calling/deepvariant/{sample}/` |
| BAM_NGSCHECKMATE (BCFTOOLS_MPILEUP + NGSCHECKMATE_NCM) | `bcftools_mpileup_ngscheckmate`, `ngscheckmate_ncm` | bcftools 1.23.1, ngscheckmate 1.0.1 | `tools_ngscheckmate = true`; per-sample mpileup `--no-version --ploidy 1 -c` with `-T` SNP bed, reheader to `{sample}-{lane}`; cohort `NCM_REF=./reference.fasta ncm.py -d . -bed <bed> -O . -N ngscheckmate -V`; outputs published to `reports/ngscheckmate/` (live-verify: ncm.py's exact output filenames, see Test) |
| Joint germline (GATK4_HAPLOTYPECALLER GVCF, GATK4_GENOMICSDBIMPORT, GATK4_GENOTYPEGVCFS, BCFTOOLS_SORT, GATK4_MERGEVCFS, GATK4_VARIANTRECALIBRATOR SNP+INDEL, GATK4_APPLYVQSR SNP+INDEL) | `gatk_haplotypecaller_gvcf`, `gatk_genomicsdbimport`, `gatk_genotypegvcfs`, `bcftools_sort_joint`, `gatk_mergevcfs_joint`, `gatk_variantrecalibrator_snp`, `gatk_variantrecalibrator_indel`, `gatk_applyvqsr_snp`, `gatk_applyvqsr_indel` | gatk4 4.6.2.0, bcftools 1.23.1 | `joint_germline = true`; VQSR resource labels from `conf/igenomes.config` GRCh38 (1000G omni2.5 SNP → `known_snps`, dbsnp; gatk+mills indels); upstream prefixes `joint_variant_calling_SNP/INDEL` (VQSR intermediates unpublished upstream — the port keeps them under `results/` for DAG handoff); final `joint_germline_recalibrated.vcf.gz`; one whole-genome interval from the fasta `.fai` by default, per-chromosome under `scatter_gatk = true` (see the scatter/gather row) |
| Joint VCF QC + VEP | `bcftools_stats_joint`, `vcftools_tstv_count_joint`, `vcftools_tstv_qual_joint`, `vcftools_filter_summary_joint`, `ensemblvep_vep_joint` | bcftools 1.23.1, vcftools 0.1.17, ensembl-vep 112.0 | upstream runs VCF_QC + VEP on the joint VCF (`vcf_all`); the per-sample QC/VEP rules are gated off in joint mode and these cohort rules take over (prefix `joint_germline_recalibrated`) |
| MULTIQC | `multiqc` | multiqc 1.35 | fan-in over all report producers (depends_on covers the gated branches — skipped rules auto-satisfy); scans the results dir with the upstream `assets/multiqc_config.yml` |
| PREPARE_INTERVALS (BUILD_INTERVALS, CREATE_INTERVALS_BED, TABIX_BGZIPTABIX) + per-interval scatter/gather of BQSR / ApplyBQSR / HaplotypeCaller / joint GenotypeGVCFs (GATK4_GATHERBQSRREPORTS, CRAM/BAM_MERGE_INDEX_SAMTOOLS, GATK4_MERGEVCFS) | `create_intervals_bed`, `tabix_interval`, `gatk_baserecalibrator_scatter`, `gatk_gatherbqsrreports`, `gatk_applybqsr_scatter`, `merge_index_samtools`, `gatk_haplotypecaller_scatter`, `gatk_mergevcfs_scatter`, `gatk_haplotypecaller_gvcf_scatter`, `gatk_genomicsdbimport_scatter`, `gatk_genotypegvcfs_scatter`, `bcftools_sort_joint_scatter`, `gatk_mergevcfs_joint_scatter` | gawk 5.3.0, samtools 1.24, gatk4 4.6.2.0 | `scatter_gatk = true` (default off). Deviation: the engine's scatter fan-out takes a **static** value list, so intervals are one per chromosome (`config.chromosomes` — keep it in sync with the fasta `.fai`) instead of upstream's duration-binned windows (`nucleotides_per_second`); every downstream gather is exact and writes the same paths as the single-job branch, so results are identical — the branch adds per-contig parallelism, it does not change outputs. Needs live verification (see Test) |
| fastp split parts beyond `0001.` (multi-part BWA_MEM + BAM_MERGE_INDEX_SAMTOOLS) | — not ported | — | **structural**: upstream's channel fan-out over split parts cannot be expressed as fixed output paths; the port caps input at one split part (~50M read pairs / 200M lines per sample). **Do not run WGS data above the cap** with this port |
| Somatic callers (Mutect2, somatic Strelka2/Manta, CNVkit, ASCAT, MSIsensor2/pro, SomaticSniper, VarDict, Control-FREEC, LoFreq, Varlociraptor) | — not ported | — | tumor/normal **pairs required**; the port's samplesheet is single-sample germline (`[[sample_groups]]`; the engine's `[[pairs]]` mechanism is a possible follow-up) |
| Sentieon / Parabricks / DRAGMAP | — not ported | — | commercial accelerators (licensed binaries); out of scope |
| VCF_QC + ENSEMBLVEP_VEP fan-out over the optional callers | `bcftools_stats_{freebayes,strelka,mpileup,deepvariant,manta,tiddit}`, `vcftools_tstv_count_{...}`, `vcftools_tstv_qual_{...}`, `vcftools_filter_summary_{...}`, `ensemblvep_vep_{...}` | bcftools 1.23.1, vcftools 0.1.17, ensembl-vep 112.0 | upstream runs VCF_QC + VEP on every caller VCF (`vcf_all`); the port now mirrors that: when `call_freebayes`/`call_strelka`/`call_mpileup`/`call_deepvariant`/`call_manta`/`call_tiddit` enables a caller, its VCF is QC'd (`reports/bcftools/<caller>/`, `reports/vcftools/<caller>/`) and annotated (`annotation/<caller>/`) with the same prefix conventions as the haplotypecaller rules; tiddit's uncompressed `.vcf` uses `--vcf` instead of `--gzvcf` (as upstream). All 30 rules are gated on their caller's flag (+ `skip_bcftools`/`skip_vcftools`/`annotate_vep`/`vep_cache_ready`) and feed the multiqc `depends_on`. Needs live verification (see Test) |

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
- **scatter/gather is opt-in and uses per-chromosome intervals**: the
  upstream per-interval fan-out (duration-binned windows) becomes
  `scatter_gatk = true` in the port — one interval per chromosome (the
  engine's scatter needs a static value list, so `config.chromosomes`
  replaces upstream's `nucleotides_per_second` binning; both split a
  whole-genome interval set, but per-chromosome granularity is coarser).
  Off by default: BQSR, ApplyBQSR, HaplotypeCaller (and joint
  GenotypeGVCFs) run as single whole-genome jobs (one interval from the
  fasta `.fai` for the joint path). Every gather (GatherBQSRReports,
  samtools merge+index, MergeVcfs) writes the same output paths in both
  branches, so the branches are exchangeable without touching downstream
  rules.
- **BAM_CONVERT_SAMTOOLS collapse (UMI path)**: at this upstream commit the
  four `samtools view` calls in `bam_convert_samtools` have no distinguishing
  `-f/-F` flags anywhere in `conf/` (grep-verified), so the view → merge →
  collate → cat machinery emits four identical BAMs and **doubles the reads**
  in the merged output. The port reproduces the functional intent with a
  single `samtools collate -O | samtools fastq` instead.
- **CNN-scored intermediate location**: upstream disables the
  CNNSCOREVARIANTS publishDir (the `{sample}.cnn.vcf.gz` stays in the task
  workdir); oxo-flow hands files between rules through `results/`, so the
  port keeps it under `results/variant_calling/cnnscorevariants/` (same for
  the joint VQSR intermediates, unpublished upstream).
- **reference-prep outputs are published**: upstream leaves built
  BWA/BWAmem2 indexes, `.dict` and `.fai` in the task workdir unless
  `save_reference`; oxo-flow requires tracked outputs, so the port publishes
  them under `results/reference/` — point `bwa_index_dir` /
  `bwa_mem2_index_dir` / `fasta_fai` / `dict` at those paths to use them.
- **single-lane model**: `meta.id` = `{sample}` from BWA-MEM onward
  (upstream: `{sample}-{lane}`); preprocessing stage files keep the upstream
  `{sample}-{lane}` prefix. Read-group IDs are `{sample}.{lane}` as upstream.
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

The two newest gated branches need live verification (container runs, not yet
performed — planned as a follow-up wave):

- `scatter_gatk = true` (per-chromosome BQSR/ApplyBQSR/Haplotypecaller/joint
  GenotypeGVCFs + the GatherBQSRReports / samtools merge+index / MergeVcfs
  gathers) — dry-run covers both `joint_germline` modes;
- the per-caller VCF_QC + VEP fan-out (30 rules over
  `call_freebayes`/`call_strelka`/`call_mpileup`/`call_deepvariant`/
  `call_manta`/`call_tiddit`) — dry-run covers all six callers at once;
  tiddit's `--vcf` vcftools invocation in particular needs a live check.

A reasonable smoke sequence: `scatter_gatk=true` on the existing live-verified
germline dataset (results must be byte-identical to the default branch), then
one optional caller with QC+VEP on a small WES slice.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md); the upstream MIT license is included verbatim at
[LICENSE.upstream](LICENSE.upstream).
