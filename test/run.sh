#!/usr/bin/env bash
# Acceptance test for oxo-flow-sarek port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow 2>&1 | grep -q '{sample}' && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "==> somatic pairs branch: dry-run with a filled pair sheet + call_tiddit"
# Upstream nf-core/sarek 3.10.0 BAM_VARIANT_CALLING_SOMATIC_TIDDIT: TIDDIT
# per pair member ({id}.tiddit.normal / {id}.tiddit.tumor, no --skip_assembly
# — the pair legs receive the bwa index) + SVDB_MERGE (args2 --output-type z).
# The port drives it through the [[pairs]] fan-out from
# config/somatic_pairs.tsv; the shipped sheet is header-only (no row = no
# somatic instance). This block adds one pair row to the sheet copy, flips
# call_tiddit, and asserts the somatic rules schedule while germline tiddit
# also schedules (its flag is the same upstream tools flag).
sed 's/^call_tiddit = false$/call_tiddit = true/' main.oxoflow > .somatic-test-tmp.oxoflow
printf 'pair_id\texperiment\tcontrol\npair1\ttest\ttest2\n' > config/somatic_pairs.tsv
trap 'rm -f .somatic-test-tmp.oxoflow; git checkout -q config/somatic_pairs.tsv' EXIT
# NOTE: no --samples here — the pair instances come from the filled
# config/somatic_pairs.tsv sheet, and --samples first:1 would prune the
# pair fan-out before the when-gate evaluates.
"$OXO" dry-run .somatic-test-tmp.oxoflow > /tmp/oxo-dryrun-somatic-$$.txt 2>&1
grep -qE "^  [0-9]+\. tiddit_sv_somatic_pair1[^ ]*  \[run" /tmp/oxo-dryrun-somatic-$$.txt \
    || { echo "somatic branch: tiddit_sv_somatic not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. svdb_merge_tiddit_pair1[^ ]*  \[run" /tmp/oxo-dryrun-somatic-$$.txt \
    || { echo "somatic branch: svdb_merge_tiddit not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. tiddit_sv_cohort_test[^ ]*  \[run" /tmp/oxo-dryrun-somatic-$$.txt \
    || { echo "somatic branch: germline tiddit_sv not scheduled"; exit 1; }
# Default config must have no somatic tiddit instance (header-only sheet).
if grep -qE "^  [0-9]+\. (tiddit_sv_somatic|svdb_merge_tiddit)[^ ]*  \[run" /tmp/oxo-dryrun-$$.txt; then
    echo "somatic branch: rule scheduled with default config"; exit 1
fi
rm -f .somatic-test-tmp.oxoflow
git checkout -q config/somatic_pairs.tsv
trap - EXIT
echo "  tiddit_sv_somatic + svdb_merge_tiddit on with a filled pair sheet; off by default"

echo "==> somatic freebayes/mpileup branch: dry-run with call_freebayes + call_mpileup"
# Upstream nf-core/sarek 3.10.0: paired somatic FreeBayes is a joint two-BAM
# call (normal first, tumor second) followed by BCFTOOLS_SORT / TABIX /
# VCFLIB_VCFFILTER / TABIX and the per-caller VCF_QC_BCFTOOLS_VCFTOOLS fan-out;
# the paired mpileup leg feeds Control-FREEC only (no published VCF), while the
# tumor-only path publishes both freebayes and mpileup VCFs into vcf_all.
# This block flips both flags and adds a pair row plus a tumor-only row
# (empty control) to the sheet copy, then asserts the paired chain and QC
# rules schedule alongside their tumor-only instances and germline
# counterparts. VEP legs stay skipped (vep_cache_ready = false by default).
sed 's/^call_freebayes = false$/call_freebayes = true/; s/^call_mpileup = false$/call_mpileup = true/' \
    main.oxoflow > .somatic-test-tmp.oxoflow
printf 'pair_id\texperiment\tcontrol\npair1\ttest\ttest2\ntumoronly\ttest\t\n' > config/somatic_pairs.tsv
trap 'rm -f .somatic-test-tmp.oxoflow; git checkout -q config/somatic_pairs.tsv' EXIT
# NOTE: no --samples here — same reason as the tiddit block above.
"$OXO" dry-run .somatic-test-tmp.oxoflow > /tmp/oxo-dryrun-fb-somatic-$$.txt 2>&1
for r in freebayes_somatic_pair1 \
         bcftools_mpileup_call_somatic_pair1 \
         freebayes_somatic_tumor_only_tumoronly \
         bcftools_mpileup_call_somatic_tumor_only_tumoronly \
         bcftools_sort_freebayes_somatic_pair1 \
         tabix_freebayes_somatic_pair1 \
         vcffilter_freebayes_somatic_pair1 \
         tabix_freebayes_somatic_filt_pair1 \
         bcftools_stats_freebayes_somatic_pair1 \
         vcftools_filter_summary_freebayes_somatic_pair1 \
         vcftools_tstv_count_freebayes_somatic_pair1 \
         vcftools_tstv_qual_freebayes_somatic_pair1; do
    grep -qE "^  [0-9]+\. ${r}[^ ]*  \[run" /tmp/oxo-dryrun-fb-somatic-$$.txt \
        || { echo "somatic freebayes/mpileup branch: ${r} not scheduled"; exit 1; }
done
grep -qE "^  [0-9]+\. freebayes_cohort_test[^ ]*  \[run" /tmp/oxo-dryrun-fb-somatic-$$.txt \
    || { echo "somatic freebayes/mpileup branch: germline freebayes not scheduled"; exit 1; }
grep -qE "^  [0-9]+\. bcftools_mpileup_call_cohort_test[^ ]*  \[run" /tmp/oxo-dryrun-fb-somatic-$$.txt \
    || { echo "somatic freebayes/mpileup branch: germline mpileup not scheduled"; exit 1; }
# Default config must have no somatic freebayes/mpileup instance.
if grep -qE "^  [0-9]+\. [a-z_]*(freebayes|mpileup)_somatic[^ ]*  \[run" /tmp/oxo-dryrun-$$.txt; then
    echo "somatic freebayes/mpileup branch: rule scheduled with default config"; exit 1
fi
rm -f .somatic-test-tmp.oxoflow
git checkout -q config/somatic_pairs.tsv
trap - EXIT
echo "  freebayes/mpileup somatic chain + tumor-only on with filled sheet; off by default"

echo "PASS"
