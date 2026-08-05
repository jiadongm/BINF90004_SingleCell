#!/usr/bin/env bash
# Assembles everything that goes to Zenodo.
#
#   bash analysis/make_zenodo_bundle.sh      # run from BINF_scRNA-seq_2026/
#
# Output: zenodo_upload/ — upload the whole contents of that folder as the
# files of a single Zenodo record.
#
# Deliberately EXCLUDED from the student bundle (build artefacts the vignettes
# never read): gene_annotation.rds, pseudobulk.rds, pseudorep_curve.rds.
#
# Deliberately NOT rehosted: heart-counts.Rds (440 MB raw counts). It is already
# public under CC-BY-4.0 at zenodo.org/records/18237749 and students do not need
# it — P1_start.rds replaces it. The record description links there instead.

set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
OUT="$ROOT/zenodo_upload"
STAGE="$OUT/.stage"

rm -rf "$OUT"; mkdir -p "$STAGE/data/checkpoints" "$STAGE/scripts"

echo "1. Staging the student bundle"
cp dist/P1_start.rds       "$STAGE/data/"
cp dist/celltype_key.rds   "$STAGE/data/"
cp dist/cluster_labels.rds "$STAGE/data/"
cp dist/markers_cache.rds  "$STAGE/data/"
cp dist/de_demo.rds        "$STAGE/data/"
cp dist/checkpoints/02_integrated_clustered.rds "$STAGE/data/checkpoints/"

echo "2. Staging provenance scripts"
cp analysis/prep_P1_data.R analysis/prep_P2_data.R analysis/make_de_demo.R \
   analysis/make_zenodo_bundle.sh "$STAGE/scripts/"

echo "3. Zipping"
( cd "$STAGE" && zip -qr "$OUT/binf90004_scrnaseq_data_2026.zip" data )
( cd "$STAGE" && zip -qr "$OUT/generation_scripts.zip" scripts )

echo "4. Checksums"
( cd "$OUT" && shasum -a 256 ./*.zip | sed 's|\./||' > CHECKSUMS.txt )

rm -rf "$STAGE"

echo
echo "=== zenodo_upload/ ==="
ls -lh "$OUT"
echo
echo "unzipped bundle contents:"
unzip -l "$OUT/binf90004_scrnaseq_data_2026.zip" | tail -n +4 | head -12
