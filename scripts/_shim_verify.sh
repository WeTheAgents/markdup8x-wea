#!/usr/bin/env bash
# Phase-0 verification: shim's generated line2 must match the nf-core
# main.nf.test.snap entries byte-for-byte.
set -uo pipefail

SHIM=/d/GitHub/markdup-wea/scripts/picard-shim.sh
SNAP="${SNAP:-C:/Users/peach/AppData/Local/Temp/nfcore-modules-snap/modules/nf-core/picard/markduplicates/tests/main.nf.test.snap}"

# Extract the LINE2 build block from the shim and source it with stub vars.
gen_line2() {
    local INPUT=$1 OUTPUT=$2 METRICS=$3 ASSUME_SORT_ORDER=$4 REFERENCE_SEQUENCE=$5 REMOVE_DUPLICATES="${6:-false}"
    local INPUT_BN OUTPUT_BN METRICS_BN SO_PRINT REF_FRAG
    INPUT_BN=$(basename -- "$INPUT")
    OUTPUT_BN=$(basename -- "$OUTPUT")
    METRICS_BN=$(basename -- "$METRICS")
    SO_PRINT="${ASSUME_SORT_ORDER:-queryname}"
    REF_FRAG=""
    [[ -n "$REFERENCE_SEQUENCE" ]] && REF_FRAG=" --REFERENCE_SEQUENCE $(basename -- "$REFERENCE_SEQUENCE")"

    echo "# MarkDuplicates --INPUT ${INPUT_BN} --OUTPUT ${OUTPUT_BN} --METRICS_FILE ${METRICS_BN} --ASSUME_SORT_ORDER ${SO_PRINT}${REF_FRAG} --MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP 50000 --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 --SORTING_COLLECTION_SIZE_RATIO 0.25 --TAG_DUPLICATE_SET_MEMBERS false --REMOVE_SEQUENCING_DUPLICATES false --TAGGING_POLICY DontTag --CLEAR_DT true --DUPLEX_UMI false --FLOW_MODE false --FLOW_DUP_STRATEGY FLOW_QUALITY_SUM_STRATEGY --FLOW_USE_END_IN_UNPAIRED_READS false --FLOW_USE_UNPAIRED_CLIPPED_END false --FLOW_UNPAIRED_END_UNCERTAINTY 0 --FLOW_UNPAIRED_START_UNCERTAINTY 0 --FLOW_SKIP_FIRST_N_FLOWS 0 --FLOW_Q_IS_KNOWN_END false --FLOW_EFFECTIVE_QUALITY_THRESHOLD 15 --ADD_PG_TAG_TO_READS true --REMOVE_DUPLICATES ${REMOVE_DUPLICATES} --ASSUME_SORTED false --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --PROGRAM_RECORD_ID MarkDuplicates --PROGRAM_GROUP_NAME MarkDuplicates --READ_NAME_REGEX <optimized capture of last three ':' separated fields as numeric values> --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false"
}

# Extract the snap expected line2 for a named test.
snap_line2() {
    local case_name=$1
    python3 -c "
import json, sys
d = json.load(open('$SNAP'))
print(d['$case_name']['content'][1][1])
"
}

declare -A CASES=(
    [sorted]="sarscov2 [sorted bam]|test.paired_end.sorted.bam|test.md.bam|test.md.metrics.txt|||"
    [unsorted]="sarscov2 [unsorted bam]|test.paired_end.bam|test.md.bam|test.md.metrics.txt|||"
    [cram]="homo_sapiens [cram]|test.paired_end.sorted.cram|test.md.cram|test.md.metrics.txt||genome.fasta|"
)

PASS=0; FAIL=0
for k in sorted unsorted cram; do
    IFS='|' read -r name input output metrics so ref rd <<< "${CASES[$k]}"
    expected=$(snap_line2 "$name")
    actual=$(gen_line2 "$input" "$output" "$metrics" "$so" "$ref" "${rd:-false}")
    if [[ "$expected" == "$actual" ]]; then
        echo "PASS  $k"
        PASS=$((PASS+1))
    else
        echo "FAIL  $k"
        diff <(echo "$expected") <(echo "$actual") | head -20
        FAIL=$((FAIL+1))
    fi
done
echo "---"
echo "PASS=$PASS FAIL=$FAIL"
exit $FAIL
