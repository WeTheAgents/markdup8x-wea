#!/usr/bin/env bash
# picard-shim.sh — Picard MarkDuplicates CLI shim over markdup-wea.
#
# Goal: be invocable as `picard [JVM-flags] MarkDuplicates --KEY VALUE ...`
# (the nf-core/modules invocation) without Java, dispatching to markdup-wea
# and rewriting the metrics file's `# MarkDuplicates ...` line so the
# nf-core/modules nf-test snapshot matches Picard 3.4.0 byte-for-byte.
#
# See: docs/nfcore-invocation.md, plan medium-path-eager-spark.md (Stage 3a).
set -uo pipefail

WEA_BIN="${WEA_BIN:-markdup-wea}"
PICARD_VERSION_STRING="3.4.0"

die() { echo "picard-shim: $*" >&2; exit 2; }

# ---- Drop JVM-style args (-Xmx*, -Xms*, -D*) before the subcommand ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        -Xmx*|-Xms*|-XX:*|-D*) shift ;;
        -jar) shift 2 ;;     # tolerate `picard -jar picard.jar`
        MarkDuplicates) shift; break ;;
        --version)
            # `picard --version` (without subcommand). nf-core uses
            # `picard MarkDuplicates --version`; handle both.
            echo "Version:${PICARD_VERSION_STRING}" >&2
            exit 0
            ;;
        --) shift; break ;;
        *) die "unsupported pre-subcommand arg: $1" ;;
    esac
done

# ---- Parse MarkDuplicates --KEY VALUE pairs ----
INPUT=""
OUTPUT=""
METRICS=""
ASSUME_SORT_ORDER=""        # if unset → emulated Picard default 'queryname' in metrics line
REFERENCE_SEQUENCE=""
BARCODE_TAG=""
READ_ONE_BARCODE_TAG=""
READ_TWO_BARCODE_TAG=""
REMOVE_DUPLICATES="false"
WANT_VERSION=0

# Args we silently swallow (JVM tunings, Picard-only knobs, defaulted bools).
swallow_arg() {
    case "$1" in
        --MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP|--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP|\
        --SORTING_COLLECTION_SIZE_RATIO|--TAG_DUPLICATE_SET_MEMBERS|\
        --REMOVE_SEQUENCING_DUPLICATES|--TAGGING_POLICY|--CLEAR_DT|--DUPLEX_UMI|\
        --ADD_PG_TAG_TO_READS|--DUPLICATE_SCORING_STRATEGY|--PROGRAM_RECORD_ID|\
        --PROGRAM_GROUP_NAME|--READ_NAME_REGEX|--OPTICAL_DUPLICATE_PIXEL_DISTANCE|\
        --MAX_OPTICAL_DUPLICATE_SET_SIZE|--VERBOSITY|--QUIET|--VALIDATION_STRINGENCY|\
        --COMPRESSION_LEVEL|--MAX_RECORDS_IN_RAM|--CREATE_INDEX|--CREATE_MD5_FILE|\
        --USE_JDK_DEFLATER|--USE_JDK_INFLATER|--ASSUME_SORTED|--TMP_DIR|\
        --MOLECULAR_IDENTIFIER_TAG|--showHidden|--help|--arguments_file)
            return 0 ;;
        --FLOW_*) return 0 ;;
        *) return 1 ;;
    esac
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --version)        WANT_VERSION=1; shift ;;
        --INPUT|-I)               INPUT="$2";              shift 2 ;;
        --OUTPUT|-O)              OUTPUT="$2";             shift 2 ;;
        --METRICS_FILE|-M)        METRICS="$2";            shift 2 ;;
        --ASSUME_SORT_ORDER)      ASSUME_SORT_ORDER="$2";  shift 2 ;;
        --REFERENCE_SEQUENCE|-R)  REFERENCE_SEQUENCE="$2"; shift 2 ;;
        --BARCODE_TAG)            BARCODE_TAG="$2";        shift 2 ;;
        --READ_ONE_BARCODE_TAG)   READ_ONE_BARCODE_TAG="$2"; shift 2 ;;
        --READ_TWO_BARCODE_TAG)   READ_TWO_BARCODE_TAG="$2"; shift 2 ;;
        --REMOVE_DUPLICATES)      REMOVE_DUPLICATES="$2";  shift 2 ;;
        *)
            if swallow_arg "$1"; then
                # Most Picard args are `--KEY VALUE`; some (--FLOW_MODE) are
                # boolean-with-value too. Always consume the next token if it
                # exists and doesn't start with `--`.
                if [[ $# -ge 2 && "$2" != --* ]]; then shift 2; else shift; fi
            else
                die "unrecognized MarkDuplicates arg: $1"
            fi
            ;;
    esac
done

if [[ $WANT_VERSION -eq 1 ]]; then
    # nf-core captures: `picard MarkDuplicates --version 2>&1 | sed -n 's/.*Version://p'`
    echo "Version:${PICARD_VERSION_STRING}" >&2
    exit 0
fi

[[ -n "$INPUT"   ]] || die "--INPUT is required"
[[ -n "$OUTPUT"  ]] || die "--OUTPUT is required"
[[ -n "$METRICS" ]] || die "--METRICS_FILE is required"

# ---- Pre-processing: CRAM → BAM via samtools ----
TMPFILES=()
cleanup() { for f in "${TMPFILES[@]:-}"; do [[ -n "$f" ]] && rm -f "$f"; done; }
trap cleanup EXIT

WEA_INPUT="$INPUT"
case "$INPUT" in
    *.cram|*.CRAM)
        [[ -n "$REFERENCE_SEQUENCE" ]] || die "CRAM input requires --REFERENCE_SEQUENCE"
        TMP_BAM=$(mktemp --suffix=.bam)
        TMPFILES+=("$TMP_BAM")
        samtools view -b --threads 4 --reference "$REFERENCE_SEQUENCE" \
            -o "$TMP_BAM" "$INPUT" \
            || die "CRAM→BAM conversion failed"
        WEA_INPUT="$TMP_BAM"
        ;;
esac

# If input BAM is not coordinate-sorted, samtools-sort into a tmp BAM.
# Picard handles unsorted/queryname input natively; wea requires coordinate.
HDR_SO=$(samtools view -H "$WEA_INPUT" 2>/dev/null | awk -F'\t' '/^@HD/{for(i=2;i<=NF;i++)if($i~/^SO:/){sub("SO:","",$i);print $i;exit}}')
if [[ "$HDR_SO" != "coordinate" ]]; then
    TMP_SORTED=$(mktemp --suffix=.bam)
    TMPFILES+=("$TMP_SORTED")
    samtools sort --threads 4 -o "$TMP_SORTED" "$WEA_INPUT" \
        || die "pre-sort (samtools sort) failed"
    WEA_INPUT="$TMP_SORTED"
fi

# CRAM output → BAM intermediate, then convert back.
WEA_OUTPUT="$OUTPUT"
NEED_OUTPUT_CRAM=0
case "$OUTPUT" in
    *.cram|*.CRAM)
        NEED_OUTPUT_CRAM=1
        TMP_OUT_BAM=$(mktemp --suffix=.bam)
        TMPFILES+=("$TMP_OUT_BAM")
        WEA_OUTPUT="$TMP_OUT_BAM"
        ;;
esac

# ---- Build wea command ----
WEA_ARGS=( "$WEA_INPUT" -o "$WEA_OUTPUT" -M "$METRICS" )
if [[ -n "$ASSUME_SORT_ORDER" ]]; then
    case "$ASSUME_SORT_ORDER" in
        coordinate) WEA_ARGS+=( --assume-sort-order coordinate ) ;;
        queryname)  : ;;  # Picard hint; wea reads @HD SO from the BAM header
        unsorted)   : ;;
        *) die "unsupported --ASSUME_SORT_ORDER: $ASSUME_SORT_ORDER" ;;
    esac
fi
[[ -n "$BARCODE_TAG"          ]] && WEA_ARGS+=( --barcode-tag "$BARCODE_TAG" )
[[ -n "$READ_ONE_BARCODE_TAG" ]] && WEA_ARGS+=( --read-one-barcode-tag "$READ_ONE_BARCODE_TAG" )
[[ -n "$READ_TWO_BARCODE_TAG" ]] && WEA_ARGS+=( --read-two-barcode-tag "$READ_TWO_BARCODE_TAG" )
[[ "$REMOVE_DUPLICATES" == "true" ]] && WEA_ARGS+=( --remove-duplicates )

"$WEA_BIN" "${WEA_ARGS[@]}"
WEA_RC=$?
[[ $WEA_RC -eq 0 ]] || exit $WEA_RC

# ---- Post: convert BAM→CRAM if requested ----
if [[ $NEED_OUTPUT_CRAM -eq 1 ]]; then
    samtools view -C --threads 4 --reference "$REFERENCE_SEQUENCE" \
        -o "$OUTPUT" "$WEA_OUTPUT" \
        || die "BAM→CRAM conversion failed"
fi

# ---- Rewrite metrics line 2 with Picard-templated provenance ----
# nf-core snapshot fingerprint: lines [0..2] of the metrics file. Wea writes
# line0/line2 already in the correct shape (`## htsjdk.samtools.metrics.StringHeader`).
# Only line 1 (`# markdup-wea ...`) needs to be replaced with the Picard
# `# MarkDuplicates --INPUT ... <full default-flags expansion> ...` template.

INPUT_BN=$(basename -- "$INPUT")
OUTPUT_BN=$(basename -- "$OUTPUT")
METRICS_BN=$(basename -- "$METRICS")
SO_PRINT="${ASSUME_SORT_ORDER:-queryname}"   # Picard default-print value
REF_FRAG=""
[[ -n "$REFERENCE_SEQUENCE" ]] && REF_FRAG=" --REFERENCE_SEQUENCE $(basename -- "$REFERENCE_SEQUENCE")"

LINE2="# MarkDuplicates --INPUT ${INPUT_BN} --OUTPUT ${OUTPUT_BN} --METRICS_FILE ${METRICS_BN} --ASSUME_SORT_ORDER ${SO_PRINT}${REF_FRAG} --MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP 50000 --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 --SORTING_COLLECTION_SIZE_RATIO 0.25 --TAG_DUPLICATE_SET_MEMBERS false --REMOVE_SEQUENCING_DUPLICATES false --TAGGING_POLICY DontTag --CLEAR_DT true --DUPLEX_UMI false --FLOW_MODE false --FLOW_DUP_STRATEGY FLOW_QUALITY_SUM_STRATEGY --FLOW_USE_END_IN_UNPAIRED_READS false --FLOW_USE_UNPAIRED_CLIPPED_END false --FLOW_UNPAIRED_END_UNCERTAINTY 0 --FLOW_UNPAIRED_START_UNCERTAINTY 0 --FLOW_SKIP_FIRST_N_FLOWS 0 --FLOW_Q_IS_KNOWN_END false --FLOW_EFFECTIVE_QUALITY_THRESHOLD 15 --ADD_PG_TAG_TO_READS true --REMOVE_DUPLICATES ${REMOVE_DUPLICATES} --ASSUME_SORTED false --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --PROGRAM_RECORD_ID MarkDuplicates --PROGRAM_GROUP_NAME MarkDuplicates --READ_NAME_REGEX <optimized capture of last three ':' separated fields as numeric values> --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 5 --MAX_RECORDS_IN_RAM 500000 --CREATE_INDEX false --CREATE_MD5_FILE false --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false"

TMP_METRICS=$(mktemp)
{
    head -n 1 "$METRICS"
    printf '%s\n' "$LINE2"
    tail -n +3 "$METRICS"
} > "$TMP_METRICS" && mv "$TMP_METRICS" "$METRICS"

exit 0
