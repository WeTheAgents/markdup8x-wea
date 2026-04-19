//! Picard-compatible barcode hash primitives.
//!
//! Replicates Java semantics (i32 wrapping arithmetic, asymmetric
//! missing-tag fallbacks) byte-for-byte. See `docs/umi_semantics.md`.

use std::fmt;

#[derive(Debug, PartialEq, Eq)]
pub enum BarcodeError {
    InvalidChar(u8),
}

impl fmt::Display for BarcodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BarcodeError::InvalidChar(c) => write!(
                f,
                "UMI contains illegal character 0x{:02x} (allowed: ATCGNatcgn-)",
                c
            ),
        }
    }
}

impl std::error::Error for BarcodeError {}

/// Java `String.hashCode()` polynomial-31 on raw bytes.
/// UMI alphabet is ASCII, so byte-wise iteration is exact.
pub fn java_string_hashcode(s: &[u8]) -> i32 {
    let mut h: i32 = 0;
    for &c in s {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h
}

/// Java `Objects.hash(umi)` — `31 + (umi.map_or(0, java_string_hashcode))`.
/// Asymmetric fallback: `None` → 31 (NOT 0). Used by BARCODE_TAG.
pub fn picard_barcode_hash(umi: Option<&[u8]>) -> i32 {
    let inner = umi.map_or(0, java_string_hashcode);
    31i32.wrapping_add(inner)
}

/// Java `String.hashCode()` with null → 0.
/// Used by READ_ONE_BARCODE_TAG / READ_TWO_BARCODE_TAG.
pub fn read_barcode_value(tag_attr: Option<&[u8]>) -> i32 {
    tag_attr.map_or(0, java_string_hashcode)
}

/// Picard's `UmiUtil.ALLOWED_UMI` regex: `^[ATCGNatcgn-]*$`.
pub fn validate_umi(s: &[u8]) -> Result<(), BarcodeError> {
    for &c in s {
        match c {
            b'A' | b'T' | b'C' | b'G' | b'N' | b'a' | b't' | b'c' | b'g' | b'n' | b'-' => {}
            other => return Err(BarcodeError::InvalidChar(other)),
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values derived from openjdk-21 on Hetzner (2026-04-19),
    // captured at /mnt/HC_Volume_105344878/tmp/umi-probe/java_hashcodes.txt.
    //
    // string        String.hashCode      Objects.hash
    // ''                          0                31
    // 'A'                        65                96
    // 'AAAAAAAA'         1094643840        1094643871
    // 'AGCTAGCT'        -2095457618       -2095457587
    // 'agctagct'         1416882862        1416882893
    // 'AGCT-TGCA'        -552756473        -552756442
    // 'AGCT-AGCT'        -553322483        -553322452
    // 'AGCTTGCA'        -2094891608       -2094891577
    // (null)                    N/A                31

    #[test]
    fn java_hashcode_empty() {
        assert_eq!(java_string_hashcode(b""), 0);
    }

    #[test]
    fn java_hashcode_single_char() {
        assert_eq!(java_string_hashcode(b"A"), 65);
    }

    #[test]
    fn java_hashcode_eight_a() {
        assert_eq!(java_string_hashcode(b"AAAAAAAA"), 1094643840);
    }

    #[test]
    fn java_hashcode_negative_wrap() {
        // Confirms i32 two's-complement wrap matches Java's int.
        assert_eq!(java_string_hashcode(b"AGCTAGCT"), -2095457618);
    }

    #[test]
    fn java_hashcode_case_distinct() {
        let upper = java_string_hashcode(b"AGCTAGCT");
        let lower = java_string_hashcode(b"agctagct");
        assert_eq!(upper, -2095457618);
        assert_eq!(lower, 1416882862);
        assert_ne!(upper, lower);
    }

    #[test]
    fn java_hashcode_dash_kept_in_input() {
        let with_dash = java_string_hashcode(b"AGCT-TGCA");
        let no_dash = java_string_hashcode(b"AGCTTGCA");
        assert_eq!(with_dash, -552756473);
        assert_eq!(no_dash, -2094891608);
        assert_ne!(with_dash, no_dash);
    }

    #[test]
    fn java_hashcode_dash_distinct_halves() {
        // AGCT-TGCA ≠ AGCT-AGCT — second half matters.
        assert_eq!(java_string_hashcode(b"AGCT-TGCA"), -552756473);
        assert_eq!(java_string_hashcode(b"AGCT-AGCT"), -553322483);
    }

    #[test]
    fn picard_barcode_hash_none_is_31() {
        assert_eq!(picard_barcode_hash(None), 31);
    }

    #[test]
    fn picard_barcode_hash_empty_is_31() {
        // `31 + "".hashCode() == 31 + 0 == 31` — same as None.
        assert_eq!(picard_barcode_hash(Some(b"")), 31);
    }

    #[test]
    fn picard_barcode_hash_eight_a() {
        assert_eq!(picard_barcode_hash(Some(b"AAAAAAAA")), 1094643871);
    }

    #[test]
    fn picard_barcode_hash_negative() {
        assert_eq!(picard_barcode_hash(Some(b"AGCTAGCT")), -2095457587);
    }

    #[test]
    fn read_barcode_value_none_is_0() {
        // Asymmetry vs picard_barcode_hash: READ_ONE/TWO null → 0, not 31.
        assert_eq!(read_barcode_value(None), 0);
    }

    #[test]
    fn read_barcode_value_empty_is_0() {
        assert_eq!(read_barcode_value(Some(b"")), 0);
    }

    #[test]
    fn read_barcode_value_eight_a() {
        // Bare String.hashCode, no +31 seed.
        assert_eq!(read_barcode_value(Some(b"AAAAAAAA")), 1094643840);
    }

    #[test]
    fn missing_tag_asymmetry() {
        // The single highest-pitfall finding from A.0:
        // BARCODE_TAG missing → 31, READ_ONE/TWO missing → 0.
        assert_ne!(picard_barcode_hash(None), read_barcode_value(None));
        assert_eq!(picard_barcode_hash(None), 31);
        assert_eq!(read_barcode_value(None), 0);
    }

    #[test]
    fn validate_umi_accepts_alphabet() {
        assert!(validate_umi(b"ATCGNatcgn-").is_ok());
    }

    #[test]
    fn validate_umi_accepts_empty() {
        // Picard's `^[ATCGNatcgn-]*$` regex permits zero-length.
        assert!(validate_umi(b"").is_ok());
    }

    #[test]
    fn validate_umi_rejects_x() {
        assert_eq!(validate_umi(b"ATCGX"), Err(BarcodeError::InvalidChar(b'X')));
    }

    #[test]
    fn validate_umi_rejects_digit() {
        assert_eq!(validate_umi(b"ATC1G"), Err(BarcodeError::InvalidChar(b'1')));
    }

    #[test]
    fn validate_umi_short_circuits_at_first_bad_char() {
        // Reports the first offender, not subsequent ones.
        assert_eq!(validate_umi(b"AXY"), Err(BarcodeError::InvalidChar(b'X')));
    }
}
