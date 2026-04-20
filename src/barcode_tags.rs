//! Runtime configuration for BARCODE_TAG / READ_ONE_BARCODE_TAG /
//! READ_TWO_BARCODE_TAG.
//!
//! `BarcodeTags` is borrowed-only — the 2-byte tag identifiers are owned by
//! `main.rs` for the duration of the run. Default value = all-None = feature
//! off, in which case every extracted hash resolves to 0 and dedup behavior
//! is byte-identical to the pre-Track-A layout (A.1 parity contract).

#[derive(Debug, Clone, Copy, Default)]
pub struct BarcodeTags<'a> {
    /// `BARCODE_TAG` — whole-pair barcode/UMI. Typical value: RX.
    pub barcode: Option<&'a [u8; 2]>,
    /// `READ_ONE_BARCODE_TAG` — barcode on the firstOfPair record.
    pub read_one: Option<&'a [u8; 2]>,
    /// `READ_TWO_BARCODE_TAG` — barcode on the !firstOfPair record.
    pub read_two: Option<&'a [u8; 2]>,
    /// `MOLECULAR_IDENTIFIER_TAG` — output-only SAM aux tag on each written
    /// record. Emitted as `{contig}:{pos_1based}/` per Picard `UmiUtil.java`.
    /// Per Picard `MarkDuplicates.java:405`, MI emission is gated on
    /// `BARCODE_TAG` being set; setting `mi` without `barcode` is a CLI
    /// validation error (mirrors Picard's `customCommandLineValidation`).
    pub mi: Option<&'a [u8; 2]>,
}

impl BarcodeTags<'_> {
    pub fn any_enabled(&self) -> bool {
        self.barcode.is_some()
            || self.read_one.is_some()
            || self.read_two.is_some()
            || self.mi.is_some()
    }
}

/// Compute the MI tag value for a single record.
///
/// Matches Picard `UmiUtil.setMolecularIdentifier` with `assignedUmi=""` and
/// `duplexUmis=false`. Format: `{contig}:{pos_1based}/`. For unmapped records
/// Picard writes the literal string `"null"` for the contig (Java
/// `StringBuilder.append(null)` semantics); we reproduce that byte-for-byte.
/// `pos_1based` is `record.alignment_start()` on reverse-strand reads and
/// `record.mate_alignment_start()` on forward-strand reads — Picard's
/// canonical "pair's representative position" convention.
pub fn format_mi_value(contig: Option<&str>, is_reverse: bool, pos_1based: i64) -> String {
    let contig_str = contig.unwrap_or("null");
    format!(
        "{}:{}/",
        contig_str,
        pos_1based_or_zero(is_reverse, pos_1based)
    )
}

fn pos_1based_or_zero(_is_reverse: bool, pos_1based: i64) -> i64 {
    if pos_1based <= 0 {
        0
    } else {
        pos_1based
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mi_format_forward_paired() {
        // Forward strand → uses mate_alignment_start.
        assert_eq!(format_mi_value(Some("chr1"), false, 1200), "chr1:1200/");
    }

    #[test]
    fn mi_format_reverse_paired() {
        // Reverse strand → uses record's own alignment_start.
        assert_eq!(format_mi_value(Some("chrX"), true, 9876), "chrX:9876/");
    }

    #[test]
    fn mi_format_unmapped_contig_is_literal_null() {
        // Picard quirk: StringBuilder.append(null) emits "null".
        assert_eq!(format_mi_value(None, false, 0), "null:0/");
    }

    #[test]
    fn mi_format_missing_pos_emits_zero() {
        assert_eq!(format_mi_value(Some("chr1"), false, -1), "chr1:0/");
    }

    #[test]
    fn any_enabled_covers_mi() {
        let mi = *b"MI";
        let t = BarcodeTags {
            mi: Some(&mi),
            ..Default::default()
        };
        assert!(t.any_enabled());
    }
}
