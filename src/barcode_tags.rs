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
}

impl BarcodeTags<'_> {
    pub fn any_enabled(&self) -> bool {
        self.barcode.is_some() || self.read_one.is_some() || self.read_two.is_some()
    }
}
