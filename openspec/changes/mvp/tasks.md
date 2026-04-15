# Tasks

## Guiding principle (revised 2026-04-15)

**Goal is NOT "better than Picard". Goal is "byte-for-byte identical FLAG 0x400 output + Picard-compatible metrics, without JVM".**

Burden of evidence for replacing a critical bioinformatics tool is enormous. The nf-core community is conservative (biologists/medics, not programmers); they care that the tool **works as expected**, not that it's faster. The Sambamba cautionary tale (silent R1/R2 reorder broke downstream methylation extraction, FelixKrueger/Bismark#170) is the failure mode we must avoid.

**Scope discipline:**
- Correctness > speed. Threading/optimization is post-MVP.
- Every deviation from Picard must be *documented*, not discovered.
- Deliverable for maintainers: "identical outputs on N BAMs + deviation table + known edge-case coverage". Let them decide if they want it.

**Out of scope for MVP (deferred to post-acceptance):**
- Threading optimization (Phase 4 — `-@` flag accepted but no-op)
- Human genome validation at full scale (Phase 5b — yeast + 1 human sample is enough for maintainer review)
- crates.io publishing, GitHub Release automation (Phase 10 — maintainers' call)
- Optical duplicate detection (not used by nf-core/rnaseq; will be documented as unsupported)

## 1. Project Skeleton (Phase 1) — DONE

- [x] 1.1 Initialize Cargo.toml (noodles chosen over rust-htslib: no C dep, simpler static build)
- [x] 1.2 Configure release profile: lto=true, codegen-units=1, strip=true, opt-level=3
- [x] 1.3 Create src/main.rs with clap CLI: INPUT, -o, -M, -@, --remove-duplicates, --assume-sort-order
- [x] 1.4 Create src/lib.rs with public module re-exports
- [x] 1.5 Create LICENSE (MIT)
- [ ] 1.6 Create README.md (project description, usage, build instructions) — deferred until Phase 5a passes
- [ ] 1.7 Create .github/workflows/ci.yml (build, test, clippy, fmt on Linux/macOS/Windows)

### GATE: `cargo build --release` compiles. CI green. — PARTIAL (build OK, CI not wired)

## 2. Core Primitives (Phase 1) — DONE

- [x] 2.1–2.8 position.rs, scoring.rs, io.rs + unit tests (45/45 green)

### GATE: All unit tests pass. — PASSED

## 3. Paired-End Duplicate Detection (Phase 2) — DONE

- [x] 3.1–3.17 pending_mates, groups, scan, markdup, library detection, flag clearing, --remove-duplicates

### GATE: All unit tests pass. — PASSED

## 4. DupSet Strategies (Phase 2) — DONE

- [x] 4.1–4.5 BitVec + HashSet strategies, identical-output test

### GATE: Both strategies compile and pass identical-output test. — PASSED

## 5. Picard Edge-Case Research (NEW — Phase 2.5)

**Blocker for Phase 6+. Sambamba-proofing.**

- [ ] 5.1 Receive `docs/picard-edge-cases.md` from research team (external deep-research task)
- [ ] 5.2 Cross-reference each edge case against current `openspec/changes/mvp/specs/` — mark covered / uncovered / divergent
- [ ] 5.3 For each **uncovered** case: decide cover-in-MVP vs document-as-deviation
- [ ] 5.4 For each **divergent** case: fix to match Picard OR document rationale in deviation table
- [ ] 5.5 Create `docs/deviations.md` — explicit table of all known behavioral differences vs Picard (even if empty, its existence signals we thought about it)
- [ ] 5.6 Backfill unit tests for every edge case we claim to cover

### GATE: Every high-risk edge case from research doc is either covered by test or listed in deviations.md.

## 6. Synthetic BAM Integration Tests (Phase 3)

**Goal: exercise edge cases identified in Phase 5 on controllable inputs before real BAMs.**

- [ ] 6.1 Create test fixture: two identical pairs at same position → one flagged
- [ ] 6.2 Create test fixture: pair with higher quality wins over lower
- [ ] 6.3 Create test fixture: quality tie → tie-break matches Picard rule (from edge-case doc)
- [ ] 6.4 Create test fixture: single-end reads mixed with paired
- [ ] 6.5 Create test fixture: mate on different chromosome (chimeric pair) — verify both flagged per Picard rule
- [ ] 6.6 Create test fixture: mate unmapped → single-end treatment
- [ ] 6.7 Create test fixture: supplementary/secondary of dup primary — verify Picard-matching behavior
- [ ] 6.8 Create test fixture: pre-existing FLAG 0x400 → cleared and re-evaluated
- [ ] 6.9 Create test fixture: orphan read (no mate in BAM) → pass through with stderr warning
- [ ] 6.10 Create test fixture: multiple libraries → separate grouping
- [ ] 6.11 Create test fixture: runtime sort-order violation → error exit
- [ ] 6.12 Create test fixture: spliced reads (N in CIGAR) — RNA-seq critical
- [ ] 6.13 Create test fixture: missing @RG / missing LB → Picard fallback behavior
- [ ] 6.14 Wire all fixtures into tests/integration.rs with FLAG-level assertions
- [ ] 6.15 Every fixture runs against BOTH dupset strategies

### GATE: All synthetic scenarios pass on both DupSet strategies. Fixture inputs/outputs committed to repo for reproducibility.

## 7. Picard-Compatible Metrics (Phase 3)

- [x] 7.1–7.6 metrics.rs skeleton (already done in Phase 2)
- [ ] 7.7 Implement metrics file writer with exact Picard header format (verified byte-exact against Picard template from edge-case research)
- [ ] 7.8 Per-library rows when multiple libraries present
- [ ] 7.9 Integration test: metrics file diff vs Picard output on synthetic BAM — zero-byte diff target
- [ ] 7.10 MultiQC parse test: `multiqc .` consumes our metrics file without warnings

### GATE: Metrics file byte-identical to Picard on at least one synthetic input. MultiQC parseable.

## 8. Record-Reuse + Minimal Logging (lightweight — no threading)

- [ ] 8.1 Ensure `reader.read(&mut record)` reuse pattern in both passes (no per-record allocation)
- [ ] 8.2 Pre-size FxHashMap from BAM file size estimate
- [ ] 8.3 Log: record count, pending_mates peak, groups resolved, wall time per pass
- [ ] 8.4 `-@` flag: accept, no-op for MVP, log "threading not implemented"

### GATE: Single-threaded run completes yeast BAM without OOM. Logs are useful for debugging.

## 9. Yeast Validation (Phase 5a — CRITICAL GATE)

**This is THE correctness gate. Maintainer decision hinges on this.**

- [ ] 9.1 Acquire 5 yeast test BAMs (run `nf-core/rnaseq -profile test,docker --skip_markduplicates`)
- [ ] 9.2 Document BAM acquisition in `tests/README.md`
- [ ] 9.3 For each sample, assert: **`samtools view -f 1024 <our>.bam | cut -f1 | sort -u`** is **byte-identical** to the same pipeline on Picard output. Not just counts — every single flagged QNAME.
- [ ] 9.4 Dup count ground truth:
  - WT_REP1: 36,810 / 180,342
  - WT_REP2: 11,688 / 90,962
  - RAP1_UNINDUCED_REP2: 78,929 / 98,201
  - RAP1_UNINDUCED_REP1: 36,294 / 48,977
  - RAP1_IAA_30M_REP1: 11,094 / 48,347
- [ ] 9.5 If any QNAME differs: isolate by `diff`, trace to edge case, fix OR add to deviations.md with rationale
- [ ] 9.6 Metrics file diff: our `.metrics.txt` vs Picard `.metrics.txt` — zero-byte diff target (or documented fields)
- [ ] 9.7 MultiQC: `multiqc .` on our outputs, verify duplication plot renders identically
- [ ] 9.8 stdin pipe test: `cat yeast.bam | markdup-wea - -o out.bam` produces same result
- [ ] 9.9 `--remove-duplicates` produces correct kept-record count
- [ ] 9.10 Static binary size check: `strip`ped release binary < 10MB

### GATE: All 5 yeast samples — byte-identical flagged QNAME set + Picard-parseable metrics. Any deviation explicitly in deviations.md. This is the go/no-go point.

## 10. Single Human Sample Smoke Test (Phase 5b-lite)

**Not full benchmark — just proof it scales. Maintainers can run the full suite if they want.**

- [ ] 10.1 One ENCODE human RNA-seq sample (~200M reads)
- [ ] 10.2 Assert flagged-QNAME set identical to Picard
- [ ] 10.3 Record wall time + peak RSS for reference (not a claim, a data point)
- [ ] 10.4 Log pending_mates peak for sizing confidence

### GATE: One human sample matches Picard. No OOM, no crash.

## 11. Handoff Package (for maintainers)

- [ ] 11.1 README: what it is, what it isn't, known deviations, how to build, how to run
- [ ] 11.2 `docs/deviations.md` — final version, all known differences vs Picard
- [ ] 11.3 `docs/picard-edge-cases.md` — research doc + our coverage status
- [ ] 11.4 `BENCHMARK.md` — observed numbers from yeast + 1 human, with methodology. Framed as "data point", not "we are better"
- [ ] 11.5 Draft nf-core Slack/GitHub post explaining scope, non-goals, and what we're asking maintainers to decide

### GATE: A maintainer can read the handoff package in 20 minutes and decide yes/no without running anything themselves.

## Deferred (not MVP)

- Threading (`-@`) actual implementation — post-acceptance optimization
- Full 8-sample human benchmark — maintainers' call
- crates.io / GitHub Release automation — maintainers' call
- Optical duplicate detection — documented as unsupported
- Windows binary — Linux x86_64 is enough for nf-core; macOS/Windows post-acceptance
