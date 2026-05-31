# Working-session prompt: roll `lagoon.hoon` over to the new SoftBLAS API

Paste everything below the line into a fresh Hoon-focused session (one that has
the Urbit/`lagoon.hoon` repo and the Vere jet sources on hand). It coordinates
the Hoon side with two planned SoftBLAS C-side changes so the jet boundary stays
bit-exact and deterministic.

---

You are working on `lagoon.hoon` (Urbit's linear-algebra library) and its Vere
jet bridge to **SoftBLAS** (the software-float BLAS that backs Lagoon's jets).
Two changes are landing on the SoftBLAS (C) side. Your job is to make the
matching Hoon + jet-glue changes so nothing drifts. Treat **determinism** as the
top constraint: the jetted result must equal the Hoon result bit-for-bit, on
every platform, forever.

## Change 1 — rounding mode becomes a typed enum (expected ABI-preserving)

SoftBLAS currently threads the rounding mode through every routine as
`uint_fast8_t rndMode` carrying an **ASCII character code**:

| code | char | meaning                              | SoftFloat mode             |
|------|------|--------------------------------------|----------------------------|
| 0x6e | `'n'`| round to nearest, ties to even       | `round_near_even`          |
| 0x7a | `'z'`| round toward zero (truncate)         | `round_minMag`             |
| 0x75 | `'u'`| round toward +∞ (up)                 | `round_max`                |
| 0x64 | `'d'`| round toward −∞ (down)               | `round_min`                |
| 0x61 | `'a'`| round to nearest, ties **away**      | `round_near_maxMag`        |

It is being changed to a `typedef enum { SB_ROUND_NEAR_EVEN = 'n', ... }`. The
enum's **underlying integer values are the same ASCII codes**, so the wire/ABI
byte the jet passes does **not** change — this is a C-side legibility change.

Your tasks:
1. Find every place in the jet glue (Vere `ext/softblas` / the jet C files) and
   in `lagoon.hoon` where the rounding mode is produced or passed across the
   boundary. Confirm it is encoded as the **byte value of the ASCII code** (e.g.
   a `@t`/cord whose first byte is `'n'`, or a literal `0x6e`), NOT an ordinal
   index `0/1/2/3/4`. If it is already the char byte, **no Hoon change is needed**
   — record that and move on.
2. If anything passes an *ordinal index* instead of the char code, that is a
   latent mismatch: fix it to pass the char byte, and add a comment pinning the
   mapping above.
3. Add/confirm a Hoon-level test that exercises at least two rounding modes
   (including one with a negative operand, the only thing that distinguishes
   `z`/`d`/`u`/`a`) and asserts the jetted output equals the Hoon output.

## Change 2 — `nrm2_B` becomes the canonical `nrm2`; the old one is retired

SoftBLAS added `snrm2_B`/`dnrm2_B`/`hnrm2_B`/`qnrm2_B`: E. Anderson's safe-scaling
2-norm (ACM TOMS Alg. 978, 2017). It is ~1.7–2× faster than the classic
running-scale `nrm2` **and strictly more correct**: the classic `snrm2`/`dnrm2`/
`hnrm2` return **NaN** for subnormal inputs (`0/0`) and for equal near-overflow
pairs (`Inf/Inf`); the `_B` variants return the correct value.

The plan: **rename the `_B` routines to the canonical names** (`snrm2_B` →
`snrm2`, etc.) and **delete the classic running-scale versions**. The jet
registration name stays `*nrm2`, so the jet hint/dampening does not change — but
**the numeric output changes for the overflow/underflow/subnormal edge cases.**

Your tasks:
1. Locate the Lagoon arm(s) that implement / jet `nrm2` (and any per-precision
   variants). Identify the jet hint and the C symbol it binds to.
2. The Hoon (unjetted) `nrm2` is the source of truth for determinism. **Decide
   and make consistent:** either update the Hoon `nrm2` to the safe-scaling
   algorithm so it matches the new C exactly, OR confirm the Hoon already
   produces the safe-scaling result. The jet and the Hoon MUST agree bit-for-bit,
   including on subnormal and near-overflow inputs. This is the crux — do not
   skip it.
3. Find every pinned expected-value test for `nrm2` in the Lagoon test suite.
   Any that encode the **old buggy** results (NaN for subnormals / equal-large
   pairs) must be updated to the new correct results. Add explicit regression
   cases for: a subnormal vector, an equal near-max-magnitude pair, and a normal
   vector (unchanged).
4. Coordinate the version bump / changelog note on the Hoon side so a downstream
   pinning the old NaN behavior is warned.

## Acceptance criteria

- Jetted vs unjetted `nrm2` agree bit-for-bit on: normal, subnormal,
  equal-near-overflow, zero-length, and single-element inputs, across all
  precisions Lagoon exposes.
- Rounding-mode round-trip (Hoon → jet → result) verified for ≥2 modes incl. a
  negative operand.
- No Lagoon test still encodes the old `nrm2` NaN-on-subnormal behavior.
- A short note in the Lagoon changelog describing the `nrm2` numeric change and
  the (ABI-stable) rounding-mode enum.

## Context you may want

- SoftBLAS rounding mapping lives in `_set_rounding` (`include/softblas.h`).
- The safe-scaling algorithm and its thresholds are in `src/blas/level1/*nrm2_B.c`
  (constants `tsml`/`tbig`/`ssml`/`sbig` are 2^k from the float exponent range).
- Determinism is the whole point of SoftBLAS: every op is canonicalized
  (`nan_unify_*`) so NaN payloads are stable. Preserve that on the Hoon side.
