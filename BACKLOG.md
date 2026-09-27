# SoftBLAS bug & test backlog

Consolidated from the five-reviewer audit, the Gnome/Angel test-coverage pass,
the Urbit-headers analysis, and the 2026-05-30 re-audit. Reconciled against
`master` after PR #25.

Severity: 🔴 critical (memory-unsafe / silently wrong) · 🟠 high · 🟡 medium · ⚪ low.

---

## Done

**Correctness (Phase 2)**
- [x] B1 🔴 `*swap`/`csrot` signed-stride OOB — params now `int64_t`, ASan-clean.
- [x] B2 🟠 `qaxpy` NaN canonicalization — per-element `nan_unify_q` in the loop.
- [x] B3 🟠 `FNEG` macro UB — width-correct unsigned masks; `c128_conj` compiles.
- [x] B4 🟡 `i*amax` N==1 returns 0 (all four precisions).
- [x] 🔴 `*rotmg` infinite loop on NaN/Inf `D1`/`D2` (srotmg/drotmg/hrotmg) —
  guarded the `_ge` rescale loops; `qrotmg` born with the guard. (+regression tests)
- [x] 🔴 7 ghost declarations implemented — `qrot`/`qrotg`/`qrotm`/`qrotmg` +
  `cscal`/`cswap`/`icamax`; added to Makefile; fixed `qrot` prototype; +tests.
- [x] 🟡 `asum` (s/d/h/q) `incX < 1` guard; `gemv` (s/d/h/q) `incX<1||incY<1`
  guard (signed `incX` OOB). `i*amax` was already guarded.
- [x] 🟡 `test_hgemm` off-by-one (`i < 4` over `M*P` arrays) → `i < M*P`.

**Architecture (Phase 3)**
- [x] A1 ⚪ `qasum` heap leak removed — stack accumulator; compiled core is
  allocation-free.
- [x] A2 ⚪ Vendored Urbit headers + `#ifdef VERE` block deleted (`include/` is
  now just `softblas.h`).
- [x] A3 ⚪ `svec/dvec/hvec/qvec` moved to `tests/blas/include/test.h`.

**Test suite honesty (Phase 1) + hygiene (Phase 5)**
- [x] T1 16 test misregistrations in `test_all.c` repointed.
- [x] T2 CI registration-honesty lint (`tools/check_test_registration.sh`, runs in CI).
- [x] T3 README rounding divergences (global→per-call, `'a'`=ties-away).
- [x] H2 `benchmarking/` kept (revived for the `snrm2`/`snrm2_B` issue-#3
  comparison; the stack-VLA footgun in `benchmark_snrm2.c` was fixed).
- [x] H3 macOS build path — OS-aware Makefile + `tools/build_softfloat_darwin.sh`
  + `macos-latest` CI job (required, green).
- [x] macOS/arm64 directional-rounding SIGBUS — was a Darwin build-config bug
  (`platform.h` set `THREAD_LOCAL _Thread_local`, routing SoftFloat state through
  macOS dyld TLV → BUS). Fixed by using the GCC `platform.h` (plain globals).

**Complex Level 2/3 (Phase 2 follow-on, PR #32)**
- [x] Unconjugated dot product `c/i/z/v`-`dotu`; complex `gemv` (`i/c/z/v`, both
  row/col layouts + the `'C'` conjugate-transpose op); complex `gemm`
  (`i/c/z/v`, incl. `'C'`). Mirrors the established stateless-`rndMode` /
  `cNN_*` / `nan_unify_*` conventions; quad (`v`) uses the value-based `c128_*`
  ops + pointer-based `nan_unify_v` like the rest of the `v` family. Closes the
  README "Complex Level 2/3" roadmap item.
- [x] Coverage: per-precision canonical vectors, `'C'` conjugate transpose, both
  gemv layouts, general α/β, `ldb≠ldc` padding, negative stride, 3×3 / 4×4, and a
  25×25 strided (`incX=incY=5`) `cgemv` whose A/x come from a deterministic
  non-periodic formula with an independently-computed expected `y` and sentinel
  values in the skipped stride slots. NaN canonicalization through `cgemm`/`cgemv`
  + a complex rounding-mode test (all four modes through the real *and* imaginary
  accumulators).

**Rotation-routine fixes from the RustFloat differential fuzzer (issues #33-#35)**
- [x] 🔴 #33 `?rotmg` hang on a negative `D1` (all four precisions). The first
  rescale loop's `d1 != 0 && d1 <= rgamsq` guard never goes false for a negative
  `d1`: the `gamsq` scaling keeps it negative down to `-inf`, which is still
  `<= rgamsq`. Reachable from valid-looking inputs — a negative `d2` whose
  `q2 = (d2*y1)*y1` underflows to `-0` slips past the `q2 < 0` singularity test,
  and the `flag = 1` branch swaps it into `d1`. Guard is now `0 < d1`, matching
  the Rust port (sigilante/RustFloat). (+4 regression tests; each hangs on the
  pre-fix sources.)
- [x] 🟠 #34 `?rotg`/`?rotmg` NaN canonicalization — the only routines that
  wrote their outputs without `nan_unify_*`, so a NaN came back carrying
  SoftFloat's propagated payload and sign instead of `SINGNAN`/`DOUBNAN`/
  `HALFNAN`/`QUADNAN`. Now unified at the stores, `?rotm`-style; `rotmg`
  unifies only the values it actually writes (the `p2 == 0` early return still
  leaves `D1`/`D2`/`X1` untouched). (+8 regression tests.) Note: under the 8086
  specialization this repo builds, a sign-only NaN is *not* a discriminating
  input — `f_abs` clears the sign before the NaN reaches the arithmetic — so the
  tests use signaling/payload-carrying NaNs.
- [x] ⚪ #35 `test_cdotu_rounding_modes` read uninitialized stack when built
  with optimization: clang at `-O1` sank the initializing stores of an in-loop
  compound literal past the loop's back-edge. Operand literals hoisted out of
  the loop as `static`. The suite now passes at `-O0`, `-O1` and `-O2` (it is
  still built without `-O` by the Makefile).

**Decision: bit-exact-only complex test oracle (PR #32)**
The complex gemm/gemv tests assert exact bit patterns against the per-op
SoftFloat-rounding oracle (bit-identical to the pure-Hoon path — the point of
SoftBLAS). This shapes which "closed-form spectra" cases are worth adding:
- *Included* — entries are Gaussian integers, so every intermediate is exact:
  DFT-4 unitarity `F·Fᴴ = 4·I` (the 4th roots of unity `{1,i,−1,−i}` are exactly
  representable; N=8 would need `1/√2`, so 4 is the ceiling), Pauli algebra
  (`σx·σy = i·σz`, `σy² = I`), and the Hermitian-Gram invariant
  `Aᴴ·A = (Aᴴ·A)ᴴ` (structural — holds for any integer A, no hand-computed
  expected matrix).
- *Excluded* — not bit-representable, so they'd require a tolerance harness that
  diverges from the bit-identical contract: general DFT/FFT round-trips,
  circulant diagonalization, arbitrary-angle rotations, eigen-spectra, Hilbert
  matrices. Deliberately out of scope for this suite.
- *Note:* an exact-cancellation component (`a−a`) is `−0` under round-toward-
  negative (`'d'`) and `+0` otherwise, so the rounding-mode test masks the sign
  bit. This is correct IEEE behaviour and doubles as proof that `rndMode` reaches
  the inner `f32` ops inside `c32_mul`.

---

## Open

**Phase 4 — coverage gaps (the main remaining work)**
- [x] **C1. NaN/Inf determinism tests.** Non-canonical-NaN→canonical now covered
  across precisions for axpy/dot/nrm2/nrm2_B/scal/asum/copy/swap/gemv/gemm, the
  quad family (qnrm2/qdot/qscal/qasum + length-2 qaxpy), and the Inf cases
  `Inf−Inf` and `0·Inf`. (Residual, folded into C4/C5: NaN through the *d/h/q*
  gemv/gemm transpose/padding paths once those tests land.)
- [x] **C2.** Negative-stride now covered for swap (`d/h/q`) and copy
  (`s/d/h/q`); dot/axpy already had it across precisions.
- [x] **C3.** `N==0` no-op: single covered per-routine; `d/h/q` L1 sweeps
  (`test_dn0`/`hn0`/`qn0`) over asum/dot/nrm2/axpy/scal/copy/swap; gemv/gemm
  `0_row`/`0_col` already cover all four precisions.
- [x] **C4.** gemv transpose for `d/h/q` (non-square 2×3); single gemm now has
  the full `{N,T}²` (added the `TT` combo). *Residual:* `d/h/q` gemm transpose
  combos — the transpose index logic is precision-independent (single validates
  the path; `d/h/q` gemm have non-square `NN` coverage already). ⚪
- [x] **C5.** Padded `lda` for gemv `d/h/q` (`test_*gemv_padlda`); gemm `_ldb`
  already covered all four precisions.
- [x] **C6.** All five rounding modes (`n`/`z`/`u`/`d`/`a`) × `saxpy`/`sdot`/
  `ddot`/`sgemv`, both signs (the negative tie distinguishes `z`/`d`/`u`/`a`).
  Reduction (L1) + matrix-vector (L2), single + double. Runs on macOS now.
- [x] **C7.** float128 `v[0]=lo / v[1]=hi` layout pinned by `test_qasum_layout`.
- [x] **C8.** Complex/rot routines build + have numeric tests; added NaN
  canonicalization (`caxpy`/`cscal`/`ccopy`) and negative-stride
  (`caxpy`/`ccopy`) coverage for the single-complex (`c`) family. *Out of scope:*
  half/double/quad complex (`i`/`z`/`v`) routines are unimplemented (README:
  "yet to implement") — a future feature, not audit cleanup.

**Phase 5 — hygiene**
- [x] **H1.** float128 `v[0]=lo / v[1]=hi` invariant documented at the top of
  `include/softblas.h`.
- [x] **B5.** `nrm2`/`nrm2_B`/`asum` (all precisions) now reject a stride where
  `(N-1)*incX` overflows `uint64_t` via the `SB_STRIDE_OVERFLOWS` guard, instead
  of walking off the buffer. ASan-verified regression test (`test_nrm2_huge_stride`).
- [x] **H4. Won't-do.** X-macro for the duplicated L1 routines: the elf and
  hobbit both argued explicit per-precision duplication is correct for an
  auditable deterministic-FP library (the `q*` pointer-API divergence is real,
  not accidental). Closing as a deliberate non-change.

**Follow-ons surfaced by the re-audit (need coordination)**
- [ ] Retire the classic `nrm2` (replace with `_B`) — gated on the `lagoon.hoon`
  rollover. Brief: `docs/lagoon-rollover-prompt.md`.
- [ ] `rndMode` → `typedef enum` (elf) — also gated on the Hoon rollover (the
  enum values equal the ASCII codes, so ABI-stable).
- [ ] Optional: split `softblas.h` into a public header + an internal header for
  the `nan_unify_*` / `_set_rounding` machinery (elf).

**Optional complex-test extensions (surfaced reviewing PR #32)**
- [ ] ⚪ Inf / signed-zero arithmetic through complex gemm — `Inf−Inf → NaN`,
  `0·Inf → NaN`, `−0` propagation. Stays bit-exact; cheap. Not yet added.
- [ ] ⚪ A tolerance-based spectral harness (separate from the bit-exact suite)
  for non-bit-exact invariants: DFT/FFT round-trips, unitary-rotation invariants,
  `‖A·x‖ = ‖x‖` for orthonormal A at arbitrary N. **Needs sign-off** — it
  introduces a non-bit-exact assertion path, changing the suite's contract; see
  the "bit-exact-only complex test oracle" decision above.
