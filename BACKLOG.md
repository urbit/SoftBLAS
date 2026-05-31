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

---

## Open

**Phase 4 — coverage gaps (the main remaining work)**
- [x] **C1. NaN/Inf determinism tests.** Non-canonical-NaN→canonical now covered
  across precisions for axpy/dot/nrm2/nrm2_B/scal/asum/copy/swap/gemv/gemm, the
  quad family (qnrm2/qdot/qscal/qasum + length-2 qaxpy), and the Inf cases
  `Inf−Inf` and `0·Inf`. (Residual, folded into C4/C5: NaN through the *d/h/q*
  gemv/gemm transpose/padding paths once those tests land.)
- [ ] **C2.** Negative-stride tests for the swap family and any strided routine
  lacking them.
- [ ] **C3.** `N==0` no-op across ALL routines (only `sswap_zero` exists today).
- [ ] **C4.** gemv/gemm transpose paths (`Trans='T'`; `transA×transB ∈ {N,T}²`),
  non-square.
- [ ] **C5.** `ldb≠ldc` / padded `lda` for `h/q` gemm and all gemv.
- [x] **C6.** All five rounding modes (`n`/`z`/`u`/`d`/`a`) × `saxpy`/`sdot`/
  `ddot`/`sgemv`, both signs (the negative tie distinguishes `z`/`d`/`u`/`a`).
  Reduction (L1) + matrix-vector (L2), single + double. Runs on macOS now.
- [x] **C7.** float128 `v[0]=lo / v[1]=hi` layout pinned by `test_qasum_layout`.
- [~] **C8.** Complex/rot routines build + have numeric tests (q-rot family,
  cscal/cswap/icamax done). **Remaining:** broaden complex coverage
  (strides, NaN payloads) for `caxpy`/`ccopy`/`cdotc`.

**Phase 5 — hygiene**
- [ ] **H1.** Document the `v[0]=lo / v[1]=hi` float128 invariant in the header.
- [ ] **B5. ⚪** `nrm2`/`asum` huge-`incX` index wrap. The `incX<1` guards stop
  zero/negative; a pathologically huge `incX` can still wrap `i*incX`. Low.
- [ ] **H4.** X-macro for the duplicated L1 routines — **recommend won't-do**: the
  elf and hobbit both argued explicit per-precision duplication is correct for an
  auditable deterministic-FP library.

**Follow-ons surfaced by the re-audit (need coordination)**
- [ ] Retire the classic `nrm2` (replace with `_B`) — gated on the `lagoon.hoon`
  rollover. Brief: `docs/lagoon-rollover-prompt.md`.
- [ ] `rndMode` → `typedef enum` (elf) — also gated on the Hoon rollover (the
  enum values equal the ASCII codes, so ABI-stable).
- [ ] Optional: split `softblas.h` into a public header + an internal header for
  the `nan_unify_*` / `_set_rounding` machinery (elf).
