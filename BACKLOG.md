# SoftBLAS bug & test backlog

Consolidated from the five-reviewer audit, the Gnome/Angel test-coverage pass, and
the Urbit-headers analysis. Ordered as a work plan. State: `master` @ PR #1 merged
(per-call `rndMode`; no `softblas_roundingMode` global; no `softblas_state.c`).

Severity: 🔴 critical (memory-unsafe / silently wrong) · 🟠 high · 🟡 medium · ⚪ low.

---

## Re-audit pass (2026-05-30) — resolved & newly found

**Resolved this pass (verified by the five-agent council + tests, suite 200/202):**
- B1 signed-stride OOB (`*swap`/`csrot`) — params now `int64_t`, ASan-clean. ✅
- B2 `qaxpy` per-element NaN canonicalization. ✅
- B3 `FNEG` macro UB — width-correct unsigned masks. ✅
- B4 `i*amax` N==1 returns 0. ✅
- A1 `qasum` heap leak — stack accumulator; compiled core is allocation-free. ✅
- A3 `svec/dvec/hvec/qvec` moved to `tests/blas/include/test.h`. ✅
- T1 test misregistrations in `test_all.c` repointed. ✅
- T3 README rounding divergences (global→per-call, `'a'`=ties-away). ✅

**Fixed this pass (new code):**
- 🔴 `*rotmg` infinite loop on NaN/Inf `D1`/`D2`. `f_ge(d,gamsq)` = `!(f_lt(...))`
  is true forever for non-finite d. Guarded the `_ge` rescale loops with a
  finiteness check in `srotmg`/`drotmg`/`hrotmg` (+ regression tests). `qrotmg`
  carries the same guard from birth (see below).
- 🔴 7 ghost declarations — public prototypes with no source/Makefile entry:
  `qrot`/`qrotg`/`qrotm`/`qrotmg` (whole quad-rotation family) and
  `cscal`/`cswap`/`icamax`. **Implemented all 7**, added to the Makefile, fixed
  `qrot`'s wrong `float16_t` prototype → `float128_t`, added munit tests.
- 🟡 `asum` (s/d/h/q) and `gemv` (s/d/h/q) lacked stride guards: `asum` now
  rejects `incX < 1` (cf. nrm2); `gemv` rejects `incX < 1 || incY < 1` (signed
  `incX` would wrap to a huge OOB index). `i*amax` was already guarded.
- 🟡 `test_hgemm` off-by-one: three loops ran `i < 4` over `M*P`-element arrays
  (one OOB crash at M*P=3, one silent under-check at M*P=15) → `i < M*P`.

**New, still open:**
- ⚪ `snrm2`/`dnrm2`/`hnrm2` (classic running-scale) return NaN for subnormal
  inputs (`0/0`) and for equal near-overflow pairs (`Inf/Inf`); the new `*nrm2_B`
  variants are correct. Slated to **retire the classic versions** once Lagoon
  rolls to `_B` (see the rndMode-enum / nrm2 rollover work). Until then the bug
  rides in the soon-to-be-removed code.
- ⚪ macOS/arm64 only: the bundled SoftFloat's directional-rounding inexact path
  SIGBUSes (reproduces with pure `softfloat.h`, no SoftBLAS code). Not a SoftBLAS
  bug; Linux x86-64 CI is green. File upstream against the SoftFloat Darwin port.
  Blocks `test_saxpy_rounding{,_modes}` on Apple Silicon only.

---

## Phase 1 — Make the test suite honest (foundation; low risk, stays green)

- [ ] **T1. 16 test misregistrations in `tests/test_all.c`.** Registered display-name
  maps to a *different* function, so "140/140" overstates coverage. The correct
  functions exist and are declared.
  - 4× `*swap_stride` slots run `*swap_two` (lines 28, 55, 81, 108).
  - 12× `idamax/ihamax/iqamax` slots (every `_0/_12345/_stride/_13542`) run the
    `isamax` functions (lines 56–59, 82–85, 109–112).
  Repoint each to its real function. (Positive-stride/existing inputs → still green.)
- [ ] **T2. CI registration-honesty lint.** Assert every `{"/test_X", fn}` has
  name-stem == `fn` identifier. Prevents recurrence of T1.
- [ ] **T3. README doc divergences** (caller-breaking after PR #1):
  - rounding is now a per-call `rndMode` arg — README still says set the global
    `softblas_roundingMode` (`README.md:33`). Update.
  - `'a'` maps to `softfloat_round_near_maxMag` = **ties-away**, not unconditional
    away-from-zero as documented. Fix doc (or mapping).
  - v1.1.0 note still says errors are `exit(-1)`; PR #9 reverted to `return`.

## Phase 2 — Correctness bugs (fix + add the regression test that catches each)

- [ ] **B1. 🔴 `*swap`/`csrot` signed-stride.** `sswap/dswap/hswap/qswap` (and
  `csrot`) declare `incX`/`incY` as `uint64_t` then test `if (incX < 0)` (always
  false) → negative stride becomes `UINT64_MAX` → OOB read/write. Change params to
  `int64_t` (header + 4 `.c` files). Test: `sswap(3,[10,20,30],-1,[1,2,3],-1,'n')`
  → reversed result, ASan-clean. (Reproduced.)
- [ ] **B2. 🟠 `qaxpy` NaN canonicalization only `QY[0]`.** `nan_unify_q(QY)` touches
  one element; non-canonical NaNs in `QY[1..]` escape → determinism break. Loop
  `nan_unify_q` over all written elements (respect stride). Test: length-2 quad with
  a non-canonical NaN payload at index 1 → expect `{0, QUADNAN}`. (Reproduced.)
- [ ] **B3. 🟠 `FNEG` macro UB.** `include/softblas.h:175` `x.v ^ (1 << (…-1))`:
  `1<<31` is UB (UBSan-confirmed), `1<<63` wrong, and `c128_conj` doesn't compile
  (`uint64_t[2] ^ int`). Use width-correct unsigned masks; handle the float128
  array case. Test: `c32_conj`/`c64_conj` imag-bit flips, built `-fsanitize=undefined`.
- [ ] **B4. 🟡 `i*amax` N==1 returns `1`, 0-based for N>1.** `isamax.c:5` (+d/h/q).
  Return `0` for N==1. Test: `idamax(1,[42.0],1)==0` for all four precisions.
- [ ] **B5. ⚪ `*nrm2` loop bound overflow.** `1+(N-1)*incX` can wrap → early exit /
  OOB. Iterate `i=0..N-1` with `ix=i*incX`, or guard. (Reproduced w/ huge stride.)

## Phase 3 — Allocation & the Urbit-headers coupling (architecture)

- [ ] **A1. ⚪ `qasum` heap leak → make the core allocation-free.** `qasum.c:5`
  `qvec`/`malloc`s a 16-byte accumulator, never freed. Replace with a stack
  `float128_t`. This is the **only** heap allocation in the compiled library.
- [ ] **A2. Remove the vendored Urbit headers + `#ifdef VERE` block.** PR #1 copied
  **20** Urbit runtime headers into `include/` (`allocate.h, noun.h, imprison.h,
  manage.h, hashtable.h, jets.h, nock.h, vortex.h, …`), reachable only via
  `#ifdef VERE → #include "noun.h"` in `softblas.h:10-14`. **`VERE` is defined
  nowhere** — not in SoftBLAS's Makefile/CI, not in Vere's `ext/softblas/build.zig`
  — so the block is inert and the headers are dead weight. Even with `VERE` set,
  Vere's *own* headers are already on the include path when SoftBLAS builds inside
  it; the vendored copies are stale duplicates that will drift. After A1 the library
  allocates nothing, so there is nothing to route through `u3a`. **Delete the 20
  headers and the block.** If the host allocator is ever needed, take it via a
  one-line host-provided shim, not a copied header tree.
- [ ] **A3. ⚪ Move `svec/dvec/hvec/qvec` out of the public header** into a test-only
  header — they `malloc`, are used only by the 40 test files (never by Vere/Lagoon),
  and are the only other allocators. Add NULL checks (or `abort`) while moving.

## Phase 4 — New coverage (the systematic gaps)

- [ ] **C1. NaN/Inf determinism tests** — *zero coverage today*, and it's the whole
  point of SoftBLAS. Feed non-canonical NaN payloads through each reduction / axpy /
  scal / gemv / gemm and assert the canonical pattern out; include length>1 quad
  (overlaps B2). Add `Inf−Inf`, `0·Inf`.
- [ ] **C2. Negative-stride tests** for the swap family (overlaps B1) and any strided
  routine lacking them.
- [ ] **C3. `N==0` no-op across ALL routines** (only `sswap_zero` exists). Assert no
  read/write; documented return (0 for asum/nrm2/dot).
- [ ] **C4. gemv/gemm transpose paths** (`Trans='T'`; `transA×transB ∈ {N,T}²`) with
  non-square matrices — entirely untested.
- [ ] **C5. `ldb≠ldc` / padded `lda` across all four precisions** — only `s/d gemm_ldb`
  exist; `h/q` gemm and all gemv lack padded-leading-dimension tests.
- [ ] **C6. All five rounding modes × ≥2 routines**, incl. a negative operand (the
  only thing separating `z`/`d`/`u`/`a`); pin `'a'` = ties-away (or fix mapping);
  test unknown-mode = no-op.
- [ ] **C7. float128 `v[0]=lo / v[1]=hi` layout** pin test — every `q*` test depends
  on it implicitly.
- [~] **C8. Complex / rot routines.** Mostly resolved: `scnrm2`, `csrot`, the
  `s/d/h rot*` family, and now the full `q*rot*` family + `cscal`/`cswap`/`icamax`
  all compile, are built, and have numeric tests. The compile errors are gone.
  **Remaining:** broaden coverage (strides, NaN payloads) for the complex routines;
  audit `caxpy`/`ccopy`/`cdotc` for the same depth.

## Phase 5 — Hygiene (lower priority)

- [ ] H1. Document the `v[0]=lo/v[1]=hi` invariant in the header (not just a `qdot.c`
  comment). H2. ~~`benchmarking/` is dead — candidate for deletion.~~ **Keep:** revived
  for the `snrm2`/`snrm2_B` comparison that backs issue #3 (`benchmark_snrm2.c`,
  `benchmark_snrm2_B.c`). The stack-VLA footgun in `benchmark_snrm2.c` was fixed
  (heap-allocated, length-capped). H3. macOS build path (Makefile is GNU-ld/Linux only). H4. Consider an
  X-macro for the 6 identical Level-1 `h/s/d` routines.
