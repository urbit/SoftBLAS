# SoftBLAS — Handoff to urbit/numerics agent

_Scratch document — not committed. Written 2026-05-30 after a five-reviewer audit
(Gnome/Dwarf/Elf/Hobbit/Angel) of `urbit/SoftBLAS` and the follow-up fixes._

## TL;DR

A correctness audit found three silently-wrong-answer bugs, a process-killing
error path, and an inert rounding-mode API. The first four are **fixed and merged
to `master`**; CI now runs on every PR. The rounding-mode work lives in **PR #1**,
which has been completed and reconciled with `master` but is **a breaking API
change** that requires a coordinated Lagoon update in Vere.

Your job picks up at: (1) finish the Lagoon side of the `rndMode` change, (2) bump
the pinned SoftBLAS commit in Vere, (3) work the remaining audit backlog below.

## What shipped this session

| PR | State | Summary |
|----|-------|---------|
| **#8** | merged | `f128_abs` masked both float128 words (corrupted negative quads); `*gemm` indexed B with `ldc` instead of `ldb`; `*gemv` folded `incX/incY` into the matrix index. All three produced silently wrong results, none caught by the old tests. + regression tests. |
| **#9** | merged | `exit(-1)` on invalid args (bad transpose/layout; `swap` N==0) → replaced with `return`. `swap` N==0 is now a no-op per BLAS spec. (Vere is the only client; it validates and bails via `u3m_bail`, so a plain `return` suffices — no error-callback machinery.) |
| **#10** | merged | GitHub Actions CI: builds SoftFloat + SoftBLAS, runs the munit suite on every push/PR. |
| **#1** | open, **mergeable**, **CI green** | Rounding-mode redesign — see below. Branch `sigilante/rounding-mode`. Merged up to `master` so it no longer reintroduces the #8 bugs. 140/140 tests pass (locally and in CI). |

`master` is green (139 tests). The bugs were invisible before because every prior
test used unit strides and `ldb == ldc`.

## PR #1 (`sigilante/rounding-mode`) — what it is and what's left

**Design:** rounding is now **per-call**. Every arithmetic routine gained a
trailing `const uint_fast8_t rndMode` parameter and calls `_set_rounding(rndMode)`
on entry (sets SoftFloat's `softfloat_roundingMode`). The old global
`softblas_roundingMode` + `softblas_state.c` are **removed**. This supersedes the
audit's "problem #5" (the old `set_rounding()` was never called, so the documented
rounding control silently did nothing).

`_set_rounding(a)` accepts `'n' 'z' 'u' 'd' 'a'`; an **unknown mode is a no-op**
(leaves SoftFloat's current mode unchanged) rather than `exit(1)`.

**Completed in this handoff's work (pushed to the branch):**
- Converted all remaining implementations to match the header's `rndMode`
  prototypes; restored the `'a'` (round-away) mode; made unknown-mode a no-op.
- Updated every test call site; added `test_saxpy_rounding`, which proves the mode
  actually changes results (`1.0 + 2⁻²⁴` → `1.0` under `'n'`, next-float-up under
  `'u'`). Would fail if `_set_rounding` were inert.
- Merged `origin/master` in: #8's fixes + `rndMode` now coexist (verified:
  `sgemm` has `ldb` + `rndMode`; `sgemv` has `A[i*lda+j]` + `rndMode`; `f128_abs`
  fixed). 140/140 pass.

**⚠️ This is a breaking API change.** Every routine signature gained a trailing
arg. All downstream callers must pass a mode.

**Deliberately NOT given `rndMode`** (they do no rounding / header omits it):
`isamax/idamax/ihamax/iqamax` (index-of-max = abs + compare), and `sdsdot/hsdot`
(their prototypes don't take it). Leave these alone.

**Converted but NOT built** (kept header-consistent for later): the complex
(`caxpy/ccopy/cdotc/csrot/cscal/scasum/scnrm2`) and rotation
(`*rot/*rotg/*rotm`) routines now carry `rndMode`, but they are **not** in the
Makefile and several still have pre-existing compile bugs (see backlog). `*rotmg`
were left unconverted (multi-line signatures; also unbuilt).

## Immediate action items (Vere / Lagoon side)

1. **Lagoon `rndMode` update — REQUIRED, do in tandem with merging PR #1.**
   `vere/pkg/noun/jets/i/lagoon.c` has **~32 BLAS call sites** (`saxpy/sdot/sscal/
   sgemm` and h/d/q variants) that call the *old* signatures. Each needs the
   trailing `rndMode` argument. Lagoon already computes the mode — see
   `_set_rounding(c3_w a)` at ~`lagoon.c:40`, which maps `%n/%z/%u/%d/%a`. Thread
   that same char (the jet's rounding-mode argument) into each BLAS call. The
   moment PR #1 merges, Vere won't compile until this is done.
   - Note: `i*amax` calls take **no** `rndMode` — don't add it there.
   - Lagoon currently sets `softfloat_roundingMode` directly (and the now-removed
     `softblas_roundingMode`); after this change, drop the
     `softblas_roundingMode = ...` writes (that global no longer exists) and rely
     on passing `rndMode` per call.

2. **Bump the pinned SoftBLAS commit in Vere.**
   `vere/ext/softblas/build.zig.zon:10` pins `cbffb33f…` (an old commit that has
   *all three* #8 bugs). After merges land, bump `.url` + `.hash` to a `master`
   commit that includes #8/#9/#10 (and PR #1 once merged). Also check
   `vere/ext/softblas/build.zig`'s source list still matches SoftBLAS's
   `BLAS_SRCS` (no `softblas_state.c` anymore).

3. **Merge PR #1** once Lagoon is ready. It's mergeable; "UNSTABLE" only means CI
   isn't wired on that branch yet (CI reached `master` via #10 after PR #1 forked;
   the merge commit already pulled `master` in, so a re-merge of latest `master`
   would attach CI).

## Remaining SoftBLAS backlog (from the audit, not yet done)

A fuller list was kept in `WORKLIST.md` (was stashed under
`stash@{0}` on the `ci/github-actions` branch: `git stash show -p stash@{0}`).
Highlights, roughly by severity:

**Correctness / determinism**
- `qaxpy`: `nan_unify_q(QY)` only canonicalizes `QY[0]`; NaNs in `QY[1..]` escape →
  breaks determinism for vectors length > 1. Unify per written element.
- `*nrm2`: loop bound `(N-1)*incX` can integer-overflow → early exit / infinite loop.
- `FNEG` macro (`1 << 31` / `1 << 63`): shifts a signed int into the sign bit = UB.
  Affects complex conjugation. Use a width-correct unsigned literal.
- `*swap`/`csrot`: declare `incX` as `uint64_t` then test `if (incX < 0)` (always
  false) → negative strides become huge positives → OOB. (May differ post-PR#1.)
- `iamax`: 0-based for N>1 but hardcodes `return 1` for N==1 — inconsistent; document
  the convention (Fortran BLAS is 1-based; this is effectively 0-based).
- `f128M_abs`: now masks the right word, but is still a structurally broken
  `(float128_t*){…}` compound literal (currently unused).

**Resource / robustness**
- `qasum` leaks 16 bytes/call (heap accumulator never freed; use a stack temp).
- `svec/dvec/hvec/qvec` don't null-check `malloc`.

**Broken / orphaned code presented as "complete"**
- `rot/rotg/rotm/rotmg` (all precisions) are declared in the header but **not built**
  → link errors if a caller uses them. Either build+fix or document as unbuilt.
- `sdsdot.c`, `hsdot.c`, `scnrm2.c` don't compile (undefined idents / complex type
  errors). `qrot` prototype uses `float16_t*` (copy-paste). `*rotmg` compute
  constants via integer multiply of hex bit-patterns (meaningless). `srot/drot/hrot`
  compare SoftFloat structs to int literals with `!=` (UB) and omit
  `#include "softblas.h"`.

**Test-harness bugs (in `tests/test_all.c`)**
- `*swap_stride` slots run `*swap_two` (the stride variant is never exercised).
- `idamax/ihamax/iqamax` stride+13542 slots point at the **`isamax`** functions —
  the d/h/q index-max tests don't test their own routine.

**Hygiene / scope**
- `v[0]=lo, v[1]=hi` float128 layout invariant should be documented in the header
  (currently an `// XXX` note in `qdot.c`).
- Consider an X-macro for the 6 genuinely-identical Level-1 routines
  (axpy/copy/dot/scal/swap/asum) across `h/s/d`; keep `q` separate (pointer API).
- Move internal helpers (`*vec`, `nan_unify_*`) to an internal header.
- `benchmarking/` is dead (benchmarking was removed from planned work, commit
  f4908d3) — candidate for deletion. Unused complex 16/128-bit (`i`/`v`) type stubs.
- Naming sovereignty is limited by downstream clients — coordinate renames.

## Build & test

- **CI / Linux (canonical):** submodules (SoftFloat is an SSH submodule — rewrite to
  HTTPS), `make -C SoftFloat/build/Linux-x86_64-GCC`, then `make tests && ./test_all`.
  See `.github/workflows/ci.yml`.
- **macOS gotcha:** the Makefile hardcodes `gcc`, `-LSoftFloat/build/Linux-x86_64-GCC`,
  and `-l:softfloat.a` (GNU-ld only) — it does **not** build out-of-the-box on
  Apple `ld`. During this work I built locally with `clang`, linking the checked-in
  (arm64) `softfloat.a` directly and compiling the `BLAS_SRCS` list by hand. A
  Darwin build path is an open backlog item.

## Landmines / things that bit me

- The whole reason the bugs survived: tests only used unit strides and `ldb==ldc`.
  Any new routine work should add non-unit-stride / `ldb≠ldc` coverage.
- PR #1's source merged cleanly with #8 because the changes were on different lines
  (signature/`_set_rounding` at top vs. index fixes in the loop body) — but always
  re-verify `f128_abs`, `*gemm` `ldb`, and `*gemv` indexing after any future merge.
- When adding `rndMode` to test calls, watch for calls with trailing comments
  (`sswap(0,…);  // …`) — a naive "ends with `);`" filter misses them.
- Merging `master` into PR #1 silently re-took master's **Makefile**, which lists
  `softblas_state.c` (deleted here) and the WIP complex routines (`scasum`/`ccopy`
  don't compile) — this broke CI. Resolved by restoring the branch's `BLAS_SRCS`
  (40 REAL routines) and re-deleting `softblas_state.c`. On any future `master`
  merge, re-check the Makefile didn't resurrect either. The branch's complex/rot
  files carry `rndMode` but are intentionally unbuilt and still WIP-broken
  (e.g. `ccopy` declares `ix/iy` but uses `iX/iY`) — fix when completing complex.
