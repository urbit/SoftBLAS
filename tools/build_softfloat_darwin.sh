#!/bin/sh
#  Build the SoftFloat static archive for macOS (Darwin) from the generic
#  FAST_INT64 template. The vendored SoftFloat ships ready-made build dirs for
#  Linux/Windows but none for Darwin, and its Makefile/platform.h carry `==>`
#  edit-markers that must be stripped to activate. This reproduces, in CI and
#  locally, the dir the Makefile expects at SoftFloat/build/Darwin-<arch>-clang.
#
#  Idempotent: re-running with an existing softfloat.a is a no-op.
set -eu

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SF_BUILD="$ROOT/SoftFloat/build"
ARCH="$(uname -m)"
DEST="$SF_BUILD/Darwin-${ARCH}-clang"
TEMPLATE="$SF_BUILD/template-FAST_INT64"

if [ -f "$DEST/softfloat.a" ]; then
    echo "SoftFloat archive already built: $DEST/softfloat.a"
    exit 0
fi

if [ ! -d "$TEMPLATE" ]; then
    echo "error: SoftFloat template not found at $TEMPLATE" >&2
    echo "       (did the SoftFloat submodule get checked out?)" >&2
    exit 1
fi

echo "Generating $DEST from $TEMPLATE ..."
rm -rf "$DEST"
cp -R "$TEMPLATE" "$DEST"

#  Strip the leading '==> ' edit-markers to activate the template Makefile.
#  BSD sed (macOS) requires an explicit empty suffix for -i.
sed -i '' 's/^==> //' "$DEST/Makefile"

#  Use the GCC/clang platform.h, NOT the bare generic template one. The generic
#  template defines `THREAD_LOCAL _Thread_local`, which makes SoftFloat's state
#  (softfloat_roundingMode/exceptionFlags/...) genuinely thread-local; on macOS
#  every access then routes through dyld's TLV machinery and BUS-faults in the
#  directional-rounding inexact path. The Linux-x86_64-GCC config leaves
#  THREAD_LOCAL empty (plain globals) and enables the GCC builtins
#  (SOFTFLOAT_BUILTIN_CLZ, SOFTFLOAT_INTRINSIC_INT128, opts-GCC.h); Apple/LLVM
#  clang is GCC-compatible, so this header is correct for arm64 macOS too.
cp "$SF_BUILD/Linux-x86_64-GCC/platform.h" "$DEST/platform.h"

make -C "$DEST"
echo "Built $DEST/softfloat.a"
