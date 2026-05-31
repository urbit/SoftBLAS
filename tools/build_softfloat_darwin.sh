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

#  Strip the leading '==> ' edit-markers to activate the template lines.
#  BSD sed (macOS) requires an explicit empty suffix for -i.
sed -i '' 's/^==> //' "$DEST/Makefile" "$DEST/platform.h"

make -C "$DEST"
echo "Built $DEST/softfloat.a"
