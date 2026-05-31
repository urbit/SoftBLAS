CC ?= gcc
#  -DSOFTFLOAT_FAST_INT64 matches the flag the SoftFloat archive is built with,
#  exposing the value-based f128_*/c128_* ops the quad-complex routines use.
CFLAGS = -DSOFTFLOAT_FAST_INT64 -Iinclude -ISoftFloat/source/include -Imunit -Itests/blas/include

#  Pick the SoftFloat build per host OS. The Linux archive links with GNU ld's
#  `-l:softfloat.a`; macOS (BSD ld) has no `-l:` form, so link the archive by
#  path. Build the macOS archive first with tools/build_softfloat_darwin.sh
#  (CI does this; see .github/workflows/ci.yml).
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SOFTFLOAT_DIR = SoftFloat/build/Darwin-$(shell uname -m)-clang
  LDFLAGS =
  LIBS = $(SOFTFLOAT_DIR)/softfloat.a
else
  SOFTFLOAT_DIR = SoftFloat/build/Linux-x86_64-GCC
  LDFLAGS = -L$(SOFTFLOAT_DIR)/
  LIBS = -l:softfloat.a
endif

MUNIT_SRC = munit/munit.c
MUNIT_OBJ = $(MUNIT_SRC:.c=.o)

BLAS_SRC_DIR_L1 = ./src/blas/level1
BLAS_SRC_DIR_L2 = ./src/blas/level2
BLAS_SRC_DIR_L3 = ./src/blas/level3
BLAS_SRCS = \
  $(BLAS_SRC_DIR_L1)/sasum.c \
  $(BLAS_SRC_DIR_L1)/dasum.c \
  $(BLAS_SRC_DIR_L1)/hasum.c \
  $(BLAS_SRC_DIR_L1)/qasum.c \
  $(BLAS_SRC_DIR_L1)/saxpy.c \
  $(BLAS_SRC_DIR_L1)/daxpy.c \
  $(BLAS_SRC_DIR_L1)/haxpy.c \
  $(BLAS_SRC_DIR_L1)/qaxpy.c \
  $(BLAS_SRC_DIR_L1)/scopy.c \
  $(BLAS_SRC_DIR_L1)/dcopy.c \
  $(BLAS_SRC_DIR_L1)/hcopy.c \
  $(BLAS_SRC_DIR_L1)/qcopy.c \
  $(BLAS_SRC_DIR_L1)/sdot.c \
  $(BLAS_SRC_DIR_L1)/ddot.c \
  $(BLAS_SRC_DIR_L1)/hdot.c \
  $(BLAS_SRC_DIR_L1)/qdot.c \
  $(BLAS_SRC_DIR_L1)/sdsdot.c \
  $(BLAS_SRC_DIR_L1)/hsdot.c \
  $(BLAS_SRC_DIR_L1)/snrm2.c \
  $(BLAS_SRC_DIR_L1)/dnrm2.c \
  $(BLAS_SRC_DIR_L1)/hnrm2.c \
  $(BLAS_SRC_DIR_L1)/qnrm2.c \
  $(BLAS_SRC_DIR_L1)/snrm2_B.c \
  $(BLAS_SRC_DIR_L1)/dnrm2_B.c \
  $(BLAS_SRC_DIR_L1)/hnrm2_B.c \
  $(BLAS_SRC_DIR_L1)/qnrm2_B.c \
  $(BLAS_SRC_DIR_L1)/sscal.c \
  $(BLAS_SRC_DIR_L1)/dscal.c \
  $(BLAS_SRC_DIR_L1)/hscal.c \
  $(BLAS_SRC_DIR_L1)/qscal.c \
  $(BLAS_SRC_DIR_L1)/sswap.c \
  $(BLAS_SRC_DIR_L1)/dswap.c \
  $(BLAS_SRC_DIR_L1)/hswap.c \
  $(BLAS_SRC_DIR_L1)/qswap.c \
  $(BLAS_SRC_DIR_L1)/isamax.c \
  $(BLAS_SRC_DIR_L1)/idamax.c \
  $(BLAS_SRC_DIR_L1)/ihamax.c \
  $(BLAS_SRC_DIR_L1)/iqamax.c \
  $(BLAS_SRC_DIR_L1)/scasum.c \
  $(BLAS_SRC_DIR_L1)/caxpy.c \
  $(BLAS_SRC_DIR_L1)/ccopy.c \
  $(BLAS_SRC_DIR_L1)/cdotc.c \
  $(BLAS_SRC_DIR_L1)/scnrm2.c \
  $(BLAS_SRC_DIR_L1)/csrot.c \
  $(BLAS_SRC_DIR_L1)/srot.c \
  $(BLAS_SRC_DIR_L1)/drot.c \
  $(BLAS_SRC_DIR_L1)/hrot.c \
  $(BLAS_SRC_DIR_L1)/srotg.c \
  $(BLAS_SRC_DIR_L1)/drotg.c \
  $(BLAS_SRC_DIR_L1)/hrotg.c \
  $(BLAS_SRC_DIR_L1)/srotm.c \
  $(BLAS_SRC_DIR_L1)/drotm.c \
  $(BLAS_SRC_DIR_L1)/hrotm.c \
  $(BLAS_SRC_DIR_L1)/srotmg.c \
  $(BLAS_SRC_DIR_L1)/drotmg.c \
  $(BLAS_SRC_DIR_L1)/hrotmg.c \
  $(BLAS_SRC_DIR_L1)/qrot.c \
  $(BLAS_SRC_DIR_L1)/qrotg.c \
  $(BLAS_SRC_DIR_L1)/qrotm.c \
  $(BLAS_SRC_DIR_L1)/qrotmg.c \
  $(BLAS_SRC_DIR_L1)/cscal.c \
  $(BLAS_SRC_DIR_L1)/cswap.c \
  $(BLAS_SRC_DIR_L1)/icamax.c \
  $(BLAS_SRC_DIR_L1)/iaxpy.c \
  $(BLAS_SRC_DIR_L1)/zaxpy.c \
  $(BLAS_SRC_DIR_L1)/vaxpy.c \
  $(BLAS_SRC_DIR_L1)/icopy.c \
  $(BLAS_SRC_DIR_L1)/zcopy.c \
  $(BLAS_SRC_DIR_L1)/vcopy.c \
  $(BLAS_SRC_DIR_L1)/idotc.c \
  $(BLAS_SRC_DIR_L1)/zdotc.c \
  $(BLAS_SRC_DIR_L1)/vdotc.c \
  $(BLAS_SRC_DIR_L1)/iscal.c \
  $(BLAS_SRC_DIR_L1)/zscal.c \
  $(BLAS_SRC_DIR_L1)/vscal.c \
  $(BLAS_SRC_DIR_L1)/iswap.c \
  $(BLAS_SRC_DIR_L1)/zswap.c \
  $(BLAS_SRC_DIR_L1)/vswap.c \
  $(BLAS_SRC_DIR_L1)/ihrot.c \
  $(BLAS_SRC_DIR_L1)/zdrot.c \
  $(BLAS_SRC_DIR_L1)/vqrot.c \
  $(BLAS_SRC_DIR_L1)/hiasum.c \
  $(BLAS_SRC_DIR_L1)/dzasum.c \
  $(BLAS_SRC_DIR_L1)/qvasum.c \
  $(BLAS_SRC_DIR_L1)/hinrm2.c \
  $(BLAS_SRC_DIR_L1)/dznrm2.c \
  $(BLAS_SRC_DIR_L1)/qvnrm2.c \
  $(BLAS_SRC_DIR_L1)/iiamax.c \
  $(BLAS_SRC_DIR_L1)/izamax.c \
  $(BLAS_SRC_DIR_L1)/ivamax.c \
  $(BLAS_SRC_DIR_L2)/sgemv.c \
  $(BLAS_SRC_DIR_L2)/dgemv.c \
  $(BLAS_SRC_DIR_L2)/hgemv.c \
  $(BLAS_SRC_DIR_L2)/qgemv.c \
  $(BLAS_SRC_DIR_L3)/sgemm.c \
  $(BLAS_SRC_DIR_L3)/dgemm.c \
  $(BLAS_SRC_DIR_L3)/hgemm.c \
  $(BLAS_SRC_DIR_L3)/qgemm.c
#   $(BLAS_SRC_DIR_L1)/srotg.c
BLAS_OBJS = $(BLAS_SRCS:.c=.o)

TEST_SRC_DIR_L1 = ./tests/blas/level1
TEST_SRC_DIR_L2 = ./tests/blas/level2
TEST_SRC_DIR_L3 = ./tests/blas/level3
TEST_SRCS = $(wildcard $(TEST_SRC_DIR_L1)/*.c) \
            $(wildcard $(TEST_SRC_DIR_L2)/*.c) \
            $(wildcard $(TEST_SRC_DIR_L3)/*.c)
TEST_OBJS = $(TEST_SRCS:.c=.o)

TEST_ALL_SRC = ./tests/test_all.c
TEST_ALL_OBJ = $(TEST_ALL_SRC:.c=.o)

# set default target to library
.DEFAULT_GOAL := library

TARGET = libsoftblas.a
TEST_TARGET = test_all

library: $(TARGET)

tests: $(TEST_TARGET)

$(TARGET): $(BLAS_OBJS)
	ar rcs $@ $^

$(TEST_TARGET): $(TEST_OBJS) $(BLAS_OBJS) $(MUNIT_OBJ) $(TEST_ALL_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIBS)

$(TEST_OBJS): $(TEST_SRC_DIR_L1)/%.o: $(TEST_SRC_DIR_L1)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_OBJS): $(TEST_SRC_DIR_L2)/%.o: $(TEST_SRC_DIR_L2)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_OBJS): $(TEST_SRC_DIR_L3)/%.o: $(TEST_SRC_DIR_L3)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_ALL_OBJ): $(TEST_ALL_SRC)
	$(CC) $(CFLAGS) -c $< -o $@

$(BLAS_OBJS): $(BLAS_SRC_DIR_L1)/%.o: $(BLAS_SRC_DIR_L1)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BLAS_OBJS): $(BLAS_SRC_DIR_L2)/%.o: $(BLAS_SRC_DIR_L2)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BLAS_OBJS): $(BLAS_SRC_DIR_L3)/%.o: $(BLAS_SRC_DIR_L3)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(MUNIT_OBJ): $(MUNIT_SRC)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TEST_OBJS) $(BLAS_OBJS) $(MUNIT_OBJ) $(TEST_ALL_OBJ)

allclean:
	rm -f $(TEST_OBJS) $(BLAS_OBJS) $(MUNIT_OBJ) $(TEST_ALL_OBJ) $(TEST_TARGET)
