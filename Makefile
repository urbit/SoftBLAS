CC = gcc
CFLAGS = -Iinclude -ISoftFloat/source/include -Imunit -Itests/blas/include
LDFLAGS = -LSoftFloat/build/Linux-x86_64-GCC/
LIBS = -l:softfloat.a

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
  $(BLAS_SRC_DIR_L1)/snrm2.c \
  $(BLAS_SRC_DIR_L1)/dnrm2.c \
  $(BLAS_SRC_DIR_L1)/hnrm2.c \
  $(BLAS_SRC_DIR_L1)/qnrm2.c \
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

TEST_SRC_DIR = ./tests/blas/level1
TEST_SRCS = $(wildcard $(TEST_SRC_DIR)/*.c)
TEST_OBJS = $(TEST_SRCS:.c=.o)

TEST_ALL_SRC = ./tests/test_all.c
TEST_ALL_OBJ = $(TEST_ALL_SRC:.c=.o)

TARGET = test_all

tests: $(TARGET)

$(TARGET): $(TEST_OBJS) $(BLAS_OBJS) $(MUNIT_OBJ) $(TEST_ALL_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIBS)

$(TEST_OBJS): $(TEST_SRC_DIR)/%.o: $(TEST_SRC_DIR)/%.c
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
	rm -f $(TEST_OBJS) $(BLAS_OBJS) $(MUNIT_OBJ) $(TEST_ALL_OBJ) $(TARGET)
