CC = gcc
CFLAGS = -Iinclude -ISoftFloat/source/include -Imunit -Itests/blas/include
LDFLAGS = -LSoftFloat/build/Linux-x86_64-GCC/
LIBS = -l:softfloat.a

MUNIT_SRC = munit/munit.c
MUNIT_OBJ = $(MUNIT_SRC:.c=.o)

BLAS_SRC_DIR = ./src/blas/level1
BLAS_SRCS = \
  $(BLAS_SRC_DIR)/sasum.c \
  $(BLAS_SRC_DIR)/dasum.c \
  $(BLAS_SRC_DIR)/hasum.c \
  $(BLAS_SRC_DIR)/qasum.c \
  $(BLAS_SRC_DIR)/saxpy.c \
  $(BLAS_SRC_DIR)/daxpy.c \
  $(BLAS_SRC_DIR)/haxpy.c \
  $(BLAS_SRC_DIR)/qaxpy.c \
  $(BLAS_SRC_DIR)/scopy.c \
  $(BLAS_SRC_DIR)/dcopy.c \
  $(BLAS_SRC_DIR)/hcopy.c \
  $(BLAS_SRC_DIR)/qcopy.c \
  $(BLAS_SRC_DIR)/sdot.c \
  $(BLAS_SRC_DIR)/ddot.c \
  $(BLAS_SRC_DIR)/hdot.c \
  $(BLAS_SRC_DIR)/qdot.c \
  $(BLAS_SRC_DIR)/snrm2.c \
  $(BLAS_SRC_DIR)/dnrm2.c \
  $(BLAS_SRC_DIR)/hnrm2.c \
  $(BLAS_SRC_DIR)/qnrm2.c \
  $(BLAS_SRC_DIR)/sscal.c \
  $(BLAS_SRC_DIR)/dscal.c \
  $(BLAS_SRC_DIR)/hscal.c \
  $(BLAS_SRC_DIR)/qscal.c \
  $(BLAS_SRC_DIR)/sswap.c \
  $(BLAS_SRC_DIR)/dswap.c \
  $(BLAS_SRC_DIR)/hswap.c \
  $(BLAS_SRC_DIR)/qswap.c \
  $(BLAS_SRC_DIR)/isamax.c \
  $(BLAS_SRC_DIR)/idamax.c \
  $(BLAS_SRC_DIR)/ihamax.c \
  $(BLAS_SRC_DIR)/iqamax.c
#   $(BLAS_SRC_DIR)/srotg.c
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

$(BLAS_OBJS): $(BLAS_SRC_DIR)/%.o: $(BLAS_SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(MUNIT_OBJ): $(MUNIT_SRC)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TEST_OBJS) $(BLAS_OBJS) $(MUNIT_OBJ) $(TEST_ALL_OBJ)

allclean:
	rm -f $(TEST_OBJS) $(BLAS_OBJS) $(MUNIT_OBJ) $(TEST_ALL_OBJ) $(TARGET)
