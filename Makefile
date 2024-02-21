CC = gcc
CFLAGS = -Iinclude -ISoftFloat/source/include -Imunit -Itests/blas/include
LDFLAGS = -LSoftFloat/build/Linux-x86_64-GCC/
LIBS = -l:softfloat.a

SRC_MUNIT = munit/munit.c
OBJ_MUNIT = $(SRC_MUNIT:.c=.o)

SRC_BLAS = src/blas/level1/sasum.c
OBJ_BLAS = $(SRC_BLAS:.c=.o)

SRC_TEST = ./tests/blas/level1/sasum.c
OBJ_TEST = $(SRC_TEST:.c=.o)

TARGET = test_sasum

all: $(TARGET)

$(TARGET): $(OBJ_TEST) $(OBJ_BLAS) $(OBJ_MUNIT)
	$(CC) $(LDFLAGS) $^ -o $(TARGET) $(LIBS)

$(OBJ_TEST): $(SRC_TEST)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_BLAS): $(SRC_BLAS)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_MUNIT): $(SRC_MUNIT)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ_TEST) $(OBJ_BLAS) $(OBJ_MUNIT) $(TARGET)

