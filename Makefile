# 0-plain C, 1-naive SSE, 2-optimized SSE
TYPE     := 2
# compiler optimization level
O        := 2
# iterations count, with ITER%4==0 A should be identity, with ITER%4==1 A==B (rotion matrix)
ITER     := 100000001

SRC      := test.c
SUFFIX   := -o$(O)-t$(TYPE)-$(CC)
OBJ      := $(SRC:.c=$(SUFFIX).o)
CPPFLAGS := $(CPPFLAGS) -DTYPE=$(TYPE) -Wall -g -O$(O)
LDFLAGS  := $(LDFLAGS) -Wall -g -O$(O)

run: test$(SUFFIX)
	./test$(SUFFIX) $(ITER)

%$(SUFFIX).o: %.c
	$(CC) $(CPPFLAGS) -c -o $@ $<

test_$(TYPE)-$(CC): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

clean:
	rm -r $(OBJ)
