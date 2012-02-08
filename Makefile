TYPE     := 2
O        := 2
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
