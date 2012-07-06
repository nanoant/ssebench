# matrix (SSE matrix multiplication), scimark2_c (C scimark)
TEST     := matrix
# 0-plain C, 1-naive SSE, 2-optimized SSE
TYPE     := 2
# compiler optimization level
O        := 2
# iterations count, with ITER%4==0 A should be identity, with ITER%4==1 A==B (rotion matrix)
ITER     := 100000001

SRC      := $(TEST).c
OBJDIR   := build.$(shell uname -s)
ARCH     := native
SUFFIX   := -o$(O)-t$(TYPE)-$(CC)-$(ARCH)
BIN      := $(OBJDIR)/$(TEST)$(SUFFIX)
OBJ      := $(addprefix $(OBJDIR)/,$(SRC:.c=$(SUFFIX).o))
CPPFLAGS := $(CPPFLAGS) -g -fopenmp -DTYPE=$(TYPE) -Wall -O$(O) -march=$(ARCH)
LDFLAGS  := $(LDFLAGS) -fopenmp -Wall -O$(O)

run: $(BIN)
	./$(BIN) $(ITER)

$(OBJDIR):; mkdir $(OBJDIR)

$(OBJDIR)/%$(SUFFIX).o: %.c $(OBJDIR)
	$(CC) $(CPPFLAGS) -c -S -o $(@:.o=.s) $<
	$(CC) $(CPPFLAGS) -c -o $@ $<

$(BIN): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

clean:
	rm -r $(OBJ)
