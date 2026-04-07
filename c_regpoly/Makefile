# ============================================================
# MinimalCode Makefile — builds POL only
# Usage: ./POL 1 example1 32poly.dat
# ============================================================

DIRPROJECT  = /home/frpan/projets/claude_projects/regpoly/MinimalCode/c_regpoly
DIRSRC      = $(DIRPROJECT)/src
DIRINCLUDE  = $(DIRPROJECT)/include
DIRLIB      = $(DIRPROJECT)/lib

# ---- ZEN library --------------------------------------------
DIRZEN      = /home/frpan/projets/claude_projects/regpoly/ZEN
DIRZENLIBZ  = $(DIRZEN)/zen/lib/linuxx8664
DIRZENLIBF  = $(DIRZEN)/zenfact/lib/linuxx8664
DIRZENINCZ  = $(DIRZEN)/zen/include
DIRZENINCF  = $(DIRZEN)/zenfact/include
ZL          = -L $(DIRZENLIBZ) -l:gcco_libzen.a -L $(DIRZENLIBF) -l:gcco_libzenfact.a
ZI          = -I $(DIRZENINCZ) -I $(DIRZENINCF)

# ---- TestU01 library ----------------------------------------
DIRTESTU01LIB = /home/frpan/projets/claude_projects/regpoly/testu01/install/lib
DIRTESTU01INC = /home/frpan/projets/claude_projects/regpoly/testu01/install/include
CLINK         = -L$(DIRTESTU01LIB) -Wl,-rpath,$(DIRTESTU01LIB) -ltestu01 -lprobdist -lmylib
CINC          = -I$(DIRTESTU01INC)

# ---- Compiler -----------------------------------------------
CC    = gcc
OPT   = -O2 -funroll-loops -ffast-math -fomit-frame-pointer
WARN  = -Wall -W -Wmissing-prototypes -Wstrict-prototypes
CFLAGS = $(OPT) $(WARN) -I$(DIRINCLUDE) $(ZI) $(CINC)

# ---- Modules ------------------------------------------------
MODULES = \
  mecf        vectorsF2   timer       combinaisons  regpoly   \
  tempering   tuplets     polylcg     MT            genf2w    \
  temperMK    permut      readfiles   tausworthe    tgfsr     \
  basisreduc  myzen       polynomials                         \
  AC1D        carry2      marsaxorshift matsumoto             \
  marsa       custom      tmsnets     GFSRWU

OBJECTS = $(addprefix $(DIRLIB)/, $(addsuffix .o, $(MODULES)))
REGPOLYLIB = $(DIRLIB)/libREGPOLY.a

# ---- Targets ------------------------------------------------
all: POL

$(REGPOLYLIB): $(OBJECTS)
	ar rcs $@ $^

$(DIRLIB)/%.o: $(DIRSRC)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

POL: $(REGPOLYLIB)
	$(CC) $(CFLAGS) $(DIRSRC)/rechercheguide.c \
	      -L$(DIRLIB) -lREGPOLY $(ZL) $(CLINK) -lm \
	      -o $(DIRPROJECT)/POL

test_transition_matrix_generator: $(REGPOLYLIB)
	$(CC) $(CFLAGS) $(DIRSRC)/test_transition_matrix_generator.c \
	      -L$(DIRLIB) -lREGPOLY $(ZL) -lm \
	      -o $(DIRPROJECT)/test_transition_matrix_generator

clean:
	rm -f $(DIRLIB)/*.o $(REGPOLYLIB) $(DIRPROJECT)/POL \
	      $(DIRPROJECT)/test_transition_matrix_generator

.PHONY: all clean POL test_transition_matrix_generator
