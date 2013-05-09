SRCDIR=src
INCDIR=include
OBJDIR=obj

HLIST   = $(wildcard $(INCDIR)/*.h)
CPPLIST = $(wildcard $(SRCDIR)/*.cxx)
OLIST = $(patsubst $(SRCDIR)/%.cxx,$(OBJDIR)/%.o,$(CPPLIST))

CC=g++
CFLAGS=-Wall `root-config --cflags` -I./$(INCDIR)
LDFLAGS=`root-config --libs` -L./$(OBJDIR)


EXECUTABLE=main

all: $(OLIST)
	@echo "Linking " $@	
	@$(CC) $(LDFLAGS) $(OLIST) -o $(EXECUTABLE)


$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CFLAGS) -c $^ -o $@


clean:
	@echo "Cleaning up..."
	@rm $(OBJDIR)/*.o 
	@rm -rf $(OBJDIR)
	@rm $(EXECUTABLE)

