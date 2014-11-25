PROGRAM_NAME = overlap
VERSION = 1.02
CC = g++
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
SDIR = code
BDIR = .

.PHONY: clean default build distclean dist 

default : build

$(BDIR)/overlap.o : $(SDIR)/overlap.cpp
	$(CC) -c $(SDIR)/overlap.cpp -o $(BDIR)/overlap.o

$(BDIR)/String.o : $(SDIR)/String.cpp
	$(CC) -c $(SDIR)/String.cpp -o $(BDIR)/String.o

$(BDIR)/common.o : $(SDIR)/common.cpp
	$(CC) -c $(SDIR)/common.cpp -o $(BDIR)/common.o

clean:
	rm -rf $(BDIR)/*.o $(BDIR)/$(PROGRAM_NAME)

distclean: clean
	rm -rf *.tar.gz

dist:
	tar -zcf $(ARCHIVE).tar.gz $(SDIR) Makefile

build : $(BDIR)/overlap.o $(BDIR)/String.o $(BDIR)/common.o
	$(CC) $(BDIR)/overlap.o $(BDIR)/String.o $(BDIR)/common.o -o $(BDIR)/overlap
