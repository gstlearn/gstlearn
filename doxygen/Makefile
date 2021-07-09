TOP = ..
include $(TOP)/config/Makefile.inc

all: doxygen

doxygen:
	doxygen Doxyfile
	$(CHMODCMD) html
	
clean:
	@echo "Cleaning doxygen ..."
	@rm -rf html
	@rm -rf xml
	@rm -rf *.tmp
