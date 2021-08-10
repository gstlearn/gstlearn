
.PHONY: all gstlearn tests testsCur total doxygen clean check-arch

ARCH_XX = XX_$(ARCH)_XX
WIN_TYPE_XX = XX_$(WIN_TYPE)_XX

all: gstlearn

gstlearn: check-arch
	@+make -sC src

tests: gstlearn
	@+make -sC tests $@
	
testsCur: gstlearn
	@+make -sC testsCur $@

total: gstlearn
	@+make -sC tests $@
	
doxygen:
	@+make -C doxygen

clean:
	@echo "Cleaning gstlearn ..."
	@+make -sC src     clean
	@+make -sC tests   clean
	@+make -sC doxygen clean

check-arch:
ifeq ($(ARCH_XX),XX__XX)
	$(error ARCH environment variable is undefined!)
endif
ifeq ($(ARCH), windows)
ifeq ($(WIN_TYPE_XX),XX__XX)
	$(error WIN_TYPE variable (32 or 64) is undefined!)
endif
endif
