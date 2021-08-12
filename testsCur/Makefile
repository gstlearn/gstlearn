
.PHONY: all testsCur total clean check-arch

ARCH_XX = XX_$(ARCH)_XX
WIN_TYPE_XX = XX_$(WIN_TYPE)_XX

all: testsCur

testsCur: check-arch
	@+make -sC cpp

total: testsCur
	@+make -sC data

clean:
	@echo "Cleaning testsCur ..."
	@+make -sC cpp  clean
	@+make -sC data clean

check-arch:
ifeq ($(ARCH_XX),XX__XX)
	$(error ARCH environment variable is undefined!)
endif
ifeq ($(ARCH), windows)
ifeq ($(WIN_TYPE_XX),XX__XX)
	$(error WIN_TYPE variable (32 or 64) is undefined!)
endif
endif
