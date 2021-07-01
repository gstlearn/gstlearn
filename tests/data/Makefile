TOP = ../..
include $(TOP)/config/Makefile.inc
HERE   = $(TOP)/tests/data
LIST = $(wildcard Jeu* PGS* PRD* Var*)

all:
	@echo "Running Non-regression tests...";
	@for VAR in $(LIST); do \
		make -s $$VAR ; \
	done

Jeu%: target_force
	@echo "---> Non-regression test for $@";
	@if [ -f $@/README ]; then cat $@/README; fi;
	@$(EXEDIR)/testim.exe $@ > $@/Result.out;
	@$(HERE)/compare.com $@ 

PGS%: target_force
	@echo "---> Non-regression test for $@";
	@if [ -f $@/README ]; then cat $@/README; fi;
	@$(EXEDIR)/testpgs.exe $@ > $@/Result.out;
	@$(HERE)/compare.com $@ 

PRD%: target_force
	@echo "---> Non-regression test for $@";
	@if [ -f $@/README ]; then cat $@/README; fi;
	@$(EXEDIR)/testprd.exe $@ > $@/Result.out;
	@$(HERE)/compare.com $@ 

Var%: target_force
	@echo "---> Non-regression test for $@";
	@if [ -f $@/README ]; then cat $@/README; fi;
	@$(EXEDIR)/testvar.exe $@ > $@/Result.out;
	@$(HERE)/compare.com $@ 

# Force target to be made
target_force:

validate:
	@for VAR in $(LIST); do \
		echo "---> Record results files as new Reference in $$VAR"; \
		$(HERE)/update.com $$VAR Result; \
		$(HERE)/update.com $$VAR Vario; \
		$(HERE)/update.com $$VAR Model; \
	done

update:
	@for VAR in $(LIST); do \
		echo "---> Update ASCII files for test for $$VAR"; \
		$(EXEDIR)/testupd.exe $$VAR; \
	done

clean:
	@echo "Cleaning tests/data ..."
	@for VAR in $(LIST); do \
		rm -rf $$VAR/*.out; \
		rm -rf __pycache__; \
		rm -rf .RData; \
		rm -rf .Rhistory; \
	done

