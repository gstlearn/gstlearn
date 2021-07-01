TOP = ../..
include $(TOP)/config/Makefile.inc

SRCPP   = $(wildcard *.cpp) $(wildcard */*.cpp)
EXECPP  = $(addprefix $(EXEDIR)/,$(SRCPP:.cpp=.exe))

.PHONY: all clean

all: $(EXECPP)

$(EXEDIR)/%.exe:%.cpp $(LA)
	@mkdir -p $(dir $@)
	@echo "File: $<"
	@$(CXX) -o $(EXEDIR)/$*.exe $(CXXFLAGS) $(LDFLAGS) $(INCLUDE) $< $(LA) $(LIBS)

$(LA):
	@make -sC $(TOP)/src $@

clean:
	@echo "Cleaning tests/cpp ..."
	rm -rf $(TOP)/$(BINDIR)
	@rm -f *~
	@rm -f Makefile.bak

