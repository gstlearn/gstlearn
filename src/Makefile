TOP = ..
include $(TOP)/config/Makefile.inc

SRCPP   = $(wildcard *.cpp) $(wildcard */*.cpp)
OBJCPP  = $(addprefix $(OBJDIR)/,$(SRCPP:.cpp=.o))

.PHONY: all clean
 
all: $(LA) $(LSO)

$(LA): $(OBJCPP)
	@echo "Making gstlearn [static $(LIB_MODE)] ..."
	@$(AR) $(LA) $(OBJCPP)

$(LSO): $(OBJCPP)
	@echo "Making gstlearn [shared $(LIB_MODE)] ..."
	@$(CXX) -shared $(LDFLAGS) $(INCLUDE) -o $(LSO) $(OBJCPP) $(LIBS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	@echo "File: $<"
	@$(CXX) -c -o $(OBJDIR)/$*.o  $(CXXFLAGS) $(INCLUDE) $<

clean:
	@echo "Cleaning src ..."
	@rm -rf $(LIBDIR)
	@rm -rf $(BINDIR)
	@rm -f *~
	@rm -f Makefile.bak

