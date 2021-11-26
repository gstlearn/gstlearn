# Minimum requirement for compiling gstlearn:
# g++/clang/make (dev), Rtools (windows only)
# git (client), Boost (dev), NetCDF (dev), doxygen, hdf5

AR       = ar r
CXX      = g++
#CXXFLAGS = -Wall -std=c++11 -Wno-array-bounds
#CXXFLAGS = -Wall -Wextra -std=c++11
CXXFLAGS = -Wall -std=c++11
LDFLAGS  =
INCLUDE  = -I$(TOP)/include -I.
LIBDIR   = $(TOP)/lib/$(ARCH)$(WIN_TYPE)

ifeq ($(DEBUG), 1)
  CXXFLAGS += -g
  TARGET    = gstlearnd
  OBJDIR    = $(LIBDIR)/$(TARGET)
  LIB_MODE  = debug
else
  CXXFLAGS += -O2
  TARGET    = gstlearn
  OBJDIR    = $(LIBDIR)/$(TARGET)
  LIB_MODE  = release
endif

BINDIR = $(TOP)/bin/$(ARCH)$(WIN_TYPE)
EXEDIR = $(BINDIR)/$(LIB_MODE)

ifeq ($(ARCH),linux)
     CXXFLAGS += -fPIC -D_REENTRANT
     LDFLAGS  += -z defs -Wl,--hash-style=both 
     LSO       = $(LIBDIR)/lib$(TARGET).so
     LA        = $(LIBDIR)/$(TARGET).a
     LIBS      = -lm -lboost_system
     CHMODCMD  = chmod 777
endif

ifeq ($(ARCH),linux64)
     CXXFLAGS += -fPIC -D_REENTRANT
     LDFLAGS  += -z defs -Wl,--hash-style=both
     LSO       = $(LIBDIR)/lib$(TARGET).so
     LA        = $(LIBDIR)/$(TARGET).a
     LIBS      = -lm -lboost_system
     CHMODCMD  = chmod 777
endif

ifeq ($(ARCH),macosx)
     CXX       = clang++
     INCLUDE  += -I$(BOOST_DIR)
     CXXFLAGS += -fPIC
     LSO       = $(LIBDIR)/lib$(TARGET).dylib
     LA        = $(LIBDIR)/$(TARGET).a
     LIBS      = -lm -framework CoreFoundation -framework IOKit
     CHMODCMD  = chmod 777
endif

ifeq ($(ARCH),windows)
     INCLUDE  += -I$(BOOST_DIR)
     CXXFLAGS += -m$(WIN_TYPE) -DNOMINMAX -DNO_TIMER
     LSO       = $(LIBDIR)/$(TARGET).dll
     LA        = $(LIBDIR)/$(TARGET).lib
     LIBS      = -liphlpapi -lrpcrt4
     CHMODCMD  = echo 
endif

ifeq ($(USE_NETCDF),1)
     NETCDF_SYM= -D_USE_NETCDF
     CXXFLAGS += -D_USE_NETCDF
# TODO : To be kept ?
     INCLUDE  += -I/usr/local/include
     LDFLAGS  += -L/usr/local/lib
     LIBS     += -lnetcdf_c++4 -lnetcdf
endif

# This is only for LINUX. To be improved for Windows
CXXFLAGS += -I/usr/include/hdf5/serial
LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_cpp
#LIBS += -lhdf5 -lhdf5_cpp
