# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
examples: $(EXTRA_PROGRAMS)
AM_CPPFLAGS=-I$(top_builddir)
AM_CXXFLAGS = @DEFAULT_CFLAGS@
AM_CPPFLAGS += -I$(top_srcdir)/src/kernel
AM_CPPFLAGS += $(GMP_CFLAGS)

LDADD = $(top_builddir)/src/libgivaro.la $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS=-static


EXTRA_PROGRAMS=recint-iterator rsa extended-int-types
CLEANFILES=$(EXTRA_PROGRAMS)

recint_iterator_SOURCES  = recint-iterator.C
extended_int_types_SOURCES  = extended-int-types.C
rsa_SOURCES       = rsa.C

#external flags
AM_CPPFLAGS += $(OPTFLAGS) 
LDADD += $(OPTLIBES)

# for compilation of new examples
GIVARO_BIN=@bindir@
define comp_new_examp
	$(AM_V_CXX)$(CXXCOMPILE) -c -o $@.$(OBJEXT) $<
	$(AM_V_CXXLD)$(CXXLINK) $@.$(OBJEXT) $(LDADD)
endef

%:%.C
	$(comp_new_examp)

%:%.cpp
	$(comp_new_examp)
