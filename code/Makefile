# Include the configuration file
include ../config/config.mk

# Lists of files to be built
objs=dopri54.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=lv_test lvs_test param_est

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(cflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $< $(lp_lflags)

lv_test: lv_test.cc dopri54.o odes.hh
	$(cxx) $(cflags) -o $@ $^

lvs_test: lvs_test.cc dopri54.o odes.hh
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

param_est: param_est.cc dopri54.o odes.hh
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)
