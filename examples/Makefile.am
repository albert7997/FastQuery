AUTOMAKE_OPTIONS=gnu
EXTRA_DIST = read.cpp writeData.cpp writeAttr.cpp queryIndex.cpp \
buildIndex.cpp histogram.cpp BPread.cpp BPbuildIndex.cpp BPqueryIndex.cpp \
MagNullBuildIndex.cpp MagNullQueryIndex.cpp MagNullQueryScan.cpp \
s3dregions.cpp s3dstream.cpp h5uc-data.h5 test.dat tos_O1_2001-2002.nc

EXAMPLES = read writeData writeAttr queryIndex buildIndex histogram
if BUILD_BP
EXAMPLES += BPread BPbuildIndex BPqueryIndex MagNullBuildIndex MagNullQueryIndex MagNullQueryScan# s3dregions s3dstream
endif

examples: $(EXAMPLES)
nompi: clean-am
	$(MAKE) FQ_MOD=_nompi CXXFLAGS="$(CXXFLAGS) -DFQ_NOMPI" $(EXAMPLES)

# Use explicit make rules to bypass libtool and avoid problems on Cray XE6 (hopper)
read$(EXEEXT): read.cpp ../src/libfastquery$(FQ_MOD).la
	$(CXX) -o $@ $(CXXFLAGS) -I../src read.cpp -L../src/$(OBJDIR) -lfastquery$(FQ_MOD) $(LIBS)

writeData$(EXEEXT): writeData.cpp ../src/libfastquery$(FQ_MOD).la
	$(CXX) -o $@ $(CXXFLAGS) -I../src writeData.cpp -L../src/$(OBJDIR) -lfastquery$(FQ_MOD) $(LIBS)

writeAttr$(EXEEXT): writeAttr.cpp ../src/libfastquery$(FQ_MOD).la
	$(CXX) -o $@ $(CXXFLAGS) -I../src writeAttr.cpp -L../src/$(OBJDIR) -lfastquery$(FQ_MOD) $(LIBS)

queryIndex$(EXEEXT): queryIndex.cpp ../src/libfastquery$(FQ_MOD).la
	$(CXX) -o $@ $(CXXFLAGS) -I../src queryIndex.cpp -L../src/$(OBJDIR) -lfastquery$(FQ_MOD) $(LIBS)

buildIndex$(EXEEXT): buildIndex.cpp ../src/libfastquery$(FQ_MOD).la
	$(CXX) -o $@ $(CXXFLAGS) -I../src buildIndex.cpp -L../src/$(OBJDIR) -lfastquery$(FQ_MOD) $(LIBS)

histogram$(EXEEXT): histogram.cpp ../src/libfastquery$(FQ_MOD).la
	$(CXX) -o $@ $(CXXFLAGS) -I../src histogram.cpp -L../src/$(OBJDIR) -lfastquery$(FQ_MOD) $(LIBS)

#$(EXEEXT): %.cpp ../src/libfastquery$(FQ_MOD).la
#	$(CXX) -o $@ $(CXXFLAGS) -I../src $< -L../src/$(OBJDIR) -lfastquery$(FQ_MOD) $(LIBS)

../src/libfastquery$(FQ_MOD).la::
	cd ../src && $(MAKE) libfastquery$(FQ_MOD).la

# has to be linked with $(CXXLINK)
#tcapi$(EXEEXT): tcapi.$(OBJEXT) $(tcapi_LDADD)
#        $(CXXLINK) tcapi.$(OBJEXT) $(tcapi_LDADD) $(LIBS)

clean-local:
	-rm -rf *.dSYM $(EXAMPLES)

distclean: distclean-am clean
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile

.PHONY: nompi clean-am

