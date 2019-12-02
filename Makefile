ifndef ROOTDIR
ROOTDIR=.
endif

REALMAKEFILE=../Makefile.2

default: mc       

all: mc                      

serialobjs: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) serialobjs)


mc: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) mc)

clean:
	cd $(ROOTDIR) && rm -f *~
	cd $(ROOTDIR) && rm -f mc.x  
	cd $(ROOTDIR) && rm -f src/*~
	@( cd $(ROOTDIR) && if [ -d src/obj ] ; \
		then cd src/obj && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf obj ; \
	fi )
	@( cd $(ROOTDIR) && if [ -d src/objpt ] ; \
		then cd src/objpt && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf objpt ; \
	fi )

objdir: 
	@( cd $(ROOTDIR) && if [ ! -d src/obj ] ; \
		then mkdir src/obj ; \
	fi ) ;

.PHONY: mc default all clean objdir objdirpt serialobjs
