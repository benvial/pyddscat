
all: fortran wrappers move

fortran:
	cd ddscat && make build

wrappers:
	cd ddscat && \
	python gen_wrapper.py ddscat && \
	python gen_wrapper.py ddpostprocess && \
	python gen_wrapper.py vtrconvert && \
	python gen_wrapper.py calltarget

move:
	mkdir -p bin
	mv ddscat/src/ddscat bin/ddscat_bin
	mv ddscat/ddscat bin/
	mv ddscat/src/ddpostprocess bin/ddpostprocess_bin 
	mv ddscat/ddpostprocess bin/
	mv ddscat/src/vtrconvert bin/vtrconvert_bin 
	mv ddscat/vtrconvert bin/
	mv ddscat/src/calltarget bin/calltarget_bin 
	mv ddscat/calltarget bin/
