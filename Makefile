bnlearn:
	R CMD build .

.PHONY: install check clean dependency

install:
	R CMD INSTALL --preclean --no-multiarch --with-keep.source .

check:
	R CMD check bnlearn_*.tar.gz

clean:
	rm -rf *.Rcheck bnlearn_*.tar.gz

dependency:
	xargs apt-get -y install < pkglist
