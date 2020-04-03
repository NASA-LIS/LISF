.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor LIS_usersguide.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf LIS_usersguide.adoc

.PHONY: clean
clean:
	rm -f LIS_usersguide.html LIS_usersguide.pdf
