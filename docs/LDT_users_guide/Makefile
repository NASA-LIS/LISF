.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor LDT_usersguide.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf LDT_usersguide.adoc

.PHONY: clean
clean:
	rm -f LDT_usersguide.html LDT_usersguide.pdf
