.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor LVT_usersguide.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf LVT_usersguide.adoc

.PHONY: clean
clean:
	rm -f LVT_usersguide.html LVT_usersguide.pdf
