.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor creating-a-custom-modulefile.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf creating-a-custom-modulefile.adoc

.PHONY: clean
clean:
	rm -f creating-a-custom-modulefile.html creating-a-custom-modulefile.pdf
