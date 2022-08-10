.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor LDT_users_guide.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf -r asciidoctor-mathematical LDT_users_guide.adoc

.PHONY: clean
clean:
	rm -f LDT_users_guide.html LDT_users_guide.pdf
