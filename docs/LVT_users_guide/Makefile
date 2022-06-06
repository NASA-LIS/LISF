.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor LVT_users_guide.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf -r asciidoctor-mathematical LVT_users_guide.adoc

.PHONY: clean
clean:
	rm -f LVT_users_guide.html LVT_users_guide.pdf
