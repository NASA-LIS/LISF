.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor LIS_users_guide.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf LIS_users_guide.adoc

.PHONY: clean
clean:
	rm -f LIS_users_guide.html LIS_users_guide.pdf
