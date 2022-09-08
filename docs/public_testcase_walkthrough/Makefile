.PHONY: all
all: html pdf

.PHONY: html
html:
	asciidoctor public_testcase_walkthrough.adoc

.PHONY: pdf
pdf:
	asciidoctor-pdf public_testcase_walkthrough.adoc

.PHONY: clean
clean:
	rm -f public_testcase_walkthrough.html public_testcase_walkthrough.pdf
