.DEFAULT_GOAL := html

RM := rm -f
make-clean := $(RM) working_with_github.pptx working_with_github.pdf working_with_github.html

.PHONY: slides
slides: working_with_github.pptx

.PHONY: pdf
pdf: working_with_github.pdf

.PHONY: html
html: working_with_github.html

.PHONY: all
all: html pdf slides

.PHONY: clean
clean:
	$(call make-clean)

%.pptx: %.xml
	pandoc -f docbook --slide-level 1 -t pptx $< -o $@

%.xml: %.adoc
	asciidoctor -b docbook5 $<

%.pdf: %.adoc
	asciidoctor-pdf $<

%.html: %.adoc
	asciidoctor $<


.PHONY: help
help:
	@echo ""
	@echo "usage: make <target>"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "clean"
	@echo "   $(make-clean)"
	@echo ""
	@echo "html"
	@echo "   generate html version of working_with_github.adoc"
	@echo "   DEFAULT target"
	@echo ""
	@echo "pdf"
	@echo "   generate pdf version of working_with_github.adoc"
	@echo ""
	@echo "slides"
	@echo "   generate pptx version of working_with_github.adoc"
	@echo ""
