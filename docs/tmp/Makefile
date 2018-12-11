.DEFAULT_GOAL := html

RM := rm -f
make-clean := $(RM) howto-create-testcases.pptx howto-create-testcases.pdf howto-create-testcases.html

.PHONY: slides
slides: howto-create-testcases.pptx

.PHONY: pdf
pdf: howto-create-testcases.pdf

.PHONY: html
html: howto-create-testcases.html

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
	@echo "   generate html version of howto-create-testcases.adoc"
	@echo "   DEFAULT target"
	@echo ""
	@echo "pdf"
	@echo "   generate pdf version of howto-create-testcases.adoc"
	@echo ""
	@echo "slides"
	@echo "   generate pptx version of howto-create-testcases.adoc"
	@echo ""
