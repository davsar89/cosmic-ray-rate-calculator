PYTHON ?= python
PDFLATEX ?= pdflatex
SOURCE_DATE_EPOCH ?= 1773333791
export SOURCE_DATE_EPOCH

.PHONY: test generate pdf artifacts clean

test:
	$(PYTHON) -m pytest

generate:
	$(PYTHON) generate_artifacts.py

pdf: generate
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error cosmic_ray_rate_formula.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error cosmic_ray_rate_formula.tex

artifacts: pdf

clean:
	rm -f \
		cosmic_ray_rate_formula.aux \
		cosmic_ray_rate_formula.log \
		cosmic_ray_rate_formula.out \
		cosmic_ray_rate_formula.synctex.gz
