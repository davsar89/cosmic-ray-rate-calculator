PYTHON ?= python
PDFLATEX ?= pdflatex

.PHONY: test generate fragments pdf artifacts clean

test:
	$(PYTHON) -m pytest

generate:
	$(PYTHON) generate_artifacts.py

fragments:
	$(PYTHON) generate_artifacts.py --skip-plot

pdf: fragments
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error cosmic_ray_rate_formula.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error cosmic_ray_rate_formula.tex

artifacts: pdf
	$(PYTHON) generate_artifacts.py

clean:
	rm -f \
		cosmic_ray_rate_formula.aux \
		cosmic_ray_rate_formula.log \
		cosmic_ray_rate_formula.out \
		cosmic_ray_rate_formula.synctex.gz
