# Usage: make help
# Everything runs with the Python of the conda env $(ENV); no `conda activate` needed.

ENV ?= sage
ifndef PYTHON
PYTHON := $(shell conda run -n $(ENV) python -c "import sys; print(sys.executable)" 2>/dev/null)
endif
FILE ?= examples/x4plus1_mod3.txt
LATEX := pdflatex -interaction=nonstopmode -halt-on-error

.DEFAULT_GOAL := help
.PHONY: help env check example run test bench notebook demo report slides docs all clean

help:
	@echo "make env       create (or update) the conda env '$(ENV)' from environment.yml"
	@echo "make check     confirm the env works: Python, Sage, pytest, Jupyter, pdflatex"
	@echo "make example   factor x^4 + 1 mod 3 with the full step-by-step trace"
	@echo "make run FILE=examples/hw3_sample2.txt   factor any HW3-format input file"
	@echo "make test      run the 66 tests (about 1 minute)"
	@echo "make bench     regenerate bench/results/ tables and plots (about 1 minute)"
	@echo "make notebook  re-execute demo.ipynb in place"
	@echo "make demo      open demo.ipynb in Jupyter Lab"
	@echo "make report    build report/report.pdf"
	@echo "make slides    build slides/slides.pdf"
	@echo "make all       test, bench, notebook, report and slides"
	@echo "make clean     remove caches and LaTeX build files"

env:
	@if conda env list | grep -qE '^$(ENV)[[:space:]]'; then \
		conda env update -n $(ENV) -f environment.yml --prune; \
	else \
		conda env create -n $(ENV) -f environment.yml; \
	fi

check: _python
	@$(PYTHON) -c "import sys, sage.version; print('Python', sys.version.split()[0], '| SageMath', sage.version.version)"
	@$(PYTHON) -m pytest --version
	@$(PYTHON) -m jupyter --version | grep -E "^(jupyterlab|nbconvert|ipywidgets)"
	@command -v pdflatex >/dev/null && echo "pdflatex found" || echo "pdflatex missing: needed only for report and slides (see README)"

example: _python
	$(PYTHON) -m kersplit --trace < $(FILE)

run: _python
	$(PYTHON) -m kersplit < $(FILE)

test: _python
	$(PYTHON) -m pytest

bench: _python
	$(PYTHON) -m bench.bench

notebook: _python
	$(PYTHON) -m jupyter nbconvert --to notebook --execute --inplace demo.ipynb

demo: _python
	$(PYTHON) -m jupyter lab demo.ipynb

report: report/report.pdf
slides: slides/slides.pdf
docs: report slides

report/report.pdf: report/report.tex bench/results/time_vs_p.png
	cd report && $(LATEX) report.tex >/dev/null && $(LATEX) report.tex >/dev/null
	@echo "built $@"

slides/slides.pdf: slides/slides.tex bench/results/time_vs_p.png
	cd slides && $(LATEX) slides.tex >/dev/null && $(LATEX) slides.tex >/dev/null
	@echo "built $@"

all: test bench notebook docs

clean:
	rm -rf .pytest_cache .ipynb_checkpoints
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
	rm -f report/*.aux report/*.log report/*.out
	rm -f slides/*.aux slides/*.log slides/*.out slides/*.nav slides/*.snm slides/*.toc

.PHONY: _python
_python:
	@test -n "$(PYTHON)" || { echo "conda env '$(ENV)' not found: run 'make env' first"; exit 1; }
