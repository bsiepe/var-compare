# Reproducibility targets for Siepe, Kloft & Heck (2024)
#
# Prerequisites:
#   - R 4.3.2 with all packages installed (see Dockerfile), or
#   - Docker: `make docker-build` then `make docker-run`
#
# Typical usage:
#   make all              # run full pipeline
#   make simulation-1     # render simulation study 1 only
#
# Note: simulation-1.Rmd and simulation-2.Rmd each contain both the
# simulation loop and the visualizations in the same script. The full
# render will execute both. For simulation study 1, intermediate results
# (output/df_eval_bggm_0705.RDS, output/df_eval_gvar_0705.RDS) are
# already provided and can be used to reproduce figures without rerunning
# the simulation. To do this, open simulation-1.Rmd (or
# revision1/simulation-1-revision.Rmd) in RStudio and run only the
# chunks after the section "IF RESULTS PRELOADED".

RSCRIPT = Rscript
RENDER  = $(RSCRIPT) -e "rmarkdown::render('$<')"

# Output sentinels (touch files to track what has been rendered)
SENTINEL_DIR = .make
$(SENTINEL_DIR):
	mkdir -p $(SENTINEL_DIR)

# ------------------------------------------------------------------ #
# Main targets                                                         #
# ------------------------------------------------------------------ #

.PHONY: all main revision clean docker-build docker-run

all: main revision

main: $(SENTINEL_DIR)/true-networks \
      $(SENTINEL_DIR)/simulation-1 \
      $(SENTINEL_DIR)/simulation-2 \
      $(SENTINEL_DIR)/empirical-example

revision: $(SENTINEL_DIR)/revision1/simulation-1-revision \
          $(SENTINEL_DIR)/revision1/simulation-2-revision

# ------------------------------------------------------------------ #
# Main scripts (order matters: true-networks must run first)          #
# ------------------------------------------------------------------ #

$(SENTINEL_DIR)/true-networks: true-networks.Rmd aux_funs.R | $(SENTINEL_DIR)
	$(RSCRIPT) -e "rmarkdown::render('true-networks.Rmd')"
	touch $@

$(SENTINEL_DIR)/simulation-1: simulation-1.Rmd aux_funs.R data/l_graphs.RDS | $(SENTINEL_DIR)
	$(RSCRIPT) -e "rmarkdown::render('simulation-1.Rmd')"
	touch $@

$(SENTINEL_DIR)/simulation-2: simulation-2.Rmd aux_funs.R data/l_graphs.RDS | $(SENTINEL_DIR)
	$(RSCRIPT) -e "rmarkdown::render('simulation-2.Rmd')"
	touch $@

$(SENTINEL_DIR)/empirical-example: empirical-example.Rmd aux_funs.R | $(SENTINEL_DIR)
	$(RSCRIPT) -e "rmarkdown::render('empirical-example.Rmd')"
	touch $@

# ------------------------------------------------------------------ #
# Revision scripts                                                     #
# ------------------------------------------------------------------ #

$(SENTINEL_DIR)/revision1/simulation-1-revision: revision1/simulation-1-revision.Rmd aux_funs.R data/l_graphs.RDS | $(SENTINEL_DIR)
	mkdir -p $(SENTINEL_DIR)/revision1
	$(RSCRIPT) -e "rmarkdown::render('revision1/simulation-1-revision.Rmd')"
	touch $@

$(SENTINEL_DIR)/revision1/simulation-2-revision: revision1/simulation-2-revision.Rmd aux_funs.R data/l_graphs.RDS | $(SENTINEL_DIR)
	mkdir -p $(SENTINEL_DIR)/revision1
	$(RSCRIPT) -e "rmarkdown::render('revision1/simulation-2-revision.Rmd')"
	touch $@

# ------------------------------------------------------------------ #
# Individual convenience targets                                       #
# ------------------------------------------------------------------ #

.PHONY: true-networks simulation-1 simulation-2 empirical-example
true-networks:    $(SENTINEL_DIR)/true-networks
simulation-1:     $(SENTINEL_DIR)/simulation-1
simulation-2:     $(SENTINEL_DIR)/simulation-2
empirical-example: $(SENTINEL_DIR)/empirical-example

# ------------------------------------------------------------------ #
# Docker                                                               #
# ------------------------------------------------------------------ #

docker-build:
	docker build -t var-compare:4.3.2 .

# Mounts the project directory so output files are written to the host.
docker-run:
	docker run --rm -v "$$(pwd)":/project var-compare:4.3.2 \
	    Rscript -e "setwd('/project'); source('Makefile')"

# Run the full pipeline inside Docker via make
docker-make:
	docker run --rm -v "$$(pwd)":/project var-compare:4.3.2 \
	    bash -c "cd /project && make all"

# ------------------------------------------------------------------ #
# Clean sentinel files (does not delete output or figures)            #
# ------------------------------------------------------------------ #

clean:
	rm -rf $(SENTINEL_DIR)
