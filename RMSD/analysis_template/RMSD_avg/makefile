include common.mk

all: figure/rmsd.png

clean:
	rm -rf data
	rm -rf figure

figure/rmsd.png: $(strip ${SCRIPT_DIR})/plot_rmsd.py data/rmsd.csv
	python "$<" \
        --path-to-rmsd_dir=data/ \
        --output-directory=figure/

data/rmsd.csv: $(strip ${SCRIPT_DIR})/average_rmsd.py
	python "$<" \
        --path-to-processing-dirs="$(strip ${PATH_TO_PROCESSING_DIRS})" \
        --output-directory=data/