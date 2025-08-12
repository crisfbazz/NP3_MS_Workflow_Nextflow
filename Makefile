.PHONY : clean run all

all: clean run

clean:
	rm -Rf .nextflow* work bin/*__pycache__* nf_output/*

run:
    # remove previous default result and execute the pre_process+run commands - default process
	rm -rf nf_output/np3_results/L754_test_nf
	nextflow run ./nf_workflow.nf -resume -c nextflow.config

init_modules:
	git submodule update --init --recursive
