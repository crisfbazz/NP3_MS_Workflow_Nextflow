#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// NP3 submodule folder
NP3_TOOL_FOLDER = "$moduleDir/bin/NP3_MS_Workflow/"

// np3 ms workflow default parms
params.np3_command = "run"
// mandatory parms
params.output_name = "L754_test_nf"
params.metadata = "$NP3_TOOL_FOLDER/test/L754_bacs/marine_bacteria_library_L754_metadata.csv"
params.raw_data_path = "$NP3_TOOL_FOLDER/test/L754_bacs/mzxml/"
params.pre_processed_data_name = "processed_data_test"
// critical parms run
params.mz_tolerance = "0.025"
params.fragment_tolerance = "0.05"
params.ion_mode = "1"
// critical parms pre_process
params.peak_width = "2,10"
params.rt_tolerance_deviation = "3.0"
params.mz_tolerance_deviation = "0.05"
params.ppm_tolerance = "15"
// Other parms for run
params.rt_tolerance = "1,2"
params.similarity_function = "np3_shifted_cosine"
params.trim_mz = "TRUE"
params.noise_cutoff = "FALSE"
// MN parms
params.similarity_mn = "0.6"
params.net_top_k = "15"
params.max_component_size = "200"
params.min_matched_peaks = "6"
// command gnps_result parms
params.result_library_search_path = "$NP3_TOOL_FOLDER/test/L754_bacs/ProteoSAFe-MOLECULAR-LIBRARYSEARCH-V2-da67f38d-download_all_identifications/MOLECULAR-LIBRARYSEARCH-V2-da67f38d-download_all_identifications-main.tsv"
params.np3_count_table = "$moduleDir/nf_output/np3_results/run/L754_test_nf/outs/L754_test_nf/count_tables/clean/L754_test_nf_spectra_clean_ann_corr_spearman.csv"

// np3 fixed params
params.output_path = "./np3_results/run/"
params.rules_ionization = "$NP3_TOOL_FOLDER/rules/np3_modifications.csv"
params.pre_processed_output_path = "./np3_results/pre_processed/"
params.gnps_result_output_path = "./np3_results/gnps_result_join/"

//This publish dir is mostly  useful when we want to import modules in other workflows, keep it here usually don't change it
params.publishdir = "$launchDir"
TOOL_FOLDER = "$moduleDir/bin"

// Augmenting with nf_output
_publishdir = "${params.publishdir}/nf_output"

process np3PreProcess {
    /* This process executes the *pre_process* command from np3*/

    publishDir "$_publishdir", mode: 'copy', overwrite: false

    conda "$TOOL_FOLDER/environment_np3_nextflow_unix.yml"

    input:
    val output_name
    val pre_processed_output_path
    val metadata
    val raw_data_path
    val pre_processed_data_name
    val mz_tolerance_deviation
    val rt_tolerance_deviation
    val ppm_tolerance
    val peak_width
    val ion_mode

    output:
    path "$pre_processed_output_path/${pre_processed_data_name}/", emit: output_pre_process
    file "$pre_processed_output_path/${pre_processed_data_name}.zip"

    script:
    """
    node $NP3_TOOL_FOLDER/np3_workflow.js pre_process \
    --metadata $metadata --data_name $output_name --raw_data_path $raw_data_path \
    --processed_data_name $pre_processed_data_name --mz_tolerance $mz_tolerance_deviation \
    --ion_mode $ion_mode --rt_tolerance $rt_tolerance_deviation \
    --ppm_tolerance $ppm_tolerance --peak_width $peak_width --processed_data_overwrite FALSE \
    --verbose 1
    # check if the pre process result exists and if yes create the pre_processed_output_path and then
    # copy the pre process result to the current work space in the pre_processed_output_path
    [ -d '$raw_data_path/$pre_processed_data_name' ] && mkdir -p '$pre_processed_output_path/' && \
        cp -r '$raw_data_path/$pre_processed_data_name' '$pre_processed_output_path/'
    [ -d '$pre_processed_output_path/${pre_processed_data_name}/' ] && cd $pre_processed_output_path/ && \
        zip -r ${pre_processed_data_name}.zip ${pre_processed_data_name}/ && cd -
    """
}

process np3Run {
    /* This process executes the *run* command from np3*/

    publishDir "$_publishdir", mode: 'copy', overwrite: false

    conda "$TOOL_FOLDER/environment_np3_nextflow_unix.yml"

    input:
    val output_name
    val output_path
    val metadata
    val raw_data_path
    val pre_processed_data_name
    val mz_tolerance
    val fragment_tolerance
    val ion_mode
    val rt_tolerance
    val similarity_function
    val trim_mz
    val noise_cutoff
    val similarity_mn
    val net_top_k
    val max_component_size
    val min_matched_peaks
    val rules_ionization

    output:
    path "$output_path/${output_name}/"
    file "$output_path/${output_name}.zip"

    script:
    """
    node $NP3_TOOL_FOLDER/np3_workflow.js run \
    --metadata $metadata --output_name $output_name \
    --output_path $output_path --raw_data_path $raw_data_path \
    --processed_data_name $pre_processed_data_name --mz_tolerance $mz_tolerance \
    --fragment_tolerance $fragment_tolerance --ion_mode $ion_mode --rt_tolerance $rt_tolerance \
    --similarity_function $similarity_function --trim_mz $trim_mz --noise_cutoff $noise_cutoff \
    --similarity_mn $similarity_mn --net_top_k $net_top_k --max_component_size $max_component_size \
    --min_matched_peaks $min_matched_peaks --rules $rules_ionization --verbose 1
    cd $output_path && zip -r ${output_name}.zip ${output_name}/ && cd -
    """
}

process np3GNPSResult {
    /* This process executes the *gnps_result* command from np3*/

    publishDir "$_publishdir", mode: 'copy', overwrite: false

    conda "$TOOL_FOLDER/environment_np3_nextflow_unix.yml"

    input:
    path np3_count_table
    val result_library_search_path
    val gnps_result_output_path

    output:
    path "$gnps_result_output_path/$np3_count_table.name"

    script:
    """
    node $NP3_TOOL_FOLDER/np3_workflow.js gnps_result \
    --result_specnets_DB_path $result_library_search_path \
    --count_file_path $np3_count_table
    mkdir -p '$gnps_result_output_path/'
    cp $np3_count_table $gnps_result_output_path/
    """
}

workflow  {
    /*
    The NP3 workflow that executes the pre_process and run command
    */

    // path inputs
    input_metadata = Channel.fromPath(params.metadata)
    input_raw_data_path = Channel.fromPath(params.raw_data_path)
    input_rules_ionization = Channel.fromPath(params.rules_ionization)
    // call the *pre_process* nf process here, then call the *run* nf process
    np3PreProcess(params.output_name, params.pre_processed_output_path, input_metadata, \
    input_raw_data_path, params.pre_processed_data_name, \
    params.mz_tolerance_deviation, params.rt_tolerance_deviation, params.ppm_tolerance, \
    params.peak_width, params.ion_mode)
    // call *run* nf process using the *pre_process* result - sequential calls
    // passing np3PreProcess.out.output_pre_process (which is equal to  params.pre_processed_output_path/params.pre_processed_data_name) instead of params.raw_data_path
    // passing "." instead of params.pre_processed_data_name because the pre_processed_data_name is already included in the output_pre_processed path and was copied to the given output in the np3PreProcess script
    np3Run(params.output_name, params.output_path, input_metadata, \
    np3PreProcess.out.output_pre_process, ".", params.mz_tolerance, \
    params.fragment_tolerance, params.ion_mode, params.rt_tolerance, params.similarity_function, \
    params.trim_mz, params.noise_cutoff, params.similarity_mn, params.net_top_k, \
    params.max_component_size, params.min_matched_peaks, input_rules_ionization)
    // TODO call library search here and then call gnps_result join
    //input_np3_count_table = Channel.fromPath(params.np3_count_table)
    //input_result_library_search_path = Channel.fromPath(params.result_library_search_path)
    // TODO call *gnps_result* process and update the resulting run zip file
    //np3GNPSResult(input_np3_count_table, input_result_library_search_path, params.gnps_result_output_path)
}
