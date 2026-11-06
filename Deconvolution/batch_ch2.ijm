dataDir="lightsheet_raw_data"


minx="516";
miny="-72";
minz="112";
maxx="1536";
maxy="1842";
maxz="1132";
excitation="642";
itenum="16";

passedArgument = call("java.lang.System.getenv", "SLURM_ARRAY_TASK_ID");

options="select_xml="+dataDir+"dataset.xml process_angle=[All angles] process_channel=[Single channel (Select from List)] process_illumination=[All illuminations] process_timepoint=[Multiple Timepoints (Select from List)] processing_channel=[channel 1] timepoint_"+passedArgument+" type_of_image_fusion=[Multi-view deconvolution] bounding_box=[Define manually] fused_image=[Save as TIFF stack] minimal_x="+minx+" minimal_y="+miny+" minimal_z="+minz+" maximal_x="+maxx+" maximal_y="+maxy+" maximal_z="+maxz+" imglib2_container=[CellImg (large images)] imglib2_container_ffts=[CellImg (large images)] save_memory type_of_iteration=[Efficient Bayesian - Optimization I (fast, precise)] image_weights=[Virtual weights (less memory, slower)] osem_acceleration=[1 (balanced)] number_of_iterations="+itenum+" use_tikhonov_regularization tikhonov_parameter=0.0060 compute=[Entire image at once] compute_on=[CPU (Java)] psf_estimation=[Provide file with PSF] psf_display=[Do not show PSFs] output_file_directory="+dataDir+"deconvolved32bit/ use_same_psf_for_all_angles/illuminations browse=/PSF_location/Mean_of_PSFs_"+excitation+".tif transform_psfs psf_file=/PSF_location/Mean_of_PSFs_"+excitation+".tif"


run("Fuse/Deconvolve Dataset", options);
eval("script", "System.exit(0);");
