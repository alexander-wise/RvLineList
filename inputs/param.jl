#params for RvLineList.jl


#module Params

#export max_spectra_to_use
#export pipeline_output_path_ebf11, pipeline_output_path_afw5465
#export inst, data_path, target_subdir, fits_target_str
#export espresso_mask_filename, orders_to_use, norm_type, ccf_mid_velocity
#export linelist_params

Params = Dict{Symbol,Any}()

push!(Params,:max_spectra_to_use => 1000000)
if Params[:max_spectra_to_use] < 200
   println("Warning: param.jl setting max_spectra_to_use to " * string(Params[:max_spectra_to_use]))
end

push!(Params,:long_output => true) #whether or not to carry all line data through to the final mask

path_params = Dict(
#:expres_data_path => "/home/awise/data/expres/"),
#:neid_data_path => "/home/awise/data/neid/"),
#:harpsn_data_path => "/home/awise/data/harps-n/",
#:harps_data_path => "/home/awise/data/harps/",

#:neid_data_path => "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pipeline/v1.1/L2/",
#:neid_data_path => "/home/awise/data/neid/solar",

#:pipeline_output_summary_path => "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.1/outputs/ebf11/Espresso/summary_espressoG2_norm=blaze&mask=1.csv",
#:pipeline_output_summary_path => "/home/awise/data/neid/solar/summary_1.csv",

:pipeline_output_summary_path => "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.1/outputs/ebf11/Espresso_cont3/summary_espressoG2_norm=cont&mask=3.csv",
#:pipeline_output_summary_path => "/home/awise/data/neid/solar/summary_espressoG2_norm=blaze&mask=1.csv",

:daily_ccfs_base_path => "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.1/outputs/ebf11/Espresso_cont3/",
#:daily_ccfs_base_path => "/home/awise/data/neid/solar/",
#:daily_ccf_fn => "daily_ccfs_1.jld2",
:daily_ccf_fn => "daily_ccfs_espressoG2_norm=cont&mask=3.jld2",

:daily_manifests_base_path => "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.1/manifest",
#:daily_manifests_base_path => "/home/awise/data/neid/solar/",
:daily_manifest_fn => "manifest.csv",

:output_dir => joinpath(pwd(),"outputs"),

)

merge!(Params,path_params)


spectra_params = Dict(

:inst => :neid,
#:target_subdir = >"101501",

#data set params
#:EXPRES_excalibur_only => true, #whether or not to only include excalibur wavelengths in expres data - this parameter has been replaced with orders_to_use
#orders_to_use = 43:72 #EXPRES excalibur range
#orders_to_use = sort(sample(43:72,10,replace=false)) #Ten random orders from the EXPRES excalibur range
#order_to_use = 12:83 #EXPRES all useable orders
#orders_to_use = orders_to_use_default(NEID.NEID2D()
#orders_to_use = orders_to_use_default(EXPRES.EXPRES2D()
#orders_to_use = 4:118 # all NEID orders with wavelength values
:orders_to_use => 17:115, # NEID orders with minimal NaNs in extracted orders from daily flux averages

:norm_type => :continuum, #normalization type. Options: :raw, :continuum, :blaze #TODO: this param might not work
)
merge!(Params,spectra_params)


linelist_params = Dict(
# VALD mask params

# lineprops = lineprops_101501 #which VALD line list to use in mask creation - not works in param file yet because mask names do not use this param as they should
#maskWavelengths = string("Reiners",target) #keyword for mask wavelengths - should be e.g. Reiners or Reiners101501
:allowBlends => 0, #number of blends allowed with each line, e.g. 0 for single lines only, 1 for doublets only, [0,1] for singlets and doublets
:overlap_cutoff => 2.0e-5, #distance between lines required for them to be categorized as a blend, expressed as a fraction of the speed of light
:depth_cutoff => 0.05,
:iron1Only => "all", # which species to use - options are "Fe1", "nonFe1", "all"
:badLineFilter => "none", #which mask to use to filter out bad line - only lines within ~ 3 km/s of one of these mask lines will be kept
:rejectTelluricSlope => 2000.0, #derivative of spectrum required for telluric line rejection - a value of 0 turns off telluric rejection
:nbin => 1, #number of mask subdivisions to make
#bin_n = 1, #which subdivision to use - one is first
:binParam => :depth, #which mask parameter to bin multiple masks by - only important for subdividing masks - can be :depth or :lambda
:depthPercentile => true, # if binning by depth, use sum of depths for bin weights


# Empirical mask params
#:discard_neg_nan => true, #whether or not lines affected by negative / nan values should be discarded #commented out as this now happens by default
:quant => "90", #quantile for stability of fit params for empirical masks, expressed as a string out of 100
:min_frac_converged => "100", #minimum fraction of the line fits that converged in the dataset for the empirical line to be used, expressed as a string out of 100
:line_width_50 => 7392.0, #NEID solar line width for ESPRESSO G2 mask, mask_scale_factor = 2.0, other ccf params default values

)

merge!(Params,linelist_params)

# RV params
#lsf_width = 2.0e3 #estimated line width for telluric avoidance
#mask_type = :gaussian #mask shape. Options: :tophat, :gaussian, :supergaussian, :halfcos
#mask_scale_factor = 4.0 #approximate number of pixels in the width scale of the mask
#fwtf = 2.0 #fraction of the FWHM of the CCF to fit
#RV_fit_type = :gaussian #type of fit to use on CCF to measure RV. Options: :quadratic, :gaussian

#global tophap_ccf_mask_scale_factor=1.6

#end


#local fits_target_str
#local espresso_mask_filename
#local ccf_mid_velocity
#local data_path

if Params[:inst] == :expres
   @assert haskey(Params,:expres_data_path)
   @assert haskey(Params,:target_subdir)
   push!(Params,:data_path => joinpath(Params[:expres_data_path],Params[:target_subdir]))
   push!(Params,:fits_target_str => Params[:target_subdir])
   if Params[:fits_target_str] == "101501"
      push!(Params,:espresso_mask_filename => "G9.espresso.mas")
      push!(Params,:ccf_mid_velocity => -5.0e3)
   end
   if Params[:fits_target_str] == "10700"
      push!(Params,:espresso_mask_filename => "G8.espresso.mas")
      push!(Params,:ccf_mid_velocity => -16640.0)
   end
   if Params[:fits_target_str] == "34411"  #G1?
      push!(Params,:espresso_mask_filename => "G2.espresso.mas")
      push!(Params,:ccf_mid_velocity => 66500.0)
   end
   if Params[:fits_target_str] == "26965" # K0?
      push!(Params,:espresso_mask_filename => "G9.espresso.mas")
      push!(Params,:ccf_mid_velocity => -40320.0)
   end
end

if Params[:inst] == :neid #we are using solar data
   push!(Params,:fits_target_str => "neid")
   #@assert haskey(Params,:neid_data_path)
   #push!(Params,:data_path => Params[:neid_data_path])
   push!(Params,:espresso_mask_filename => "G2.espresso.mas")
   push!(Params,:ccf_mid_velocity => 0.0)
end

if Params[:inst] == :harps #the target is Alpha Centauri B
   @assert haskey(Params,:harps_data_path)
   @assert haskey(Params,:target_subdir)
   push!(Params,:data_path => joinpath(Params[:harps_data_path],Params[:target_subdir]))
   push!(Params,:espresso_mask_filename => "K2.espresso.mas")
   push!(Params,:ccf_mid_velocity => -22000.0)
end

if Params[:inst] == :harpsn #we are using solar data
   @assert haskey(Params,:harpsn_data_path)
   @assert haskey(Params,:target_subdir)
   push!(Params,:data_path => joinpath(Params[:harpsn_data_path],Params[:target_subdir]))
   push!(Params,:espresso_mask_filename => "G2.espresso.mas")
   push!(Params,:ccf_mid_velocity => 0.0)
end
