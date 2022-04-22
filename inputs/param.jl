#params for RvLineList.jl

#module params
#export max_spectra_to_use
#export inst, data_path, target_subdir, fits_target_str
#export espresso_mask_filename, orders_to_use, norm_type, ccf_mid_velocity


max_spectra_to_use = 1000000
if max_spectra_to_use < 200
   @warn "param.in setting max_spectra_to_use to " * max_spectra_to_use
end

#expres_data_path = "/home/awise/data/expres/"
#neid_data_path="/home/awise/data/neid/"
#harpsn_data_path="/home/awise/data/harps-n/"
#harps_data_path="/home/awise/data/harps/"

pipeline_output_path_ebf11 = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pipeline/v1.1/outputs/ebf11/test1/"
#pipeline_output_path_ebf11 = "/home/awise/data/neid/solar/"

pipeline_output_path_afw5465 = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pipeline/v1.1/outputs/afw5465/test1/"
#pipeline_output_path_afw5465 = "/home/awise/data/neid/solar/"

neid_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pipeline/v1.1/L2/"

output_dir = joinpath(pwd(),"outputs")

inst = NEID
#target_subdir="101501"

local fits_target_str
local espresso_mask_filename
local ccf_mid_velocity
local data_path

if inst==EXPRES
   @assert @isdefined expres_data_path
   data_path = joinpath(expres_data_path,target_subdir)
   fits_target_str = target_subdir
   if fits_target_str == "101501"
      espresso_mask_filename = "G9.espresso.mas"
      ccf_mid_velocity = -5.0e3
   end
   if fits_target_str == "10700"
      espresso_mask_filename = "G8.espresso.mas"
      ccf_mid_velocity = -16640.0
   end
   if fits_target_str == "34411"  #G1?
      espresso_mask_filename = "G2.espresso.mas"
      ccf_mid_velocity = 66500.0
   end
   if fits_target_str == "26965" # K0?
      espresso_mask_filename = "G9.espresso.mas"
      ccf_mid_velocity = -40320.0
   end
end

if inst==NEID #we are using solar data
   fits_target_str="neid"
   @assert @isdefined neid_data_path
   data_path = neid_data_path
   espresso_mask_filename = "G2.espresso.mas"
   ccf_mid_velocity = 0.0
end

if inst==HARPS #the target is Alpha Centauri B
   @assert @isdefined harps_data_path
   data_path = joinpath(harps_data_path,target_subdir)
   espresso_mask_filename = "K2.espresso.mas"
   ccf_mid_velocity = -22000.0
end

if inst==HARPSN #we are using solar data
   @assert @isdefined harpsn_data_path
   data_path = joinpath(harpsn_data_path,target_subdir)
   espresso_mask_filename = "G2.espresso.mas"
   ccf_mid_velocity = 0.0
end



#data set params
#EXPRES_excalibur_only = true #whether or not to only include excalibur wavelengths in expres data - this parameter has been replaced with orders_to_use
#orders_to_use = 43:72 #EXPRES excalibur range
#orders_to_use = sort(sample(43:72,10,replace=false)) #Ten random orders from the EXPRES excalibur range
#order_to_use = 12:83 #EXPRES all useable orders
#orders_to_use = orders_to_use_default(NEID2D())
#orders_to_use = 4:118 # all NEID orders with wavelength values
orders_to_use = 17:115 # NEID orders with minimal NaNs in extracted orders from daily flux averages

norm_type = :continuum #normalization type. Options: :raw, :continuum, :blaze



# VALD mask params

# lineprops = lineprops_101501 #which VALD line list to use in mask creation - not works in param file yet because mask names do not use this param as they should
#maskWavelengths = string("Reiners",target) #keyword for mask wavelengths - should be e.g. Reiners or Reiners101501
allowBlends = 0 #number of blends allowed with each line, e.g. 0 for single lines only, 1 for doublets only, [0,1] for singlets and doublets
overlap_cutoff = 1e-5 #distance between lines required for them to be categorized as a blend, expressed as a fraction of the speed of light
depth_cutoff=0.05
iron1Only = "all" # which species to use - options are "Fe1", "nonFe1", "all"
badLineFilter = "none" #which mask to use to filter out bad line - only lines within ~ 3 km/s of one of these mask lines will be kept
rejectTelluricSlope = 0.0 #derivative of spectrum required for telluric line rejection - a value of 0 turns off telluric rejection
nbin = 1 #number of mask subdivisions to make
#bin_n = 1 #which subdivision to use - one is first
binParam = :depth #which mask parameter to bin multiple masks by - only important for subdividing masks - can be :depth or :wavelength
depthPercentile = true # if binning by depth, use sum of depths for bin weights


# Empirical mask params
quant = "90" #quantile for stability of fit params for empirical masks
#min_frac_converged = "90" #minimum fraction of the line fits that converged in the dataset for the empirical line to be used
line_width_50 = 7392.0 #NEID solar line width for ESPRESSO G2 mask, mask_scale_factor = 2.0, other ccf params default values


# RV params
#lsf_width = 2.0e3 #estimated line width for telluric avoidance
#mask_type = :gaussian #mask shape. Options: :tophat, :gaussian, :supergaussian, :halfcos
#mask_scale_factor = 4.0 #approximate number of pixels in the width scale of the mask
#fwtf = 2.0 #fraction of the FWHM of the CCF to fit
#RV_fit_type = :gaussian #type of fit to use on CCF to measure RV. Options: :quadratic, :gaussian

#global tophap_ccf_mask_scale_factor=1.6

#end