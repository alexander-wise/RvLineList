#local_dir = "/home/awise/Desktop/RvSpectML/RvLineList"
#local_dir = "/storage/work/afw5465/RvSpectML"

#using Pkg
# Pkg.activate(".")
using RvLineList, DataFrames, CSV

local_dir = pkgdir(RvLineList)
cd(local_dir)

#load param file
include(joinpath(local_dir,"inputs","param.jl"))

#generate empirical NEID mask

pipeline_plan = RvLineList.PipelinePlan()
RvLineList.RvSpectMLBase.Pipeline.save_data!(pipeline_plan, :fit_lines) #tell the pipeline plan to save the line fits

#####
#begin generateEmpiricalMask() function code
#####

params = Params
output_dir = params[:output_dir]
verbose = true

RvLineList.assert_params_exist(params, [:output_dir, :norm_type, :inst, :daily_ccfs_base_path, :daily_manifests_base_path, :pipeline_output_summary_path, :daily_ccf_fn, :daily_manifest_fn, :orders_to_use, :fits_target_str, :line_width_50, :min_frac_converged, :quant])

RvLineList.reset_all_needs!(pipeline_plan) #TODO: revise module PipelinePlan or my usage of the module so this line is not needed.

#if verbose println("# Reading in customized parameters from param.jl.")  end
#paths_to_search_for_param = [pwd(),joinpath(pwd(),"inputs"),joinpath(pkgdir(RvLineList),"inputs")]
#eval(read_data_paths(paths_to_search=paths_to_search_for_param, filename="param.jl"))

if params[:norm_type] != :continuum
   println("Warning: norm_type=:continuum is required by generateEmpiricalMask() but found norm_type=:", String(params[:norm_type]))
end

if RvLineList.need_to(pipeline_plan,:read_spectra)
   if params[:inst] == :expres
      all_spectra = RvLineList.EXPRES_read_spectra(data_path) #note this is in expres_old.jl and needs to be loaded manually for now. Also this uses data_path and max_spectra_to_use param (not included in params_to_check since this is deprecated)
   elseif params[:inst] == :neid
      all_spectra = RvLineList.combine_NEID_daily_obs(params[:daily_ccfs_base_path], params[:daily_ccf_fn], params[:daily_manifests_base_path], params[:daily_manifest_fn], RvLineList.get_NEID_best_days(params[:pipeline_output_summary_path],startDate=RvLineList.Date(2021,01,01), endDate=RvLineList.Date(2021,09,30), nBest=100))
   else
      print("Error: spectra failed to load; inst not supported.")
   end
   RvLineList.dont_need_to!(pipeline_plan,:read_spectra)
end

#make sure EXPRES spectra are normalized
if typeof(RvLineList.get_inst(all_spectra)) <: RvLineList.AnyEXPRES   RvLineList.continuum_normalize_spectra!(all_spectra)   end

#extract orders_to_use from all_spectra. We use remove_bad_chunks=false because neid has some nans due to a bad detector column
order_list_timeseries = RvLineList.extract_orders(all_spectra, pipeline_plan, orders_to_use=params[:orders_to_use], remove_bad_chunks=false, recalc=true);

if verbose println("# Removing and tracking negative and nan values from spectra.") end
@time nan_wavelength_intervals, neg_wavelength_intervals = RvLineList.remove_and_track_nans_negatives!(order_list_timeseries)


if RvLineList.need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
   if verbose println("# Making template spectrum.")  end
   @assert !RvLineList.need_to(pipeline_plan,:extract_orders)
   GC.gc()   # run garbage collector for deallocated memory
   #map(i->order_list_timeseries.metadata[i][:rv_est] = 0.0, 1:length(order_list_timeseries) )
   @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvLineList.RvSpectML.make_template_spectra(order_list_timeseries, smooth_factor=2.0) #, alg=:Sinc)
   #@time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, alg=:Sinc)
   if RvLineList.save_data(pipeline_plan, :template)
      #using JLD2, FileIO
      RvLineList.save(joinpath(output_dir,params[:fits_target_str] * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   end
   RvLineList.dont_need_to!(pipeline_plan, :template);
end


if RvLineList.need_to(pipeline_plan,:fit_lines)
   if verbose println("# Initializing search for lines in template spectrum.")  end
   cl = RvLineList.ChunkList(map(grid->RvLineList.ChunkOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
   clt = RvLineList.ChunkListTimeseries([0.0], [cl], inst=first(all_spectra).inst, metadata=[Dict{Symbol,Any}()] )
   if verbose println("# Removing and tracking negative and nan values from template.") end
   #add neg/nan wavelength intervals from template to the original DataFrames
   nan_wavelength_intervals, neg_wavelength_intervals = RvLineList.remove_and_track_nans_negatives!(clt, nan_wavelength_intervals=nan_wavelength_intervals, neg_wavelength_intervals=neg_wavelength_intervals, use_inst_bad_col_ranges=false)
   """ #setup code for measuring RVs using template matching technique -- currently not used
   cl_derivs = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ, deriv, deriv2, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
   template_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl[chid].λ,cl[chid].flux,extrapolation_bc=Flat()), 1:length(cl))
   template_deriv_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl_derivs[chid].λ,cl_derivs[chid].flux,extrapolation_bc=Flat()), 1:length(cl_derivs))
   """
   # We're done with the spectral_orders_matrix, so we can release the memory now
   spectral_orders_matrix = nothing
   GC.gc()
   RvLineList.need_to!(pipeline_plan,:template)
   #lines_in_template_logy = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=true), verbose=false)  # TODO: Automate threshold for finding a line
   if verbose println("# Performing a fresh search for lines in template spectra.")  end
   @time lines_in_template = RvLineList.RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvLineList.RvSpectML.LineFinder.LineFinderPlan(line_width=params[:line_width_50],min_deriv2=0.5, use_logλ=true, use_logflux=false), verbose=false)  # TODO: Automate threshold for finding a line
   RvLineList.match_bad_wavelength_intervals_with_lines!(lines_in_template, neg_wavelength_intervals, cl, col_name = :neg_bad_line)
   RvLineList.match_bad_wavelength_intervals_with_lines!(lines_in_template, nan_wavelength_intervals, cl, col_name = :nan_bad_line)

   too_close_to_telluric = RvLineList.getTelluricIndices(lines_in_template, true, params[:overlap_cutoff], vel_slope_threshold = params[:rejectTelluricSlope], RV_offset = 0.0, RV_range = 1e-4)
   lines_in_template[!, :rejectTelluricSlope] = too_close_to_telluric

   @assert (eltype(lines_in_template[!, :neg_bad_line]) == Bool) && (eltype(lines_in_template[!, :nan_bad_line]) == Bool) && (eltype(too_close_to_telluric) == Bool) #make sure these are Booleans so the next line is valid
   good_line_idx = findall(.~(lines_in_template[:, :neg_bad_line] .| lines_in_template[:, :nan_bad_line] .| too_close_to_telluric))
   bad_line_idx = findall((lines_in_template[:, :neg_bad_line] .| lines_in_template[:, :nan_bad_line] .| too_close_to_telluric))
   if verbose println("# Found " * string(nrow(lines_in_template) - length(good_line_idx)) * " lines contaminated by nan or negative values or tellurics.") end

   lines_to_fit = lines_in_template[good_line_idx,:]

   if verbose println("# Fitting uncontaminated lines in all spectra.")  end
   @time fits_to_lines = RvLineList.RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_to_fit )
   #@time line_RVs = fit_all_line_RVs_in_chunklist_timeseries(order_list_timeseries, template_linear_interp, template_deriv_linear_interp, lines_in_template )

   RvLineList.dont_need_to!(pipeline_plan,:fit_lines);
end

#select lines with good frac converged, median depth, std velocity width, and local continuum slope
fit_distrib_df = RvLineList.select_line_fits_with_good_depth_width_slope(params, fits_to_lines, parse(Float64,params[:quant]) / 100, verbose=verbose, output_dir=output_dir, fits_target_str=params[:fits_target_str])

#create a new Dataframe with full line_in_template size to track neg/nan lines
df_total = similar(fit_distrib_df,size(lines_in_template,1))
df_total[good_line_idx,:] = fit_distrib_df
#fix the line_ids since the fits didn't know about the neg/nan lines
df_total[good_line_idx,:line_id] = good_line_idx
df_total[bad_line_idx,:line_id] = bad_line_idx
#add neg/nan/telluric columns to the Dataframe
df_total[!,:bool_filter_neg_bad_line] = .~lines_in_template[:, :neg_bad_line]
df_total[!,:bool_filter_nan_bad_line] = .~lines_in_template[:, :nan_bad_line]
df_total[!,:bool_filter_rejectTelluricSlope] = .~lines_in_template[:, :rejectTelluricSlope]

empirical_mask = df_total

#####
#end generateEmpiricalMask() function code
#####

rename!(empirical_mask, :median_λc => :lambda)
rename!(empirical_mask, :median_depth => :depth)

using PyCall
try
   pyimport("pandas") #make sure pandas is installed
catch e
   pyimport_conda("pandas", "pandas") #if the module fails to load, install it to Conda.jl environment
   pyimport("pandas") #make sure pandas is installed
end
import Pandas.DataFrame as pd_df

@pyinclude("src/make_VALD_line_list.py")
@time VALD_mask, VALD_mask_long = py"getVALDmasks"(overlap_cutoff=Params[:overlap_cutoff], depth_cutoff=Params[:depth_cutoff], iron1Only=Params[:iron1Only], badLineFilter=Params[:badLineFilter], allowBlends=Params[:allowBlends])

function pd_df_to_df(df_pd)
   df = DataFrame()
   for col in df_pd.columns
      df[!,col] = getproperty(df_pd,col).values
   end
   df
end

VALD_mask_long_df = pd_df_to_df(VALD_mask_long)
VALD_mask_long_df[!,:species] .= convert.(String,VALD_mask_long_df[!,:species])

combined_mask_df = mask_intersection(empirical_mask, VALD_mask_long_df, threshold=500.0)

combined_mask_pd = pd_df(combined_mask_df)

combined_mask_df[!,:passed_all_bool_filters] = (combined_mask_df[:,:bool_filter_min_frac_converged]
.&& combined_mask_df[:,:bool_filter_median_depth_between_5percent_and_1]
.&& combined_mask_df[:,:bool_filter_std_velocity_width_quant]
.&& combined_mask_df[:,:bool_filter_std_local_continuum_slope_quant]
.&& combined_mask_df[:,:bool_filter_neg_bad_line]
.&& combined_mask_df[:,:bool_filter_nan_bad_line]
.&& combined_mask_df[:,:bool_filter_depth_cutoff]
.&& combined_mask_df[:,:bool_filter_allowBlends]
.&& combined_mask_df[:,:bool_filter_iron1Only]
.&& combined_mask_df[:,:bool_filter_rejectTelluricSlope]
.&& combined_mask_df[:,:bool_filter_badLineFilter])

combined_mask_df_filtered = combined_mask_df[ combined_mask_df[!,:passed_all_bool_filters], :]

mask = combined_mask_df_filtered



frac_depths = [0.2,0.3,0.4,0.5,0.6,0.625,0.65,0.675,0.7,0.75,0.8,0.9]

#frac_depths = 0.6:0.001:0.75

line_RVs = zeros(size(mask,1),length(order_list_timeseries),length(frac_depths))
times = zeros(size(line_RVs))
for j in 1:length(order_list_timeseries)
   for i in 1:size(mask,1)
      for k in eachindex(frac_depths)
         line_id = Int(mask.line_id[i])
         pixels = RvLineList.RvSpectMLBase.find_pixels_for_line_in_chunklist(order_list_timeseries.chunk_list[j], lines_in_template.fit_min_λ[line_id], lines_in_template.fit_max_λ[line_id], lines_in_template.chunk_id[line_id]).pixels
         pixels = max(pixels[1]-15,1) : min(pixels[end]+15,6501)
         line_RVs[i,j,k] = (RvLineList.RvSpectMLBase.calc_line_bisector_at_frac_depth(order_list_timeseries[j][lines_in_template.chunk_id[line_id]].λ[pixels], 
         order_list_timeseries[j][lines_in_template.chunk_id[line_id]].flux[pixels]; frac_depth=frac_depths[k]) - lines_in_template.fit_λc[line_id]) / lines_in_template.fit_λc[line_id] * RvLineList.C_m_s
         times[i,j,k] = order_list_timeseries.times[j]
      end
   end
end 


line_RVs_template = zeros(size(mask,1),length(frac_depths))
for i in 1:size(mask,1)
   for k in eachindex(frac_depths)
      line_id = Int(mask.line_id[i])
      pixels = RvLineList.RvSpectMLBase.find_pixels_for_line_in_chunklist(cl, lines_in_template.fit_min_λ[line_id], lines_in_template.fit_max_λ[line_id], lines_in_template.chunk_id[line_id]).pixels
      pixels = max(pixels[1]-15,1) : min(pixels[end]+15,6499)
      line_RVs_template[i,k] = (RvLineList.RvSpectMLBase.calc_line_bisector_at_frac_depth(cl[lines_in_template.chunk_id[line_id]].λ[pixels], 
      cl[lines_in_template.chunk_id[line_id]].flux[pixels]; frac_depth=frac_depths[k]) - lines_in_template.fit_λc[line_id]) / lines_in_template.fit_λc[line_id] * RvLineList.C_m_s
   end
end


deep_lines = findall(mask.depth .> 0.7)
middle_lines = findall(0.5 .< mask.depth .<= 0.7)
shallow_lines = findall(mask.depth .<= 0.5)

using Plots

#scatter(times[deep_lines,:,1]',line_RVs[deep_lines,:,8]')
#histogram(line_RVs[deep_lines,:,8]')

line_RVs_nonan = deepcopy(line_RVs)
line_RVs_nonan[findall(isnan.(line_RVs_nonan))] .= 1e5


stds_deep = [RvLineList.median(RvLineList.std(line_RVs_nonan[deep_lines,:,i],dims=2)) for i in 1:length(frac_depths)]
stds_middle = [RvLineList.median(RvLineList.std(line_RVs_nonan[middle_lines,:,i],dims=2)) for i in 1:length(frac_depths)]
stds_shallow = [RvLineList.median(RvLineList.std(line_RVs_nonan[shallow_lines,:,i],dims=2)) for i in 1:length(frac_depths)]

scatter(frac_depths,stds_deep ./ RvLineList.median(stds_deep))
scatter!(frac_depths,stds_middle ./ RvLineList.median(stds_middle))
scatter!(frac_depths,stds_shallow ./ RvLineList.median(stds_shallow))



#RvLineList.std(line_RVs_template[deep_lines,8])
#RvLineList.std(RvLineList.median(line_RVs[deep_lines,:,8],dims=2))
#line_RVs_template[deep_lines,8] - RvLineList.median(line_RVs[deep_lines,:,8],dims=2)
