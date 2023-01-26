#local_dir = "/home/awise/Desktop/RvSpectML/RvLineList"
#local_dir = "/storage/work/afw5465/RvSpectML"

#using Pkg
# Pkg.activate(".")
using RvLineList, DataFrames, CSV, ArgParse

local_dir = pkgdir(RvLineList)
cd(local_dir)

#load param file
include(joinpath(local_dir,"inputs","param.jl"))

s = ArgParseSettings()

@add_arg_table s begin
   "--allowBlends"
      help = "number of blends allowed with each line, e.g. 0 for single lines only, 1 for doublets only, [0,1] for singlets and doublets"
      arg_type = Int
   "--overlap_cutoff"
      help = "distance between lines required for them to be categorized as a blend, expressed as a fraction of the speed of light"
      arg_type = Float64
   "--rejectTelluricSlope"
      help = "derivative of spectrum required for telluric line rejection, in units of normalized depth per fraction of speed of light - a value of 0 turns off telluric rejection, 2000-10000 is recommended range in values"
      arg_type = Float64
   "--badLineFilter"
      help = "which mask to use to filter out bad lines. Tested Options: none, ESPRESSOG2"
      arg_type = String
   "--quant"
      help = "quantile for stability of fit params for empirical masks - BROKEN in command line - works in param file" # BROKEN in command line - works in param file
      arg_type = String
   "--nbin"
      help = "number of mask subdivisions to make"
      arg_type = Int
   "--output_dir"
      help = "directory to save mask outputs in"
      arg_type = String
   "--long_output"
      help = "whether or not to carry all available line data through to the final mask"
      arg_type = Bool
end

   parsed_args = parse_args(s)


if parsed_args["allowBlends"] !== nothing
   println("Found command-line arg allowBlends. Overwriting param file definition of this arg.")
   Params[:allowBlends] = parsed_args["allowBlends"]
end

if parsed_args["overlap_cutoff"] !== nothing
   println("Found command-line arg overlap_cutoff. Overwriting param file definition of this arg.")
   Params[:overlap_cutoff] = parsed_args["overlap_cutoff"]
end

if parsed_args["rejectTelluricSlope"] !== nothing
   println("Found command-line arg rejectTelluricSlope. Overwriting param file definition of this arg.")
   Params[:rejectTelluricSlope] = parsed_args["rejectTelluricSlope"]
end

if parsed_args["badLineFilter"] !== nothing
   println("Found command-line arg badLineFilter. Overwriting param file definition of this arg.")
   Params[:badLineFilter] = parsed_args["badLineFilter"]
end

if parsed_args["quant"] !== nothing
   println("Found command-line arg quant. Overwriting param file definition of this arg.")
   Params[:quant] = parsed_args["quant"]
end

if parsed_args["nbin"] !== nothing
   println("Found command-line arg nbin. Overwriting param file definition of this arg.")
   Params[:nbin] = parsed_args["nbin"]
end

if parsed_args["output_dir"] !== nothing
   println("Found command-line arg output_dir. Overwriting param file definition of this arg.")
   Params[:output_dir] = parsed_args["output_dir"]
end

if parsed_args["long_output"] !== nothing
   println("Found command-line arg long_output. Overwriting param file definition of this arg.")
   Params[:long_output] = parsed_args["long_output"]
end


#generate empirical NEID mask

#set up output directories
if !isdir(Params[:output_dir])
   mkdir(Params[:output_dir])
end
if !isdir(joinpath(Params[:output_dir],"clean_masks"))
   mkdir(joinpath(Params[:output_dir],"clean_masks"))
end
if !isdir(joinpath(Params[:output_dir],"VALD_masks"))
   mkdir(joinpath(Params[:output_dir],"VALD_masks"))
end
if !isdir(joinpath(Params[:output_dir],"mask_bins"))
   mkdir(joinpath(Params[:output_dir],"mask_bins"))
end
if !isdir(joinpath(Params[:output_dir],"linefinder"))
   mkdir(joinpath(Params[:output_dir],"linefinder"))
end

pipeline_plan = RvLineList.PipelinePlan()
RvLineList.RvSpectMLBase.Pipeline.save_data!(pipeline_plan, :fit_lines) #tell the pipeline plan to save the line fits

#####
#begin generateEmpiricalMask() function code
#####

params = Params
output_dir = params[:output_dir]
verbose = true

assert_params_exist(params, [:output_dir, :norm_type, :inst, :daily_ccfs_base_path, :daily_manifests_base_path, :pipeline_output_summary_path, :daily_ccf_fn, :daily_manifest_fn, :orders_to_use, :fits_target_str, :line_width_50, :min_frac_converged, :quant])

reset_all_needs!(pipeline_plan) #TODO: revise module PipelinePlan or my usage of the module so this line is not needed.

#if verbose println("# Reading in customized parameters from param.jl.")  end
#paths_to_search_for_param = [pwd(),joinpath(pwd(),"inputs"),joinpath(pkgdir(RvLineList),"inputs")]
#eval(read_data_paths(paths_to_search=paths_to_search_for_param, filename="param.jl"))

if params[:norm_type] != :continuum
   println("Warning: norm_type=:continuum is required by generateEmpiricalMask() but found norm_type=:", String(params[:norm_type]))
end

if need_to(pipeline_plan,:read_spectra)
   if params[:inst] == :expres
      all_spectra = EXPRES_read_spectra(data_path) #note this is in expres_old.jl and needs to be loaded manually for now. Also this uses data_path and max_spectra_to_use param (not included in params_to_check since this is deprecated)
   elseif params[:inst] == :neid
      all_spectra = combine_NEID_daily_obs(params[:daily_ccfs_base_path], params[:daily_ccf_fn], params[:daily_manifests_base_path], params[:daily_manifest_fn], get_NEID_best_days(params[:pipeline_output_summary_path],startDate=Date(2021,01,01), endDate=Date(2021,09,30), nBest=100))
   else
      print("Error: spectra failed to load; inst not supported.")
   end
   dont_need_to!(pipeline_plan,:read_spectra)
end

#make sure EXPRES spectra are normalized
if typeof(get_inst(all_spectra)) <: AnyEXPRES   RvLineList.continuum_normalize_spectra!(all_spectra)   end

#extract orders_to_use from all_spectra. We use remove_bad_chunks=false because neid has some nans due to a bad detector column
order_list_timeseries = extract_orders(all_spectra, pipeline_plan, orders_to_use=params[:orders_to_use], remove_bad_chunks=false, recalc=true);

if verbose println("# Removing and tracking negative and nan values from spectra.") end
@time nan_wavelength_intervals, neg_wavelength_intervals = remove_and_track_nans_negatives!(order_list_timeseries)


if need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
   if verbose println("# Making template spectrum.")  end
   @assert !need_to(pipeline_plan,:extract_orders)
   GC.gc()   # run garbage collector for deallocated memory
   #map(i->order_list_timeseries.metadata[i][:rv_est] = 0.0, 1:length(order_list_timeseries) )
   @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, smooth_factor=2.0) #, alg=:Sinc)
   #@time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, alg=:Sinc)
   if save_data(pipeline_plan, :template)
      #using JLD2, FileIO
      save(joinpath(output_dir,params[:fits_target_str] * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   end
   dont_need_to!(pipeline_plan, :template);
end


if need_to(pipeline_plan,:fit_lines)
   if verbose println("# Initializing search for lines in template spectrum.")  end
   cl = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
   clt = ChunkListTimeseries([0.0], [cl], inst=first(all_spectra).inst, metadata=[Dict{Symbol,Any}()] )
   if verbose println("# Removing and tracking negative and nan values from template.") end
   #add neg/nan wavelength intervals from template to the original DataFrames
   nan_wavelength_intervals, neg_wavelength_intervals = remove_and_track_nans_negatives!(clt, nan_wavelength_intervals=nan_wavelength_intervals, neg_wavelength_intervals=neg_wavelength_intervals, use_inst_bad_col_ranges=false)
   """ #setup code for measuring RVs using template matching technique -- currently not used
   cl_derivs = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ, deriv, deriv2, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
   template_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl[chid].λ,cl[chid].flux,extrapolation_bc=Flat()), 1:length(cl))
   template_deriv_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl_derivs[chid].λ,cl_derivs[chid].flux,extrapolation_bc=Flat()), 1:length(cl_derivs))
   """
   # We're done with the spectral_orders_matrix, so we can release the memory now
   spectral_orders_matrix = nothing
   GC.gc()
   need_to!(pipeline_plan,:template)
   #lines_in_template_logy = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=true), verbose=false)  # TODO: Automate threshold for finding a line
   if verbose println("# Performing a fresh search for lines in template spectra.")  end
   @time lines_in_template = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=params[:line_width_50],min_deriv2=0.5, use_logλ=true, use_logflux=false), verbose=false)  # TODO: Automate threshold for finding a line
   match_bad_wavelength_intervals_with_lines!(lines_in_template, neg_wavelength_intervals, cl, col_name = :neg_bad_line)
   match_bad_wavelength_intervals_with_lines!(lines_in_template, nan_wavelength_intervals, cl, col_name = :nan_bad_line)
   @assert (eltype(lines_in_template[!, :neg_bad_line]) == Bool) && (eltype(lines_in_template[!, :nan_bad_line]) == Bool) #make sure these are Booleans so the next line is valid
   good_line_idx = findall(.~(lines_in_template[:, :neg_bad_line] .| lines_in_template[:, :nan_bad_line]))
   if verbose println("# Found " * string(nrow(lines_in_template) - length(good_line_idx)) * " lines contaminated by nan or negative values.") end

   """
   if params[:discard_neg_nan]
      if verbose println("# Removing " * string(nrow(lines_in_template) - length(good_line_idx)) * " lines due to matching nan or negative values.") end
      lines_to_fit = lines_in_template[good_line_idx,:] #TODO save "accepted" and "rejected" line lists instead of "lines" and "good_lines"
   else
      lines_to_fit = lines_in_template
   end
   """

   lines_to_fit = lines_in_template[good_line_idx,:]


   frac_depths = [0.2,0.3,0.4,0.5,0.6,0.625,0.65,0.675,0.7,0.75,0.8,0.9]

   frac_depths = 0.6:0.001:0.75

   line_RVs = zeros(size(mask,1),length(order_list_timeseries),length(frac_depths))
   times = zeros(size(line_RVs))
   for j in 1:length(order_list_timeseries)
      for i in 1:size(mask,1)
         for k in eachindex(frac_depths)
            line_id = Int(mask.line_id[i])
            pixels = RvSpectMLBase.find_pixels_for_line_in_chunklist(order_list_timeseries.chunk_list[j], lines_to_fit.fit_min_λ[line_id], lines_to_fit.fit_max_λ[line_id], lines_to_fit.chunk_id[line_id]).pixels
            pixels = max(pixels[1]-15,1) : min(pixels[end]+15,6501)
            line_RVs[i,j,k] = (RvSpectMLBase.calc_line_bisector_at_frac_depth(order_list_timeseries[j][lines_to_fit.chunk_id[line_id]].λ[pixels], 
            order_list_timeseries[j][lines_to_fit.chunk_id[line_id]].flux[pixels]; frac_depth=frac_depths[k]) - lines_to_fit.fit_λc[line_id]) / lines_to_fit.fit_λc[line_id] * C_m_s
            times[i,j,k] = order_list_timeseries.times[j]
         end
      end
   end 


   line_RVs_template = zeros(size(mask,1),length(frac_depths))
   for i in 1:size(mask,1)
      for k in eachindex(frac_depths)
         line_id = Int(mask.line_id[i])
         pixels = RvSpectMLBase.find_pixels_for_line_in_chunklist(cl, lines_to_fit.fit_min_λ[line_id], lines_to_fit.fit_max_λ[line_id], lines_to_fit.chunk_id[line_id]).pixels
         pixels = max(pixels[1]-15,1) : min(pixels[end]+15,6499)
         line_RVs_template[i,k] = (RvSpectMLBase.calc_line_bisector_at_frac_depth(cl[lines_to_fit.chunk_id[line_id]].λ[pixels], 
         cl[lines_to_fit.chunk_id[line_id]].flux[pixels]; frac_depth=frac_depths[k]) - lines_to_fit.fit_λc[line_id]) / lines_to_fit.fit_λc[line_id] * C_m_s
      end
   end


   deep_lines = findall(mask.depth .> 0.7)
   middle_lines = findall(0.5 .< mask.depth .<= 0.7)
   shallow_lines = findall(mask.depth .<= 0.5)

   scatter(times[deep_lines,:,1]',line_RVs[deep_lines,:,8]')
   histogram(line_RVs[deep_lines,:,8]')

   line_RVs_nonan = deepcopy(line_RVs)
   line_RVs_nonan[findall(isnan.(line_RVs_nonan))] .= 1e5
   stds_deep = [median(std(line_RVs_nonan[deep_lines,:,i],dims=2)) for i in 1:length(frac_depths)]
   stds_middle = [median(std(line_RVs_nonan[middle_lines,:,i],dims=2)) for i in 1:length(frac_depths)]
   stds_shallow = [median(std(line_RVs_nonan[shallow_lines,:,i],dims=2)) for i in 1:length(frac_depths)]

   scatter(frac_depths,stds_deep)
   scatter(frac_depths,stds_middle)
   scatter(frac_depths,stds_shallow)




   std(line_RVs_template[deep_lines,8])
   std(median(line_RVs[deep_lines,:,8],dims=2))
   line_RVs_template[deep_lines,8] - median(line_RVs[deep_lines,:,8],dims=2)
