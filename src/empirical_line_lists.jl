# empirical_line_lists.jl
# generate empirical mask, based on Eric Ford's RvSpectML.jl/examples/expres_analyze_line_by_line.jl, adapted to NEID data by Alex Wise

#TODO: Make sure barycentric shifting is working in template making (i thought its only 20 m/s but its actualy closer to 1 km/s)

""" Generate an empirical mask by finding spectral lines, fiting them with gaussians with a linear offset, and filtering the fitted line params

Inputs (must be defined in param.jl):
   inst::Module - instrument used to collect the data. Works with: :expres, :neid

   data_path::String - string specifying the path to the data

Optional inputs:

   times_to_use::UnitRange{Int64} - indices of the data to use out of the result of inst.make_manifest(). #TODO: change this to take a list of datetimes in a universal format to mitigate risk of index mismatch

   isSolarData::Bool - boolean specifying whether or not the spectra are taken of the Sun

   output_dir::String - path from the current directory to the directory where the empricial mask should be saved (Default: "outputs/masks/")

   pipeline_plan::PipelinePlan - pipeline tracker for RvSpectML

   verbose:Bool - determines whether descriptions of each step are printed

Steps:
   1) add functions/scripts to load data: EXPRES, NEID, HARPS-N, solarDataOrNighttimeDataBoolean param, times_to_use param, make sure data are normalized
   1.5) make a function to generate times_to_use for NEID solar data, using +/- X hours of solar noon (see 1.1 in to-do list.txt)
   2) make a template spectrum using the loaded data (see :template section below)
   3) find/fit lines in template spectrum (:fit lines section below) - FUTURE - add option to fit doublets / preset list of lines
   4) select lines based on fits (already done below probs)

Outputs:
   Empirical Mask::DataFrame: empirically-fitted and filtered RV mask including (specify fields)
"""
function generateEmpiricalMask( ; output_dir::String="", pipeline_plan::PipelinePlan = PipelinePlan(), verbose::Bool=true)

   #isSolarData=false
   #output_dir = "outputs/masks/"
   #pipeline_plan = PipelinePlan()
   #verbose = true
   reset_all_needs!(pipeline_plan) #TODO: revise module PipelinePlan or my usage of the module so this line is not needed.


   if verbose println("# Reading in customized parameters from param.jl.")  end
   paths_to_search_for_param = [pwd(),joinpath(pwd(),"inputs"),joinpath(pkgdir(RvLineList),"inputs")]
   eval(read_data_paths(paths_to_search=paths_to_search_for_param, filename="param.jl"))

   global norm_type
   if norm_type != :continuum
      println("Warning: norm_type=:continuum is required by generateEmpiricalMask() but found norm_type=:", String(norm_type), ". Setting global norm_type=:continuum in order to continue.")
      norm_type = :continuum
   end

   if need_to(pipeline_plan,:read_spectra)
      if inst == :expres
         all_spectra = EXPRES_read_spectra(data_path) #note this is in expres_old.jl and needs to be loaded manually for now
      elseif inst == :neid
         all_spectra = combine_NEID_daily_obs(get_NEID_best_days(startDate=Date(2021,01,01), endDate=Date(2021,09,30), nBest=100))
      else
         print("Error: spectra failed to load; inst not supported.")
      end
      dont_need_to!(pipeline_plan,:read_spectra)
   end

   #make sure EXPRES spectra are normalized
   if typeof(get_inst(all_spectra)) <: AnyEXPRES   RvLineList.continuum_normalize_spectra!(all_spectra)   end

   #extract orders_to_use from all_spectra. We use remove_bad_chunks=false because neid has some nans due to a bad detector column
   order_list_timeseries = extract_orders(all_spectra, pipeline_plan, orders_to_use=orders_to_use, remove_bad_chunks=false, recalc=true);

   #remove nans from data before making template spectra - TODO: turn this into an interpolation or modify template construction to allow nans
   for i in eachindex(order_list_timeseries)
     for j in eachindex(order_list_timeseries[1])
       order_list_timeseries[i][j].flux[isnan.(order_list_timeseries[i][j].flux)] .= 1.0
     end
   end

   if need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
      if verbose println("# Making template spectra.")  end
      @assert !need_to(pipeline_plan,:extract_orders)
      GC.gc()   # run garbage collector for deallocated memory
      #map(i->order_list_timeseries.metadata[i][:rv_est] = 0.0, 1:length(order_list_timeseries) )
      # Smothing is broken with GP interpolation.  Need to fix.  In mean time, here's a Sinc interpolation workaround
      @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, smooth_factor=2.0) #, alg=:Sinc)
      #@time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, alg=:Sinc)
      if save_data(pipeline_plan, :template)
         #using JLD2, FileIO
         save(joinpath(output_dir,fits_target_str * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
      end
      dont_need_to!(pipeline_plan, :template);
   end
   
   if need_to(pipeline_plan,:fit_lines)
      if verbose println("# Performing fresh search for lines in template spectra.")  end
      cl = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
      cl_derivs = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ, deriv, deriv2, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
      template_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl[chid].λ,cl[chid].flux,extrapolation_bc=Flat()), 1:length(cl))
      template_deriv_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl_derivs[chid].λ,cl_derivs[chid].flux,extrapolation_bc=Flat()), 1:length(cl_derivs))
      # We're done with the spectral_orders_matrix, so we can release the memory now
      spectral_orders_matrix = nothing
      GC.gc()
      need_to!(pipeline_plan,:template)
      #lines_in_template_logy = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=true), verbose=false)  # TODO: Automate threshold for finding a line
      lines_in_template = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=false), verbose=false)  # TODO: Automate threshold for finding a line
   
      if verbose println("# Finding above lines in all spectra.")  end
      @time fits_to_lines = RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_in_template )
      @time line_RVs = fit_all_line_RVs_in_chunklist_timeseries(order_list_timeseries, template_linear_interp, template_deriv_linear_interp, lines_in_template )
   
      if save_data(pipeline_plan,:fit_lines)
         CSV.write(joinpath(output_dir,"linefinder",fits_target_str * "_linefinder_lines.csv"), lines_in_template )
         CSV.write(joinpath(output_dir,"linefinder",fits_target_str * "_linefinder_line_fits.csv"), fits_to_lines )
         CSV.write(joinpath(output_dir,"linefinder",fits_target_str * "_linefinder_line_RVs.csv"), line_RVs )
      end
   
      dont_need_to!(pipeline_plan,:fit_lines);
   end
   
   select_line_fits_with_good_depth_width_slope(fits_to_lines, parse(Float64,quant) / 100, verbose=verbose, output_dir=output_dir)

end #end function generateEmpiricalMask()



#author: Eric Ford
#adapted from: RvSpectML.jl/examples/expres_analyze_line_by_line.jl
function select_line_fits_with_good_depth_width_slope(line_fits_df::DataFrame, quantile_threshold::Real; verbose::Bool = false, output_dir::String = "" )
   fit_distrib = line_fits_df |> @groupby(_.line_id) |>
            @map( { median_a=median(_.fit_a), median_b=median(_.fit_b), median_depth=median(_.fit_depth), median_σ²=median(_.fit_σ²), median_λc=median(_.fit_λc),
                   std_a=std(_.fit_a), std_b=std(_.fit_b), std_depth=std(_.fit_depth), std_σ²=std(_.fit_σ²), std_λc=std(_.fit_λc),
                   line_id=first(_.line_id),  frac_converged=mean(_.fit_converged)  } ) |>
            @filter(_.frac_converged == 1.0 ) |> DataFrame

   std_depth_treshold = quantile(fit_distrib.std_depth,quantile_threshold)
   median_σ_width_treshold = 2000 # quantile(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,1-quantile_threshold)
   std_σ_width_treshold = quantile(sqrt.(fit_distrib.std_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,quantile_threshold)
   std_b_treshold = quantile(fit_distrib.std_b,quantile_threshold)
   std_a_treshold = quantile(fit_distrib.std_a,quantile_threshold)
   good_lines_alt = fit_distrib |>
      @filter( 0.05 <= _.median_depth <= 0.95 ) |>
      #@filter( sqrt.(_.median_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps >= median_σ_width_treshold ) |>
      @filter( sqrt.(_.std_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps <= std_σ_width_treshold ) |>
      @filter( _.std_b < std_b_treshold) |>

      DataFrame
   if verbose
      println("# Found ", size(good_lines_alt,1), " good lines (std_depth_width_slope), rejected ", size(fit_distrib,1)-size(good_lines_alt,1), " lines.")
   end
   if length(output_dir) > 0
      val_str = Printf.@sprintf("%1.2f",quantile_threshold)
      airVacString = get_airVacString(inst)
      CSV.write(joinpath(output_dir,fits_target_str * "_good_lines_fit_quant=" * val_str * airVacString * ".csv"), good_lines_alt )
   end
   return good_lines_alt
end


#generate function that finds 100 days between jan 1, 2021 and sept 30, 2021, with highest number of num_rvs.good according to summary_1.csv
function get_NEID_best_days(; pipeline_output_summary_path::String=joinpath(pipeline_output_path_ebf11,"summary_1.csv"), startDate::Date=Date(2021,01,01), endDate::Date=Date(2021,09,30), nBest::Int64=100)
   @assert endDate > startDate
   #summary_1 = CSV.read("/home/awise/data/neid/solar/summary_1.csv", DataFrame)
   summary_1 = CSV.read(pipeline_output_summary_path, DataFrame)
   summary_1_in_date_range = summary_1[(summary_1[!,"obs_date.string"] .>= startDate) .& (summary_1[!,"obs_date.string"] .<= endDate),:]
   indices = sort(reverse(sortperm(summary_1_in_date_range[!,"num_rvs.good"]))[1:nBest])
   NEID_best_days = summary_1_in_date_range[indices,:]
   println(string("Taking all NEID days with at least ",minimum(NEID_best_days."num_rvs.good")," good RVs."))
   return NEID_best_days."obs_date.string"
end

#given a vector of dates, as well as filenames for ccf and manifest files, get a vector of all relevant file paths for these dates.
function get_NEID_daily_obs_list_and_manifest_list(dates::Vector{Date}; pipeline_output_path=pipeline_output_path_ebf11, ccf_fn::String = "daily_ccfs_1.jld2", manifest_fn = "manifest.csv")
   date_dirs = reduce.(joinpath,split.(Dates.format.(dates,ISODateFormat),"-"))
   obs_list = joinpath.(pipeline_output_path,date_dirs,ccf_fn)
   manifest_list = joinpath.(pipeline_output_path,date_dirs,manifest_fn)
   return obs_list, manifest_list
end

#combine NEID clean daily mean spectra from daily_ccfs files into a vector of type Spectra2DBasic
function combine_NEID_daily_obs(dates::Vector{Date})
   obs_list, manifest_list = get_NEID_daily_obs_list_and_manifest_list(dates)
   @time obs_lambda_flux_var = [getindex.(Ref(load(fn)),["mean_lambda", "mean_clean_flux_continuum_normalized", "mean_clean_var"]) for fn in obs_list]
   daily_rvs = zeros(length(manifest_list))
   for i in eachindex(manifest_list)
      manifest_i = CSV.read(manifest_list[i], DataFrame)
      index_minimum_hour_angle = argmin(abs.(manifest_i[:,"hour_angle"]))
      daily_rvs[i] = manifest_i[index_minimum_hour_angle,"ssbz"] * speed_of_light_mps
   end
   obs_lambda = [obs[1] for obs in obs_lambda_flux_var]
   obs_flux = [obs[2] for obs in obs_lambda_flux_var]
   obs_var = [obs[3] for obs in obs_lambda_flux_var]
   obs_metadata = [Dict{Symbol,Any}(:date=>dates[i], :rv_est=>daily_rvs[i], :bjd=>datetime2julian(DateTime(Dates.year(dates[i]), Dates.month(dates[i]), Dates.day(dates[i]), 12))) for i in 1:length(dates)]
   all_spectra = map(obs -> Spectra2DBasic(obs[1], obs[2], obs[3], NEID2D(), metadata=obs[4]), eachrow(hcat(obs_lambda,obs_flux,obs_var,obs_metadata))) #I think there's a better way to do this, perhaps by using obs_lambda_flux_var direclty
   return all_spectra
end







function fit_all_line_RVs_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, template, template_deriv, lines::DataFrame)
   @assert size(lines,1) >= 2
   line_idx::Int64 = 1
   λmin = lines[line_idx,:fit_min_λ]
   λmax = lines[line_idx,:fit_max_λ]
   chid = lines[line_idx,:chunk_id]
   df = fit_line_RVs_in_chunklist_timeseries(clt, template, template_deriv, λmin, λmax, chid)
   df[!,:line_id] .= line_idx
   for line_idx in 2:size(lines,1)
     λmin = lines[line_idx,:fit_min_λ]
     λmax = lines[line_idx,:fit_max_λ]
     chid = lines[line_idx,:chunk_id]
     df_tmp = fit_line_RVs_in_chunklist_timeseries(clt, template, template_deriv, λmin, λmax, chid)
     df_tmp[!,:line_id] .= line_idx
     append!(df,df_tmp)
   end
   return df
 end

 

""" `fit_line_in_chunklist_timeseries( chunk_list_timeseries, λmin, λmax, chunk_index)`
Return DataFrame with results of fits to each line in a given chunk of chunk_list_timeseries (including each observation time)
Inputs:
- chunk_list_timeseries: Data to fit
- template: high SNR data template
- λmin
- λmax
- chunk_index:  Restricts fitting to specified chunk
Outputs a DataFrame with keys:
- fit_z: doppler factor for line fit
- fit_covar: Covariance matrix for fit parameters
- χ²_per_dof: quality of fit to chunk
- fit_converged: Bool indicating if `LsqFit` converged
- obs_idx: index of observation in chunk_list_timeseries
- chunk_id: index of chunk in chunk_list_timeseries
- pixels: range of pixels that was fit
"""
function fit_line_RVs_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, template, template_deriv, λmin::Real, λmax::Real, chid::Integer)
  nobs = length(clt.chunk_list)
  fit_z = Vector{Float64}(undef,nobs)
  pixels = Vector{UnitRange{Int64}}(undef,nobs)
  global count_msgs

  for t in 1:nobs
    pixels[t] = RvSpectMLBase.find_pixels_for_line_in_chunklist(clt.chunk_list[t], λmin, λmax, chid).pixels

    mean_flux = mean(clt.chunk_list[t][chid].flux[pixels[t]])
    flux = clt.chunk_list[t][chid].flux ./ mean_flux
    var = clt.chunk_list[t][chid].var ./ mean_flux^2

    mean_template_flux = mean(template[chid](clt.chunk_list[t][chid].λ)[pixels[t]])
    template_flux = template[chid](clt.chunk_list[t][chid].λ) ./ mean_template_flux
    deriv = template_deriv[chid](clt.chunk_list[t][chid].λ) ./ mean_template_flux #is this correct? to divide the derivative by the mean flux to normalize?

    fit_z[t] = fit_line_RV(flux, var, template_flux, deriv, pixels[t])
  end
  return DataFrame(:fit_RV=>fit_z * C_m_s, :obs_idx =>collect(1:nobs), :chunk_id=>chid, :pixels=>pixels )
end

C_m_s = 2.99792458e8 #speed of light in m/s


#propagate variance in this calculation to get z_err
#have option to assume variance across line is constant (median across pixels)
function fit_line_RV(flux::AbstractArray{T1,1}, var::AbstractArray{T2,1}, template_flux, deriv, idx::UnitRange) where { T1<:Real, T2<:Real }
   @assert length(flux) == length(var) == length(template_flux) == length(deriv)

   #z = deriv[idx] \ (flux[idx] - template_flux[idx]) #unweighted z
   z = ((1/sqrt.(var[idx]))' .* deriv[idx]) \ ((1/sqrt.(var[idx]))' .* (flux[idx] - template_flux[idx])) #weighted z
   """
   scatter(flux[idx]-template_flux[idx],z*deriv[idx]) #visual check of the fit - should be a strong correlation here
   xlabel!("data - linear_interp(template)")
   ylabel!("z * d_template_d_z")

   z1 = 1 / (deriv[idx]' * deriv[idx]) * (deriv[idx]' * (flux[idx] - template_flux[idx])) #this should be equal to unweighted z

   z1 = dot(deriv[idx] .* (flux[idx] - template_flux[idx]), 1/var[idx]) / dot(deriv[idx] .* deriv[idx], 1/var[idx]) #this should be equal to weighted z

   """

   return z
 end

"""


fits_to_lines = CSV.read(joinpath(output_dir,"masks/linefinder",fits_target_str * "_linefinder_line_fits.csv"), DataFrame )
line_RVs = CSV.read(joinpath(output_dir,"masks/linefinder",fits_target_str * "_linefinder_line_RVs.csv"), DataFrame )
lines_in_template = CSV.read(joinpath(output_dir,"masks/linefinder",fits_target_str * "_linefinder_lines.csv"), DataFrame )

C_m_s = 2.99792458e8 #speed of light in m/s


using RvLineList
using Dates
using Query
using Plots
using StatsBase

days = Dates.dayofyear.(RvLineList.get_NEID_best_days())

imax = size(lines_in_template)[1]

idx = sample(mask[:,:line_id], 100, replace = false)

for i in idx


j=0



j+=1
i = idx[j]
gaussian_fits = fits_to_lines |> @filter(_.line_id==i) |> DataFrame
gaussian_fit_avg_center = mean(gaussian_fits[:,:fit_λc])
gaussian_fit_RVs = (gaussian_fits[:,:fit_λc] .- gaussian_fit_avg_center) ./ gaussian_fit_avg_center .* C_m_s

template_fits = line_RVs |> @filter(_.line_id==i) |> DataFrame
template_fit_RVs = template_fits[:,:fit_RV] .- mean(template_fits[:,:fit_RV])
line_id = template_fits[1,:line_id]
mask_id = findfirst(i-> i==line_id,mask[:,:line_id])
species1 = mask[mask_id,:species]

scatter(days,gaussian_fit_RVs, label="gaussian fit")
scatter!(days,template_fit_RVs,label="template fit")
title!("species = " * string(species1) * ", wavelength = " * string(mask[mask_id,:lambda]))
xlabel!("day of year in 2021")
display(ylabel!("RV (m/s)"))

chunk_id = gaussian_fits[1,:chunk_id]
obs_ids = 75:85
pixels = eval(Meta.parse(gaussian_fits[1,:pixels]))
obs_id = obs_ids[1]

λ_to_use = [order_list_timeseries[obs_id].data[chunk_id].λ[pixels] for obs_id in obs_ids]
flux_to_use = [order_list_timeseries[obs_id].data[chunk_id].flux[pixels] for obs_id in obs_ids]
mean_λ = [mean([λ_to_use[x][y] for x in 1:length(obs_ids)]) for y in 1:length(pixels)]
mean_flux = [mean([flux_to_use[x][y] for x in 1:length(obs_ids)]) for y in 1:length(pixels)]


plot(λ_to_use[1], flux_to_use[1] - mean_flux)
for k in 2:length(obs_ids)
   plot!(λ_to_use[k], flux_to_use[k] -  mean_flux)
end
display(ylabel!("flux"))

#line id = 4960, check why its not symmetric, 
#change output to file to not be PyObjects

end

i+=1

scatter(template_fit_RVs,gaussian_fit_RVs)
xlabel("template fit RV (m/s)")
ylabel("gaussian fit RV (m/s)")

i+=1
"""