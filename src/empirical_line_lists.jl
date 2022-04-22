# empirical_line_lists.jl
# generate empirical mask, based on Eric Ford's RvSpectML.jl/examples/expres_analyze_line_by_line.jl, adapted to NEID data by Alex Wise

#TODO: Make sure barycentric shifting is working in template making (i thought its only 20 m/s but its actualy closer to 1 km/s)

""" Generate an empirical mask by finding spectral lines, fiting them with gaussians with a linear offset, and filtering the fitted line params

Inputs (must be defined in param.jl):
   inst::Module - instrument used to collect the data. Works with: EXPRES

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
      if inst == EXPRES
         all_spectra = EXPRES_read_spectra(data_path)
      elseif inst == NEID
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
   for i in 1:length(order_list_timeseries)
     for j in 1:length(order_list_timeseries[1])
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
      # We're done with the spectral_orders_matrix, so we can release the memory now
      #spectral_orders_matrix = nothing
      #GC.gc()
      need_to!(pipeline_plan,:template)
      #lines_in_template_logy = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=true), verbose=false)  # TODO: Automate threshold for finding a line
      lines_in_template = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=false), verbose=false)  # TODO: Automate threshold for finding a line
   
      if verbose println("# Finding above lines in all spectra.")  end
      @time fits_to_lines = RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_in_template )
   
      if save_data(pipeline_plan,:fit_lines)
         CSV.write(joinpath(output_dir,"linefinder",fits_target_str * "_linefinder_lines.csv"), lines_in_template )
         CSV.write(joinpath(output_dir,"linefinder",fits_target_str * "_linefinder_line_fits.csv"), fits_to_lines )
      end
   
      dont_need_to!(pipeline_plan,:fit_lines);
   end
   
   select_line_fits_with_good_depth_width_slope(fits_to_lines, parse(Float64,quant) / 100, verbose=verbose, output_dir=output_dir)

end #end function generateEmpiricalMask()

#read spectra code for EXPRES data
#TODO: move this somewhere else?
function EXPRES_read_spectra(data_path::String; verbose::Bool=false)
   if verbose println("# Finding what data files are avaliable.")  end
   if isfile(joinpath(pipeline_output_path_afw5465,"manifest.csv"))
      if verbose println("# Reading in manifest from manifest.csv") end
      df_files  = CSV.read(joinpath(pipeline_output_path_afw5465,"manifest.csv"), DataFrame)
      @assert size(df_files,1) >= 1
      @assert hasproperty(df_files,:Filename)
      @assert hasproperty(df_files,:bjd)
   else
      if verbose println("# Generating manifest file manifest.csv") end
      df_files = inst.make_manifest(data_path, inst)
      CSV.write(joinpath(pipeline_output_path_afw5465,"manifest.csv"), df_files)
   end
   df_files_use = df_files |> @take(max_spectra_to_use) |> DataFrame
   if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
   @time all_spectra = map(inst.read_data,eachrow(df_files_use))
end


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
function get_NEID_best_days(; startDate::Date=Date(2021,01,01), endDate::Date=Date(2021,09,30), nBest::Int64=100)
   @assert endDate > startDate
   #summary_1 = CSV.read("/home/awise/data/neid/solar/summary_1.csv", DataFrame)
   summary_1 = CSV.read(joinpath(pipeline_output_path_ebf11,"summary_1.csv"), DataFrame)
   summary_1_in_date_range = summary_1[(summary_1[!,"obs_date.string"] .>= startDate) .& (summary_1[!,"obs_date.string"] .<= endDate),:]
   indices = sort(reverse(sortperm(summary_1_in_date_range[!,"num_rvs.good"]))[1:nBest])
   NEID_best_days = summary_1_in_date_range[indices,:]
   println(string("Taking all NEID days with at least ",minimum(NEID_best_days."num_rvs.good")," good RVs."))
   return NEID_best_days."obs_date.string"
end

#combine NEID clean daily mean spectra from daily_ccfs files into a vector of type Spectra2DBasic
function combine_NEID_daily_obs(dates::Vector{Date})
   date_dirs = reduce.(joinpath,split.(Dates.format.(dates,ISODateFormat),"-"))
   obs_list = joinpath.(pipeline_output_path_ebf11,date_dirs,"daily_ccfs_1.jld2")
   manifest_list = joinpath.(pipeline_output_path_ebf11,date_dirs,"manifest.csv")
   #manifest_list = joinpath.(pipeline_output_path_ebf11,date_dirs,"manifest.csv")
   #reduce(vcat,CSV.read.(manifest_list,DataFrame))
   @time obs_lambda_flux_var = [getindex.(Ref(load(fn)),["mean_lambda", "mean_clean_flux_continuum_normalized", "mean_clean_var"]) for fn in obs_list]
   daily_rvs = zeros(length(manifest_list))
   for i in 1:length(manifest_list)
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




