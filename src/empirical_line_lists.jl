# empirical_line_lists.jl
# generate empirical mask, based on Eric Ford's RvSpectML.jl/examples/expres_analyze_line_by_line.jl, adapted to NEID data by Alex Wise

#TODO: Make sure barycentric shifting is working in template making (i thought its only 20 m/s but its actualy closer to 1 km/s)

""" Check that params module includes required params for a function """
function assert_params_exist(params::Dict{Symbol,Any}, params_to_check::Vector{Symbol})
   for p in params_to_check
      @assert haskey(params,p) string(p," not defined in params")
   end
end

""" Get wavelength of pixel left edge for a given time t, order o, and pixel p """
function get_pixel_left_edge(order_list_timeseries::ACLT, t::Int64, o::Int64, p::Int64) where { ACLT<:AbstractChunkListTimeseries }
   λ = order_list_timeseries[t][o].λ
   if p == 1
      return λ[p] - ((λ[p+1] - λ[p]) / 2)
   else
      return (λ[p-1] + λ[p]) / 2
   end
end

""" Get wavelength of pixel right edge for a given time t, order o, and pixel p """
function get_pixel_right_edge(order_list_timeseries::ACLT, t::Int64, o::Int64, p::Int64) where { ACLT<:AbstractChunkListTimeseries }
   λ = order_list_timeseries[t][o].λ
   if p == length(λ)
      return λ[p] + ((λ[p] - λ[p-1]) / 2)
   else      
      return (λ[p] + λ[p+1]) / 2
   end
end

" Remove nans and negative values from chunk_list_timeseries and return a list of intervals for the wavelengths of affected pixels"
function remove_and_track_nans_negatives!(chunk_list_timeseries::ACLT;
   nan_wavelength_intervals::DataFrame = DataFrame(:chunk => Int[], :wavelength_intervals => Interval{Float64}[]), 
   neg_wavelength_intervals::DataFrame = DataFrame(:chunk => Int[], :wavelength_intervals => Interval{Float64}[]), 
   use_inst_bad_col_ranges::Bool = true, track_only_dont_remove::Bool = false) where { ACLT<:AbstractChunkListTimeseries }
   """
   chunk_list_timeseries = order_list_timeseries
   nan_wavelength_intervals = DataFrame(:chunk => Int[], :wavelength_intervals => Interval{Float64}[])
   neg_wavelength_intervals = DataFrame(:chunk => Int[], :wavelength_intervals => Interval{Float64}[])
   use_inst_bad_col_ranges = true
   track_only_dont_remove = true
   """

   @assert hasproperty(nan_wavelength_intervals,:chunk) && hasproperty(nan_wavelength_intervals,:wavelength_intervals) # make sure input DataFrame has the expected columns for wavelength intervals
   @assert hasproperty(neg_wavelength_intervals,:chunk) && hasproperty(neg_wavelength_intervals,:wavelength_intervals) # make sure input DataFrame has the expected columns for wavelength intervals

   #get inst and max pixels in order
   inst = chunk_list_timeseries.inst
   bad_idx_ij_max = falses(max_pixel_in_order(inst) - min_pixel_in_order(inst) + 1)

   n_times = length(chunk_list_timeseries)
   n_chunks = length(chunk_list_timeseries[1])
   
   current_interval::Interval{Float64} = 0..1 #initialze to reduce allocations in loops

   for j in 1:n_chunks
      #find the existing DataFrame rows for this chunk
      nan_indices_j = findall(nan_wavelength_intervals[:,:chunk] .== j)
      neg_indices_j = findall(neg_wavelength_intervals[:,:chunk] .== j)
      #initialize a vector of intervals for nans and negs for this chunk
      chunk_j_nan_wavelength_intervals = nan_wavelength_intervals[:,:wavelength_intervals][nan_indices_j]
      chunk_j_neg_wavelength_intervals = neg_wavelength_intervals[:,:wavelength_intervals][neg_indices_j]
      #remove the old Dataframe rows for this chunk, these will be replaced at the end of this for loop
      #we do this to make sure the wavelength intervals for each chunk remain non_overlapping
      deleteat!(nan_wavelength_intervals,nan_indices_j) 
      deleteat!(neg_wavelength_intervals,neg_indices_j)

      for i in 1:n_times
         if use_inst_bad_col_ranges
            #get indices of chunk in parent data
            (pixels, ordr) = parentindices(chunk_list_timeseries[i][j].flux)
            #get known bad pixel ranges
            bad_ranges = bad_col_ranges(inst,ordr)
            #convert parent pixel indices to chunk_list_timeseries pixel indices
            for idx in eachindex(bad_ranges)
               bad_ranges[idx] = intersect(bad_ranges[idx],pixels) .- (pixels[1] - 1)
            end
            #add bad pixels ranges to chunk_j_nan_wavelength_intervals
            for r in bad_ranges
               if length(r) > 0
                  current_interval = get_pixel_left_edge(chunk_list_timeseries,i,j,r[1])..get_pixel_right_edge(chunk_list_timeseries,i,j,r[end])
                  chunk_j_nan_wavelength_intervals = union_with_non_overlapping_intervals!(chunk_j_nan_wavelength_intervals,current_interval)
               end
            end
         end
         #check how many pixels are in this chunk, so we know how much of the preallocated bad_idx_ij_max to use
         n_pixels = length(chunk_list_timeseries[i][j].flux)
         bad_idx_ij = view(bad_idx_ij_max,1:n_pixels)
         #check for nans in flux
         bad_idx_ij .= isnan.(chunk_list_timeseries[i][j].flux)
         for k in findall(bad_idx_ij)
            current_interval = get_pixel_left_edge(chunk_list_timeseries,i,j,k)..get_pixel_right_edge(chunk_list_timeseries,i,j,k)
            chunk_j_nan_wavelength_intervals = union_with_non_overlapping_intervals!(chunk_j_nan_wavelength_intervals,current_interval)
         end
         if !track_only_dont_remove
            chunk_list_timeseries[i][j].flux[bad_idx_ij] .= 1.0
         end
         #check for negative flux
         bad_idx_ij .= chunk_list_timeseries[i][j].flux .< 0.0
         for k in findall(bad_idx_ij)
            current_interval = get_pixel_left_edge(chunk_list_timeseries,i,j,k)..get_pixel_right_edge(chunk_list_timeseries,i,j,k)
            chunk_j_neg_wavelength_intervals = union_with_non_overlapping_intervals!(chunk_j_neg_wavelength_intervals,current_interval)
         end
         if !track_only_dont_remove
            chunk_list_timeseries[i][j].flux[bad_idx_ij] .= 0.01
            chunk_list_timeseries[i][j].var[bad_idx_ij] .= 1.0
         end
         #check for nans in var
         bad_idx_ij .= isnan.(chunk_list_timeseries[i][j].var)
         for k in findall(bad_idx_ij)
            current_interval = get_pixel_left_edge(chunk_list_timeseries,i,j,k)..get_pixel_right_edge(chunk_list_timeseries,i,j,k)
            chunk_j_nan_wavelength_intervals = union_with_non_overlapping_intervals!(chunk_j_nan_wavelength_intervals,current_interval)
         end
         if !track_only_dont_remove
            chunk_list_timeseries[i][j].var[bad_idx_ij] .= 1.0
         end
         #check for negative var
         bad_idx_ij .= chunk_list_timeseries[i][j].var .< 0.0
         for k in findall(bad_idx_ij)
            current_interval = get_pixel_left_edge(chunk_list_timeseries,i,j,k)..get_pixel_right_edge(chunk_list_timeseries,i,j,k)
            chunk_j_neg_wavelength_intervals = union_with_non_overlapping_intervals!(chunk_j_neg_wavelength_intervals,current_interval)
         end
         if !track_only_dont_remove
         chunk_list_timeseries[i][j].var[bad_idx_ij] .= 1.0
         end
      end
      chunk_j_nan_wavelength_intervals = discard_emptysets(chunk_j_nan_wavelength_intervals) #in a test run, 59% of neg wavelength intervals were empty sets
      chunk_j_neg_wavelength_intervals = discard_emptysets(chunk_j_neg_wavelength_intervals) #in a test run, 6% of nan wavelength intervals were empty sets
      chunk_j_nan_df = DataFrame(:chunk => fill(j,length(chunk_j_nan_wavelength_intervals)), :wavelength_intervals => chunk_j_nan_wavelength_intervals)
      chunk_j_neg_df = DataFrame(:chunk => fill(j,length(chunk_j_neg_wavelength_intervals)), :wavelength_intervals => chunk_j_neg_wavelength_intervals)
      append!(nan_wavelength_intervals,chunk_j_nan_df)
      append!(neg_wavelength_intervals,chunk_j_neg_df)
   end

   return nan_wavelength_intervals, neg_wavelength_intervals
end


""" take the union of an interval with a vector of non-overlapping intervals """
function union_with_non_overlapping_intervals!(non_overlapping_intervals::Vector{Interval{Float64}}, new_interval::Interval{Float64})
   idx_overlap_first = zero(1)
   for idx in eachindex(non_overlapping_intervals)
      #loop through non_overlapping_intervals and check for overlap with each one
      if new_interval ∩ non_overlapping_intervals[idx] != ∅
         #keep track of the first index of overlap
         if iszero(idx_overlap_first)
            idx_overlap_first = idx
            #combine the first overlapping interval with the new_interval
            non_overlapping_intervals[idx_overlap_first] = non_overlapping_intervals[idx_overlap_first] ∪ new_interval
         else
            #combine with the interval at idx_overlap_first, which already includes new_interval
            non_overlapping_intervals[idx_overlap_first] = non_overlapping_intervals[idx_overlap_first] ∪ non_overlapping_intervals[idx]
            non_overlapping_intervals[idx] = ∅ #set the merged interval to \emptyset to preserve non_overlapping property
         end
      end
   end
   if iszero(idx_overlap_first)
      return push!(non_overlapping_intervals,new_interval)
   else
      return non_overlapping_intervals
   end
end

function discard_emptysets(intervals::Vector{Interval{Float64}})
   return intervals[.!isequal.(intervals,∅)]
end

""" take a vector of intervals and combine the ones that overlap """
function non_overlapping_intervals(possibly_overlapping_intervals::Vector{Interval{Float64}})
   intervals = deepcopy(possibly_overlapping_intervals)
   if length(intervals) > 1
      for idx in eachindex(intervals)[2:end]
         #loop through all previous intervals and check for overlap.
         for idx2 in 1:(idx-1)
            if intervals[idx] ∩ intervals[idx2] != ∅
               #if the intervals overlap, combine them at the current index and set the earlier index to \emptyset
               intervals[idx] = intervals[idx] ∪ intervals[idx2]
               intervals[idx2] = ∅
            end
         end #for idx2 in 1:(idx-1)
      end #for idx in eachindex(intervals)[2:end]
   end #if length(intervals) > 1
   return intervals[.!isequal.(intervals,∅)]
end #function non_overlapping_intervals(intervals::Interval{Float64}[])


""" debugging function used to show where negative and nan values are in flux and var """
function make_save_bad_pixel_plots(order_list_timeseries::ACLT; outdir::String = "/home/awise/Desktop/bad_pixel_plots/") where { ACLT<:AbstractChunkListTimeseries }
   #remove nans from data and track them
   n_times = length(order_list_timeseries)
   n_orders = length(order_list_timeseries[1])
   n_pixels = length(order_list_timeseries[1][1].flux)
   flux_nan_idx = falses(n_times,n_orders,n_pixels) #note these 3-D boolean array sizes assume all orders in order_list_timeseries have the same number of pixels
   flux_neg_idx = falses(n_times,n_orders,n_pixels) 
   var_nan_idx = falses(n_times,n_orders,n_pixels) 
   var_neg_idx = falses(n_times,n_orders,n_pixels)
   flux_nan_idx_ij = falses(n_pixels)
   flux_neg_idx_ij = falses(n_pixels)
   var_nan_idx_ij = falses(n_pixels)
   var_neg_idx_ij = falses(n_pixels)
   for i in 1:n_times
      for j in 1:n_orders
         #check for nans in flux
         flux_nan_idx_ij .= isnan.(order_list_timeseries[i][j].flux)
         flux_nan_idx[i,j,flux_nan_idx_ij] .= true
         #check for negative flux
         flux_neg_idx_ij .= order_list_timeseries[i][j].flux .< 0.0
         flux_neg_idx[i,j,flux_neg_idx_ij] .= true
         #check for nans in var
         var_nan_idx_ij .= isnan.(order_list_timeseries[i][j].var)
         var_nan_idx[i,j,var_nan_idx_ij] .= true
         #check for negative var
         var_neg_idx_ij .= order_list_timeseries[i][j].var .< 0.0
         var_neg_idx[i,j,var_neg_idx_ij] .= true
      end
   end
   neg_idx = flux_neg_idx .| var_neg_idx
   #neg_idx = .|([neg_idx[i,:,:] for i in 1:n_times]...)
   nan_idx = flux_nan_idx .| var_nan_idx
   #nan_idx = .|([nan_idx[i,:,:] for i in 1:n_times]...)

   flux_neg_idx_sum = dropdims(sum(flux_neg_idx,dims=1),dims=1)
   flux_nan_idx_sum = dropdims(sum(flux_nan_idx,dims=1),dims=1)
   var_neg_idx_sum = dropdims(sum(var_neg_idx,dims=1),dims=1)
   var_nan_idx_sum = dropdims(sum(var_nan_idx,dims=1),dims=1)


   flux_neg_idx_bool = flux_neg_idx_sum .!= 0
   flux_nan_idx_bool = flux_nan_idx_sum .!= 0
   var_neg_idx_bool = var_neg_idx_sum .!= 0
   var_nan_idx_bool = var_nan_idx_sum .!= 0

   flux_neg_idx2 = Float64.(flux_neg_idx_sum)
   flux_nan_idx2 = Float64.(flux_nan_idx_sum)
   var_neg_idx2 = Float64.(var_neg_idx_sum)
   var_nan_idx2 = Float64.(var_nan_idx_sum)

   flux_neg_idx2[flux_neg_idx2 .== 0] .= NaN
   flux_nan_idx2[flux_nan_idx2 .== 0] .= NaN
   var_neg_idx2[var_neg_idx2 .== 0] .= NaN
   var_nan_idx2[var_nan_idx2 .== 0] .= NaN

   #using Colors
   #using Plots
   #using Plots.PlotMeasures

   heatmap(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), .~(flux_neg_idx_bool), c="black", xlabel="order #", ylabel="pixel #", title="flux negative values", size=(1600,800), left_margin=7mm, bottom_margin=5mm)
   heatmap!(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), flux_neg_idx2, c=colormap("Greens"), label="flux negative")
   savefig(joinpath(outdir,"flux_neg.png"))

   heatmap(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), .~(flux_nan_idx_bool), c="black", xlabel="order #", ylabel="pixel #", title="flux nan values", size=(1600,800), left_margin=7mm, bottom_margin=5mm)
   heatmap!(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), flux_nan_idx2, c=colormap("Blues"), label="flux nan")
   savefig(joinpath(outdir,"flux_nan.png"))

   heatmap(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), .~(var_neg_idx_bool), c="black", xlabel="order #", ylabel="pixel #", title="var neg values", size=(1600,800), left_margin=7mm, bottom_margin=5mm)
   heatmap!(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), var_neg_idx2, c=colormap("Reds"), label="var negative")
   savefig(joinpath(outdir,"var_neg.png"))

   heatmap(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), .~(var_nan_idx_bool), c="black", xlabel="order #", ylabel="pixel #", title="var nan values", size=(1600,800), left_margin=7mm, bottom_margin=5mm)
   heatmap!(1:size(flux_neg_idx_sum,2), 1:size(flux_neg_idx_sum,1), var_nan_idx2, c=colormap("Purples"), label="var nan")
   savefig(joinpath(outdir,"var_nan.png"))
end


""" check if intervals a and b have any overlap using only their endpoints as input """
function intervals_overlap(a_low::Float64, a_high::Float64, b_low::Float64, b_high::Float64)
   a_low <= b_high && b_low <= a_high
end


""" Take in bad wavelength intervals and match with line list to flag lines as affected by a bad pixel. Output is a bit vector where 1 indicates a match"""
function match_bad_wavelength_intervals_with_lines!(line_list::DataFrame, bad_waves::DataFrame, cl::AbstractChunkList; col_name::Symbol = :bad_line)
   @assert ("pixels" in names(line_list)) && ("chunk_id" in names(line_list)) #make sure line_list has the expected columns for line wavelengths
   @assert hasproperty(bad_waves,:chunk) && hasproperty(bad_waves,:wavelength_intervals) # make sure bad_waves has the expected columns for wavelength intervals
   n_lines = nrow(line_list)
   bad_lines = falses(n_lines)

   for line in eachrow(line_list)
      λs = cl[line[:chunk_id]].λ[line[:pixels]]
      #check one extra pixel on each side of the line
      λlow = λs[1] - (λs[2]-λs[1])
      λhigh = λs[end] + (λs[end]-λs[end-1])
      λ_interval = λlow..λhigh
      #find bad lines in the same chunk
      bad_waves_idx_in_chunk = findall(bad_waves[:,:chunk] .== line[:chunk_id])
      #look for overlap between line and any bad lines in the same chunk
      for bad_wavelength_interval in bad_waves[bad_waves_idx_in_chunk,:wavelength_intervals]
         if bad_wavelength_interval ∩ λ_interval != ∅
            bad_lines[rownumber(line)] = true
         end
      end
   end
   line_list[!, col_name] = bad_lines
end

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
function generateEmpiricalMask(params::Dict{Symbol,Any} ; output_dir::String=params[:output_dir], pipeline_plan::PipelinePlan = PipelinePlan(), verbose::Bool=true)

   """
   params = Params
   output_dir = params[:output_dir]
   pipeline_plan = PipelinePlan()
   verbose = true
   """

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

      cl_derivs = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ, deriv, deriv2, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
      template_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl[chid].λ,cl[chid].flux,extrapolation_bc=Flat()), 1:length(cl))
      template_deriv_linear_interp = map(chid -> Interpolations.LinearInterpolation(cl_derivs[chid].λ,cl_derivs[chid].flux,extrapolation_bc=Flat()), 1:length(cl_derivs))
      # We're done with the spectral_orders_matrix, so we can release the memory now
      spectral_orders_matrix = nothing
      GC.gc()
      need_to!(pipeline_plan,:template)
      #lines_in_template_logy = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=true), verbose=false)  # TODO: Automate threshold for finding a line
      if verbose println("# Performing a fresh search for lines in template spectra.")  end
      @time lines_in_template = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=params[:line_width_50],min_deriv2=0.5, use_logλ=true, use_logflux=false), verbose=false)  # TODO: Automate threshold for finding a line
      match_bad_wavelength_intervals_with_lines!(lines_in_template, neg_wavelength_intervals, cl, col_name = :neg_bad_line)
      match_bad_wavelength_intervals_with_lines!(lines_in_template, nan_wavelength_intervals, cl, col_name = :nan_bad_line)
      @assert (eltype(lines_in_template[!, :neg_bad_line]) == Bool) && (eltype(lines_in_template[!, :nan_bad_line]) == Bool) #make sure the next statement is comparing Booleans
      good_line_idx = findall(.~(lines_in_template[:, :neg_bad_line] .| lines_in_template[:, :nan_bad_line]))

      if params[:discard_neg_nan]
         if verbose println("# Removing " * string(nrow(lines_in_template) - length(good_line_idx)) * " lines due to matching nan or negative values.") end
         lines_to_fit = lines_in_template[good_line_idx,:]
      else
         lines_to_fit = lines_in_template
      end

      if verbose println("# Fitting above lines in all spectra.")  end
      @time fits_to_lines = RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(order_list_timeseries, lines_to_fit )
      @time line_RVs = fit_all_line_RVs_in_chunklist_timeseries(order_list_timeseries, template_linear_interp, template_deriv_linear_interp, lines_to_fit )
   
      if save_data(pipeline_plan,:fit_lines)
         CSV.write(joinpath(output_dir,"linefinder",params[:fits_target_str] * "_linefinder_lines.csv"), lines_in_template )
         CSV.write(joinpath(output_dir,"linefinder",params[:fits_target_str] * "_linefinder_good_lines.csv"), lines_to_fit )
         CSV.write(joinpath(output_dir,"linefinder",params[:fits_target_str] * "_linefinder_line_fits.csv"), fits_to_lines )
         CSV.write(joinpath(output_dir,"linefinder",params[:fits_target_str] * "_linefinder_line_RVs.csv"), line_RVs )
      end
   
      dont_need_to!(pipeline_plan,:fit_lines);
   end
   
   #select lines with good frac converged, median depth, std velocity width, and local continuum slope
   fit_distrib_df = select_line_fits_with_good_depth_width_slope(params, fits_to_lines, parse(Float64,params[:quant]) / 100, verbose=verbose, output_dir=output_dir, fits_target_str=params[:fits_target_str])

   #add neg/nan filters to the Dataframe
   fit_distrib_df[!,:bool_filter_neg_bad_line] = .~lines_in_template[:, :neg_bad_line]
   fit_distrib_df[!,:bool_filter_nan_bad_line] = .~lines_in_template[:, :nan_bad_line]

   return fit_distrib_df

end #end function generateEmpiricalMask()



#original author: Eric Ford
#adapted from: RvSpectML.jl/examples/expres_analyze_line_by_line.jl
function select_line_fits_with_good_depth_width_slope(params::Dict{Symbol,Any}, line_fits_df::DataFrame, quantile_threshold::Real; verbose::Bool = false, output_dir::String = "", fits_target_str::String = "" )

   assert_params_exist(params, [:inst, :min_frac_converged])

   min_frac_converged = parse(Float64,params[:min_frac_converged]) / 100

   fit_distrib = line_fits_df |> @groupby(_.line_id) |>
            @map( { median_a=median(_.fit_a), median_b=median(_.fit_b), median_depth=median(_.fit_depth), median_σ²=median(_.fit_σ²), median_λc=median(_.fit_λc),
                   std_a=std(_.fit_a), std_b=std(_.fit_b), std_depth=std(_.fit_depth), std_σ²=std(_.fit_σ²), std_λc=std(_.fit_λc),
                   line_id=first(_.line_id),  frac_converged=mean(_.fit_converged)  } ) |> DataFrame

   std_depth_threshold = quantile(fit_distrib.std_depth,quantile_threshold)
   median_σ_width_threshold = 2000 # quantile(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,1-quantile_threshold)
   std_σ_width_threshold = quantile(sqrt.(fit_distrib.std_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,quantile_threshold)
   std_b_threshold = quantile(fit_distrib.std_b,quantile_threshold)
   std_a_threshold = quantile(fit_distrib.std_a,quantile_threshold)

   fit_distrib[!,:bool_filter_min_frac_converged] = fit_distrib.frac_converged .>= min_frac_converged
   fit_distrib[!,:bool_filter_median_depth_between_5percent_and_1] = 0.05 .<= fit_distrib.median_depth .<= 1.0
   fit_distrib[!,:bool_filter_std_velocity_width_quant] = sqrt.(fit_distrib.std_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps .<= std_σ_width_threshold 
   fit_distrib[!,:bool_filter_std_local_continuum_slope_quant] = fit_distrib.std_b .<= std_b_threshold

   #good_lines_alt = fit_distrib |>
   #   @filter( 0.05 <= _.median_depth <= 1.0 ) |>
   #   #@filter( sqrt.(_.median_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps >= median_σ_width_threshold ) |>
   #   @filter( sqrt.(_.std_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps <= std_σ_width_threshold ) |>
   #   @filter( _.std_b < std_b_threshold) |>
   #   DataFrame
   
   good_lines_bool = fit_distrib[:,:bool_filter_min_frac_converged] .&& fit_distrib[:,:bool_filter_median_depth_between_5percent_and_1] .&& fit_distrib[:,:bool_filter_std_velocity_width_quant] .&& fit_distrib[:,:bool_filter_std_local_continuum_slope_quant]

   good_lines_alt = fit_distrib[good_lines_bool,:]

   if verbose
      println("# Found ", size(good_lines_alt,1), " good lines (std_depth_width_slope), rejected ", size(fit_distrib,1)-size(good_lines_alt,1), " lines.")
   end
   if length(output_dir) > 0
      quant_str = Printf.@sprintf("%1.2f",quantile_threshold)
      frac_str = Printf.@sprintf("%1.2f",min_frac_converged)
      airVacString = get_airVacString(params[:inst])
      CSV.write(joinpath(output_dir,fits_target_str * "_all_lines_fit_distrib_quant=" * quant_str * "min_frac_converged=" * frac_str * airVacString * ".csv"), fit_distrib )
   end
   return fit_distrib
end


#generate function that finds 100 days between jan 1, 2021 and sept 30, 2021, with highest number of num_rvs.good according to summary_1.csv
function get_NEID_best_days(pipeline_output_summary_path::String; startDate::Date=Date(2021,01,01), endDate::Date=Date(2021,09,30), nBest::Int64=100)
   @assert endDate > startDate
   #summary_1 = CSV.read("/home/awise/data/neid/solar/summary_1.csv", DataFrame)
   summary_1 = CSV.read(pipeline_output_summary_path, DataFrame)
   summary_1_in_date_range = summary_1[(summary_1[!,"obs_date.string"] .>= startDate) .& (summary_1[!,"obs_date.string"] .<= endDate),:]
   indices = sort(reverse(sortperm(summary_1_in_date_range[!,"num_rvs.good"]))[1:nBest])
   NEID_best_days = summary_1_in_date_range[indices,:]
   println(string("# Taking all NEID days with at least ",minimum(NEID_best_days."num_rvs.good")," good RVs."))
   return NEID_best_days."obs_date.string"
end

#given filenames for daily files and a vector of dates, return a vector of all relevant file paths for these dates.
function get_NEID_daily_paths(daily_base_path::String, daily_fn::String, date::Date)
   date_dir = reduce(joinpath,split(Dates.format(date,ISODateFormat),"-"))
   return joinpath(daily_base_path,date_dir,daily_fn)
end

#combine NEID clean daily mean spectra from daily_ccfs files into a vector of type Spectra2DBasic
function combine_NEID_daily_obs(daily_ccfs_base_path::String, daily_ccf_fn::String, daily_manifests_base_path::String, daily_manifest_fn::String, dates::Vector{Date})
   obs_list = get_NEID_daily_paths.(daily_ccfs_base_path, daily_ccf_fn, dates)
   manifest_list = get_NEID_daily_paths.(daily_manifests_base_path, daily_manifest_fn, dates)
   @time obs_lambda_flux_var = [getindex.(Ref(load(fn)),["mean_lambda", "mean_clean_flux_continuum_normalized", "mean_clean_var_continuum_normalized"]) for fn in obs_list] #Ref() is used here so we don't load the file 3 times per observation (once for each of λ, flux, and var fields)
   daily_rvs = zeros(length(manifest_list))
   bjds = zeros(length(manifest_list))
   for i in eachindex(manifest_list)
      manifest_i = CSV.read(manifest_list[i], DataFrame)
      index_minimum_hour_angle = argmin(abs.(manifest_i[:,"hour_angle"]))
      daily_rvs[i] = manifest_i[index_minimum_hour_angle,"ssbz"] * RvSpectML.speed_of_light_mps
      bjds[i] = manifest_i[index_minimum_hour_angle,"bjd"]
   end
   obs_lambda = getindex.(obs_lambda_flux_var,1)
   obs_flux = getindex.(obs_lambda_flux_var,2)
   obs_var = getindex.(obs_lambda_flux_var,3)
   obs_metadata = [Dict{Symbol,Any}(:date=>dates[i], :rv_est=>daily_rvs[i], :bjd=>bjds[i]) for i in 1:length(dates)]
   all_spectra = map(obs -> Spectra2DBasic(obs[1], obs[2], obs[3], NEID2D(), metadata=obs[4]), zip(obs_lambda,obs_flux,obs_var,obs_metadata))
   return all_spectra
end





function fit_one_line_RVs_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, template, template_deriv, lines::DataFrame, idx::Int64)
   @assert 1 <= idx <= size(lines,1)
   λmin = lines[idx,:fit_min_λ]
   λmax = lines[idx,:fit_max_λ]
   chid = lines[idx,:chunk_id]
   df = fit_line_RVs_in_chunklist_timeseries(clt, template, template_deriv, λmin, λmax, chid)
   df[!,:line_id] .= idx
   return df
end

function fit_all_line_RVs_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, template, template_deriv, lines::DataFrame)
   df = DataFrame()
   map(idx -> append!(df,fit_one_line_RVs_in_chunklist_timeseries(clt, template, template_deriv, lines, idx)) , 1:size(lines,1))
   return df
 end


"""
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
"""
 

""" `fit_line_RVs_in_chunklist_timeseries( clt, template, template_deriv, λmin, λmax, chid)`
Return DataFrame with results of fits to each line in a given chunk of chunk_list_timeseries (including each observation time)
Inputs:
- clt: Data to fit
- template: high SNR data template
- template_deriv: derivative of template
- λmin
- λmax
- chunk_index:  Restricts fitting to specified chunk
Outputs a DataFrame with keys:
- fit_z: doppler factor for line fit
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
    deriv = template_deriv[chid](clt.chunk_list[t][chid].λ) ./ mean_template_flux

    fit_z[t] = fit_line_RV(flux, var, template_flux, deriv, pixels[t])
  end
  return DataFrame(:fit_RV=>fit_z * C_m_s, :obs_idx =>collect(1:nobs), :chunk_id=>chid, :pixels=>pixels )
end

C_m_s = 2.99792458e8 #speed of light in m/s


#propagate variance in this calculation to get z_err
#have option to assume variance across line is constant (median across pixels)
function fit_line_RV(flux::AbstractArray{T1,1}, var::AbstractArray{T2,1}, template_flux, deriv, idx::UnitRange) where { T1<:Real, T2<:Real }
   @assert length(flux) == length(var) == length(template_flux) == length(deriv)

   #z = 1 / (deriv[idx]' * deriv[idx]) * (deriv[idx]' * (flux[idx] - template_flux[idx])) #unweighted z

   z = dot(deriv[idx] .* (flux[idx] - template_flux[idx]), 1/var[idx]) / dot(deriv[idx] .* deriv[idx], 1/var[idx]) #weighted z

   #z1 = deriv[idx] \ (flux[idx] - template_flux[idx]) #this should be equal to unweighted z
   #z1 = ((1/sqrt.(var[idx]))' .* deriv[idx]) \ ((1/sqrt.(var[idx]))' .* (flux[idx] - template_flux[idx])) #this should be equal to weighted z
   """
   scatter(flux[idx]-template_flux[idx],z*deriv[idx]) #visual check of the fit - should be a strong correlation here
   xlabel!("data - linear_interp(template)")
   ylabel!("z * d_template_d_z")
   """

   return z
 end

"""


fits_to_lines = CSV.read(joinpath(params[:output_dir],"masks/linefinder",params[:fits_target_str] * "_linefinder_line_fits.csv"), DataFrame )
line_RVs = CSV.read(joinpath(params[:output_dir],"masks/linefinder",params[:fits_target_str] * "_linefinder_line_RVs.csv"), DataFrame )
lines_in_template = CSV.read(joinpath(params[:output_dir],"masks/linefinder",params[:fits_target_str] * "_linefinder_lines.csv"), DataFrame )

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