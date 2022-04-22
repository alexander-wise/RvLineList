# functions for neid solar spectra

#function to read a set of NEID solar spectra and filter based on hour angle
#deprecated because data set is too large to load into memory at once - using daily averages instead (see: combine_NEID_daily_obs)
function read_and_filter_all_spectra(data_path::String; verbose::Bool=false, max_solar_hour_angle::Float=0.0)
   if verbose println("# Finding what data files are avaliable.")  end
   if isfile(joinpath(pipeline_output_path_afw5465,"manifest.csv"))
      if verbose println("# Reading in manifest from manifest.csv") end
      df_files  = CSV.read(joinpath(pipeline_output_path_afw5465,"manifest.csv"), DataFrame)
      @assert size(df_files,1) >= 1
      @assert hasproperty(df_files,:Filename)
      @assert hasproperty(df_files,:bjd)
   else
      if verbose println("# Generating manifest file manifest.csv") end
      if inst==NEID
         df_files = make_manifest_NEID(data_path, max_spectra_to_use=max_spectra_to_use)
      else
         df_files = inst.make_manifest(data_path, inst)
      end
      CSV.write(joinpath(pipeline_output_path_afw5465,"manifest.csv"), df_files)
   end
   if @isdefined fits_target_str
      df_files_use = df_files |>
      @filter( _.target == parse(Int64,fits_target_str) ) |>
      @filter( max_solar_hour_angle == 0.0 || abs(_.hour_angle) <= max_solar_hour_angle ) |>
      @take(max_spectra_to_use) |> DataFrame
   else
      df_files_use = df_files |>
      @filter( max_solar_hour_angle == 0.0 || abs(_.hour_angle) <= max_solar_hour_angle ) |>
      @take(max_spectra_to_use) |> DataFrame
   end
   if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
   @time all_spectra = map(NEID.read_data,eachrow(df_files_use))
end



##########################################################################################
#scripts that were aimed towards making a daily template spectrum - keeping this here for example code but should not assume it is working
##########################################################################################

pipeline_plan = PipelinePlan()
verbose = true
reset_all_needs!(pipeline_plan)

if verbose println("# Reading in customized parameters from param.jl.")  end
paths_to_search_for_param = [pwd(),joinpath(pwd(),"inputs"),joinpath(pkgdir(RvSpectMLBase),"..","RvLineList","inputs")]
eval(read_data_paths(paths_to_search=paths_to_search_for_param, filename="param.jl"))

"""
global norm_type
if norm_type != :continuum
   println("Warning: norm_type=:continuum is required by generateEmpiricalMask() but found norm_type=:", String(norm_type), ". Setting global norm_type=:continuum in order to continue.")
   norm_type = :continuum
end
"""

if need_to(pipeline_plan,:read_spectra)
   all_spectra = read_spectra(joinpath(data_path,date_subdir), max_solar_hour_angle=1.0)
   dont_need_to!(pipeline_plan,:read_spectra)
end

#make sure spectra are normalized
#all_spectra = include(joinpath(pkgdir(EchelleInstruments),"examples/read_neid_solar_data_20190918.jl"))
#if typeof(get_inst(all_spectra)) <: AnyEXPRES   RvLineList.continuum_normalize_spectra!(all_spectra)   end

order_list_timeseries = extract_orders(all_spectra, pipeline_plan, orders_to_use=orders_to_use, recalc=true);

# First, make line list using default orders, where laser comb calibration is available
#linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",espresso_mask_filename)
#line_list_df = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = RvSpectMLBase.max_bc, recalc=true, verbose=true)
#(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_scale_factor=2.0, ccf_mid_velocity=ccf_mid_velocity, recalc=true)
#line_width_50 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.5)
#line_width_05 = RvSpectMLBase.calc_line_width(v_grid,view(ccfs,:,1),frac_depth=0.05)
#line_list_excalibur_df = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan,  orders_to_use=43:72, v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = line_width_05, recalc=true, verbose=true)
#((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries, line_list_excalibur_df, pipeline_plan,  mask_scale_factor=8.0, range_no_mask_change=line_width_05, ccf_mid_velocity=ccf_mid_velocity, v_step=100, calc_ccf_var=true, recalc=true)
#alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=2)
#rvs_ccf = calc_rvs_from_ccf_total(ccfs, ccf_vars, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)

global f_mean
if need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
   if verbose println("# Making template spectra.")  end
   @assert !need_to(pipeline_plan,:extract_orders)
   GC.gc()   # run garbage collector for deallocated memory
   map(i->order_list_timeseries.metadata[i][:rv_est] = 0.0, 1:length(order_list_timeseries) )
   # Smothing is broken with GP interpolation.  Need to fix.  In mean time, here's a Sinc interpolation workaround
   @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, smooth_factor=2.0) #, alg=:Sinc)
   #@time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, alg=:Sinc)
   #if save_data(pipeline_plan, :template)
   using JLD2, FileIO
   save(joinpath(pipeline_output_path_afw5465,date_subdir * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
   #end
   dont_need_to!(pipeline_plan, :template);
end



using RCall
afs_src = joinpath(pwd(),"AFS.R")
R"source($afs_src)"

function calc_continuum_model(spectrum::AbstractSpectra2D; order_idx::Integer )
    possible_pix = get_pixel_range(get_inst(spectrum),order_idx)
    bad_pix = bad_col_ranges(get_inst(spectrum),order_idx)
    pix_rng = EchelleInstruments.calc_complement_index_ranges(possible_pix,bad_pix)
    pix = mapreduce(p->collect(p),vcat,pix_rng)
    afs_inputs = zeros(Float64,length(pix),2)
    afs_inputs[:,1] .= spectrum.λ[pix,order_idx]
    afs_inputs[:,2] .= spectrum.flux[pix,order_idx]
    @assert !any(isnan.(afs_inputs))
    #=
    wv = mapreduce(p->spec.λ[p,order_idx],vcat,pix_rng)
    @assert !any(isnan.(wv))
    inten = mapreduce(p->convert(Vector{Float64},spec.flux[p,order_idx]),vcat,pix_rng)
    @assert !any(isnan.(inten))
    afs_inputs = hcat(wv,inten)
    =#
    #df = DataFrame("wv"=>wv,"intes"=>inten)
    afs_output_R = R"AFS($afs_inputs,0.95,0.25)"
    afs_output = rcopy(afs_output_R)
    continuum = zeros(eltype(spectrum.flux),size(spectrum.flux[:,order_idx]))
    continuum = fill(NaN,size(spectrum.flux[:,order_idx]))
    continuum[pix] .= afs_output
    return continuum
end



using EchelleInstruments.NEID
function calc_continuum_model(spectrum::AbstractSpectra2D )
    vec_of_orders = pmap(ord->calc_continuum_model(spectrum,order_idx=ord), orders_to_use )
    output = fill(NaN, size(spectrum.flux))
    for (i,ord) in enumerate(min_order(get_inst(spectrum)):max_order(get_inst(spectrum)))
        output[:,ord] .= vec_of_orders[i]
    end
    return output
end


spec = all_spectra[1]
spec.flux = f_mean
continuum = calc_continuum_model(spec)
save(joinpath(pipeline_output_path_afw5465,date_subdir * "_continuum.jld2"), Dict("continuum"=>continuum,"orders_to_use"=>orders_to_use, "AFS_params"=>[0.95,0.25]))
