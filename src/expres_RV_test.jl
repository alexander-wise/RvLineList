local_dir = "/home/awise/Desktop/RvSpectML"
#local_dir = "/storage/work/afw5465/RvSpectML"

using Pkg
 Pkg.activate(joinpath(local_dir,"alex_sandbox/"))
 using Query, DataFrames, DataStructures, CSV, JLD, Statistics, StatsBase, Plots, JLD2, FileIO, ArgParse
 using RvSpectMLBase, EchelleInstruments, EchelleCCFs, RvSpectML, Scalpels

 include(joinpath(local_dir,"data_paths.jl"))
 include(joinpath(local_dir,"RvLineList/inputs/param.jl"))
 include(joinpath(local_dir,"RvLineList/src/masks.jl"))

 s = ArgParseSettings()

 @add_arg_table s begin
    "--norm_type"
       help = "normalization type. Options: :raw, :continuum, :blaze"
       arg_type = Symbol
    "--EXPRES_excalibur_only"
       help = "whether or not to only include excalibur wavelengths in expres data - could be replaced with orders_to_use param here"
       arg_type = Bool
   "--mask_type"
      help = "mask shape. Options: :tophat, :gaussian, :supergaussian, :halfcos"
      arg_type = Symbol
   "--mask_scale_factor"
      help = "approximate number of pixels in the width scale of the mask"
      arg_type = Float64
   "--overlap_cutoff"
      help = "distance between lines required for them to be categorized as a blend, expressed as a fraction of the speed of light x 10^-5"
      arg_type = Int
   "--quant"
      help = "quantile for stability of fit params for empirical masks"
      arg_type = String
   "--nbin"
      help = "number of mask subdivisions to make"
      arg_type = Int
   "--bin_n"
      help = "which subdivision from nbin to use - one is first as this is a julia index"
      arg_type = Int
 end
 
 parsed_args = parse_args(s)
 
if parsed_args["norm_type"] != nothing
   println("Found command-line arg norm_type. Overwriting param file definition of this arg.")
   global norm_type = parsed_args["norm_type"]
end

if parsed_args["EXPRES_excalibur_only"] != nothing
   println("Found command-line arg EXPRES_excalibur_only. Overwriting param file definition of this arg.")
   global EXPRES_excalibur_only = parsed_args["EXPRES_excalibur_only"]
end

if parsed_args["mask_type"] != nothing
   println("Found command-line arg mask_type. Overwriting param file definition of this arg.")
   global mask_type = parsed_args["mask_type"]
end

if parsed_args["mask_scale_factor"] != nothing
   println("Found command-line arg mask_scale_factor. Overwriting param file definition of this arg.")
   global mask_scale_factor = parsed_args["mask_scale_factor"]
end

if parsed_args["overlap_cutoff"] != nothing
   println("Found command-line arg overlap_cutoff. Overwriting param file definition of this arg.")
   global overlap_cutoff = parsed_args["overlap_cutoff"]
end

if parsed_args["quant"] != nothing
   println("Found command-line arg quant. Overwriting param file definition of this arg.")
   global quant = parsed_args["quant"]
end

if parsed_args["nbin"] != nothing
   println("Found command-line arg nbin. Overwriting param file definition of this arg.")
   global nbin = parsed_args["nbin"]
end

if parsed_args["bin_n"] != nothing
   println("Found command-line arg bin_n. Overwriting param file definition of this arg.")
   global bin_n = parsed_args["bin_n"]
end


function hasOneToOneMatch(X::A1,Y::A2) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}}
   xTrue = zeros(Bool,size(X,1))
   for y in 1:length(Y)
      n_matches = 0
      match_i = -1
      for x in 1:length(X)
         if (abs(X[x]-Y[y])/X[x])<1e-5
            n_matches += 1
            match_i = x
         end
      end
      if n_matches == 1
         xTrue[match_i]=true
      end
   end
   return xTrue
end

#norm_type_results = Dict()

ccfs = Dict()

v_grids = Dict()

pipeline_plans = Dict()

rvs = Dict()

times = Dict()

targs = ["101501", "26965", "34411", "10700"]

for targ in targs

   GC.gc()

   global target = targ

   maskWavelengths = string("Reiners",target) #keyword for mask wavelengths - should be e.g. Reiners or Reiners101501
   badLineFilter = string("HD",target)

   linelist_names = SortedDict(Dict())
 
   push!(linelist_names, "empirical_mask"=>string(target,"_good_lines_fit_quant=0.",quant,"_frac=",min_frac_converged,"_fmt=espresso.mas"))
 
   VALDmaskStr = string("VALD", depthPercentile ? "" : "_nonDP", "_species=", iron1Only, "_depthcutoff=", depth_cutoff, "_overlapcutoff=", overlap_cutoff, "e-05_allowBlends=", allowBlends, "_badLineFilter=", badLineFilter, "_rejectTelluricSlope=", rejectTelluricSlope, "_waves=", maskWavelengths, "_nbin=", 1, "_binParam=", binParam, "_n=", 0, ".mas")
   push!(linelist_names, "vald_mask"=>VALDmaskStr)

   push!(linelist_names, "espresso"=>"G2.espresso.mas")

   global fits_target_str = target
   verbose=true
   #path_to_EXPRES_data = joinpath("/home/awise/data/",target)
   path_to_EXPRES_data = joinpath(expres_data_path,target)
   default_paths_to_search = [local_dir]
   pipeline_plan = PipelinePlan()
   dont_make_plot!(pipeline_plan, :movie)

   reset_all_needs!(pipeline_plan)
   if need_to(pipeline_plan,:read_spectra)
      if verbose println("# Finding what data files are avaliable.")  end
      global df_files = EXPRES.make_manifest(path_to_EXPRES_data, fits_target_str)

      if verbose println("# Reading in customized parameters from param.jl.")  end
      include(joinpath(local_dir,"param.jl"))

      if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
      @time all_spectra = map(row->EXPRES.read_data(row,store_min_data=true, store_tellurics=true, normalization=norm_type), eachrow(df_files_use))
      #for i in 1:length(all_spectra)
      #  all_spectra[i].metadata[:bjd] = all_spectra[i].metadata[:BARYMJD]
      #end
      #RvSpectMLBase.discard_blaze(all_spectra)
      RvSpectMLBase.discard_continuum(all_spectra)
      dont_need_to!(pipeline_plan,:read_spectra)
   end

   # ccf_mid_velocity = 0.0 #since the empirical mask is at v=0, and we shifted the VALD masks to be at v=0 too, we must reset this global variable to zero.


   pixels_to_use = map(ord->min_col_default(all_spectra[1].inst,ord):max_col_default(all_spectra[1].inst,ord),orders_to_use)
   order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=orders_to_use);
   n_spectra = length(order_list_timeseries)
   n_orders = length(orders_to_use)

   ### NEW RUN: empirical mask
   # ccf_mid_velocity = 0.0 #since the empirical mask is at v=0, and we shifted the VALD masks to be at v=0 too, we must reset this global variable to zero.

   ll = string("empirical_mask")
   #ll = "espresso"
   ll_v = string("vald_mask")

   line_width_50 = 7000.

   need_to!(pipeline_plan,:read_line_list)
   need_to!(pipeline_plan,:clean_line_list_tellurics)
   empty!(pipeline_plan.cache)
   println(linelist_names[ll])
   linelist_path = joinpath(local_dir,"RvLineList","inputs","empirical_masks",linelist_names[ll])
   mask_ll = prepare_line_list(linelist_path, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=0.0,
      Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()),  orders_to_use=orders_to_use, verbose=true,
      convert_air_to_vacuum=false)

   need_to!(pipeline_plan,:read_line_list)
   need_to!(pipeline_plan,:clean_line_list_tellurics)
   empty!(pipeline_plan.cache)
   println(linelist_names[ll_v])
   linelist_path = joinpath(local_dir,"RvLineList","outputs","masks",linelist_names[ll_v])
   mask_ll_v = prepare_line_list(linelist_path, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=0.0,
      Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()),  orders_to_use=orders_to_use, verbose=true,
      convert_air_to_vacuum=false)

   nf = length(mask_ll.lambda)
   # hm = zeros(Bool,nf)
   # for i in 1:nf
   #   hm[i] = hasOneMatch(mask_ll.lambda[i],mask_ll_v.lambda)
      #hm[i] = hasOneMatch(linelists[ll].lambda[i]*(1.0-5e3/3e8),linelists[ll_v].lambda)
   hm = hasOneToOneMatch(mask_ll.lambda,mask_ll_v.lambda)
   #end
   mask_empirical_total = mask_ll[hm,:]

   masks_empirical = binMask(mask_empirical_total,nbin)
   mask_empirical = masks_empirical[bin_n]

   println("# mask_scale_factor = ", mask_scale_factor, ", frac_width_to_fit =  ", fwtf, ". ")
   @time (ccfs_empirical, v_grid) = ccf_total(order_list_timeseries, mask_empirical, pipeline_plan,  mask_scale_factor=mask_scale_factor, range_no_mask_change=5*line_width_50, ccf_mid_velocity=0.0, v_step=100,
      #v_max=RvSpectMLBase.max_bc,
      v_max= 0*RvSpectMLBase.max_bc+5*line_width_50, #+2*line_width_50,
      calc_ccf_var=false, recalc=true, mask_type = mask_type)
   push!(ccfs, string("ccfs_empirical_",target) => ccfs_empirical)
   push!(v_grids, string("v_grid_empirical_",target) => v_grid)
   push!(pipeline_plans, string("pipeline_plan_empirical_",target) => pipeline_plan)
   #if nbin == 1 #if nbin > 1, rv fits sometimes fail
   alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
   rvs_ccf_empirical = RvSpectML.calc_rvs_from_ccf_total(ccfs_empirical, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
   push!(rvs, string("rvs_empirical_",target) => rvs_ccf_empirical)
   #end



   k="espresso"
   need_to!(pipeline_plan,:read_line_list)
   need_to!(pipeline_plan,:clean_line_list_tellurics)
   empty!(pipeline_plan.cache)
   println(linelist_names[k])
   linelist_path = joinpath("alex_sandbox","masks",linelist_names[k])
   line_list_df = prepare_line_list(linelist_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
   Δv_to_avoid_tellurics = 2*RvSpectMLBase.max_bc+3*lsf_width, verbose=true)

   masks_espresso = binMask(line_list_df,nbin)
   mask_espresso = masks_espresso[bin_n]

   println("# mask_scale_factor = ", mask_scale_factor, ", frac_width_to_fit =  ", fwtf, ". ")
   @time (ccfs_espresso, v_grid) = ccf_total(order_list_timeseries, mask_espresso, pipeline_plan,  mask_scale_factor=mask_scale_factor, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=100,
      #v_max=RvSpectMLBase.max_bc,
      v_max= 0*RvSpectMLBase.max_bc+5*line_width_50, #+2*line_width_50,
      calc_ccf_var=false, recalc=true, mask_type = mask_type)
   push!(ccfs, string("ccfs_espresso_",target) => ccfs_espresso)
   push!(v_grids, string("v_grid_espresso_",target) => v_grid)
   push!(pipeline_plans, string("pipeline_plan_espresso_",target) => pipeline_plan)
   #if nbin == 1 #if nbin > 1, rv fits sometimes fail
   alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
   rvs_ccf_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
   push!(rvs, string("rvs_espresso_",target) => rvs_ccf_espresso)
   #end


   times[target] = order_list_timeseries.times

end

GC.gc()

save(joinpath(local_dir,"RvLineList","outputs",string("param_study_results","_norm_type=",norm_type,
   "_EXPRES_excalibur_only=",EXPRES_excalibur_only,"_mask_type=",mask_type,"_mask_scale_factor=",mask_scale_factor,
   "_overlap_cutoff=",overlap_cutoff,"_quant=",quant,"_nbin=",nbin,"_bin_n=",bin_n,"_orders_to_use=",join(orders_to_use),".jld")),
   "ccfs",ccfs,"v_grids",v_grids,"pipeline_plans",pipeline_plans,"rvs",rvs,"times",times,"orders_to_use",orders_to_use)

"""
sample job submission scripts:

for norm_type in raw blaze continuum; do for EXPRES_excalibur_only in true false; do for mask_type in gaussian tophat; do for mask_scale_factor in 1 2 3 4 5 6 7 8; do for overlap_cutoff in 2; do qsub -v norm_type=$norm_type,EXPRES_excalibur_only=$EXPRES_excalibur_only,mask_type=$mask_type,mask_scale_factor=$mask_scale_factor,overlap_cutoff=$overlap_cutoff expres_mask_param_study.pbs; done; done; done; done; done

for norm_type in raw continuum; do for EXPRES_excalibur_only in true; do for mask_type in gaussian; do for mask_scale_factor in 0.5 1 1.5; do for overlap_cutoff in 1 2 3 4 6; do for quant in 50 70 90; do for nbin in 1; do for bin_n in 0; do qsub -v norm_type=$norm_type,EXPRES_excalibur_only=$EXPRES_excalibur_only,mask_type=$mask_type,mask_scale_factor=$mask_scale_factor,overlap_cutoff=$overlap_cutoff,quant=$quant,nbin=$nbin,n=$n expres_mask_param_study.pbs; done; done; done; done; done; done; done; done

for norm_type in raw continuum; do for EXPRES_excalibur_only in true; do for mask_type in gaussian; do for mask_scale_factor in 1 4 8; do for overlap_cutoff in 1 6; do for quant in 90; do for nbin in 8; do for bin_n in 1 2 3 4 5 6 7 8; do qsub -v norm_type=$norm_type,EXPRES_excalibur_only=$EXPRES_excalibur_only,mask_type=$mask_type,mask_scale_factor=$mask_scale_factor,overlap_cutoff=$overlap_cutoff,quant=$quant,nbin=$nbin,n=$n expres_mask_param_study.pbs; done; done; done; done; done; done; done; done
   
for i in $(seq 1 50); do for norm_type in raw continuum; do for EXPRES_excalibur_only in true; do for mask_type in gaussian; do for mask_scale_factor in 1; do for overlap_cutoff in 1; do for quant in 90; do for nbin in 1; do for bin_n in 1; do qsub -v norm_type=$norm_type,EXPRES_excalibur_only=$EXPRES_excalibur_only,mask_type=$mask_type,mask_scale_factor=$mask_scale_factor,overlap_cutoff=$overlap_cutoff,quant=$quant,nbin=$nbin,n=$n expres_mask_param_study.pbs; done; done; done; done; done; done; done; done; done

"""