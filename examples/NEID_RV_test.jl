local_dir = "/home/awise/.julia/dev"
#local_dir = "/storage/work/afw5465/RvSpectML"

using Pkg
 Pkg.activate(joinpath(local_dir,"alex_sandbox/"))
 using Query, DataFrames, DataStructures, CSV, JLD, Statistics, StatsBase, Plots, JLD2, FileIO, ArgParse
 using RvSpectMLBase, EchelleInstruments, EchelleCCFs, RvSpectML, Scalpels, RvLineList

 include(joinpath(local_dir,"RvLineList/inputs/param.jl"))



linelist_names = SortedDict(Dict())

push!(linelist_names, "new_mask"=>joinpath("/home/awise/Desktop/neid_masks/clean_masks","RvLineList_allowBlends=0_overlapcutoff=2.0e-5_rejectTelluricSlope=2000.0_badLineFilter=none,_quant=90_nbin=1_DP=true_binParam=depth_n=1_VALD_output=false_VACUUM.csv"))
push!(linelist_names, "espresso"=>joinpath(pkgdir(RvLineList),"inputs","ESPRESSO_masks","G2.espresso.mas"))

mask1 = RvLineList.read_mask_air(linelist_names["espresso"])
mask2 = RvLineList.read_mask(linelist_names["new_mask"])
 


pipeline_plan = PipelinePlan()
reset_all_needs!(pipeline_plan)

verbose=true
if need_to(pipeline_plan,:read_spectra)
   if verbose println("# Finding what data files are avaliable.")  end
   global df_files = RvLineList.get_NEID_best_days(Params[:pipeline_output_summary_path],startDate=RvLineList.Date(2021,01,01), endDate=RvLineList.Date(2022,06,14), nBest=300)

   #if verbose println("# Reading in customized parameters from param.jl.")  end
   #include(joinpath("/home/awise/Desktop/RvSpectML","param.jl"))
   global ccf_mid_velocity = 0.0
   global df_files_use = df_files

   if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
   @time all_spectra = RvLineList.combine_NEID_daily_obs(Params[:daily_ccfs_base_path], Params[:daily_ccf_fn], Params[:daily_manifests_base_path], Params[:daily_manifest_fn], df_files_use)
   #for i in 1:length(all_spectra)
   #  all_spectra[i].metadata[:bjd] = all_spectra[i].metadata[:BARYMJD]
   #end
   #RvSpectMLBase.discard_blaze(all_spectra)
   RvSpectMLBase.discard_continuum(all_spectra)
   dont_need_to!(pipeline_plan,:read_spectra)
end

# ccf_mid_velocity = 0.0 #since the empirical mask is at v=0, and we shifted the VALD masks to be at v=0 too, we must reset this global variable to zero.

orders_to_use = Params[:orders_to_use]
pixels_to_use = map(ord->min_col_default(all_spectra[1].inst,ord):max_col_default(all_spectra[1].inst,ord),orders_to_use)
order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=orders_to_use, verbose=true)
n_spectra = length(order_list_timeseries)
n_orders = length(orders_to_use)

### NEW RUN: empirical mask
# ccf_mid_velocity = 0.0 #since the empirical mask is at v=0, and we shifted the VALD masks to be at v=0 too, we must reset this global variable to zero.

ll = string("espresso")
#ll = "espresso"
ll_v = string("new_mask")

line_width_50 = 7000.

need_to!(pipeline_plan,:read_line_list)
need_to!(pipeline_plan,:clean_line_list_tellurics)
empty!(pipeline_plan.cache)
println(linelist_names[ll])
linelist_path = joinpath(local_dir,"RvLineList","inputs","empirical_masks",linelist_names[ll])
mask_ll = prepare_line_list(linelist_path, all_spectra, pipeline_plan,  v_center_to_avoid_tellurics=0.0,
   Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()),  orders_to_use=orders_to_use, verbose=true,
   convert_air_to_vacuum=true)

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


dont_need_to!(pipeline_plan,:clean_line_list_tellurics)

lsf_width = 2.0e3 #estimated line width for telluric avoidance
mask_type = :gaussian #mask shape. Options: :tophat, :gaussian, :supergaussian, :halfcos
mask_scale_factor = 4.0 #approximate number of pixels in the width scale of the mask
fwtf = 2.0 #fraction of the FWHM of the CCF to fit
RV_fit_type = :gaussian #type of fit to use on CCF to measure RV. Options: :quadratic, :gaussian
tophap_ccf_mask_scale_factor=1.6


println("# mask_scale_factor = ", mask_scale_factor, ", frac_width_to_fit =  ", fwtf, ". ")
#@time (ccfs_empirical, v_grid) = ccf_total(order_list_timeseries, mask_ll, pipeline_plan,  mask_scale_factor=mask_scale_factor,
#   range_no_mask_change=5*line_width_50, ccf_mid_velocity=0.0, v_step=100,
#   v_max= 0*RvSpectMLBase.max_bc+5*line_width_50, #+2*line_width_50,
#   calc_ccf_var=false, recalc=true, mask_type = mask_type)

@time (order_ccfs, order_ccf_vars, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, mask_ll, pipeline_plan, mask_type=:gaussian,
   #Δfwhm=Δfwhm,
   mask_scale_factor=mask_scale_factor, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=100,
   #v_max=max(range_no_mask_change*line_width_50,2*max_bc),
   v_max= 0*RvSpectMLBase.max_bc+5*line_width_50, #+2*line_width_50
   orders_to_use=orders_to_use, allow_nans=true, calc_ccf_var=true, recalc=true)

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