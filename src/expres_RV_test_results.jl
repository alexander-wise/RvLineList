local_dir = "/home/awise/Desktop/RvSpectML"
#local_dir = "/storage/work/afw5465/RvSpectML"

using Pkg
 Pkg.activate(joinpath(local_dir,"alex_sandbox/"))
 using Query, DataFrames, DataStructures, CSV, JLD, Statistics, StatsBase, Plots, JLD2, FileIO, ArgParse
 using RvSpectMLBase, EchelleInstruments, EchelleCCFs, RvSpectML, Scalpels

 include(joinpath(local_dir,"data_paths.jl"))
 include(joinpath(local_dir,"RvLineList/inputs/param.jl"))

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

"""
ccfs = Dict()
rvs = Dict()
times = Dict()
"""

targs = ["101501", "26965", "10700", "34411"]
norm_types = [:raw, :continuum, :blaze]
EXPRES_excalibur_onlys = [true, false]
mask_types = [:tophat, :gaussian]
mask_scale_factors = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]
overlap_cutoffs = [2]
nbin=1
n=0

rms_rvs_espresso = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[], rv_rms=Float64[])
rms_rvs_espresso_nightly = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[], rv_rms=Float64[])
rms_rvs_espresso_within_night = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[], rv_rms=Float64[])
rms_rvs_empirical = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[], rv_rms=Float64[])
rms_rvs_empirical_nightly = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[], rv_rms=Float64[])
rms_rvs_empirical_within_night = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[], rv_rms=Float64[])


for norm_type in norm_types
   for EXPRES_excalibur_only in EXPRES_excalibur_onlys
      for mask_type in mask_types
         for mask_scale_factor in mask_scale_factors
            for overlap_cutoff in overlap_cutoffs
               results = load(joinpath(local_dir,"RvLineList","outputs",string("param_study_results","_norm_type=",norm_type,
               "_EXPRES_excalibur_only=",EXPRES_excalibur_only,"_mask_type=",mask_type,"_mask_scale_factor=",mask_scale_factor,
               "_overlap_cutoff=",overlap_cutoff,"_quant=",quant,"_nbin=",nbin,"_n=",n,".jld")))
               for targ in targs
                  #times = results["times"][targ]
                  #ccfs = results["ccfs"][targ]
                  rvs_espresso = results["rvs"][string("rvs_espresso_",targ)]
                  #rv_rms_espresso = std(rvs_espresso)
                  rv_rms_espresso = std(rvs_espresso .- mean(rvs_espresso))
                  rv_rms_espresso_nightly = std(bin_rvs_nightly(times=results["times"][targ],rvs=rvs_espresso .- mean(rvs_espresso)))
                  rv_rms_espresso_within_night = rms_rvs_within_night(times=results["times"][targ],rvs=rvs_espresso .- mean(rvs_espresso))
                  push!(rms_rvs_espresso,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,rv_rms_espresso])
                  push!(rms_rvs_espresso_nightly,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,rv_rms_espresso_nightly])
                  push!(rms_rvs_espresso_within_night,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,rv_rms_espresso_within_night])
                  rvs_empirical = results["rvs"][string("rvs_empirical_",targ)]
                  #rv_rms_empirical = std(rvs_empirical)
                  rv_rms_empirical = std(rvs_empirical .- mean(rvs_empirical))
                  rv_rms_empirical_nightly = std(bin_rvs_nightly(times=results["times"][targ],rvs=rvs_empirical .- mean(rvs_empirical)))
                  rv_rms_empirical_within_night = rms_rvs_within_night(times=results["times"][targ],rvs=rvs_empirical .- mean(rvs_empirical))
                  push!(rms_rvs_empirical,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,rv_rms_empirical])
                  push!(rms_rvs_empirical_nightly,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,rv_rms_empirical_nightly])
                  push!(rms_rvs_empirical_within_night,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,rv_rms_empirical_within_night])
               end
            end
         end
      end
   end
end
#commands used to visualize best RV params:
"""
i=1

sort(filter(row -> (row.targ == targs[i]) && (row.overlap_cutoff == 2), rms_rvs_espresso),:rv_rms)
sort(filter(row -> (row.targ == targs[i]) && (row.overlap_cutoff == 2), rms_rvs_espresso_nightly),:rv_rms)
sort(filter(row -> (row.targ == targs[i]) && (row.overlap_cutoff == 2), rms_rvs_espresso_within_night),:rv_rms)

sort(filter(row -> (row.targ == targs[i]) && (row.overlap_cutoff == 2), rms_rvs_empirical),:rv_rms)
sort(filter(row -> (row.targ == targs[i]) && (row.overlap_cutoff == 2), rms_rvs_empirical_nightly),:rv_rms)
sort(filter(row -> (row.targ == targs[i]) && (row.overlap_cutoff == 2), rms_rvs_empirical_within_night),:rv_rms)
"""
#index                     1        2        3        4
#best RMSs target:         101501,  26965,   10700,   34411
#best espresso RMSs:       3.44,    3.09,    1.94,    1.64
#espresso nightly:         3.36,    2.94,    1.19,    1.40
#espresso within-night:    0.51,    0.68,    1.32,    1.04

#best empirical RMSs:      4.14,    2.94,    1.69,    1.81
#empirical nightly:        4.08,    2.82,    0.91,    1.41
#empirical within-night:   0.60,    0.66,    1.33,    1.14

#HD101501 best ESPRESSO params: norm_type=:continuum, mask_type=:gaussian, mask_scale_factor=8.0
#     within-night, half the mask_scale_factor and norm_type=:raw

#all 3 other stars best ESPRESSO params: norm_type=:raw, mask_type=either, mask_scale_factor=1.0

#HD101501 empirical: continuum + excalibur_only is best, continum + not only excalibur is worst. gaussian > tophat, mask_scale_factor = max
#     within night, quarter the mask_scale_factor and norm_type=:raw
#HD26962 empirical: continuum + excalibur_only is best, mask_shape and mask_scale_factor have little effect
#     within night, raw is better than continuum
#HD10700 empirical: raw + excalibur_only is best, mask_shape has little effect and mask_scale_Factor = min
#     nightly, continuum now is better than raw. 
#HD34411 empirical: continuum + excalibur_only is best, mask_shape has little effect and mask_scale_Factor has little effect but should be min
#     within night, raw is better than continuum

#usually continuum is better, but for HD10700 empirical and for ESPRESSO masks raw is better, need to keep testing both
#excalibur_only is best, mask type: gaussian >= tophat, mask_scale_factor =1.0
#within night, raw is better (also overall for 10700 since it is dominated by within-night), while nightly, continuum is better

#neid worst model in their grid for which we would still use the observations: air mass at 2 hrs from noon on dec 23rd, highest barometric pressure
#correlate water column (measured from spectra) with RVs to test effect


targs = ["101501", "26965", "10700", "34411"]
norm_types = [:raw, :continuum]
EXPRES_excalibur_onlys = [true]
mask_types = [:gaussian]
mask_scale_factors = [0.5,1.0,1.5]
overlap_cutoffs = [1,2,3,4,6]
quants = ["50", "70", "90"]
nbins = [1]
ns = [0]
fwtf = 2.0

rms_rvs_espresso = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[],quant=String[],nbin=Int[],n=Int[], rv_rms=Float64[])
rms_rvs_espresso_nightly = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[],quant=String[],nbin=Int[],n=Int[], rv_rms=Float64[])
rms_rvs_espresso_within_night = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[],quant=String[],nbin=Int[],n=Int[], rv_rms=Float64[])
rms_rvs_empirical = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[],quant=String[],nbin=Int[],n=Int[], rv_rms=Float64[])
rms_rvs_empirical_nightly = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[],quant=String[],nbin=Int[],n=Int[], rv_rms=Float64[])
rms_rvs_empirical_within_night = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[],quant=String[],nbin=Int[],n=Int[], rv_rms=Float64[])


for norm_type in norm_types
   for EXPRES_excalibur_only in EXPRES_excalibur_onlys
      for mask_type in mask_types
         for mask_scale_factor in mask_scale_factors
            for overlap_cutoff in overlap_cutoffs
               for quant in quants
                  for nbin in nbins
                     for n in ns
                        try
                          global results = load(joinpath(local_dir,"RvLineList","outputs",string("param_study_results","_norm_type=",norm_type,
                           "_EXPRES_excalibur_only=",EXPRES_excalibur_only,"_mask_type=",mask_type,"_mask_scale_factor=",mask_scale_factor,
                           "_overlap_cutoff=",overlap_cutoff,"_quant=",quant,"_nbin=",nbin,"_n=",n,".jld")))
                        catch ArgumentError
                           global results = NaN
                        end
                        for targ in targs
                           #if typeof(results)==Float64
                           #   rv_rms_espresso = NaN
                           #   rv_rms_empirical = NaN
                           #else
                           times = results["times"][targ]
                           ccfs_espresso = results["ccfs"][string("ccfs_espresso_",targ)]
                           v_grid_espresso = results["v_grids"][string("v_grid_espresso_",targ)]
                           pipeline_plan_espresso = results["pipeline_plans"][string("pipeline_plan_espresso_",targ)]
                           alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
                           #try
                           rvs_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, pipeline_plan_espresso, v_grid=v_grid_espresso, times = times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
                           #rv_rms_espresso = std(rvs_espresso)
                           rv_rms_espresso = std(rvs_espresso .- mean(rvs_espresso))
                           rv_rms_espresso_nightly = std(bin_rvs_nightly(times=results["times"][targ],rvs=rvs_espresso .- mean(rvs_espresso)))
                           rv_rms_espresso_within_night = rms_rvs_within_night(times=results["times"][targ],rvs=rvs_espresso .- mean(rvs_espresso))
                           #catch
                           #   rv_rms_espresso = NaN
                           #end
                           ccfs_empirical = results["ccfs"][string("ccfs_empirical_",targ)]
                           v_grid_empirical = results["v_grids"][string("v_grid_empirical_",targ)]
                           pipeline_plan_empirical = results["pipeline_plans"][string("pipeline_plan_empirical_",targ)]
                           alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
                           #try
                           rvs_empirical = RvSpectML.calc_rvs_from_ccf_total(ccfs_empirical, pipeline_plan_empirical, v_grid=v_grid_empirical, times = times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
                           rv_rms_empirical = std(rvs_empirical .- mean(rvs_empirical))
                           rv_rms_empirical_nightly = std(bin_rvs_nightly(times=results["times"][targ],rvs=rvs_empirical .- mean(rvs_empirical)))
                           rv_rms_empirical_within_night = rms_rvs_within_night(times=results["times"][targ],rvs=rvs_empirical .- mean(rvs_empirical))
                           #catch
                           #   rv_rms_empirical = NaN
                           #end
                           #end
                           push!(rms_rvs_espresso,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,quant,nbin,n,rv_rms_espresso])
                           push!(rms_rvs_espresso_nightly,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,quant,nbin,n,rv_rms_espresso_nightly])
                           push!(rms_rvs_espresso_within_night,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,quant,nbin,n,rv_rms_espresso_within_night])
                           push!(rms_rvs_empirical,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,quant,nbin,n,rv_rms_empirical])
                           push!(rms_rvs_empirical_nightly,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,quant,nbin,n,rv_rms_empirical_nightly])
                           push!(rms_rvs_empirical_within_night,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,quant,nbin,n,rv_rms_empirical_within_night])
                        end
                     end
                  end
               end
            end
         end
      end
   end
end


#My hypothesis is that the best RVs to be in overlap_cutoff=1 and quant=90.
#Result: This is not generally true. Often overlap_cutoff=2 and/or quant=70 work better.
"""
i=1

sort(filter(row -> (row.targ == targs[i]), rms_rvs_empirical),:rv_rms)
sort(filter(row -> (row.targ == targs[i]), rms_rvs_empirical_nightly),:rv_rms)
sort(filter(row -> (row.targ == targs[i]), rms_rvs_empirical_within_night),:rv_rms)

"""

#index                     1        2        3        4
#best RMSs target:         101501,  26965,   10700,   34411
#best espresso RMSs:       3.44,    3.09,    1.94,    1.64
#espresso nightly:         3.36,    2.94,    1.19,    1.40
#espresso within-night:    0.51,    0.68,    1.32,    1.04

#round 1:
#best empirical RMSs:      4.14,    2.94,    1.69,    1.81
#empirical nightly:        4.08,    2.82,    0.91,    1.41
#empirical within-night:   0.60,    0.66,    1.33,    1.14

#round 2:
#best empirical RMSs:      4.2,     2.95,    1.60,    1.74
#empirical nightly:        3.95,    2.84,    0.86,    1.43
#empirical within-night:   0.50,    0.66,    1.33,    1.21

#HD101501 empirical: continuum, quant=50 is best; overlap_cutoff=6 helps a bit. Within-night, quant=70,overlap_cutoff=1 is best; raw helps a bit
#HD26962 empirical: continuum, quant=90, overlap_cutoff=1 is best. Within_night, same but raw is best
#HD10700 empirical: raw is best normally and within_night. overlap_cutoff = 2 to 5 normally and nightly, but 1 is okay for within_night. Quant always 90 or 70
#HD34411 empirical: overlap_cutoff=1 seems most important always. continuum is a bit better nightly, raw bit better within night





# Next, using overlap_cutoff=1 and quant=90, I compare summed empirical mask bins RV RMS to the same thing with
# 0, or 2 bins rejected
# problem: need to move the mask binning to AFTER mapping the mask onto the spectrum, so the bins are equally represented on the spectrum.


targs = ["101501", "26965", "10700", "34411"]
norm_types = [:raw, :continuum]
EXPRES_excalibur_onlys = [true]
mask_types = [:gaussian]
mask_scale_factors = [1.0, 4.0, 8.0]
overlap_cutoffs = [1, 6]
quants = ["90"]
nbins = [8]
ns = [1,2,3,4,5,6,7,8]
fwtf = 1.5

times = Dict()
ccfs = Dict()
v_grids = Dict()
pipeline_plans = Dict()

rms_rvs_empirical = DataFrame(targ=String[],norm_type=Symbol[],excalibur_only=Bool[],mask_type=Symbol[],mask_scale_factor=Float64[],overlap_cutoff=Int[],quant=String[],nbin=Int[],n=Int[], rv_rms=Float64[])

for norm_type in norm_types
   for EXPRES_excalibur_only in EXPRES_excalibur_onlys
      for mask_type in mask_types
         for mask_scale_factor in mask_scale_factors
            for overlap_cutoff in overlap_cutoffs
               for quant in quants
                  for nbin in nbins
                     for targ in targs
                        times[targ] = Dict()
                        ccfs[targ] = Dict()
                        v_grids[targ] = Dict()
                        pipeline_plans[targ] = Dict()
                        for n in ns
                           try
                              global results = load(joinpath(local_dir,"RvLineList","outputs",string("param_study_results","_norm_type=",norm_type,
                              "_EXPRES_excalibur_only=",EXPRES_excalibur_only,"_mask_type=",mask_type,"_mask_scale_factor=",mask_scale_factor,
                              "_overlap_cutoff=",overlap_cutoff,"_quant=",quant,"_nbin=",nbin,"_n=",n,".jld")))
                           catch ArgumentError
                              global results = NaN
                           end
                           if typeof(results)==Float64
                           else
                              times[targ][n] = results["times"][targ]
                              ccfs[targ][n] = results["ccfs"][string("ccfs_empirical_",targ)]
                              v_grids[targ][n] = results["v_grids"][string("v_grid_empirical_",targ)]
                              pipeline_plans[targ][n] = results["pipeline_plans"][string("pipeline_plan_empirical_",targ)]

                              alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
                              rvs_empirical = RvSpectML.calc_rvs_from_ccf_total(ccfs[targ][n], pipeline_plans[targ][n], v_grid=v_grids[targ][n], times = times[targ][n], recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
                              rv_rms_empirical = std(rvs_empirical .- mean(rvs_empirical))
                              push!(rms_rvs_empirical,[targ,norm_type,EXPRES_excalibur_only,mask_type,mask_scale_factor,overlap_cutoff,quant,nbin,n,rv_rms_empirical])
                           end
                        end
                     end
                  end
               end
            end
         end
      end
   end
end


rms_rvs = DataFrame(targ=String[], number_of_bins_rejected=Int[],rv_rms=Float64[])
rms_rvs_nightly = DataFrame(targ=String[], number_of_bins_rejected=Int[],rv_rms=Float64[])
rms_rvs_within_night1 = DataFrame(targ=String[], number_of_bins_rejected=Int[],rv_rms=Float64[])

for targ in targs
   n_sorted = sort(filter(row -> (row.targ == targ) && (row.overlap_cutoff == 6) && (row.norm_type == :continuum) && (row.mask_scale_factor == 1.0) && (row.quant == "90"), rms_rvs_empirical),:rv_rms).n
   for n_reject in 0:(length(ns) -1 )
      sum_ccfs = ccfs[targ][n_sorted[1]]
      for i in 2:(length(ns)-n_reject)
         try
            sum_ccfs += ccfs[targ][i]
         catch end
      end
      v_grid = v_grids[targ][n_sorted[1]]
      pipeline_plan = pipeline_plans[targ][n_sorted[1]]
      alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
      try
         rvs = RvSpectML.calc_rvs_from_ccf_total(sum_ccfs, pipeline_plan, v_grid=v_grid, times = times[targ][n_sorted[1]], recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
         rv_rms = std(rvs .- mean(rvs))
         rv_rms_nightly = std(bin_rvs_nightly(times=times[targ][1],rvs=rvs .- mean(rvs)))
         rv_rms_within_night = rms_rvs_within_night(times=times[targ][1],rvs=rvs .- mean(rvs))

         push!(rms_rvs,[targ,n_reject,rv_rms])      
         push!(rms_rvs_nightly,[targ,n_reject,rv_rms_nightly])      
         push!(rms_rvs_within_night1,[targ,n_reject,rv_rms_within_night])      
      catch
         rv_rms = NaN
      end
   end
end

 
"""
i=1

sort(filter(row -> (row.targ == targs[i]), rms_rvs),:rv_rms)

"""




#best RMSs target:    101501,  26965,   10700,   34411
#best espresso RMSs:  3.43,    3.09,    1.94,    1.66
#best empirical RMSs: 4.48,    3.59,    4.03,    2.69

#HD101501 best ESPRESSO params: norm_type=:continuum, mask_type=:gaussian, mask_scale_factor=8.0

#all 3 other stars best ESPRESSO params: norm_type=:raw, mask_type=either, mask_scale_factor=1.0

outputdir = joinpath(local_dir,"RvLineList","outputs")
filelist = readdir(outputdir)
filelist = filter(fn -> occursin("orders_to_use",fn), filelist)

rv_rms = DataFrame(targ=String[], norm_type=Symbol[], orders_to_use = Vector{Int64}[],
   rv_rms=Float64[], rv_rms_nightly=Float64[], rv_rms_within_night=Float64[])

for i in 1:length(filelist)
   for targ in targs
      fn = filelist[i]
      results = load(joinpath(outputdir,fn))
      norm_type = occursin("norm_type=raw",fn) ? :raw : :continuum
      orders_to_use = results["orders_to_use"]

      rvs_empirical = results["rvs"][string("rvs_empirical_",targ)]
      #rv_rms_empirical = std(rvs_empirical)
      rv_rms_empirical = std(rvs_empirical .- mean(rvs_empirical))
      rv_rms_empirical_nightly = std(bin_rvs_nightly(times=results["times"][targ],rvs=rvs_empirical .- mean(rvs_empirical)))
      rv_rms_empirical_within_night = rms_rvs_within_night(times=results["times"][targ],rvs=rvs_empirical .- mean(rvs_empirical))
      push!(rv_rms,[targ,norm_type,orders_to_use,rv_rms_empirical,rv_rms_empirical_nightly,rv_rms_empirical_within_night])
   end
end

rv_rms_continuum = filter(row -> (row.norm_type == :continuum), rv_rms)
rv_rms_raw = filter(row -> (row.norm_type == :raw), rv_rms)

"""
i=1

sort(filter(row -> (row.targ == targs[i]) & (row.norm_type == :continuum), rv_rms), :rv_rms_nightly)

"""

# OLD RESULTS:
#index                     1        2        3        4
#best RMSs target:         101501,  26965,   10700,   34411
#best empirical RMSs:      4.2,     2.95,    1.60,    1.74
#empirical nightly:        3.95,    2.84,    0.86,    1.43
#empirical within-night:   0.50,    0.66,    1.33,    1.21


#for each order in 43:72, find its average RMS, averaged over all of runs it appeared in

order_rms = DataFrame(targ=String[], order=Int64[], count=Int64[], rv_rms=Float64[], rv_rms_nightly=Float64[], rv_rms_within_night=Float64[])
for targ in targs
   global rv_rms_targ = filter(row -> (row.targ == targ), rv_rms_continuum)
   for o in 43:72
      for i in 1:48
         if o in rv_rms_targ[i,:orders_to_use]
            #if length(filter(row -> (row.order == o) && (row.targ == targ), order_rms)) >= 1
            #   order_rms[]
            push!(order_rms,[targ,o,1,rv_rms_targ[i,:rv_rms],rv_rms_targ[i,:rv_rms_nightly],rv_rms_targ[i,:rv_rms_within_night]])
         end
      end
   end
end
order_rms_avg_targ = DataFrame(targ=String[], order=Int64[], count=Int64[], rv_rms=Float64[], rv_rms_nightly=Float64[], rv_rms_within_night=Float64[])
for targ in targs
   for o in 43:72
      df_local = filter(row -> (row.targ == targ) && (row.order==o), order_rms)
      push!(order_rms_avg_targ,[targ,o,sum(df_local[:,:count]),mean(df_local[:,:rv_rms]),mean(df_local[:,:rv_rms_nightly]),mean(df_local[:,:rv_rms_within_night])])
   end
end

order_rms_avg = DataFrame(targ=String[], order=Int64[], count=Int64[], rv_rms=Float64[], rv_rms_nightly=Float64[], rv_rms_within_night=Float64[])
for o in 43:72
   df_local = filter(row ->  (row.order==o), order_rms)
   push!(order_rms_avg,["all",o,sum(df_local[:,:count]),mean(df_local[:,:rv_rms]),mean(df_local[:,:rv_rms_nightly]),mean(df_local[:,:rv_rms_within_night])])
end



"""
i=1

println(sort(filter(row -> (row.targ == targs[i]), order_rms_avg_targ), :rv_rms_within_night))

i=2 

println(sort(filter(row -> (row.targ == targs[i]), order_rms_avg_targ), :rv_rms_within_night))

i=3

println(sort(filter(row -> (row.targ == targs[i]), order_rms_avg_targ), :rv_rms_within_night))

i=4

println(sort(filter(row -> (row.targ == targs[i]), order_rms_avg_targ), :rv_rms_within_night))

"""


"""
i=1

println(sort( order_rms_avg, :rv_rms_within_night))

"""
