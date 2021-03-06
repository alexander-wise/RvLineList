
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
   "--VALD_output"
      help = "whether or not to carry VALD line data through to the final mask"
      arg_type = Bool
end

   parsed_args = parse_args(s)


if parsed_args["allowBlends"] !== nothing
   println("Found command-line arg allowBlends. Overwriting param file definition of this arg.")
   global allowBlends = parsed_args["allowBlends"]
end

if parsed_args["overlap_cutoff"] !== nothing
   println("Found command-line arg overlap_cutoff. Overwriting param file definition of this arg.")
   global overlap_cutoff = parsed_args["overlap_cutoff"]
end

if parsed_args["rejectTelluricSlope"] !== nothing
   println("Found command-line arg rejectTelluricSlope. Overwriting param file definition of this arg.")
   global rejectTelluricSlope = parsed_args["rejectTelluricSlope"]
end

if parsed_args["badLineFilter"] !== nothing
   println("Found command-line arg badLineFilter. Overwriting param file definition of this arg.")
   global badLineFilter = parsed_args["badLineFilter"]
end

if parsed_args["quant"] !== nothing
   println("Found command-line arg quant. Overwriting param file definition of this arg.")
   global quant = parsed_args["quant"]
end

if parsed_args["nbin"] !== nothing
   println("Found command-line arg nbin. Overwriting param file definition of this arg.")
   global nbin = parsed_args["nbin"]
end

if parsed_args["output_dir"] !== nothing
   println("Found command-line arg output_dir. Overwriting param file definition of this arg.")
   global output_dir = parsed_args["output_dir"]
end

if parsed_args["VALD_output"] !== nothing
   println("Found command-line arg VALD_output. Overwriting param file definition of this arg.")
   global VALD_output = parsed_args["VALD_output"]
end


#generate empirical NEID mask

#set up output directories
if !isdir(joinpath(output_dir,"clean_masks"))
   mkdir(joinpath(output_dir,"clean_masks"))
end
if !isdir(joinpath(output_dir,"VALD_masks"))
   mkdir(joinpath(output_dir,"VALD_masks"))
end
if !isdir(joinpath(output_dir,"mask_bins"))
   mkdir(joinpath(output_dir,"mask_bins"))
end
if !isdir(joinpath(output_dir,"linefinder"))
   mkdir(joinpath(output_dir,"linefinder"))
end

pipeline_plan = RvLineList.PipelinePlan()
RvLineList.RvSpectMLBase.Pipeline.save_data!(pipeline_plan, :fit_lines) #tell the pipeline plan to save the line fits

@time empirical_mask = generateEmpiricalMask(output_dir=output_dir, pipeline_plan=pipeline_plan) #TODO: figure out why I can't run this with eval(param.jl) statement commented out in empriical_line_lists.jl (note param file was loaded earlier into global scope) - this eval statement resets quant

using PyCall
try
   pyimport("pandas") #make sure pandas is installed
catch e
   pyimport_conda("pandas", "pandas") #if the module fails to load, install it to Conda.jl environment
   pyimport("pandas") #make sure pandas is installed
end
import Pandas.DataFrame as pd_df

@pyinclude("src/make_VALD_line_list.py")
@time VALD_masks, VALD_masks_long = py"getVALDmasks"(overlap_cutoff=overlap_cutoff, depth_cutoff=depth_cutoff, iron1Only=iron1Only, badLineFilter=badLineFilter, allowBlends=allowBlends)
VALD_mask = VALD_masks[0]
VALD_mask_long = VALD_masks_long[0]


empirical_mask_3col = select(empirical_mask,[:median_??c, :median_depth, :line_id])
rename!(empirical_mask_3col, [:lambda, :depth, :line_id])
#debug empirical_mask_2col = CSV.read(joinpath(output_dir,"RvLineList_species=all_depthcutoff=0.05_overlapcutoff=1.0e-5_allowBlends=0_badLineFilter=none_rejectTelluricSlope=0.0_nbin=1_DP=true_binParam=depth_n=0_VACUUM.csv"), DataFrame)
#debug empirical_mask_2col = Dict(pairs(eachcol(empirical_mask_2col)))


empirical_mask_pd = pd_df(empirical_mask_3col)

combined_mask = py"mask_intersection"(empirical_mask_pd, VALD_output ? VALD_mask_long : VALD_mask, threshold=500.0)

function pd_df_to_df(df_pd)
   df = DataFrame()
   for col in df_pd.columns
      df[!,col] = getproperty(df_pd,col).values
   end
   df
end

combined_mask_df = pd_df_to_df(combined_mask)

telluric_indices = py"getTelluricIndices"(combined_mask, true, overlap_cutoff, vel_slope_threshold=rejectTelluricSlope, RV_offset = 0.0, RV_range = 1e-4)

mask = combined_mask_df[map(!,telluric_indices),:]

mask.weight = mask.depth

#binned_masks = binMask(rename(mask,[:lambda,:weight]), nbin, binParam=:weight)
binned_masks = binMask(mask, nbin, binParam=binParam)
#rename!(binned_masks[1], [:lambda, :depth])

for bin_n in 1:nbin
   saveStr = "RvLineList" * "_allowBlends="*string(allowBlends) * "_overlapcutoff="*string(overlap_cutoff) * "_rejectTelluricSlope="*string(rejectTelluricSlope) * "_badLineFilter="*badLineFilter* "_quant="*quant * "_nbin="*string(nbin) * "_DP="*string(depthPercentile) * "_binParam="*string(binParam) * "_n="*string(bin_n) * "_VALD_output="*string(VALD_output) * "_VACUUM" * ".csv"
   CSV.write(joinpath(output_dir,"clean_masks",saveStr),binned_masks[bin_n])
end

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=1e-5 --rejectTelluricSlope=2000 --badLineFilter="none", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=2e-5 --rejectTelluricSlope=2000 --badLineFilter="none", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=1e-5 --rejectTelluricSlope=2000 --badLineFilter="ESPRESSOG2", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=2e-5 --rejectTelluricSlope=2000 --badLineFilter="ESPRESSOG2", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"


"""

mask1 = RvLineList.read_mask_air(joinpath(pkgdir(RvLineList),"inputs","ESPRESSO_masks","G2.espresso.mas"))
mask2 = RvLineList.read_mask(joinpath("/home/awise/Desktop/neid_masks/clean_masks","RvLineList_allowBlends=0_overlapcutoff=2.0e-5_rejectTelluricSlope=2000.0_badLineFilter=none,_quant=90_nbin=1_DP=true_binParam=depth_n=1_VALD_output=false_VACUUM.csv"))

using Plots, JLD2, FileIO

neidSolar = load("/home/awise/data/neid/solar/2021/02/21/daily_ccfs_1.jld2")
lam = neidSolar["mean_lambda"]
flux = neidSolar["mean_clean_flux_continuum_normalized"]

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot()
for i in 1:122
   plot!(lam[:,i],flux[:,i], linecolor="yellow", legend=false)
end
xlims!(3500,10000)
ylims!(-0.1,1.1)
xlabel!("Wavelength (Angstroms)")
ylabel!("Normalized Flux")

for i in 1:(size(mask1,1)-1)
   plot!(rectangle(mask1[i,:lambda]*1e-5,mask1[i,:depth],mask1[i,:lambda],1-mask1[i,:depth]), opacity=.5, linecolor="red",fillcolor="red")
end
i=size(mask1,1)
display(plot!(rectangle(mask1[i,:lambda]*1e-5,mask1[i,:depth],mask1[i,:lambda],1-mask1[i,:depth]), opacity=.5, linecolor="red",fillcolor="red"))


for i in 1:(size(mask2,1)-1)
   plot!(rectangle(mask2[i,:lambda]*1e-5,mask2[i,:depth],mask2[i,:lambda],1-mask2[i,:depth]), opacity=.5, linecolor="blue",fillcolor="blue")
end
i=size(mask2,1)
display(plot!(rectangle(mask2[i,:lambda]*1e-5,mask2[i,:depth],mask2[i,:lambda],1-mask2[i,:depth]), opacity=.5, linecolor="blue",fillcolor="blue"))

savefig("/home/awise/Desktop/masks.png")

xlims!(5137,5145)
savefig("/home/awise/Desktop/masks_zoomed1.png")

xlims!(8600,8800)
savefig("/home/awise/Desktop/masks_zoomed2.png")

"""