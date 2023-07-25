
#using Pkg
# Pkg.activate(".")
using RvLineList, DataFrames, CSV, ArgParse

local_dir = pkgdir(RvLineList)
cd(local_dir)

#load param file
include(joinpath(local_dir,"inputs","param.jl"))

s = ArgParseSettings()

@add_arg_table s begin
   "--blend_RV_factors_filename"
      help = "filename for calculated blend RV factors. If the file is present, then these factors will be used instead of allowBlends and overlap_cutoff."
      arg_type = String
   "--blend_RV_cutoff"
      help = "cutoff value, in m/s per K, below which lines are considered not significantly contaminated by blends. Only used if :blend_RV_factors_filename is specified."
      arg_type = Float64
   "--allowBlends"
      help = "number of blends allowed with each line, e.g. 0 for single lines only, 1 for doublets only, [0,1] for singlets and doublets"
      arg_type = Int
   "--overlap_cutoff"
      help = "distance between lines required for them to be categorized as a blend, expressed as a fraction of the speed of light"
      arg_type = Float64
   "--depth_cutoff"
      help = "line depth, as a fraction of continuum, to be considered negligible. Note this parameter does not affect calculated blend RV factors."
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



if parsed_args["blend_RV_factors_filename"] !== nothing
   println("Found command-line arg blend_RV_factors_filename. Overwriting param file definition of this arg.")
   Params[:blend_RV_factors_filename] = parsed_args["blend_RV_factors_filename"]
end

if parsed_args["blend_RV_cutoff"] !== nothing
   println("Found command-line arg blend_RV_cutoff. Overwriting param file definition of this arg.")
   Params[:blend_RV_cutoff] = parsed_args["blend_RV_cutoff"]
end

if parsed_args["allowBlends"] !== nothing
   println("Found command-line arg allowBlends. Overwriting param file definition of this arg.")
   Params[:allowBlends] = parsed_args["allowBlends"]
end

if parsed_args["overlap_cutoff"] !== nothing
   println("Found command-line arg overlap_cutoff. Overwriting param file definition of this arg.")
   Params[:overlap_cutoff] = parsed_args["overlap_cutoff"]
end

if parsed_args["depth_cutoff"] !== nothing
   println("Found command-line arg depth_cutoff. Overwriting param file definition of this arg.")
   Params[:depth_cutoff] = parsed_args["depth_cutoff"]
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

@time empirical_mask = generateEmpiricalMask(Params, output_dir=Params[:output_dir], pipeline_plan=pipeline_plan)

#rename columns median_λc and median_depth to lambda and depth
rename!(empirical_mask, :median_λc => :lambda)
rename!(empirical_mask, :median_depth => :depth)

#filter empirical mask to only contain good lines according to empirical mask flags
#empirical_mask_filtered = empirical_mask[ (empirical_mask[:,:bool_filter_min_frac_converged]
#         .&& empirical_mask[:,:bool_filter_median_depth_between_5percent_and_1]
#         .&& empirical_mask[:,:bool_filter_std_velocity_width_quant]
#         .&& empirical_mask[:,:bool_filter_std_local_continuum_slope_quant]
#         .&& empirical_mask[:,:bool_filter_neg_bad_line] .&& empirical_mask[:,:bool_filter_nan_bad_line]), :]

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

#empirical_mask_3col = select(empirical_mask_filtered,[:lambda, :depth, :line_id])
#rename!(empirical_mask_3col, [:lambda, :depth, :line_id])
#debug empirical_mask_2col = CSV.read(joinpath(output_dir,"RvLineList_species=all_depthcutoff=0.05_overlapcutoff=1.0e-5_allowBlends=0_badLineFilter=none_rejectTelluricSlope=0.0_nbin=1_DP=true_binParam=depth_n=0_VACUUM.csv"), DataFrame)
#debug empirical_mask_2col = Dict(pairs(eachcol(empirical_mask_2col)))

#empirical_mask_pd = Params[:long_output] ? pd_df(empirical_mask) : pd_df(empirical_mask_3col) 
#empirical_mask_pd = pd_df(empirical_mask)

#combined_mask = py"mask_intersection"(empirical_mask_pd, Params[:long_output] ? VALD_mask_long : VALD_mask, threshold=500.0)
#combined_mask = py"mask_intersection"(empirical_mask_pd, VALD_mask_long, threshold=500.0)

"""
thresholds = [100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500,5000,5500,6000]
lengths = zeros(length(thresholds))
for i in 1:length(thresholds)
   @time lengths[i] = length(py"mask_intersection"(empirical_mask_pd, VALD_mask_long, threshold=thresholds[i]))
end
using Plots
plot(thresholds,lengths, legend=false)
xlims!(0,2000)
xlabel!("velocity threshold")
ylabel!("number of total mask matches found")
savefig("/home/awise/Desktop/mask_matching_cutoff.png")
"""

function pd_df_to_df(df_pd)
   df = DataFrame()
   for col in df_pd.columns
      df[!,col] = getproperty(df_pd,col).values
   end
   df
end

VALD_mask_long_df = pd_df_to_df(VALD_mask_long)
VALD_mask_long_df[!,:species] .= convert.(String,VALD_mask_long_df[!,:species])

#filter VALD mask by depth before merging with empirical mask to make sure the merging process is not affected by insignificant depth lines
VALD_mask_depth_filtered_df = VALD_mask_long_df[ VALD_mask_long_df[!,:bool_filter_depth_cutoff], :]

#merge VALD and empirical masks
combined_mask_df = mask_intersection(empirical_mask, VALD_mask_depth_filtered_df, threshold=Params[:matching_threshold])

combined_mask_pd = pd_df(combined_mask_df)

#combined_mask = py"mask_intersection"(empirical_mask_pd, VALD_mask_long, threshold=500.0)
#combined_mask_df = pd_df_to_df(combined_mask)

#combined_mask_df[!,:VALD_index] = convert.(Int,combined_mask_df[!,:VALD_index])

#the following 2 lines are commented out since telluric rejection was added to generateEmpiricalMask() so they are no longer needed
#telluric_indices = py"getTelluricIndices"(combined_mask_pd, true, Params[:overlap_cutoff], vel_slope_threshold=Params[:rejectTelluricSlope], RV_offset = 0.0, RV_range = 1e-4)
#combined_mask_df[!,:bool_filter_rejectTelluricSlope] = map(!,telluric_indices)


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


mask.weight = mask.depth

#binned_masks = binMask(rename(mask,[:lambda,:weight]), nbin, binParam=:weight)
binned_masks = binMask(mask, Params[:nbin], binParam=Params[:binParam])
#rename!(binned_masks[1], [:lambda, :depth])

for bin_n in 1:Params[:nbin]
   saveStr = "RvLineList" * "_allowBlends="*string(Params[:allowBlends]) * "_overlapCutoff="*string(Params[:overlap_cutoff]) * "_depthCutoff="*string(Params[:depth_cutoff]) * "_rejectTelluricSlope="*string(Params[:rejectTelluricSlope]) * "_badLineFilter="*Params[:badLineFilter]* "_quant="*Params[:quant] * "_nbin="*string(Params[:nbin]) * "_DP="*string(Params[:depthPercentile]) * "_binParam="*string(Params[:binParam]) * "_n="*string(bin_n) * "_long_output="*string(Params[:long_output]) * "_VACUUM" * ".csv"
   CSV.write(joinpath(Params[:output_dir],"clean_masks",saveStr),binned_masks[bin_n])
end

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=1e-5 --rejectTelluricSlope=2000 --badLineFilter="none", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=2e-5 --rejectTelluricSlope=2000 --badLineFilter="none", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=1e-5 --rejectTelluricSlope=2000 --badLineFilter="ESPRESSOG2", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"

#julia --project=RvSpectML/RvLineList RvSpectML/RvLineList/examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=2e-5 --rejectTelluricSlope=2000 --badLineFilter="ESPRESSOG2", --nbin=1 --output_dir="/home/awise/Desktop/neid_masks"


"""

mask1 = RvLineList.read_mask_air(joinpath(pkgdir(RvLineList),"inputs","ESPRESSO_masks","G2.espresso.mas"))
mask2 = RvLineList.read_mask(joinpath("/home/awise/Desktop/neid_masks/clean_masks","RvLineList_allowBlends=0_overlapcutoff=2.0e-5_rejectTelluricSlope=2000.0_badLineFilter=none,_quant=90_nbin=1_DP=true_binParam=depth_n=1_long_output=false_VACUUM.csv"))

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