
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

combined_mask_df = mask_intersection(empirical_mask, VALD_mask_long_df, threshold=500.0)

combined_mask_pd = pd_df(combined_mask_df)

#combined_mask = py"mask_intersection"(empirical_mask_pd, VALD_mask_long, threshold=500.0)
#combined_mask_df = pd_df_to_df(combined_mask)

#combined_mask_df[!,:VALD_index] = convert.(Int,combined_mask_df[!,:VALD_index])

#the following 2 lines are commented out since telluric rejection was added to generateEmpiricalMask() so they are no longer needed
#telluric_indices = py"getTelluricIndices"(combined_mask_pd, true, Params[:overlap_cutoff], vel_slope_threshold=Params[:rejectTelluricSlope], RV_offset = 0.0, RV_range = 1e-4)
#combined_mask_df[!,:bool_filter_rejectTelluricSlope] = map(!,telluric_indices)


ns = String[]
ns2 = String[]
for n in names(combined_mask_df)
   if occursin("bool_filter_",n)
      push!(ns,n)
      push!(ns2,n[13:end])
   end
end

ns_short = ["min_frac_converged", "depth_5%_to_1", "std_width_quant", "std_slope_quant", "neg_value", "nan_value", "depth_cut", "allow_blends", "iron1", "telluric", "BLF"]

idx_used = [1,2,3,4,8,10]
ns = ns[idx_used]
ns2 = ns2[idx_used]
ns_short = ns_short[idx_used]

f_overlaps = zeros(Float64,length(ns),length(ns))
n_overlaps = zeros(Int64,length(ns),length(ns))

   
for (i,n) in enumerate(ns)
   idx_true = findall(skipmissing((combined_mask_df[:,n]) .== false))
   for (j,n2) in enumerate(ns)
      n_overlaps[j,i] = sum(skipmissing(combined_mask_df[idx_true,n2]) .== false)
      f_overlaps[j,i] = sum(skipmissing(combined_mask_df[idx_true,n2]) .== false) / (length(idx_true)+0.01)
   end
end

using Plots

plt = heatmap(ns_short,ns_short,f_overlaps,size=(1500,1000))
for i in 1:length(ns)
   for j in 1:length(ns)
      annotate!(i-0.5,j-0.5,n_overlaps[i,j])
   end
end
display(plt)
