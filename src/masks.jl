# functions to generate and manipulate masks


""" Read mask containing air wavelengths in csv format.
   input: filename for mask with columns: `lambda` in air and `depth`
   output: DataFrame containing the line wavlengths in vacuum and line depths
"""
function read_mask_air(fn::String)
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","depth"])
    @assert issorted(df[!,:lambda])
    df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda]) #convert air to vacuum wavelength
    return df
end

""" Read mask containing vacuum wavelengths in csv format.
   input: filename for mask with columns: `lambda` in vacuum and `depth`
   output: DataFrame containing the line wavlengths in vacuum and line depths
"""
function read_mask_vacuum(fn::String)
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","depth"])
    @assert issorted(df[!,:lambda])
    return df
end

""" Read mask in csv format. Filename must include the substring _AIR" or "_VACUUM" denoting wavelengths medium
   input: filename for mask with columns: `lambda` in vacuum and `depth`
   output: DataFrame containing the line wavlengths in vacuum and line depths
"""
function read_mask(fn::String)
    @assert occursin("_VACUUM",fn) || occursin("_AIR",fn)
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","depth", "weight"],skipto=2)
    @assert issorted(df[!,:lambda])
    if occursin("_AIR",fn)
      df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda]) #convert air to vacuum wavelength
    end
    return df
end

function get_airVacString(inst::Symbol)
   local airVacString
   if inst in (:neid, :expres)
      airVacString="_VACUUM"
   elseif inst in (:harps, :harpsn)
      airVacString="_AIR"
   end
   return airVacString
end

function binMask(mask::DataFrame, nbin::Int; binParam::Symbol = :weight, depthPercentile::Bool = true)
    bins = Dict()
    if depthPercentile #this should be updated to include blaze/normalization/some sort of SNR
        pSorted = sortperm(mask[:,binParam])
        dCumSum = cumsum(mask[pSorted,:weight])
        for i in 1:nbin
            dLower = dCumSum[end]/nbin * (i-1)
            dUpper = dCumSum[end]/nbin * i
            push!(bins, i=>pSorted[(dCumSum .>= dLower) .& (dCumSum .<= dUpper)])
        end
    else
        for i in 1:nbin
            pLower = percentile(mask[:,binParam], 100/nbin * (i-1))
            pUpper = percentile(mask[:,binParam], 100/nbin * i)
            push!(bins, i=>findall((mask[:,binParam] .>= pLower) .& (mask[:,binParam] .< pUpper)))
            if i==nbin
                append!(bins[i],argmax(mask[:,binParam]))
                sort!(bins[8])
            end
        end
    end
    return [sort(mask[bins[i],:], :lambda) for i in 1:nbin]
end

function get_param_range(mask::DataFrame, param::Symbol, param_min::Number, param_max::Number)
   @assert String(param) in names(mask)
   @assert param_min < param_max
   indices = sortperm(mask[!,param])
   index_min = findfirst(mask[indices,param] .> param_min)
   index_max = findlast(mask[indices,param] .< param_max)
   view(mask,sort(indices[index_min:index_max]), :)
end

function get_depth_range(mask::DataFrame, depth_min::Number, depth_max::Number)
   get_param_range(mask, :depth, depth_min, depth_max)
end

function get_lambda_range(mask::DataFrame, lambda_min::Number, lambda_max::Number)
   get_param_range(mask, :lambda, lambda_min, lambda_max)
end




#check if two mask wavelengths are within 50 m/s (approx 0.001 angstroms)
#VALD numerical precision is 0.0001 angstroms, or about 5 m/s. This number (50 m/s) was chosen to be a factor of a few greater than VALD precision at all wavelengths we consider.
function wave_equal(w1, w2; threshold::Number=50.0)
   if isequal(w1,missing) || isequal(w2,missing)
      return false
   else
      return (abs(w1-w2)/((w1+w2)/2)) < (threshold / C_m_s)
   end
end



#take the intersection of two masks (pandas 2d dataframes) to generate one sub_mask
function mask_intersection(mask1::DataFrame, mask2::DataFrame; default_data::String="first", add_label_column::Bool=false, combine_mask_data::Bool=true, threshold::Number=50.0)
   """
   Take the intersection of two masks to generate one submask. Two mask entries are treated as equal if they have wavelengths (lambdas) within threshold (50 m/s by default) of each other, or if they are both equal to the same third mask entry.

   Parameters:
      mask1: mask #1 including columns for lambda (required) and depth (required if default_data="max_depth")

      mask2: mask #2 including columns for lambda (required) and depth (required if default_data="max_depth")

      default_data: keyword for which mask data to preserve in the case of matching (intersecting) lines - all but one are discarded.
         Posible values:
            "first" defaults to mask1's data for matches across masks. Will default to max_depth for matches within masks.
            "max_depth" picks the line with max depth in the case of a match. If multiple lines with the same max depth, defaults to the lower lambda line within the match.

      add_label_column: whether or not to add a column, "mask_df_name", to the output supermask containing labels "mask1" and "mask2" tracking which input mask each output mask entry originated from.

      combine_mask_data: whether or not to add extra data in mask2 to mask1 (only works when default_data=="first")

      threshold: velocity (in m/s) separation between adjacent mask entries for them to be considered equal.

   Returns:
      sub_mask: intersection of mask1 and mask2, sorted by lambda
   """
   @assert (isequal(default_data,"first") || isequal(default_data,"max_depth")) "ERROR: invalid default_data selection."

   @assert (("lambda" in names(mask1)) && ("lambda" in names(mask2))) "ERROR: lambda column not found."
   if ("depth" in names(mask1)) && ("depth" in names(mask2))
      hasDepths=true
   else
      @assert (default_data != "max_depth") "ERROR: max_depth selected but depth keyword missing from mask1.columns or mask2.columns."
      hasDepths=false
   end
   mask1sorted = sort(mask1,[:lambda])
   mask2sorted = sort(mask2,[:lambda])
   mask1sorted[!,"mask_df_name"] .= "mask1"
   mask2sorted[!,"mask_df_name"] .= "mask2"
   sub_mask = DataFrame()
   mask12combined = vcat(mask1sorted,mask2sorted, cols=:union)
   mask12sorted = sort(mask12combined, [:lambda])
   m12 = mask12sorted[!,:lambda]
   i=1
   while i < (length(m12)) #loop through combined mask index i
      j=1
      match_i = [i]
      while wave_equal(m12[i+j-1],m12[i+j],threshold=threshold)
         append!(match_i,i+j)
         j+=1
      end
      if j > 1 #if there are any matches to mask index i
         matches = mask12sorted[match_i,:]
         mask1matches = matches[:,"mask_df_name"] .== "mask1"
         mask2matches = matches[:,"mask_df_name"] .== "mask2"
         if (any(mask1matches) && any(mask2matches))
            if default_data == "first"
               matches1 = matches[mask1matches,:]
               matches2 = matches[mask2matches,:]
               max_depth_match1 = DataFrame(matches1[findmax(matches1[:,:depth])[2],:])
               max_depth_match2 = DataFrame(matches2[findmax(matches2[:,:depth])[2],:])
               if combine_mask_data
                 for k in names(mask2)
                     if ismissing(max_depth_match1[1,k])
                        max_depth_match1[!,k] .= max_depth_match2[1,k]
                     end
                  end
               end
               append!(sub_mask,max_depth_match1)
            else
               max_depth_match = DataFrame(matches[findmax(matches[:,:depth])[2],:])
               append!(sub_mask,max_depth_match)
            end
         end
      end
      i += j
   end
   if add_label_column
      sub_mask_out = sub_mask
   else
      sub_mask_out = select!(sub_mask, Not([:mask_df_name]))
   end
   return sort(sub_mask_out,[:lambda])
end

#take the intersection of two masks (pandas 2d dataframes) to generate one sub_mask
function match_VALD_to_empirical(empirical_mask::DataFrame, VALD_mask::DataFrame; add_label_column::Bool=false, combine_mask_data::Bool=true, threshold::Number=500.0)
   """
   Take an empirical mask, and find nearby VALD matches. They are considered a match wavelengths (lambdas) within threshold (500 m/s by default) of each other.

   Parameters:
      empirical_mask: empirical mask including column for lambda

      VALD_mask: VALD mask including column for lambda and depth

      add_label_column: whether or not to add a column, "mask_df_name", to the output supermask containing labels "empirical_mask" and "VALD_mask" tracking which input mask each output mask entry originated from.

      combine_mask_data: whether or not to add extra data in VALD_mask to empirical_mask

      threshold: velocity (in m/s) separation between adjacent mask entries for them to be considered a match.

   Returns:
      sub_mask: intersection of empirical_mask and VALD_mask, sorted by lambda
   """

   @assert (("lambda" in names(empirical_mask)) && ("lambda" in names(VALD_mask))) "ERROR: lambda column not found."
   empirical_mask_sorted = sort(empirical_mask,[:lambda])
   VALD_mask_sorted = sort(VALD_mask,[:lambda])
   empirical_mask_sorted[!,"mask_df_name"] .= "empirical_mask"
   VALD_mask_sorted[!,"mask_df_name"] .= "VALD_mask"
   sub_mask = DataFrame()
   #masks_combined = vcat(empirical_mask_sorted,VALD_mask_sorted, cols=:union)
   #mask12sorted = sort(mask12combined, [:lambda])
   #m12 = mask12sorted[!,:lambda]
   j0=1
   j=1
   for i in 1:size(empirical_mask)[1]
      j=findfirst(wave_equal.(VALD_mask_sorted[j0:end,:lambda],empirical_mask_sorted[i,:lambda],threshold=threshold))
      if ~isnothing(j)
         j += j0 - 1
         j0 = j
         match_j = [j]
         while wave_equal(VALD_mask_sorted[j+1,:lambda],empirical_mask_sorted[i,:lambda],threshold=threshold)
            append!(match_j,j+1)
            j+=1
         end
         matches = VALD_mask_sorted[match_j,:]
         empirical_line = DataFrame(empirical_mask_sorted[[i],:])
         max_depth_VALD_match = DataFrame(matches[findmax(matches[:,:depth])[2],:])
         if combine_mask_data
            for k in names(VALD_mask)
               if ~(k in names(empirical_line))
                  empirical_line[!,k] .= max_depth_VALD_match[1,k]
               end
            end
         end
         append!(sub_mask,empirical_line)
      end
   end
   if add_label_column
      sub_mask_out = sub_mask
   else
      sub_mask_out = select!(sub_mask, Not([:mask_df_name]))
   end
   sub_mask_out = sort(sub_mask_out,[:lambda])
   #remove duplicate mask entries due to the same empirical line showing up on multiple orders
   mask_out = DataFrame()
   i0=1
   i1=i0+1
   while i0 <= size(sub_mask_out)[1]
      while (i1 <= size(sub_mask_out)[1]) && wave_equal(sub_mask_out[i0,:lambda], sub_mask_out[i1,:lambda], threshold=threshold)
         i1+=1
      end
      matches = DataFrame(sub_mask_out[i0:i1-1,:])
      min_var_match = DataFrame(matches[findmin(matches[:,:mean_template_var])[2],:])
      append!(mask_out,min_var_match)
      i0=i1
   end
   return sort(mask_out,[:lambda])
end



#this function was translated by ChatGPT from the function of the same name in make_VALD_line_list.py, and then edited to resolve errors
function airVacuumConversion(w; toAir=true)
   ss0 = 10^4 ./ w
   n0 = 1.0 .+ 0.0000834254 .+ 0.02406147 ./ (130.0 .- ss0.^2) .+ 0.00015998 ./ (38.9 .- ss0.^2)
   if toAir
       return w ./ n0
   else
       return w .* n0
   end
end
