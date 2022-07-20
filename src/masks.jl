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
    @assert occursin("_VACUUM",fn) | occursin("_AIR",fn)
    local df = CSV.read(fn,DataFrame,threaded=false,header=["lambda","depth", "weight"],skipto=2)
    @assert issorted(df[!,:lambda])
    if occursin("_AIR",fn)
      df[!,:lambda] .= λ_air_to_vac.(df[!,:lambda]) #convert air to vacuum wavelength
    end
    return df
end

function get_airVacString(inst::Module)
   local airVacString
   if inst in (NEID, EXPRES)
      airVacString="_VACUUM"
   elseif inst in (HARPS, HARPSN)
      airVacString="_AIR"
   end
   return airVacString
end

function binMask(mask::DataFrame, nbin::Int; Inst::Module = EXPRES, orders_to_use::UnitRange{Int64} = 43:72, binParam::Symbol = :weight, depthPercentile::Bool = true)
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