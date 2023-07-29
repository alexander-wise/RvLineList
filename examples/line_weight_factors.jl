#line weight factors.jl

using CSV, DataFrames, QuadGK

C_m_s = 2.99792458e8 #speed of light in m/s
line_width_sigma = C_m_s * 1e-5 #guassian sigma factor for line width in m/s


function line_gauss(λ::Float64, λ0::Float64)
   return exp( (-C_m_s^2 / (2*line_width_sigma^2)) * ((λ - λ0) / λ0)^2 )
end

function line_gauss_deriv(λ::Float64, λ0::Float64)
   return exp( (-C_m_s^2 / (2*line_width_sigma^2)) * ((λ - λ0) / λ0)^2 ) * (-C_m_s^2 / (line_width_sigma^2)) * (λ - λ0) / λ0^2
end

function Flux(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64})
   result = 1.0
   for line in 1:length(λ0)
      result *= (1 - D0[line] * line_gauss(λ, λ0[line]))
   end
   return result
end

function dFdλ(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64})
   result = 0.0
   for line in 1:length(λ0)
      result += Flux(λ, λ0, D0) / (1 - D0[line] * line_gauss(λ, λ0[line])) * (-D0[line] * line_gauss_deriv(λ,λ0[line]))
   end
   return result
end

function dFdlogλ(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64})
   return λ * dFdλ(λ, λ0, D0)
end

function dFdT(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64}, ϵ::Vector{Float64}, T_form::Vector{Union{Missing, Float64}})
   result = 0.0
   for line in 1:length(λ0)
      result += Flux(λ, λ0, D0) / (1 - D0[line] * line_gauss(λ, λ0[line])) * (-line_gauss(λ, λ0[line]) * dDdT(λ0[line], D0[line], ϵ[line], T_form[line]))
   end
   return result
end

#same as above function, but simple means it uses a constant for dDdT (and the central line has dDdT=0)
function dFdT_simple(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64}, ϵ::Vector{Float64}, T_form::Vector{Union{Missing, Float64}}, index_center=0)
   result = 0.0
   for line in 1:length(λ0)
      if line != index_center
         result += Flux(λ, λ0, D0) / (1 - D0[line] * line_gauss(λ, λ0[line])) * (-line_gauss(λ, λ0[line]) * dDdT_simple(λ0[line], D0[line], ϵ[line], T_form[line]))
      end
   end
   return result
end

h = 6.62607015e-34 #planck's constant in SI units
k = 1.380649e-23 #boltzmann's constant in SI units
eV = 1.60218e-19 #conversion from electron volts to joules
Ω = 0.0015 #a constant from equation 4 of Wise et al. (2022): Spectral line Depth Variability in Radial Velocity Spectra
Tsun = 5772.0
T_form_old = 0.85 * Tsun #Using 15% less than effective temperature for temperature of line-forming region, following Wise et al. (2022)

#taken from equation 9 of Wise et al. (2022): Spectral line Depth Variability in Radial Velocity Spectra
function dDdT( λ0::Float64, D0::Float64, ϵ::Float64, T_form::Float64)
   return D0*(2.5/T_form + eV*(ϵ+0.75)/(k*T_form*T_form) - Ω) - (h*C_m_s/(λ0*k*T_form*T_form)) / (1 - exp(-h*C_m_s/(λ0*k*T_form)))
end

#for this function, we simply assume all lines have the same constant fractional depth variation
function dDdT_simple( λ0::Float64, D0::Float64, ϵ::Float64, T_form::Float64)
   return D0*0.01
end


#note, we are missing a factor of T here, so we are actualy calculating delta RV / T
function numerator0(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64}, ϵ::Vector{Float64}, T_form::Vector{Union{Missing, Float64}})
   return dFdT(λ, λ0, D0, ϵ, T_form) * dFdlogλ(λ, λ0, D0)
end

#same as above function, but simple means it uses a constant for dDdT (and the central line has dDdT=0)
function numerator0_simple(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64}, ϵ::Vector{Float64}, T_form::Vector{Union{Missing, Float64}}, index_center=0)
   return dFdT_simple(λ, λ0, D0, ϵ, T_form, index_center) * dFdlogλ(λ, λ0, D0)
end


function denominator0(λ::Float64, λ0::Vector{Float64}, D0::Vector{Float64})
   return dFdlogλ(λ, λ0, D0) * dFdlogλ(λ, λ0, D0)
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

λ_NEID_low = 3884.0 - 5.0
λ_NEID_high = 9892.0 + 5.0

#load VALD data
VALD_data = CSV.read(joinpath("inputs", "VALD_extract_stellar", "VALD-Solar-0.001-merged.txt"), delim=',', skipto=4, footerskip=108, header = 3, DataFrame)
select!(VALD_data, [2,3,10])
rename!(VALD_data, ["wavelength", "EE", "depth"])
VALD_indices_NEID = findall(λ_NEID_low .< VALD_data[!,:wavelength] .< λ_NEID_high)
filter!(row -> (λ_NEID_low < row.wavelength < λ_NEID_high), VALD_data)
VALD_data[:,:wavelength] .= airVacuumConversion.(VALD_data[:,:wavelength], toAir=false)

function getRelevantLines(index_center::Int64, VALD_data::DataFrame)
   @assert 1 <= index_center <= size(VALD_data,1)
   λ = VALD_data[index_center, :wavelength]
   Δλ = λ * 6 * line_width_sigma / C_m_s #consider lines within 6 sigma of center line, so the tails at 3 sigma of blended lines are considered
   λ_low = λ - Δλ
   λ_high = λ + Δλ
   return findall(λ_low .< VALD_data[:,:wavelength] .< λ_high) #filter(row -> (λ_low < row.wavelength < λ_high), VALD_data)
end


#formation temperature data provided by private email communcation with Khaled Al Moulla @ University of Geneva
T_form_data = CSV.read(joinpath("inputs", "Formation_Temperature", "T1o2_spec.csv"), DataFrame)

function interpolate_missing_values(T_form_data::DataFrame, T_form_indices::Vector{Int64})
   T_form = T_form_data[T_form_indices,"T1o2"]
   valid_indices = .!ismissing.(T_form)
   missing_indices = ismissing.(T_form)
   #in the case of a missing value, find the closest non-missing value and use that
   for i in findall(missing_indices)
      T_form[i] = T_form[findmin(abs.(T_form_indices[valid_indices] .- T_form_indices[i]))[2]]
   end
   return T_form
end

#calculate the blend factor for a given VALD index
function ΔRV(index_center::Int64, VALD_data::DataFrame, T_form_data::DataFrame)
   @assert 1 <= index_center <= size(VALD_data,1)
   VALD_indices = getRelevantLines(index_center, VALD_data)
   nn = length(VALD_indices)
   λ0 = VALD_data[VALD_indices, :wavelength]
   D0 = VALD_data[VALD_indices, :depth]
   nd = sum(D0)
   ϵ = VALD_data[VALD_indices, :EE]
   T_form_indices = zeros(Int,length(λ0))
   T_form_indices[1] = findmin(abs.(T_form_data[!,"wave"] .- λ0[1]))[2]
   #use the fact that T_form_data wavelengths are in increments of 0.01 to speed up finding other indices
   λ_differences = λ0 .- T_form_data[!,"wave"][T_form_indices[1]]
   for i in 2:length(λ0)
      T_form_indices[i] = T_form_indices[1] + Int.(round.(λ_differences[i]*100))
   end
   #these next 3 commented-out lines are the equivalent but slow way to calculate T_form_indices
   #for i in eachindex(λ0)
   #   T_form_indices[i] = findmin(abs.(T_form_data[!,"wave"] .- λ0[i]))[2]
   #end
   T_form = T_form_data[T_form_indices,"T1o2"]
   λ_center = VALD_data[index_center, :wavelength]
   Δλ = λ_center * 3 * line_width_sigma / C_m_s #integrate +/- 3 sigma of each line
   if any(ismissing.(T_form))
      return missing, nn, nd
   end
   integral_numerator, err_numerator = quadgk(λ -> numerator0(λ, λ0, D0, ϵ, T_form), λ_center - Δλ, λ_center + Δλ)
   if err_numerator > 0.01 * abs(integral_numerator)
      printn("Warning: numerical integral error is greater than 1% for index ", index_center)
   end
   #println("numerator integral result: ", integral_numerator, ", ", err_numerator)
   integral_denominator, err_denominator = quadgk(λ -> denominator0(λ, λ0, D0), λ_center - Δλ, λ_center + Δλ)
   if err_denominator > 0.01 * abs(integral_denominator)
      printn("Warning: numerical integral error is greater than 1% for index ", index_center)
   end
   #println("denominator integral result: ", integral_denominator, ", ", err_denominator)
   return C_m_s * integral_numerator / integral_denominator, nn, nd
end

#same as above function, but simple means it uses a constant for dDdT (and the central line has dDdT=0)
function ΔRV_simple(index_center::Int64, VALD_data::DataFrame, T_form_data::DataFrame)
   @assert 1 <= index_center <= size(VALD_data,1)
   VALD_indices = getRelevantLines(index_center, VALD_data)
   new_index_center = findall(VALD_indices .== index_center)[1]
   nn = length(VALD_indices)
   λ0 = VALD_data[VALD_indices, :wavelength]
   D0 = VALD_data[VALD_indices, :depth]
   nd = sum(D0)
   ϵ = VALD_data[VALD_indices, :EE]
   T_form_indices = zeros(Int,length(λ0))
   T_form_indices[1] = findmin(abs.(T_form_data[!,"wave"] .- λ0[1]))[2]
   #use the fact that T_form_data wavelengths are in increments of 0.01 to speed up finding other indices
   λ_differences = λ0 .- T_form_data[!,"wave"][T_form_indices[1]]
   for i in 2:length(λ0)
      T_form_indices[i] = T_form_indices[1] + Int.(round.(λ_differences[i]*100))
   end
   #these next 3 commented-out lines are the equivalent but slow way to calculate T_form_indices
   #for i in eachindex(λ0)
   #   T_form_indices[i] = findmin(abs.(T_form_data[!,"wave"] .- λ0[i]))[2]
   #end
   T_form = T_form_data[T_form_indices,"T1o2"]
   λ_center = VALD_data[index_center, :wavelength]
   Δλ = λ_center * 3 * line_width_sigma / C_m_s #integrate +/- 3 sigma of each line
   if any(ismissing.(T_form))
      return missing, nn, nd
   end
   integral_numerator, err_numerator = quadgk(λ -> numerator0_simple(λ, λ0, D0, ϵ, T_form, new_index_center), λ_center - Δλ, λ_center + Δλ)
   if err_numerator > 0.01 * abs(integral_numerator)
      printn("Warning: numerical integral error is greater than 1% for index ", index_center)
   end
   #println("numerator integral result: ", integral_numerator, ", ", err_numerator)
   integral_denominator, err_denominator = quadgk(λ -> denominator0(λ, λ0, D0), λ_center - Δλ, λ_center + Δλ)
   if err_denominator > 0.01 * abs(integral_denominator)
      printn("Warning: numerical integral error is greater than 1% for index ", index_center)
   end
   #println("denominator integral result: ", integral_denominator, ", ", err_denominator)
   return C_m_s * integral_numerator / integral_denominator, nn, nd
end


#calculate blend factors for every line in a VALD line list
function make_VALD_data_blend_factors_file(VALD_data::DataFrame, T_form_data::DataFrame)
   n_VALD = size(VALD_data)[1]
   blend_factors = zeros(Union{Missing, Float64},n_VALD)
   for i in 1:n_VALD
      if i/100 == i ÷ 100
         print(i)
      end
      blend_factors[i] = ΔRV(i, VALD_data, T_form_data)[1]
   end
   df = DataFrame(VALD_index = VALD_indices_NEID, blend_RV_factor = blend_factors)
   CSV.write(joinpath("inputs", "blend_factors","blend_factors_0.001.csv"), df)
   return df
end

#same as above function, but simple means it uses a constant for dDdT (and the central line has dDdT=0)
function make_VALD_data_blend_factors_file_simple(VALD_data::DataFrame, T_form_data::DataFrame)
   n_VALD = size(VALD_data)[1]
   blend_factors = zeros(Union{Missing, Float64},n_VALD)
   for i in 1:n_VALD
      if i/100 == i ÷ 100
         print(i)
      end
      blend_factors[i] = ΔRV_simple(i, VALD_data, T_form_data)[1]
   end
   df = DataFrame(VALD_index = VALD_indices_NEID, blend_RV_factor = blend_factors)
   CSV.write(joinpath("inputs", "blend_factors","blend_factors_simple_0.001.csv"), df)
   return df
end



make_VALD_data_blend_factors_file(VALD_data, T_form_data)

#make_VALD_data_blend_factors_file_simple(VALD_data, T_form_data)



###ABOVE THIS COMMENT: MAKE BLEND FACTORS FILE FOR VALD LINE LIST
###BELOW THIS COMMENT: EXPLORE BLEND FACTORS RESULTS


"""
#clean_line_list = CSV.read(joinpath("outputs", "clean_masks", "RvLineList_allowBlends=[0, 1, 2, 3, 4, 5, 6, 7, 8]_overlapCutoff=2.0e-5_depthCutoff=0.05_rejectTelluricSlope=2000.0_badLineFilter=none_quant=100_nbin=1_DP=true_binParam=depth_n=1_long_output=true_VACUUM.csv"), DataFrame)
#clean_line_list = CSV.read(joinpath("outputs", "clean_masks", "RvLineList_allowBlends=[0, 1, 2, 3, 4, 5, 6, 7, 8]_overlapCutoff=2.0e-5_depthCutoff=0.05_rejectTelluricSlope=2000.0_badLineFilter=none_quant=100_nbin=1_DP=true_binParam=depth_n=1_long_output=true_VACUUM_nov_dec_2021.csv"), DataFrame)
clean_line_list = CSV.read(joinpath("outputs", "clean_masks", "RvLineList_allowBlends=[0, 1, 2, 3, 4, 5, 6, 7, 8]_overlapCutoff=3.0e-5_depthCutoff=0.01_rejectTelluricSlope=2000.0_badLineFilter=none_quant=100_nbin=1_DP=true_binParam=depth_n=1_long_output=true_VACUUM.csv"), DataFrame)

select!(clean_line_list,[:lambda, :depth, :std_λc, :blend_number, :species])

function getVALDindex(λ::Float64, VALD_data::DataFrame)
   Δλ = λ * 2.0 * line_width_sigma / C_m_s
   λ_low = λ - Δλ
   λ_high = λ + Δλ
   VALD_candidates = findall(λ_low .< VALD_data[:,:wavelength] .< λ_high)
   if length(VALD_candidates)==0
      println("error: ",λ)
      return 1
   end
   return VALD_candidates[findmax(VALD_data[VALD_candidates,:depth])[2]]
end


ΔRV_per_K = zeros(size(clean_line_list,1))
NN = zeros(size(clean_line_list,1))
ND = zeros(size(clean_line_list,1))

@time(
for (i,line) in enumerate(eachrow(clean_line_list))
   VALD_index = getVALDindex(line[:lambda], VALD_data)
   ΔRV_per_K[i], NN[i], ND[i] = ΔRV(VALD_index,VALD_data, T_form_data)
end)



using Plots


deep_lines = findall(clean_line_list[:,:depth] .> 0.3)
#blends = findall(clean_line_list[:,:blend_number] .>= 2)
blends = findall(NN .>= 25)
deep_blends = blends[findall(in(deep_lines),blends)]

iron = findall(isequal.(clean_line_list[:,:species],"'Fe 1'"))
iron_blends = iron[findall(in(blends),iron)]
deep_iron = iron[findall(in(deep_lines),iron)]
deep_iron_blends = iron[findall(in(deep_blends),iron)]

ii = deep_iron_blends

x = clean_line_list[:, :std_λc][ii] ./ clean_line_list[:, :lambda][ii] .* C_m_s
#x = clean_line_list[:, :depth][ii]
y = abs.(ΔRV_per_K)[ii]
#z = clean_line_list[:, :blend_number][ii]
z = ND[ii]


plot(x, y, zcolor = z, seriestype = :scatter, markersize=5, color = palette(:heat, length(unique(z))) )

xlims!(3,40)
ylims!(0,3.49)
#ylims!(0,12)
xlabel!("line std(λ) (m/s)")
#xlabel!("line depth")
ylabel!("ΔRV per K (m/s/K)")

"""