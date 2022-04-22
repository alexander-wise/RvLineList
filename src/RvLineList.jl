"""
Author: Alex Wise and collaborators (see header of each file)
Created: October 2021
Contact: https://github.com/RvSpectML/RvLineList.jl
"""

module RvLineList

using JLD2, FileIO #objects used: save()

using CSV #objects used: CSV

using DataFrames #objects used: DataFrame

using Query #objects used: @filter

using RvSpectMLBase #objects used: read_data_paths

using EchelleInstruments #objects used: NEID, EXPRES, HARPSN
export NEID, EXPRES, HARPS, HARPSN

import EchelleCCFs:λ_air_to_vac #objects used: EchelleCCFs.λ_air_to_vac
#λ_air_to_vac = EchelleCCFs.λ_air_to_vac
#export λ_air_to_vac

using RvSpectML #objects used: extract_orders

using Dates #objects used: Dates.date
#filter datetime objects to get datetime list to use
#solar noon should use shift by longitude/15 to get from UTC to the thing we want

import Printf #objects used: @sprintf

using Statistics #objects used: quantile

using PyCall #makes it possible to call python functions from make_VALD_line_list.py in a julia environment
import Pandas.DataFrame as pd_df #used to convert julia DataFrame to python pandas.DataFrame

#mask manipulation
include("masks.jl")
export read_mask_air, read_mask_vacuum, binMask
export get_param_range, get_lambda_range, get_depth_range

#normaliation of EXPRES spectra - this file will not be needed once these functions are replaced with working versions from EchelleInstruments
include("expres_norm.jl")

#empirical mask generation
include("empirical_line_lists.jl")
export generateEmpiricalMask


# Packages we'll use in many places, but try to keep these minimal.
#using LinearAlgebra
#using NaNMath

end



