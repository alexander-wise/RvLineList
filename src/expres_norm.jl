
# spectra normalization functions, copied from an old version of RvSpectML.EchelleInstruments to here due to the most up-to-date version giving errors
# TODO: push changes to EchelleInstruments to fix these errors, and update RvLineList to use the official version of EchelleInstruments

""" Normalize spectrum based on continuum model from FITS file. """
function continuum_normalize_spectrum!(spectrum::ST) where { ST<:AbstractSpectra }
    @assert haskey(spectrum.metadata,:continuum)
    if spectrum.metadata[:normalization] == :raw
        spectrum.flux ./= spectrum.metadata[:continuum] .* spectrum.metadata[:blaze]
        spectrum.var ./= (spectrum.metadata[:continuum].* spectrum.metadata[:blaze] ) .^2
        spectrum.metadata[:normalization] = :continuum
    elseif spectrum.metadata[:normalization] == :blaze
        spectrum.flux ./= spectrum.metadata[:continuum]
        spectrum.var ./= (spectrum.metadata[:continuum] ) .^2
        spectrum.metadata[:normalization] = :continuum
    elseif spectrum.metadata[:normalization] == :continuum
        # do nothing
    else
        @error "Normalizing from " * string(spectrum.metadata.normalization) * " to continuum isn't implemented yet"
    end
    return spectrum
end

""" Normalize each spectrum based on continuum model from FITS files. """
function continuum_normalize_spectra!(spectra::AS) where { ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    for spectrum in spectra
        continuum_normalize_spectrum!(spectrum)
    end
    return spectra
end