#old expres code used in RvLineList scripts


#read spectra code for EXPRES data
function EXPRES_read_spectra(data_path::String; manifest_path = joinpath(pipeline_output_path_afw5465,"manifest.csv"), verbose::Bool=false)
   if verbose println("# Finding what data files are avaliable.")  end
   if isfile(manifest_path)
      if verbose println("# Reading in manifest from manifest.csv") end
      df_files  = CSV.read(manifest_path, DataFrame)
      @assert size(df_files,1) >= 1
      @assert hasproperty(df_files,:Filename)
      @assert hasproperty(df_files,:bjd)
   else
      if verbose println("# Generating manifest file manifest.csv") end
      df_files = EXPRES.make_manifest(data_path)
      CSV.write(manifest_path, df_files)
   end
   df_files_use = df_files |> @take(max_spectra_to_use) |> DataFrame
   if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
   @time all_spectra = map(EXPRES.read_data,eachrow(df_files_use))
end