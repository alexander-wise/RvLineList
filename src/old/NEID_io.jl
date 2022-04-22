# io functions for neid solar spectra


#make_manifest copied from EchelleInstruments/src/files.jl and modified to walk through neid file structure
function get_filenames_NEID(data_path::String, r::Regex; verbose::Bool = true)
   if !isdir(data_path)
      @error "Can't access data directory ", data_path, "."
   end
   if verbose  println("# Creating manifest of files to process.")    end

   #dir_filelist = readdir(data_path,join=true)
   dir_filelist = Vector{String}()
   for (root, dirs, files) in walkdir(data_path, follow_symlinks=true)
      for file in files
          append!(dir_filelist,[joinpath(root, file)]) # path to files
      end
   end

   idx_spectra = map(fn->occursin(r, last(split(fn,'/')) ),dir_filelist)
   #spectra_filelist = dir_filelist[idx_spectra]

   df_files = DataFrame(:Filename => dir_filelist[idx_spectra])

   if size(df_files,1) < 1
      @error("Did not find any files in ", data_path ,".")
   end

   if match(r"(\d{8})[T_](\d{6})", last(split(df_files.Filename[1],'/')) ) != nothing
      results = map(fn->match(r"(\d{8})[T_](\d{6})", last(split(fn,'/')) ), df_files.Filename)
      @assert all(map(r->length(r.captures),results).==2)
      df_files.date_str = map(r->String(r.captures[1]),results)
      df_files.time_str = map(r->String(r.captures[2]),results)
   end

   return df_files
end

"""Create Dataframe containing filenames and key data for all files neid*.fits in directory"""
function make_manifest_NEID(data_path::String ; max_spectra_to_use::Int = 1000 )
    df_filenames = get_filenames_NEID(data_path,r"^neidL[1,2]_\d+[T\.]\d+\.fits$")
    df_files = DataFrame(NEID.read_metadata(df_filenames.Filename[1]))
    keys = propertynames(df_files)
    allowmissing!(df_files, keys[map(k->k∉[:Filename, :bjd, :target, :airmass],keys)] )
    #println("df_files = ", df_files)


    #if length(spectra_filelist) >= 2
    if length(df_filenames.Filename) >= 2
        max_idx = min(max_spectra_to_use,length(df_filenames.Filename))
        map(fn->NEID.add_metadata_from_fits!(df_files,fn),df_filenames.Filename[2:max_idx])
    end
    df_files

end