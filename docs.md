# RVLineList quick start guide

0. Note: RvLineList is currently configured for use on the Roar Collab cluster. To use on another system, various reduced data and metadata files would need to be generated, which is beyond the scope of this guide.

1. Download RvLinelist from the [github page](https://github.com/alexander-wise/RvLineList).

2. Add the package to julia.
```
using Pkg
Pkg.add(path="/path/where/you/downloaded/RvLineList")
```

3. Install python dependencies.
```
Pkg.add("PyCall")
using PyCall
try
   pyimport("pandas") #make sure pandas is installed
catch e
   pyimport_conda("pandas", "pandas") #if the module fails to load, install it to Conda.jl environment
   pyimport("pandas") #make sure pandas is installed
end
```

4. Edit desired params in param.jl. The following params are also available as command line arguments: allowBlends, overlap_cutoff, rejectTelluricSlope, badLineFilter, quant, nbin, output_dir, and long_output. See comments in param.jl for param descriptions.

5. Run RvLineList in the terminal.
```
cd /path/where/you/downloaded/RvLineList/

julia --project=. examples/NEID_test_script.jl --allowBlends=0 --overlap_cutoff=1e-5 --rejectTelluricSlope=10000 --badLineFilter="none", --nbin=1 --output_dir=outputs/test1
```

# RvLineList Output Files

RvLineList creates four subdirectories in your **output_dir**. These are **clean_masks**, **linefinder**, **VALD_masks**, and **mask_bins**. The final results are stored in **output_dir/clean_masks**, and intermediate products are stored in the other directories.

## clean_masks

The final output mask of RvLineList, in either **default** or **long_output** format. Long_output contains additional columns with line parameters from VALD or computed by RvLineList to determine mask membership.

The columns are:

1. 1
2. 2
3. etc

## linefinder

### neid_linefinder_lines.csv

This file contains the results of automatically searching for and fitting lines in a high-SNR template spectrum constructed from the 100 "best" days NEID data in 2021. 

### neid_linefinder_line_fits.csv

To get the line wavelengths for the mask, we fit our model to the lines in neid_linefinder_lines.csv for each of the 100 "best" days of 2021 NEID solar data individually. Therefore, there are 100 x (the number of lines in neid_linefinder_lines.csv) rows in the neid_linefinder_line_fits.csv file. To match this file with the mask, the line_id columns can be matched up. E.g. if the first line in the clean_mask has line_id = 5, then its 100 model fits are the rows in neid_linefinder_line_fits.csv that have line_id = 5.

### neid_linefinder_line_RVs.csv

This file contains line-by-line RVs measured on 100 NEID daily average spectra from 2021. The are not currently used elsewhere in RvLineList.

## VALD_masks

The default behavior of RvLineList is to leave this empty, but if modified, the VALD-branch outputs can be output here.

## mask_bines

This folder used to be used to track which mask entry ended in which bin when nbins > 1, but it is not longer output, so the folder will be empty unless the user puts something there.

