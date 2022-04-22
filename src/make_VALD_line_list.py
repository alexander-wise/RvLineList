# Generate spectral line lists based on the VALD database.

import os
import numpy as np
import pandas as pd
import time


#VALD "extract stellar" parameters used for VALD-Solar.txt: Starting Wavelength: 3000, Ending Wavelength: 10000, Detection Threshold: 0.05, Microturbulence: 1.0, Teff: 5778, log g: 4.44, Chemical composition: Fe: -4.54, Extraction format: short format, all other params default

C_m_s = 2.99792458e8 #speed of light in m/s
C_km_s = C_m_s / 1000.

#set data directories

VALD_dir = 'inputs/VALD_extract_stellar' # VALD parameters
ESPRESSO_dir = 'inputs/ESPRESSO_masks'
empirical_dir = 'inputs/empirical_fits'
outdir_masks = 'outputs/masks' #stores .mas mask files that mirror the espresso mask format
outdir_bins = 'outputs/mask_bins' #stores useful data files when masks are split into bins (e.g. 10 bins of increasing wavelength)



#load AlphaCenB VALD data:
#load line species, wavelengths, excitation energies, oscillator strengths, and depths from VALD
lineprops_K1 = pd.read_table(os.path.join(VALD_dir, 'VALD-AlphaCenB.txt'),delimiter=',',skiprows=3, skipfooter=98,usecols=(0,1,2,4,9), names=("species", "lambda", "excitation_energy", "oscillator_strength", "depth"), engine="python")

#load solar VALD data:
lineprops_G2 = pd.read_table(os.path.join(VALD_dir, 'VALD-Solar.txt'),delimiter=',',skiprows=3, skipfooter=106,usecols=(0,1,2,4,9), names=("species", "lambda", "excitation_energy", "oscillator_strength", "depth"), engine="python")

#load VALD data for HD101501:
lineprops_G8 = pd.read_table(os.path.join(VALD_dir, 'VALD-HD101501.txt'),delimiter=',',skiprows=3, skipfooter=109,usecols=(0,1,2,4,9), names=("species", "lambda", "excitation_energy", "oscillator_strength", "depth"), engine="python")

#################### BASIC FUNCTIONS AND ARRAY MANIPULATIONS ####################

#convert the wavelength domain to air wavelengths using the equations from http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
def airVacuumConversion(w, toAir=True):
   ss0 = 10.**4 / w
   n0 = 1.0 + 0.0000834254 + 0.02406147 / (130. - ss0**2) + 0.00015998 / (38.9 - ss0**2)
   if toAir:
      return w / n0
   else:
      return w * n0


#function to transform espresso masks from fits files to .mas formatted masks
def genMaskFileFromMaskFitsFile(mask_label):
   dir_in = '/Users/aww/Desktop/espdr-kit-2.2.1/espdr-calib-2.2.1/cal/'
   dir_out = '/Users/aww/Desktop/activity/masks/'
   fits_data = fits.open(dir_in+'ESPRESSO_'+mask_label+'.fits')[1].data
   new_mask = np.zeros((len(fits_data),2))
   for i in range(len(fits_data)):
      new_mask[i] = fits_data[i]
   savetxt(dir_out+'ESPRESSO_'+mask_label+'_mask.mas', new_mask)


#check if two mask wavelengths are within 50 m/s (approx 0.001 angstroms)
#VALD numerical precision is 0.0001 angstroms, or about 5 m/s. This number (50 m/s) was chosen to be a factor of a few greater than VALD precision at all wavelengths we consider.
def wave_equal(w1, w2, threshold=50.0):
   return (abs(w1-w2)/((w1+w2)/2)) < (threshold / C_m_s)


#take the union of two masks (pandas 2d dataframes) to generate one super_mask
def mask_union(mask1, mask2, default_data="first", keep_all_mask1=False, add_label_column=False, threshold=50.0):
   """
   Take the union of two masks to generate one supermask. By default, two mask entries are treated as equal if they have wavelengths (lambdas) within threshold (50 m/s by default) of each other, or if they are both equal to the same third mask entry.

   Parameters:
      mask1 (2-D pandas.Datafrmae): mask #1 including columns for lambda (required) and depth (required if default_data="max_depth")

      mask2 (2-D pandas.Datafrmae): mask #2 including columns for lambda (required) and depth (required if default_data="max_depth")

      default_data (str): keyword for which mask data to preserve in the case of discarding duplicate lines.
         Posible values:
            "first" defaults to mask1's data for matches across masks. Will default to max_depth for matches within masks.
            "max_depth" picks the line with max depth in the case of a match. If multiple lines with the same max depth, defaults to the lower lambda line within the match.

      keep_all_mask1 (bool): whether or not to preserve all mask1 entries in the output. This would mean only mask2 entries have the possibility to be thrown away because they match with another mask entry.

      add_label_column (bool): whether or not to add a column, "mask_df_name", to the output supermask containing labels "mask1" and "mask2" tracking which input mask each output mask entry originated from.

      threshold (float): velocity (in m/s) separation between adjacent mask entries for them to be considered equal.

   Returns:
      super_mask (2-D pandas.Dataframe): union of mask1 and mask2, sorted by lambda
   """
   mask1 = pd.DataFrame(mask1)
   mask2 = pd.DataFrame(mask2)
   assert (default_data=="first" or default_data=="max_depth"), "ERROR: invalid default_data selection."
   s1 = len(mask1.shape)
   s2 = len(mask2.shape)
   assert s1 == s2, "ERROR: mask inputs have different shapes."
   assert s1 == 2, "ERROR: mask inputs must be 2-D pandas.Dataframes."
   assert ("lambda" in mask1.columns) and ("lambda" in mask2.columns), "ERROR: lambda column not found."
   if ("depth" in mask1.columns) and ("depth" in mask2.columns):
      hasDepths=True
   else:
      assert default_data != "max_depth", "ERROR: max_depth selected but depth keyword missing from mask1.columns or mask2.columns."
      #assert union_within_masks == False, "ERROR: union_within_masks=True but depth keyword missing from mask1.columns or mask2.columns. Unable to decide which entry to pick for combined matches within mask."
      hasDepths=False
   mask1sorted = mask1.sort_values(by="lambda", ignore_index=True)
   mask2sorted = mask2.sort_values(by="lambda", ignore_index=True)
   mask1sorted["mask_df_name"] = "mask1"
   mask2sorted["mask_df_name"] = "mask2"
   super_mask = []
   mask12combined = pd.concat([mask1sorted,mask2sorted], ignore_index=True)
   mask12sorted = mask12combined.sort_values(by="lambda", ignore_index=True)
   m12 = mask12sorted["lambda"]
   i=0
   while i < (len(m12)-1): #loop through combined mask index i
      j=1
      match_i = [i]
      while(wave_equal(m12[i+j-1],m12[i+j],threshold=threshold)):
         match_i.append(i+j)
         j+=1
      if j > 1: #if there are any matches to mask index i
         matches = mask12sorted.iloc[match_i]
         mask1matches = matches["mask_df_name"] == "mask1"
         if any(mask1matches):
            if default_data == "first":
               matches = matches.loc[matches.index[mask1matches]]
            max_depth_match = matches.loc[matches["depth"].idxmax()]
            if keep_all_mask1:
               for k,m in matches.loc[matches.index[matches["mask_df_name"] == "mask1"]].iterrows():
                  super_mask.append(m)
               if max_depth_match["mask_df_name"] != 'mask1':
                  super_mask.append(max_depth_match)
            else:
               super_mask.append(max_depth_match)
         else:
            max_depth_match = matches.loc[matches["depth"].idxmax()]
            super_mask.append(max_depth_match)
      else: #there are no matches to mask index i
         super_mask.append(mask12sorted.iloc[i])				
      i += j
   super_mask_out = pd.DataFrame(data=super_mask)#,columns=mask12sorted.columns)
   if not add_label_column:
      super_mask_out = super_mask_out.drop(columns="mask_df_name")
   return super_mask_out.sort_values(by="lambda", ignore_index=True)	


#take the intersection of two masks (pandas 2d dataframes) to generate one sub_mask
def mask_intersection(mask1, mask2, default_data="first", add_label_column=False, combine_mask_data=True, threshold=50.0):
   """
   Take the intersection of two masks to generate one submask. Two mask entries are treated as equal if they have wavelengths (lambdas) within threshold (50 m/s by default) of each other, or if they are both equal to the same third mask entry.

   Parameters:
      mask1 (2-D pandas.Datafrmae): mask #1 including columns for lambda (required) and depth (required if default_data="max_depth")

      mask2 (2-D pandas.Datafrmae): mask #2 including columns for lambda (required) and depth (required if default_data="max_depth")

      default_data (str): keyword for which mask data to preserve in the case of matching (intersecting) lines - all but one are discarded.
         Posible values:
            "first" defaults to mask1's data for matches across masks. Will default to max_depth for matches within masks.
            "max_depth" picks the line with max depth in the case of a match. If multiple lines with the same max depth, defaults to the lower lambda line within the match.

      add_label_column (bool): whether or not to add a column, "mask_df_name", to the output supermask containing labels "mask1" and "mask2" tracking which input mask each output mask entry originated from.

      combine_mask_data (bool): whether or not to add extra data in mask2 to mask1 (only works when default_data=="first")

      threshold (float): velocity (in m/s) separation between adjacent mask entries for them to be considered equal.

   Returns:
      sub_mask (2-D pandas.Dataframe): intersection of mask1 and mask2, sorted by lambda
   """
   mask1 = pd.DataFrame(mask1)
   mask2 = pd.DataFrame(mask2)
   assert (default_data=="first" or default_data=="max_depth"), "ERROR: invalid default_data selection."
   s1 = len(mask1.shape)
   s2 = len(mask2.shape)
   assert s1 == s2, "ERROR: mask inputs have different shapes."
   assert s1 == 2, "ERROR: mask inputs must be 2-D pandas.Dataframes."
   assert ("lambda" in mask1.columns) and ("lambda" in mask2.columns), "ERROR: lambda column not found."
   if ("depth" in mask1.columns) and ("depth" in mask2.columns):
      hasDepths=True
   else:
      assert default_data != "max_depth", "ERROR: max_depth selected but depth keyword missing from mask1.columns or mask2.columns."
      #assert union_within_masks == False, "ERROR: union_within_masks=True but depth keyword missing from mask1.columns or mask2.columns. Unable to decide which entry to pick for combined matches within mask."
      hasDepths=False
   mask1sorted = mask1.sort_values(by="lambda", ignore_index=True)
   mask2sorted = mask2.sort_values(by="lambda", ignore_index=True)
   mask1sorted["mask_df_name"] = "mask1"
   mask2sorted["mask_df_name"] = "mask2"
   sub_mask = []
   mask12combined = pd.concat([mask1sorted,mask2sorted], ignore_index=True)
   mask12sorted = mask12combined.sort_values(by="lambda", ignore_index=True)
   m12 = mask12sorted["lambda"]
   i=0
   while i < (len(m12)-1): #loop through combined mask index i
      j=1
      match_i = [i]
      while(wave_equal(m12[i+j-1],m12[i+j],threshold=threshold)):
         match_i.append(i+j)
         j+=1
      if j > 1: #if there are any matches to mask index i
         matches = mask12sorted.iloc[match_i]
         mask1matches = matches["mask_df_name"] == "mask1"
         mask2matches = matches["mask_df_name"] == "mask2"
         if (any(mask1matches) and any(mask2matches)):
            if default_data == "first":
               matches1 = matches.loc[matches.index[mask1matches]]
               matches2 = matches.loc[matches.index[mask2matches]]
               max_depth_match1 = matches1.loc[matches1["depth"].idxmax()]
               max_depth_match2 = matches2.loc[matches2["depth"].idxmax()]
               if combine_mask_data:
                 for k in mask2.keys():
                     if k not in mask1.keys():
                        max_depth_match1[k] = max_depth_match2[k]
               sub_mask.append(max_depth_match1)
            else:
               max_depth_match = matches.loc[matches["depth"].idxmax()]
               sub_mask.append(max_depth_match)
      i += j
   sub_mask_out = pd.DataFrame(data=sub_mask)#,columns=mask12sorted.columns)
   if not add_label_column:
      sub_mask_out = sub_mask_out.drop(columns="mask_df_name")
   return sub_mask_out.sort_values(by="lambda", ignore_index=True)	



#################### MASK GENERATION AND PROJECTION ####################



#blend rejection function
def getBlends(mask0, overlap_cutoff, allowBlends):
   """
   Find the number of mask entries which overlap with each mask entry, and construct a boolean array telling is which lines to keep based on the specified criteria in allowBlends.

   Parameters:
      mask0 (pandas.DataFrame): line list which must contain fields "lambda" and "depth"

      overlap_cutoff (float): distance from a line center at which another line's center would be considered within that line, expressed as a fraction fo the speed of light

      allowBlends (int or list of int): number of allowed blends which each line. Use 0 to allow no blends, 1 for doublets, [1,2] for doublets and triplets, etc.

   """
   assert isinstance(mask0, pd.DataFrame)
   assert isinstance(overlap_cutoff, float)
   assert (isinstance(allowBlends, list) or isinstance(allowBlends, int) or isinstance(allowBlends, float))
   maskCenters = mask0["lambda"]
   nm = len(mask0)
   depths = mask0["depth"]
   widths = np.array([overlap_cutoff if (depths[i]<0.65) else (1.+(depths[i]-0.65)*4.)*overlap_cutoff for i in range(len(mask0))]) #note: this function is based on visual analysis of the plot of measured FWHM vs line depth from VALD made in measureBISandLW()
   #widths = np.array([overlap_cutoff if (depths[i]<0.6) else (1.+(depths[i]-0.6)*5.)*overlap_cutoff for i in range(len(mask0))]) #note: this function was derived from visual analysis of HWHM vs line depth in HD101501, using a "fit" mask, i.e. data-derived. this mask may contain blends, and blends seem to cause the differerence between these formulae
   rights = maskCenters*(1.+widths)
   lefts = maskCenters*(1.-widths)
   nOverlap = np.zeros(nm,dtype=int)
   olisBoolean = np.zeros(nm,dtype=bool)
   for i in range(nm):
      nOverlap[i] = sum(rights[:i]>maskCenters[i])+sum(lefts[i+1:]<maskCenters[i])
   
   if isinstance(allowBlends, float):
      print("Warning: allowBlends is type float but required to be int or list of int. Converting to int automatically.")
      allowBlends = int(allowBlends)

   if isinstance(allowBlends,int):
      olisBoolean = (allowBlends==nOverlap)

   if isinstance(allowBlends,list):
      for i in allowBlends:
         olisBoolean += (i==nOverlap)
   return olisBoolean, nOverlap


def getTelluricIndices(mask, maskWavelengthsAreVacuum, overlap_cutoff, vel_slope_threshold=2000.0, RV_offset = 0.0, RV_range = 1e-4):
   # telluric_RV_threshold is the Earth's barycentric velocity as a fraction of speed of the light - to be added to overlap_cutoff for stars other than the sun, can be set to 0 for the sun
   fname = 'inputs/Telluric_Rejection/telluric_waves_'+str(vel_slope_threshold)+'.npy'
   if os.path.isfile(fname):
      telluric_waves = np.load(fname)*10. #the *10. converts from nm to Angstroms
   else:
      #get telluric wavelengths using telluric spectrum slope threshold criteria
      #vel_slope_threshold=2000.0
      X = np.genfromtxt('inputs/Telluric_Rejection/tapas_000002.ipac', skip_header=26)
      dD = X[:,1][1:] - X[:,1][:-1] #change in depth between adjacent pixels
      dV = (X[:,0][1:] - X[:,0][:-1]) / (X[:,0][:-1]) #change in velocity between adjacent pixels
      dD_dV = dD / dV #change in depth with respect to velocity
      telluric_indices = np.where((abs(dD_dV)>vel_slope_threshold) | (X[:,1][1:] < 0.1))[0] #where slope is significant or telluric line is saturated
      """
      #plot tellurics - commented out as this segement is only used in debugging
      clf()
      plot(X[:,0]*10.,X[:,1], color='brown', label='Telluric Spectrum')
      plot(X[:,0][Y]*10.,X[:,1][Y],'b.', label='Slope > '+str(+slope_threshold))
      xlabel(r'Wavelength ($\AA$)', size=20)
      ylabel('Normalized Flux', size=20)
      title('Telluric Rejection based on Slope Criteria', size=15)
      subplots_adjust(0.15,0.15,0.95,0.9)
      legend(loc=3)
      ylim(0,1.05)

      clf()
      scatter(10*abs(X[:,0][telluric_indices]),abs(X[:,1][telluric_indices]), c=abs(dD_dV[telluric_indices]), cmap="Oranges", s=0.1, vmin=2000, vmax=20000)
      ylabel("Normalized Flux", size=15)
      xlabel("Wavelength (Angstroms)", size=15)
      cbar = colorbar()
      cbar.set_label("d_Depth_d_Velocity")
      savefig("inputs/Telluric_Rejection/slope_colorscale_plot.png")

      dD_dV[telluric_indices[abs(dD_dV[telluric_indices])<vel_slope_threshold]] = 1e6 #make sure deep telluric lines have a big weight
      np.save("inputs/Telluric_Rejection/telluric_weights_"+str(vel_slope_threshold)+".npy", abs(dD_dV[telluric_indices])) #save telluric weights to desktop
      """
      np.save(fname, X[:,0][telluric_indices])
      telluric_waves = np.load(fname)*10. #the *10. converts from nm to Angstroms
   #telluric_waves are read in as air wavelengths, so we make sure they are both air or both vacuum
   if maskWavelengthsAreVacuum:
      telluric_waves = airVacuumConversion(telluric_waves, toAir=False)
   #shift to RV_offset
   telluric_waves = telluric_waves * (1+RV_offset/C_m_s)
   #compare mask line centers with telluric lines
   maskCenters = mask["lambda"]
   too_close_to_telluric = np.zeros(len(mask),dtype=bool)
   for i in range(len(mask)):
      if (np.amin(abs(maskCenters[i]-telluric_waves)) / maskCenters[i]) < (overlap_cutoff+RV_range):
         too_close_to_telluric[i] = True
   return too_close_to_telluric






def hasMaskMatch(mask, maskWavelengthsAreVacuum, line_width=1e-5, mask_name='G2.espresso.mas', allowMultipleMatches=True):
   maskCenters = mask["lambda"]
   if mask_name[-4:] == '.csv':
      print("Reading empirical line list for mask matching...")
      filterCenters = pd.read_csv(os.path.join(empirical_dir,mask_name))["median_λc"]
      #telluric_waves are read in as air wavelengths, so we make sure they are both air or both vacuum
      filterWavelengthsAreVacuum = False if "_AIR" in mask_name else True
      if filterWavelengthsAreVacuum:
            if "_VACUUM" not in mask_name:
               print("Warning: could not find string \"_AIR\" or \"_VACUUM\" in kwarg mask_name. Assuming filter contains vacuum wavelengths.")
   elif mask_name[-4:] == ".mas":
      print("Reading espresso line list for mask matching...")
      maskFilter = np.genfromtxt(os.path.join(ESPRESSO_dir,mask_name))
      filterCenters = maskFilter[:,0]
      filterWavelengthsAreVacuum = False
   else:
      print("Error: file extension not supported.")
      raise ValueError
   #make sure the mask and filter have the same wavelength medium
   if maskWavelengthsAreVacuum != filterWavelengthsAreVacuum:
      if maskWavelengthsAreVacuum:
         filterCenters = airVacuumConversion(filterCenters, toAir=False)
      else:
         filterCenters = airVacuumConversion(filterCenters, toAir=True)
   #find matches
   hasMatch = np.zeros(len(mask),dtype=bool)
   for i in range(len(filterCenters)):
      matches = np.where((filterCenters[i]>maskCenters-maskCenters*line_width)&(filterCenters[i]<maskCenters+maskCenters*line_width))[0]
      nMatches = len(matches)
      if nMatches>0:
         if nMatches>1:
            if allowMultipleMatches:
               hasMatch[matches]=True
         else:
            hasMatch[matches]=True
   return hasMatch



#Use a VALD line list to make an RV line list. Note the output is converted to vacuum wavelengths by default.
def getVALDmasks(nbin=1, binParam = "depth", lineprops=lineprops_G2, iron1Only='all', depthPercentile=True, overlap_cutoff=1e-5, depth_cutoff=0.05, maskWavelengths = 'Reiners', allowBlends=0, rejectTelluricSlope=0.0, badLineFilter='none', outputVacuumWavelengths=True, saveMasks=False):
   #lineprops are lambdas, excitation energies, oscillator strengths, and depths
   """
   nbin=1
   binParam = "depth"
   lineprops=lineprops_G2
   iron1Only='all'
   depthPercentile=True
   overlap_cutoff=1e-5
   depth_cutoff=0.05
   maskWavelengths='Reiners'
   allowBlends=0
   rejectTelluricSlope=0.0
   badLineFilter='ESPRESSOG2'
   """
   #local variable to track whether wavelengths are air or vacuum
   maskWavelengthsAreVacuum = False #initialized to false because VALD lineprops variables are read in as air wavelengths

   #add a column of random numbers to lineprops to test whether the effects we see are real - to test this use binParam = "random"
   randnums = np.random.rand(len(lineprops))
   lineprops["random"] = randnums

   #apply line depth cutoff
   #in older version, strictly_no_overlap=True meant overlap_cutoff = 4e-5 and depth_cutoff=0.05
   keep_indices = np.where(lineprops["depth"] >= depth_cutoff)[0]
   lineprops = lineprops.iloc[keep_indices]

   #start the mask from lineprops
   mask0 = lineprops[["lambda","depth"]].copy(deep=True)

   #trim mask to only include non-overlapping lines - this happens before any other filtering, since other filters could remove part of a blend
   olisBoolean = getBlends(mask0, overlap_cutoff, allowBlends)[0]
   mask0 = mask0[olisBoolean]
   lineprops = lineprops[olisBoolean]

   if 'Reiners' in maskWavelengths:
      coef_factors = {'':1.0, '101501':0.48, '10700':0.47, '26965':0.36, '34411':0.76} #from alex_sandbox.jl / get_line_shapes.jl code to generate plots of line RV vs depth
      baseline_RVs = {'':0.0, '101501':-5.0e3, '10700':-16640.0, '26965':-40320.0, '34411':66500.0}
      coef_factor = coef_factors[maskWavelengths[7:]]
      reinersEQ2 = lambda d: coef_factor*(-504.891 - 43.7963*d - 145.56*d*d + 884.308*d*d*d) #Multiple of all the coefficients in reiners et al. 2016 for HD101501
      dVs = reinersEQ2(mask0["depth"]) + baseline_RVs[maskWavelengths[7:]]
      mask0["lambda"] = mask0["lambda"] * (1.0+dVs/C_m_s)
      lineprops["lambda"] = lineprops["lambda"] * (1.0+dVs/C_m_s)

   #only keep Fe 1 lines or all but Fe 1 lines:
   isFe1 = lineprops["species"]=="'Fe 1'"
   if iron1Only=='Fe1':
      mask0 = mask0[isFe1]
      lineprops = lineprops[isFe1]
   elif iron1Only=='nonFe1':
      mask0 = mask0[~isFe1]
      lineprops = lineprops[~isFe1]
   elif iron1Only=='all':
      pass
   else:
      print('ERROR: invalid iron1Only code!')

   #reject mask entries too close to a telluric line
   if rejectTelluricSlope != 0.:
      too_close_to_telluric = getTelluricIndices(mask0, maskWavelengthsAreVacuum, overlap_cutoff, vel_slope_threshold=rejectTelluricSlope)
      mask0 = mask0[~too_close_to_telluric]
      lineprops = lineprops[~too_close_to_telluric]

   #badLineFilter contains a few different options that attempt to remove bad lines from the mask
   if 'ESPRESSO' in badLineFilter: #only keep mask entries with a corresponding ESPRESSO mask entry
      mask_label = badLineFilter[8:10]
      hem = hasMaskMatch(mask0, maskWavelengthsAreVacuum, mask_name=mask_label+".espresso.mas", allowMultipleMatches=False if ('strict' in badLineFilter) else True)
      mask0 = mask0[hem]
      lineprops = lineprops[hem]

   if badLineFilter[:2] == "HD":
      hasMatch = hasMaskMatch(mask0, maskWavelengthsAreVacuum, mask_name=badLineFilter[2:]+"_good_lines_fit_quant=0.90.csv", allowMultipleMatches=True)
      mask0 = mask0[hasMatch]
      lineprops = lineprops[hasMatch]

   if badLineFilter[:2] == "20":
      hasMatch = hasMaskMatch(mask0, maskWavelengthsAreVacuum, mask_name=badLineFilter[:8]+"_good_lines_fit_quant=1.00.csv", allowMultipleMatches=True)
      mask0 = mask0[hasMatch]
      lineprops = lineprops[hasMatch]

   if badLineFilter=='ExcitationEnergy': #only keep mask entries with EE > 0.2 eV and EE < 8 eV. These thresholds are chosen to exclude the balmer series (EE>8 eV), and the line most sensitive to stellar activity (these tend to have EE < 0.2).
      okayExcitationEnergy = np.where((lineprops["excitation_energy"]>0.2)&(lineprops["excitation_energy"]<8.))[0]
      mask0 = mask0[okayExcitationEnergy]
      lineprops = lineprops[okayExcitationEnergy]

   #convert air to vacuum if necessary.
   if outputVacuumWavelengths != maskWavelengthsAreVacuum:
      if maskWavelengthsAreVacuum:
            mask0["lambda"] = airVacuumConversion(mask0["lambda"], toAir=True)
            lineprops["lambda"] = airVacuumConversion(lineprops["lambda"], toAir=True)
      else:
            mask0["lambda"] = airVacuumConversion(mask0["lambda"], toAir=False)
            lineprops["lambda"] = airVacuumConversion(lineprops["lambda"], toAir=False)

   #get line indices for each bin
   bins = {}
   if depthPercentile: #this should be updated to include Arvind's normalization functions
      pSorted = np.argsort(lineprops[binParam])
      dCumSum = np.cumsum(lineprops["depth"].iloc[pSorted])
      for i in range(nbin):
         dLower = dCumSum.iloc[-1]/nbin * i
         dUpper = dCumSum.iloc[-1]/nbin * (i+1)
         bins[i] = pSorted.iloc[np.where((dCumSum >= dLower) & (dCumSum <= dUpper))[0]]
   else:
      pmin = np.amin(lineprops[binParam])
      pmax = np.amax(lineprops[binParam])
      prange = pmax-pmin
      for i in range(nbin):
         #plower = pmin + (prange / nbin * i)
         plower = np.percentile(lineprops[binParam], 100./nbin * i)
         #pupper = pmin + (prange / nbin * (i+1))
         pupper = np.percentile(lineprops[binParam], 100./nbin * (i+1))
         bins[i] = np.where((lineprops[binParam] >= plower) & (lineprops[binParam] < pupper))[0]

   #construct the masks using either the original or fitted lambdas
   masks = {}
   masks_long = {}
   for i in range(nbin):
      masks[i] = mask0.iloc[bins[i]]
      masks_long[i] = lineprops.iloc[bins[i]]
      masks[i] = masks[i].iloc[np.argsort(masks[i]["lambda"])] #sort masks by wavelength
      masks_long[i] = masks_long[i].iloc[np.argsort(masks_long[i]["lambda"])] #sort masks by wavelength
      if saveMasks:
         saveStr = 'VALD'+('_nonDP' if not depthPercentile else '')+'_species='+iron1Only+'_depthcutoff='+str(depth_cutoff)+'_overlapcutoff='+str(overlap_cutoff)+'_allowBlends='+(','.join(str(j) for j in allowBlends) if isinstance(allowBlends,list) else str(allowBlends))+'_badLineFilter='+badLineFilter+'_rejectTelluricSlope='+str(rejectTelluricSlope)+'_waves='+maskWavelengths+'_nbin='+str(nbin) + '_binParam='+binParam +'_n='+ str(i) + ("_VACUUM" if outputVacuumWavelengths else "_AIR")
         np.savetxt(os.path.join(outdir_masks, saveStr + '.mas'), masks[i])
         np.savetxt(os.path.join(outdir_masks, saveStr + '_long.mas'), masks_long[i])
         np.save(os.path.join(outdir_bins, saveStr + '.npy'), lineprops[binParam].iloc[bins[i]])
   return masks, masks_long




"""

#starting masks for EXPRES

for binParam in ["depth"]:
   for nbin in [1,8]:
      for maskWavelengths in ['Reiners101501']:
         for overlap_cutoff in [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5]:
            for iron1Only in ['all']:
               for allowBlends in [0]: #, [0,1], [0,1,2]]:
                  for rejectTelluricSlope in [0.]:
                     for badLineFilter in ['HD101501']:
                        print(time.time())
                        getVALDmasks(nbin=nbin, binParam = binParam, lineprops=lineprops_G8, iron1Only=iron1Only, overlap_cutoff=overlap_cutoff, maskWavelengths = maskWavelengths, allowBlends=allowBlends, rejectTelluricSlope=rejectTelluricSlope, badLineFilter=badLineFilter)

print(time.time())


for binParam in ["depth"]:
   for nbin in [1,8]:
      for maskWavelengths in ['Reiners26965']:
         for overlap_cutoff in [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5]:
            for iron1Only in ['all']:
               for allowBlends in [0]: #, [0,1], [0,1,2]]:
                  for rejectTelluricSlope in [0.]:
                     for badLineFilter in ['HD26965']:
                        print(time.time())
                        getVALDmasks(nbin=nbin, binParam = binParam, lineprops=lineprops_K1, iron1Only=iron1Only, overlap_cutoff=overlap_cutoff, maskWavelengths = maskWavelengths, allowBlends=allowBlends, rejectTelluricSlope=rejectTelluricSlope, badLineFilter=badLineFilter)

print(time.time())


for binParam in ["depth"]:
   for nbin in [1,8]:
      for maskWavelengths in ['Reiners34411']:
         for overlap_cutoff in [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5]:
            for iron1Only in ['all']:
               for allowBlends in [0]: #, [0,1], [0,1,2]]:
                  for rejectTelluricSlope in [0.]:
                     for badLineFilter in ['HD34411']:
                        print(time.time())
                        getVALDmasks(nbin=nbin, binParam = binParam, lineprops=lineprops_G2, iron1Only=iron1Only, overlap_cutoff=overlap_cutoff, maskWavelengths = maskWavelengths, allowBlends=allowBlends, rejectTelluricSlope=rejectTelluricSlope, badLineFilter=badLineFilter)

print(time.time())


for binParam in ["depth"]:
   for nbin in [1,8]:
      for maskWavelengths in ['Reiners10700']:
         for overlap_cutoff in [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5]:
            for iron1Only in ['all']:
               for allowBlends in [0]: #, [0,1], [0,1,2]]:
                  for rejectTelluricSlope in [0.]:
                     for badLineFilter in ['HD10700']:
                        print(time.time())
                        getVALDmasks(nbin=nbin, binParam = binParam, lineprops=lineprops_G8, iron1Only=iron1Only, overlap_cutoff=overlap_cutoff, maskWavelengths = maskWavelengths, allowBlends=allowBlends, rejectTelluricSlope=rejectTelluricSlope, badLineFilter=badLineFilter)

print(time.time())

"""


