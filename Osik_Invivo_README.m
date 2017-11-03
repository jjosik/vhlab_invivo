%Osik In-vivo analysis README file
%
%%%%General info
% Datasets and results for VF and all analyses other than STRF can be found
% under Documents/intracellular_data folder.  In general, all analysis code must
% be run under this home directory.  It is important to know that when
% running batch-level code, the selection window the pops up is a little
% clunky in it's rationale:  first a directory window will open, if the
% window that you want to select is on the current level then directly click
% 'okay' and select the folder from the following dialog/selection box that pops up.
% If the folder or multiple folders (as in the case of multiple cell
% folders) you want are one level down, click the appropriate top
% level folder in the directory window (e.g. the 'old' folder), then click
% 'okay' and select the subfolders from the resulting pop-up window (you
% can select multiple cell folders here).  NOTE: some of the newer batch level
% codes can handle multiple selections at the top level (e.g. both 'old'
% and 'young') and as it runs through the sets will ask you to select the
% appropriate cell folders.  Older batch level codes don't do this.  
% In those cases, you'll have to run separate iterations of the
% analysis for 'old' and 'young' folders. 

% Recent batch analysis variables and results have been saved under .mat
% files directly under the intracellular_data folder.  For example,
% collected_params.mat contain the VF tanh fit parameters for all cells
% used in the most recent analyses (Including:  all cells listed under the
% 'old' folder, and all cells under the 'young' folder EXCEPT cell 121).
% No_flankers_alltrial_variables.mat contains the results of the most
% recent prediction analysis, etc.

%%%%*******VF analysis
% The most straightforward way to run VF analysis is to use
% *****VF_allstim_rerunX.m - 

%       NOTE - this is one example of a batch analysis that will only 
%       run 'old' and 'young' folders separately.
%
%       
%The intention behind the code is to run the VF fits for all individual cells, so it 
% accesses invivo_VF_rerunX.m for each.  If you then want to take these results 
% and get a grand averaged fit, use compare_collVF_byage.m

%%%%% VF_allstim_rerunX.m analysis parameter notes
%      -- It shouldn't be necessary to change any basic input parameter
%       settings at the top of this file.  But it might be useful to know
%       what they do.  Save_it, overwrite, and fit_it are exactly what they
%       sound like: save everything, make all fits.  Overwrite has been overrided, so
%       currently it doesn't matter what the setting is, everything gets overwritten if the
%       whole analysis is rerun.  I've been dealing with this by using the
%       low-tech solution of changing the mat-file naming conventions
%       everytime I run a new iteration. This can be done under the 'code'
%       setting; any change in this name will be incorporated into the
%       saved mat files in each cell. (See below for information about the
%       current naming convention).  
%       display_spikes - will show the full experiment Vm trace with
%           threshold locations as the analysis is run; this can be nice, but
%           the plots don't close automatically.  So if you're running over
%           a large number of cells the analysis gets really slowed down
%           due to a large number of open plots.  Recommend '0' setting.
%       reset2phys - allows you to keep the baseline corrected Vm at
%           the absolute values rather than the relative value with
%           respect to the pre-trial baseline (i.e. 5 mV above a baseline
%           of -75 after the correction isn't 5 mV, it's -70)
%       auto_detect - detects spike through coarse threshold of -20 mV
%           rather than a user supplied cutoff, when set to '0' will prompt the
%           user.
%       use_detrend - is basically defunct; it's baseline correction using the
%           LOWESS smoothing procedure;  has been replaced by
%           use_trial_baseline
%       use_trial_baseline - the baseline correction method we want to use
%       default - basically a switch to determine VF fitting options; we've
%           been using option '2'
%       use_global_filter - only applies to the optimized filter for assessing FR, so
%           setting shouldn't matter when 50 ms flat filter
%       filter_type - option '1' sets up "flat" filter size of 50 ms
%       adapt_bins - adapts bin size in the VF heatmap
%       anchor_Vth - gets intercept by finding mean of Vm distribution for
%           which FR = 0, then is used as constraint in fit
%       model - setting no longer matters, as all three fit options rectified linear,
%            power law, and tanh are made and saved
%       code - see below
%       get_VmSTA - set to '1' to make additional Vm STA computations
%       stack_rasters - see notes
%       smooth_Vm - sets smoothing parameters; options are misleading as
%           they no longer match '0/1' to 'no/yes'; long story short: '0' is
%           the smoothing option we want
%       rebin_trials - see notes
%       supp_singletrial_FR - see notes
%       
%
% *****invivo_VF_rerunX.m - 
%       This function runs all the nitty-gritty analysis for each cell.
%       It contains all the different sub-functions that also may be of
%       interest.  For example separate functions for Vm extraction, spike
%       identification, spike thresholds, fits are all found under this
%       function.  Most have names like find_spikes.m or get_thresholds.m,
%       and so should be fairly self-descriptive.
%      --if you use this function to run cell-by-cell analysis from
%       scratch, be aware that it is set up for some user involvement.  For
%       example, coarse thresholds for spike detection are requested for
%       each dataset.  The code holds for user input by clicking on the 
%       appropriate location in the Vm trace (in my experience, choosing as high
%       a cutoff as possible gives the best results) and pressing 'enter'.  The
%       code also prompts the user if they want to edit the full datatrace.
%       This is useful in cases were the cell was lost before the
%       record ended, so the baseline is very high and not consistent
%       with the rest of the record, but otherwise the rest of the record
%       is good.  If you choose to use this option, it's useful for
%       chopping off a part of the record that will mess up the rest of the
%       analysis.  Answer is '1' for yes, '0' for no.  If the answer is '1', it will
%       open another prompt box asking for the number of cuts to the made
%       to the record.  This is a little misleading since, at the moment, the code is
%       only set up to handle one cut, and that cut simply excises
%       everything in the record to the right (i.e. later in time).  One
%       cut has been sufficient however in the analyses to date, as there are
%       only one or two cells where this is needed, and in those cases the
%       cuts where at the very end of the recording.
%       
%
% *****compare_collVF_byage.m -
%       This function generates the grand averaged VF plots from saved cell fits.
%       key facts:
%      -- this function requires three inputs: at this point you are
%       unlikely to want any input values other than plot_cells = 0,
%       bin_cellX = 1, and use_trial_baseline = 1
%      -- note: as this function runs through all cells desired in the
%       analysis, it has a slightly labor-intensive element in that it
%       requires the user to specify which .mat file is needed from the
%       cell folder.  So as the analysis runs, you'll be prompted to select 
%       from a large list of files for each individual cell.  These are the 
%       MAT files generated by invivo_VF_rerunX.m.  Each of these files 
%       have a naming systems that's probably too complex to detail here, 
%       but it's sufficient to know that the most recent analysis for every 
%       cell is saved under the same mat-file name:
%******VERY IMPORTANT - the naming convention used in the final analysis
%is:    GsmoothM1_F_BL_X_collstim_Xbinfitdata.mat
%This is the saved analysis mat-file that you want to access for every cell.
%Other variants saved include:  na10M1_F_BL_X... etc. for example, which is
%the analysis done on Vm data extracted using a 10 ms median filter.
%
%      --note regarding grand-averaged VF plot output:  this fit is NOT
%       rectified.  Applying rectification at the fit level would prevent the fit 
%       from finding the intercept in an unbiased manner.  So, in all cases where
%       a tanh fit has been made, it's necessary to apply rectification only at 
%       the plot-generating level, or when using these fits to run switched model 
%       predictions.  In fact, it was noted in the defense discussion that one 
%       of the prediction datapoints yields a DSI measure higher than 1, this is 
%       likely due to the negative values from an unrectified tanh fit. 
%      --all variables from analysis under compare_collVF_byage.m are saved
%      under the intracellular_data folder as collected_params_vars.mat.
%      -- ALSO NOTE - this function needed some minor troubleshooting/revision that
%       I was able to just fix today (6/17/17), so you'll need to access these
%       updated m-files rather than the versions you obtained when copying
%       the Documents folder.  

% *****get_ResponseErrorBars.m
%       This is the prediction analysis applied to time-averaged trial
%       response data
% *****get_ResponseErrorBars_trials2.m
%       This is the timestep-by-timestep prediction analysis
%      --The code as it stands right now is hard-coded to analyze the
%       pref and null responses only (by referencing the second entry
%       in an array that contains the peak response (i.e. on pref) and it's
%       flanking stimulus directions where the order of the array is
%       [pref-30 pref pref+30]).  Getting it to collate the three response
%       sets as pref or null again is simply a matter of reinstating the index for
%       the appropriate loop  (i.e. index 'i' was replaced with '2').
 
%Final General note
%   Additional analysis options were written into earlier code invivo_cellmain3.m.  
%   Results concerning in-trial CV and CV2 by stimulus, phase portraits, 
%   single cell spike parameters, cycle-averaged results, just to name a few things, 
%   can be found there.  

