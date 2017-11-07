function [ ] = VF_allstim_rerunX( varargin )
%VF_ALLSTIM_RERUN - runs an isolated analysis of VF plots over all
%cells that collates all stim responses, disregarding stimulus direction.
%Works at the same level of analysis as batch_rerun.m.  Needs to access 
%invivo_VF_rerunX.m

save_it = 1;
overwrite = 0;
collate_stims = 1;  
display_spikes = 0;
reset2phys = 0;
auto_detect = 0; 
use_detrend = 0;
use_trial_baseline = 1; %probably but not strictly mutually exclusive with use_detrend
default = 2;
use_global_filter = 0;
filter_type = 1;
adapt_bins = 1;
anchor_Vth = 1;
fit_it = 1;
model = 1; %set to '0' if rectified power law, '1' if rectilinear
%code = 'na10M1_F_BL_X';
code = 'GsmoothM1_F_BL_X';
get_VmSTA = 0;
stack_rasters = 1;  %use '1' if stack and average trials, '0' if find FR from single trial rasters
smooth_Vm = 0;
rebin_trials = 0;   %use this to create alternate response estimates from bins over 
                    %individual trials; 
supp_singletrial_FR = 0;    %for trial rebinning and analysis, set to 1 if rebinning should be 
                            %applied to single-trial FR estimate (note:
                            %overrides stack_rasters setting)

if nargin < 2,
    use_pre_analyzed = 0;
else
    use_pre_analyzed = varargin{1};
end

sub_folder = uigetdir();
d = dir(sub_folder);
isub = [d(:).isdir];
sel_Snames = {d(isub).name}';
sel_Snames(ismember(sel_Snames,{'.','..','Codes'})) = [];
[s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',sel_Snames);

directories = char(sel_Snames(s));
cd(sub_folder);

[~,subfile,~] = fileparts(sub_folder);
for i = 1:size(directories,1),
    cur_dir = directories(i,:);
    cd(directories(i,:));
    load([sub_folder filesep directories(i,:) filesep 'stims.mat']);
    load([sub_folder filesep directories(i,:) filesep 'data.mat']);
    load([sub_folder filesep directories(i,:) filesep 'stimorder.mat']);
    load([sub_folder filesep directories(i,:) filesep 'stimvalues.mat']);
    load([sub_folder filesep directories(i,:) filesep 'analyzed_spike.mat']);
    [Vth_est_coll,sta_vm_pref,sta_vm_null,rebin_sub_trial_array] = ...
        invivo_VF_rerunX(save_it,rebin_trials,overwrite,model,collate_stims,...
        use_pre_analyzed,display_spikes,reset2phys,auto_detect,use_detrend,...
        use_trial_baseline,sub_folder,cur_dir,default,use_global_filter,...
        filter_type,adapt_bins,anchor_Vth,fit_it,code,get_VmSTA,stack_rasters,...
        smooth_Vm,supp_singletrial_FR);

end



