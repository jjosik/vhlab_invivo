function [ ] = get_rebin_trials( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   This function was created to get better resolution within-trial binned 
%estimates of Vm and FR, for use in bootstrapped DSI predictions
%(previously done with full trial time-averages; yielded poor predictions)
%rebin_size - in ms, has been run using 50 ms

[folder,directories] = walk_directories();
cd(folder);

opref_cell_rp_rasters = {};
opref_cell_rs_vmtrials = {};
opref_cell_rp_mISIrate = {};
opref_cell_rebin_trialFR = {};
onull_cell_rp_rasters = {};
onull_cell_rs_vmtrials = {};
onull_cell_rp_mISIrate = {};
onull_cell_rebin_trialFR = {};
ypref_cell_rp_rasters = {};
ypref_cell_rs_vmtrials = {};
ypref_cell_rp_mISIrate = {};
ypref_cell_rebin_trialFR = {};
ynull_cell_rp_rasters = {};
ynull_cell_rs_vmtrials = {};
ynull_cell_rp_mISIrate = {};
ynull_cell_rebin_trialFR = {};

for i = 1:size(directories,1),
    cd(strtrim(directories(i,:)));
    sub_folder = uigetdir();
    d = dir(sub_folder);
    isub = [d(:).isdir];
    sel_Snames = {d(isub).name}';
    sel_Snames(ismember(sel_Snames,{'.','..','Codes'})) = [];
    [s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',sel_Snames);
    sub_dir = char(sel_Snames(s))
    [~,subfile,~] = fileparts(sub_folder);
    for j = 1:length(sub_dir),
        cur_dir = sub_dir(j,:);
        cd(cur_dir);
        load('analyzed_Vm.mat','total_vm_byAngle','stimvalues');
        load('analyzed_spike.mat','total_spike_byAngle');
        %load('na10M1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        %load('GsmoothM1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        load('GsmoothM1_F_BL_X_REBINNED_trial_responses.mat','pref_ind','null_ind',...
            'rebin_sub_trial_array','rebin_part_rasters','rebin_part_mISI',...
            'rebin_trial_FR');
        all_dirind = 1:length(stimvalues)-1;
        down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
        up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
        down_i = circshift(all_dirind,1,2);
        up_i = circshift(all_dirind,-1,2);
        pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
        null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
        if strcmp(subfile,'old'),
            opref_cell_rp_rasters{end+1,1} = rebin_part_rasters(:,pref_set_index(:,1));
            opref_cell_rs_vmtrials{end+1,1} = rebin_sub_trial_array(pref_set_index(:,1),:);
            opref_cell_rp_mISIrate{end+1,1} = rebin_part_mISI(:,pref_set_index(:,1));
            opref_cell_rebin_trialFR{end+1,1} = rebin_trial_FR(pref_set_index(:,1),:);
            onull_cell_rp_rasters{end+1,1} = rebin_part_rasters(:,null_set_index(:,1));
            onull_cell_rs_vmtrials{end+1,1} = rebin_sub_trial_array(null_set_index(:,1),:);
            onull_cell_rp_mISIrate{end+1,1} = rebin_part_mISI(:,null_set_index(:,1));
            onull_cell_rebin_trialFR{end+1,1} = rebin_trial_FR(null_set_index(:,1),:);
        else
            ypref_cell_rp_rasters{end+1,1} = rebin_part_rasters(:,pref_set_index(:,1));
            ypref_cell_rs_vmtrials{end+1,1} = rebin_sub_trial_array(pref_set_index(:,1),:);
            ypref_cell_rp_mISIrate{end+1,1} = rebin_part_mISI(:,pref_set_index(:,1));
            ypref_cell_rebin_trialFR{end+1,1} = rebin_trial_FR(pref_set_index(:,1),:);
            ynull_cell_rp_rasters{end+1,1} = rebin_part_rasters(:,null_set_index(:,1));
            ynull_cell_rs_vmtrials{end+1,1} = rebin_sub_trial_array(null_set_index(:,1),:);
            ynull_cell_rp_mISIrate{end+1,1} = rebin_part_mISI(:,null_set_index(:,1));
            ynull_cell_rebin_trialFR{end+1,1} = rebin_trial_FR(null_set_index(:,1),:);
        end
        clear rebin_part_rasters rebin_sub_trial_array;
        clear rebin_part_mISI rebin_trial_FR;
        clear pref_set_index null_set_index;
        cd ../
    end
    if strcmp(subfile,'old')
        old_cellnames = sub_dir;
    else
        young_cellnames = sub_dir;
    end
    cd ../
end

save('rebin_variables.mat');
        
        
end

