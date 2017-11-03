function [ output_args ] = get_all_trials( )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[folder,directories] = walk_directories();
cd(folder);

opref_cell_FR_alltrials = {};
onull_cell_FR_alltrials = {};
ypref_cell_FR_alltrials = {};
ynull_cell_FR_alltrials = {};
opref_cell_Vm_alltrials = {};
onull_cell_Vm_alltrials = {};
ypref_cell_Vm_alltrials = {};
ynull_cell_Vm_alltrials = {};
o_samprate_bycell = [];
y_samprate_bycell = [];


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
            'rebin_trial_FR','FR_inst','sub_trial_array','sampling_rate');
        for a = 1:size(sub_trial_array,1),
            for b = 1:size(sub_trial_array,2),
                corr_sub_trial_array{a,b} = sub_trial_array{a,b}(1:length(sub_trial_array{a,b})-...
                    (round(0.5*sampling_rate)));
            end
        end
        all_dirind = 1:length(stimvalues)-1;
        down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
        up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
        down_i = circshift(all_dirind,1,2);
        up_i = circshift(all_dirind,-1,2);
        pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
        null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
        if strcmp(subfile,'old'),
            opref_cell_FR_alltrials{end+1,1} = FR_inst(pref_set_index(2,1),:);
            onull_cell_FR_alltrials{end+1,1} = FR_inst(null_set_index(2,1),:);
            opref_cell_Vm_alltrials{end+1,1} = corr_sub_trial_array(pref_set_index(2,1),:);
            onull_cell_Vm_alltrials{end+1,1} = corr_sub_trial_array(null_set_index(2,1),:);
            o_samprate_bycell(end+1,1) = sampling_rate;
        else
            ypref_cell_FR_alltrials{end+1,1} = FR_inst(pref_set_index(2,1),:);
            ynull_cell_FR_alltrials{end+1,1} = FR_inst(null_set_index(2,1),:);
            ypref_cell_Vm_alltrials{end+1,1} = corr_sub_trial_array(pref_set_index(2,1),:);
            ynull_cell_Vm_alltrials{end+1,1} = corr_sub_trial_array(null_set_index(2,1),:);
            y_samprate_bycell(end+1,1) = sampling_rate;
        end
        clear rebin_part_rasters rebin_sub_trial_array;
        clear rebin_part_mISI rebin_trial_FR;
        clear pref_set_index null_set_index;
        clear trial_FR_inst sub_trial_array corr_sub_trial_array;
        cd ../
    end
    if strcmp(subfile,'old')
        old_cellnames = sub_dir;
    else
        young_cellnames = sub_dir;
    end
    cd ../
end

%save('alltrial_variables.mat');
save('noflankers_alltrial_variables.mat');

end

