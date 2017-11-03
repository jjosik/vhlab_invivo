function [ output_args ] = get_ResponseErrorBarsBinned( allow_flip )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%THIS IS WRONG - DO NOT USE - ALSO DO NOT USE CONTENTS OF
%REBIN_VARIABLES.MAT -- NOT WHAT WAS NEEDED
load('collected_params.mat');

[folder,directories] = walk_directories();
cd(folder);

se_opv = [];
se_ops = [];
se_onv = [];
se_ons = [];
se_ypv = [];
se_yps = [];
se_ynv = [];
se_yns = [];
se_odsi = [];
se_ydsi = [];
o_m_mps = [];
o_m_mns = [];
o_se_mps = [];
o_se_mns = [];
y_m_mps = [];
y_m_mns = [];
y_se_mps = [];
y_se_mns = [];
m_model_odsi = [];
m_model_ydsi = [];
se_model_odsi = [];
se_model_ydsi = [];

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
        load('GsmoothM1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        load('rebin_variables.mat','opref_cell_rs_vmtrials','onull_cell_rs_vmtrials',...
            'ypref_cell_rs_vmtrials','ynull_cell_rs_vmtrials','opref_cell_rebin_trialFR',...
            'onull_cell_rebin_trialFR','ypref_cell_rebin_trialFR','ynull_cell_rebin_trialFR',...
            'opref_cell_rp_mISIrate','onull_cell_rp_mISIrate','ypref_cell_rp_mISIrate',...
            'ynull_cell_rp_mISIrate','old_cellnames','young_cellnames');
        if strcmp(subfile,'old'),
            pref_Vm = opref_cell_rs_vmtrials(j,1);
            pref_spike = opref_cell_rebin_trialFR(j,1);
            pref_spike_ = opref_cell_rp_mISIrate(j,1);
            null_Vm = onull_cell_rs_vmtrials(j,1);
            null_spike = onull_cell_rebin_trialFR(j,1);
            null_spike_ = onull_cell_rp_mISIrate(j,1);
        else
            pref_Vm = ypref_cell_rs_vmtrials(j,1);
            pref_spike = ypref_cell_rebin_trialFR(j,1);
            pref_spike_ = ypref_cell_rp_mISIrate(j,1);
            null_Vm = ynull_ell_rs_vmtrials(j,1);
            null_spike = ynull_cell_rebin_trialFR(j,1);
            null_spike_ = ynull_cell_rp_mISIrate(j,1);
            
        end
        nsims = 1000;
        %first generate errorbars for actual and predicted mean cell responses
        samp_dist_p_meas = [];
        samp_dist_p_pred = [];
        samp_dist_n_meas = [];
        samp_dist_n_pred = [];
        for n = 1:nsims,
            if allow_flip == 1,
                samp_pref_vm = datasample(pref_Vm,size(pref_Vm,1)*size(pref_Vm,2));
        
        

end

