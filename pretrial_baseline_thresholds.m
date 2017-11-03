function [ AG_Vm_th,trial_baseline,trial_epochbase ] = ...
    pretrial_baseline_thresholds( spike_trace,t_vec,...
    spike_locations,th_loc,Vm_th,sample_interval,stimvalues,stimorder,...
    nStim_ON,nStim_OFF)
%PRETRIAL_BASELINE_THRESHOLDS - calculates spike thresholds in experiment
%block relative to mean Vm baseline from 1 second immediately preceding stimulus
%onset for each trial 
%
%******INPUTS:      SPIKE_TRACE - full trial, preprocessed waveform
%                       containing all spikes
%                   T_VEC - time vector at acquisition sampling rate
%                   SPIKE_LOCATIONS - index values of individual spike
%                       times
%                   TH_LOC - index values of individual spike threshold
%                       times
%                   VM_TH - value in mV of threshold membrane potential
%                   SAMPLE_INTERVAL - time in seconds of single
%                       timestep
%                   STIMVALUES - monotonically increasing sequence of
%                       stimulus orientation angles used in experiment
%                   STIMORDER - sequence indicating the order of
%                       stimvalues used over the course of the full
%                       experiment
%                   NSTIM_ON - time in seconds into full experiment
%                       block that a new stimulus in the stimorder
%                       sequence turns ON.
%                   NSTIM_OFF - corresponding time in seconds into
%                       sequence that stimulus turns OFF.
%*****OUTPUTS:      AG_VM_TH - (AG denotes Azouz & Gray) n x 1 cell
%                       array where each entry contains full vector of
%                       relative Vm thresholds for each n stimulus trials
%                   TRIAL_BASELINE - n x 1 array of pre-stim mean Vm baseline
%                       values for each n stimulus trials
%                   TRIAL_EPOCHBASE - n x 1 cell array where each entry
%                       contains full vector of time indices over which
%                       each trial_baseline is evaluated

trial_epochbase = cell(length(stimorder),1);
for i = 1:length(stimorder),
    trial_epochbase{i,1} = round((nStim_ON(i,1)-1)/sample_interval):...
        round(nStim_ON(i,1)/sample_interval);
    trial_baseline(i,1) = mean(spike_trace(trial_epochbase{i,1},1));
end

intrial_epoch = cell(length(stimorder),1);
th_loc_subset = cell(length(stimorder),1);
Vm_th_subset = cell(length(stimorder),1);
for j = 1:length(stimorder),
    intrial_epoch{j,1} = round(nStim_ON(j,1)/sample_interval):...
        round(nStim_OFF(j,1)/sample_interval);
    th_loc_subset{j,1} = intersect(th_loc,intrial_epoch{j,1});
    Vm_th_subset{j,1} = spike_trace(th_loc_subset{j,1}(1:end));
end

AG_Vm_th = cell(length(stimorder),1);
for k = 1:length(stimorder),
    temp = [];
    for n = 1:length(Vm_th_subset{k,1}),
        temp = [temp;(Vm_th_subset{k,1}(n,1)-trial_baseline(k,1))];
    end
    AG_Vm_th{k,1} = temp;
end


end

