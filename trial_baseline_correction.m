function [ trial_baseline,trial_baseline_array,epoch_loc_array,corrected_Vmth_array ] =...
    trial_baseline_correction( Vm_waveform,t_vec,...
    sample_interval,stimvalues,stimorder,nStim_ON,nStim_OFF,th_loc,Vm_th)
%TRIAL_BASELINE_CORRECTION - Returns baseline pretrial Vm-mean values (over
%a 1 second epoch) in two ways:  1. in chronological order of experiment
%stim. presentations, and 2. by stim. value (M) and rep. (N) in an MxN array. 

trial_epochbase = cell(length(stimorder),1);
for i = 1:length(stimorder),
    trial_epochbase{i,1} = round((nStim_ON(i,1)-1)/sample_interval):...
        round(nStim_ON(i,1)/sample_interval);
    trial_baseline(i,1) = mean(Vm_waveform(trial_epochbase{i,1},1));
    trial_epoch{i,1} = round(nStim_ON(i,1)/sample_interval):...
        round(nStim_OFF(i,1)/sample_interval);
    epoch_loc{i,1} = intersect(th_loc,trial_epoch{i,1});
    corrected_Vmth{i,1} = Vm_waveform(epoch_loc{i,1}(1:end,1),1)-trial_baseline(i,1);
end

for j = 1:length(stimvalues),
    trial_nums = find(stimorder==j);
    for k = 1:length(trial_nums),
        trial_baseline_array(j,k) = trial_baseline(trial_nums(1,k),1);
        epoch_loc_array{j,k} = epoch_loc{trial_nums(1,k),1};
        corrected_Vmth_array{j,k} = corrected_Vmth{trial_nums(1,k),1};
    end
end


