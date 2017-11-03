function [ Vm_trialAve_smooth ] = Vm_SGsmooth( process_Vm_trialAve,sample_interva )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

win = round(0.05/sample_interval);
if mod(win,2)==0,
    win = win + 1;
else
end
Vm_trialAve_smooth = sgolayfilt(process_Vm_trialAve,5,win);

end

