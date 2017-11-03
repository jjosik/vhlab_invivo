function [ Vm_trialAve_smooth ] = Vm_GaussianSmooth( process_Vm_trialAve,sample_interval )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

win = 0.05;
t_vec = sample_interval:sample_interval:length(process_Vm_trialAve);
Xn = [];
%for qq = 2:length(t_vec),
%    Xn(end+1) = mean([t_vec(1,qq) t_vec(1,qq-1)]);
%end
Xn = (t_vec(2:end)+t_vec(1:end-1))/2;
t_index = t_vec/sample_interval;

    if ~isempty(process_Vm_trialAve),
        pad_len = round(win/sample_interval/2);
        padded_Vm = [fliplr(process_Vm_trialAve(1:pad_len));process_Vm_trialAve;...
            fliplr(process_Vm_trialAve(length(process_Vm_trialAve)-pad_len:end))];
        Yn = padded_Vm;
        g_kernel = gausswin(round(win/sample_interval));
        norm_g_kernel = g_kernel./sum(g_kernel);
        temp_smooth = conv(Yn,norm_g_kernel,'same');
        Vm_trialAve_smooth = temp_smooth(pad_len+1:length(temp_smooth)-pad_len-1);
    else
        Vm_trialAve_smooth = zeros(length(t_vec)-1,1);
    end
   

end

