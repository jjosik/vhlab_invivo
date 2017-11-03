function [ waveform_dnlp,sample_interval,stimvalues,stimorder,nStim_ON,...
    nStim_OFF,otpref_total,stim_duration,reps_stimmatch,p ] = prep_rawspikes( reps,...
    spike_prefile,stims,stimsv,stimso,data)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

load(spike_prefile);
load(stims);
load(stimsv);
load(stimso);
load(data);

n = numStims(saveScript);
stimmatch_reps = length(stimorder)/n;
fullLeg_order = stimorder;  %replacement here and in line 48 bypasses myriad replacements below
reps_stimmatch = length(reps);
stimorder = [];
keep_indices = [];
if max(reps) > stimmatch_reps,
    reps = reps-(max(reps)-stimmatch_reps);
    repCheck_flag = 1;
else
    repCheck_flag = 0;
end
for p = 1:length(reps),
    stimorder = [stimorder,fullLeg_order(1,((reps(p)-1)*n)+1:reps(p)*n)];
    keep_indices = [keep_indices;(((reps(p)-1)*n)+1:reps(p)*n)'];
end
 
if length(stimorder) < length(fullLeg_order),
    nStim_ON = spike2data_Ch22.times(keep_indices,1);  
else
    nStim_ON = spike2data_Ch22.times;
end
nStim_OFF = [];
%**>>>review code below for conflicts with new code, lines 43-56
if mod(length(nStim_ON),n) > 0,
    nStim_ON = nStim_ON(1:length(nStim_ON)-mod(length(nStim_ON),n));
else
end
if length(nStim_ON) > length(stimorder),
    nStim_ON = nStim_ON(1:length(stimorder),1);
else
end
%**<<<
for j=1:length(nStim_ON),
    p = getparameters(get(saveScript,stimorder(j)));
    if isfield(p,'nCycles'),
        nCycles(j,1) = p.nCycles;
        stim_duration(j,1) = p.nCycles/p.tFrequency;
        temp_freq(j,1) = p.tFrequency;
        stim_angle(j,1) = p.angle;
        nStim_OFF(j,1) = nStim_ON(j,1)+stim_duration(j,1);
    else
        nCycles(j,1) = NaN;
        stim_duration(j,1) = NaN;  %results for blank stim are always NaN
        temp_freq(j,1) = NaN;
        stim_angle(j,1) = NaN;
        nStim_OFF(j,1) = NaN;
    end
end
stim_duration = stim_duration(~isnan(stim_duration));
nCycles = nCycles(~isnan(nCycles));
try
    sampling_rate = round(1/spike2data_Ch21.interval);
    sample_interval = spike2data_Ch21.interval;
    waveform = spike2data_Ch21.values.*100; %in mV
catch
    if ~exist('sampling_rate','var'),
        warning(['Problem with requested channel structure.',...
            '   Sampling rate not found under expected channel.',...
        '  Checking alternate source.']);
        sampling_rate = round(1/spike2data_Ch1.interval);
        sample_interval = spike2data_Ch1.interval;
        waveform = spike2data_Ch1.values.*100;
    end
end
            
[b,a] = butter(3,[(60-1.0)/floor(sampling_rate/2),...
    (60+0.5)/floor(sampling_rate/2)],'stop');
waveform_denoised = filter(b,a,waveform);
[bb,aa] = butter(3,(3000/floor(sampling_rate/2)),'low');
waveform_dnlp = filter(bb,aa,waveform_denoised);
clear a aa b bb;
clear waveform_denoised;

t_vec = sample_interval:sample_interval:(length(waveform_dnlp)*sample_interval);
f = figure;
plot(t_vec,waveform_dnlp,'k');
hold on;
prompt = 'Does this spike trace require editing?  Enter ''1'' or ''0'' ';
truncate_trial = input(prompt);
if (truncate_trial ~= 0) && (truncate_trial ~= 1),
    truncate_trial = 0;
else
end

if truncate_trial == 1,
    %***IMPORTANT - current implementation only supports excision to the
    %right of cut point, i.e. one cut; would NOT recommend working around 
    %this as there's no reasonable basis for putting detached parts back
    %together.
    prompt = {'Enter the number of cuts to be made: '};
    dlg_title = 'Input';
    num_lines = 1;
    defaultans = {'1'};
    ncuts = inputdlg(prompt,dlg_title,num_lines,defaultans);
    ncuts = eval(cell2mat(ncuts));
    display(['Use cursor to make ',num2str(ncuts),' cut(s):']);
    [x_break,~] = ginput(ncuts);
    waveform_dnlp = waveform_dnlp(1:round(x_break/sample_interval));
    %waveform = waveform(1:round(x_break/sample_interval));
    close(f);
else
end

if ~exist('otpref_total'),
    otpref_total = otpref_1;
else
end

end

