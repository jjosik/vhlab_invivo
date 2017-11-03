function invivo_cellmain3(save_it,use_windowing,detrend,adapt_bins,fit_it,...
    anchor_Vth,use_flanktest,sub_folder,directories,default,reps,...
    use_global_filter,filter_type,overwrite)

%Alectryon path:
%top = ('/Users/vhlab/Documents/intracellular_data/');
%osik_mac path:
%top = ('Users/osik_mac/Documents/MATLAB/data');
%p = addpath(genpath(top));

%select cell folder
%sub_folder = uigetdir();
%d = dir(sub_folder);
%isub = [d(:).isdir];
%sel_Snames = {d(isub).name}';
%sel_Snames(ismember(sel_Snames,{'.','..','Codes'})) = [];
%[s,v] = listdlg('PromptString','Select folders:',...
%    'SelectionMode','multiple',...
%    'ListString',sel_Snames);

%directories = char(sel_Snames(s));
%cd(sub_folder);

%%%%%////////////////////////////////////////////////////
%%%%%consider routing above code to higher level function - DONE 12/7
%%%%%////////////////////////////////////////////////////
%%
load([sub_folder filesep directories filesep 'stims.mat']);
load([sub_folder filesep directories filesep 'data.mat']);
load([sub_folder filesep directories filesep 'stimorder.mat']);
load([sub_folder filesep directories filesep 'stimvalues.mat']);
load([sub_folder filesep directories filesep 'analyzed_spike.mat']);

if overwrite == 0,
    d = dir(pwd);
    isub = [d(:).isdir];
    sub_d = {d(isub).name}';
    target = char(datetime('today'));
    if ~any(strcmp(target,sub_d)),
        mkdir(char(datetime('today')));
        cd(char(datetime('today')));
    else
        cd(char(datetime('today')));
    end
else
end

%shouldn't be necessary, as all data are in spike2 times - but using to
    %confirm
[mti2_,starttime_] = tpcorrectmti(MTI2,[sub_folder...
    filesep directories filesep 'stimtimes.txt'],1);
    
%rel_stims = find(stimorder==pref_stim_idx);
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


%%
 %*************PLOTTING RAW TRACES BY TRIAL*****************************
    %*************************************************************
    %plot the trace of each stim rep. and save plots under cell directory
margins = round(3.0/sample_interval);  %pad each stim trace with 3 s pre- 
                                       %and post-stimulus
t_vec = sample_interval:sample_interval:...
    (length(waveform_dnlp)*sample_interval);
rep_counts = zeros(length(stimvalues),1);
try 
    pref = otpref_total;
catch
    try
        pref = pref_stim_idx;
    catch
        pref = otpref_1;
    end
end
for k=1:length(nStim_ON),
    rep_counts(stimorder(1,k)) = rep_counts(stimorder(1,k))+1;    
    f = figure;
    a = round(nStim_ON(k,1)/sample_interval);
    if ~isnan(nStim_OFF(k,1)),
        b = round(nStim_OFF(k,1)/sample_interval);
        if ((a-margins)<0||(b+margins)>length(t_vec))&&(k==length(nStim_ON)),
            margins = min([a;(length(t_vec)-b)]);
        else
            margins = round(3.0/sample_interval);
        end
        plot(t_vec(1,(a-margins):(b+margins)),...
            waveform_dnlp(a-margins:b+margins,1),'k');
        hold on;
        xlabel('Time');
        ylabel('Membrane potential (mV)');
        he = line([t_vec(1,b) t_vec(1,b)],[-120 80]);
        he.LineStyle = '--';
    else
        plot(t_vec(1,(a-margins):(a+round(5.0/sample_interval)+margins)),...
            waveform_dnlp((a-margins):(a+round(5.0/sample_interval)...
            +margins),1),'k');
        hold on;
        xlabel('Time');
        ylabel('Membrane potential (mV)');
        he = line([t_vec(1,a+round(5.0/sample_interval)) t_vec(1,...
            (a+round(5.0/sample_interval)))],[-120 80]);
        he.LineStyle = '--';
    end
        hs = line([t_vec(1,a) t_vec(1,a)],[-120 80]);
        hs.LineStyle = '--';
        hold off;
        if save_it == 1,
            if stim_angle(k,1) ~= pref,          % old variable name:...
                                                        %pref_stim_idx,
                saveas(gcf,[pwd filesep 'raw_trace' '_angle',...
                    num2str(stimvalues(stimorder(1,k))), '_rep',...
                    num2str(rep_counts(stimorder(1,k))), '.fig']);
            else
                saveas(gcf,[pwd filesep 'raw_trace' 'PREF_angle',...
                    num2str(stimvalues(stimorder(1,k))), '_rep',...
                    num2str(rep_counts(stimorder(1,k))), '.fig']);
            end
        else
        end
        close(f);
end


%%
%**************Generate the firing rate histograms*****************
%******************************************************************
auto_detect = 0;  %add to inputs
bins = 0.1;
h_margins = 0.5;
if auto_detect == 1,
    th_ = -20.0;  %in mV
else
    f = figure;
    plot(t_vec,waveform_dnlp,'k');
    hold on;
    display(['Use cursor to select coarse detection threshold for the experiment trace.  Press ''Return'' when finished']);
    [~,th_] = ginput;
    close(f);
end

spike_locations = [];
spike_stepLocations = find(waveform_dnlp>th_);
if isempty(spike_stepLocations),
    msgID = 'MATLAB:th_test';
    msgtext = 'No thresholds were detected in this dataset';
    ME = MException(msgID,msgtext);
    throw(ME)
else
end
onsets = diff(spike_stepLocations);
sp_edge = find(onsets>1);
ind_edge = [spike_stepLocations(1,1);spike_stepLocations(sp_edge(:,1)+1,1)];
%if ~exist('modulation_ratio_spike','var'),
%    modulation_ratio_spike = 1.1;
%else
%end                  %*this block corrected from process_spiketrace_FFT 3/5/2017
%if modulation_ratio_spike > 1,
%    for ii = 1:length(ind_edge),
%        if ii<length(ind_edge),
%            subset = find((spike_stepLocations>ind_edge(ii,1))&...
%                (spike_stepLocations<ind_edge(ii+1,1)));
%            %subset = spike_stepLocations(find((spike_stepLocations>ind_edge(ii,1))&...
%            %    (spike_stepLocations<ind_edge(ii+1,1))));
%        else
%            subset = find(spike_stepLocations>ind_edge(ii,1));
%            %subset = spike_stepLocations(find(spike_stepLocations>ind_edge(ii,1)));
%        end
%        total_subset = [ind_edge(ii,1);spike_stepLocations(subset,1)];
%        %total_subset = [ind_edge(ii,1);subset];
%        [~,sub_loc] = max(waveform_dnlp(total_subset(1:end),1));
%        spike_locations(ii,1) = ind_edge(ii,1)+(sub_loc-1);
%    end
%else
%    for ii = 1:length(ind_edge),
%        if ii<length(ind_edge),
%            subset = spike_stepLocations(find((spike_stepLocations>ind_edge(ii,1))&...
%                (spike_stepLocations<ind_edge(ii+1,1))));
%        else
%            subset = spike_stepLocations(find(spike_stepLocations>ind_edge(ii,1)));
%        end
%        total_subset = [ind_edge(ii,1);subset];
%        [~,sub_loc] = max(waveform_dnlp(total_subset(1:end),1));
%        spike_locations(ii,1) = ind_edge(ii,1)+(sub_loc-1);
%    end
%end
for ii = 1:length(ind_edge),
    if ii<length(ind_edge),
        subset = find((spike_stepLocations>ind_edge(ii,1))&...
            (spike_stepLocations<ind_edge(ii+1,1)));
    else
        subset = find(spike_stepLocations>ind_edge(ii,1));
    end
    total_subset = [ind_edge(ii,1);spike_stepLocations(subset,1)];
    [~,sub_loc] = max(waveform_dnlp(total_subset(1:end),1));
    spike_locations(ii,1) = ind_edge(ii,1)+(sub_loc-1);
end
    
        
elim_stim = [];
bc_array = cell(length(nStim_ON),1);
be_array = cell(length(nStim_ON),1);
psth_array = cell(length(nStim_ON),1);
for m=1:length(nStim_ON),
    if ~isnan(nStim_OFF(m,1)),
        if ismember(stimorder(1,m),elim_stim),
            continue;
        else
        end
        rep_locs = find(stimorder==stimorder(1,m));
        elim_stim = [elim_stim, stimorder(1,m)];
        xRange = (nStim_OFF(m,1)-nStim_ON(m,1)+h_margins);
        newxRange = (floor(xRange*(1/bins)))/(1/bins);
        disp_edges = [0:bins:newxRange];
        psth = zeros(length(disp_edges),1);
        for n = 1:length(rep_locs),
            count_edges = (nStim_ON(rep_locs(n))/sample_interval):...
                (bins/sample_interval):...
                ((nStim_ON(rep_locs(n))/sample_interval)+...
                (newxRange/sample_interval));
            %new_t_vec = 0:sample_interval:newxRange;
            %count_range = [nStim_ON(rep_locs(n))/sample_interval nStim_OFF(rep_locs(n))/sample_interval];
            %new_spikes = find(spike_locations>count_range(1,1)&spike_locations<count_range(1,2));
            psth =  psth + reshape(histc(spike_locations,count_edges),...
                length(psth),1);
        end
        bin_centers = (disp_edges(1:end-1)+disp_edges(2:end))/2;
        alt_bin_centers = (count_edges(1:end-1)+count_edges(2:end))/2;
        bc_array{stimorder(1,m),1} = alt_bin_centers;
        be_array{stimorder(1,m),1} = disp_edges;
        psth = psth(1:end-1)/bins;
        psth_array{stimorder(1,m),1} = psth;
        f1 = figure;
        bar(bin_centers,psth,1.0);
        hold on;
        xlim([-1 newxRange+1]);
        ylim([-0.5 max(psth)+0.1*(max(psth))]);
        xlabel('Time');
        ylabel('Spike rate');
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'psth_angle',...
                num2str(stimvalues(stimorder(rep_locs(1)))), '.fig']);
        else
        end
        close(f1);
    else
    end
end
%%
        %*******Assess biophysical spike threshold using Azouz & Gray, 1999
        %*******methods****************************************************
        %******************************************************************
        reset2phys = 0;         %add to function inputs
        use_buffer = 0;
        
        alt_t_slope = gradient(waveform_dnlp,sample_interval); 
        slope_criterion = 0.03; 
        %fraction of peak dV/dt calibrated by visual inspection to match threshold
        ex_flag = [];
        if use_buffer == 1,
            th_sort = diff(spike_locations);
            sub_locations = spike_locations;
            sub_locations(find(th_sort>round(0.02/sample_interval))+1) = [];
        else
            sub_locations = spike_locations;
        end
        peak_slope = zeros(length(sub_locations),1);
        slope_Vm = zeros(length(sub_locations),1);
        threshold = zeros(length(sub_locations),1);
        Vm_th = zeros(length(sub_locations),1);
        for q = 1:length(sub_locations),
            search_size = ceil(0.002/sample_interval);
            if search_size > sub_locations(q,1),
                ex_flag = [ex_flag;q];
                continue;
            else
            end
            search_pad = (sub_locations(q,1)-search_size):...
                sub_locations(q,1)-1;
            [max_slope,peak_ind] = max(alt_t_slope(reshape(search_pad,...
                length(search_pad),1),1));
            peak_slope(q,1) = search_pad(1,1)+(peak_ind-1);
            new_search = search_pad(1,1):peak_slope(q,1);
            slope_Vm(q,1) = waveform_dnlp((search_pad(1,1)+(peak_ind-1)),1);
            %[th_slope,th_ind] = min(abs((slope_criterion*max_slope)-...
            %    alt_t_slope(search_pad,1)));
            [th_slope,th_ind] = min(abs((slope_criterion*max_slope)-...
                alt_t_slope(new_search,1)));
            threshold(q,1) = search_pad(1,1)+(th_ind-1);
            Vm_th(q,1) = waveform_dnlp((search_pad(1,1)+(th_ind-1)),1);
        end
        if any(ex_flag),
            threshold(ex_flag) = [];
            Vm_th(ex_flag) = [];
        else
        end
        f4 = figure;
        plot(t_vec(1,:),waveform_dnlp(:,1),'k');
        hold on;
        scatter(t_vec(1,threshold(:,1)),Vm_th(:,1),'x','r');
        %***Compute low-pass trend-line fit to thresholds for this cell
        %***recording.  Use result to subtract baseline for trend
        %***correction that is sensitive to electrical drift without altering
        %***state-dependent changes.
        if detrend == 1,
            %***alt. 1 (these methods require fixing edge effects)
            %xq = 20:20:max(t_vec);
            %tr_line = interp1(t_vec(1,threshold(:,1)),Vm_th(:,1),xq,'linear');
            %fc = 0.9;
            %[b,a] = butter(4,fc,'low');
            %lp_trline = filter(b,a,tr_line);
            %***alt. 2
            %ord = 3;
            %frame = 501;
            %sg_filt = sgolayfilt(Vm_th,ord,frame);
            %***"chosen" method: lowess smoothing
            default_span = 51;
            check_interval = diff(t_vec(1,threshold(:,1)));
            if max(check_interval) < 5000,
                span = default_span;
            else
                span = 501;  %may not need to change this; need more test cases
            end
            ytr_line = smooth(t_vec(1,threshold(:,1)),Vm_th(:,1),span,'lowess');
            Vq = interp1(t_vec(1,threshold(:,1)),ytr_line,t_vec(1,:),'linear','extrap');
            %alt. nonuniform interp method (not optimal due to edge effects)
            %Vq = resample(ytr_line,t_vec(1,threshold(:,1)),sampling_rate);
            sub_wv = waveform_dnlp-Vq';
            if reset2phys == 1,
                initial_Vm = median(spike_trace(1:default_span));
                new_trace = sub_wv+initial_Vm;
            else
                new_trace = sub_wv;
            end
            
            plot(t_vec,Vq,'k--');
            plot(t_vec(1,threshold(:,1)),ytr_line,'b--');
        else
        end
        if save_it == 1,
            saveas(gcf,[pwd filesep 'Sp_threshold_locations.fig']);
            close(f4);
        else
        end
        if detrend == 1,
            f5 = figure;
            plot(t_vec,sub_wv,'k');
            hold on;
            scatter(t_vec(1,threshold(:,1)),Vm_th-(Vq(1,threshold(:,1))'),'x','r');
            if save_it == 1,
                saveas(gcf,[pwd filesep 'Sp_threshold_locations_DETRENDED.fig']);
                close(f5);
            else
            end
        else
        end
        
%%
%*****Remove spikes to get raw Vm traces and trial-averaged Vm trace**
%*********************************************************************
win_d = ceil(0.004/sample_interval);
%alt: win_d = 41;  %corresponds to the 4-ms used in Priebe et al. 2004
if mod(win_d,2) == 0,
    win_d = win_d+1;
else
end
if detrend == 1,
    Vm_waveform = medfilt1(sub_wv,win_d);
else
    Vm_waveform = medfilt1(waveform_dnlp,win_d);
end
trial_Range = cell(3,1);
Vm_array = cell(length(stimvalues),1);
%Vm_trial_array = cell(length(stimvalues),stimmatch_reps);
Vm_trial_array = cell(length(stimvalues),length(reps));
dRange_array = cell(length(stimvalues),1);
for jj = 1:length(stimvalues),
    %if ~isnan(nStim_OFF(jj,1)),
    trial_nums = find(stimorder==jj);       
    for ii = 1:length(trial_nums),
        if isnan(nStim_OFF(trial_nums(ii),1)),
            nan_flag = 1;
            continue;
        else
        end
        trial_Range{ii} = [(nStim_ON(trial_nums(ii),1)...
            /sample_interval) (nStim_OFF(trial_nums(ii),1)...
            /sample_interval)+(h_margins/sample_interval)];
        %ave_Range = nStim_OFF(trial_nums(ii),1)-nStim_ON(trial_nums(ii),1)+h_margins;
        %newAve_Range = (floor(ave_Range*(1/bins)))/(1/bins);
        disp_Range = sample_interval:sample_interval:...
            (nStim_OFF(trial_nums(ii),1)+...
            h_margins-nStim_ON(trial_nums(ii),1));
        Vm_trial_array{jj,ii} = Vm_waveform(trial_Range{ii}(1,1):...
            trial_Range{ii}(1,2)-1);
    end
    if exist('nan_flag'),
        continue;
    else
    end
    Vm_trialAve = zeros(length(disp_Range),1);
    for ii = 1:length(trial_nums),
        Vm_trialAve = Vm_trialAve + Vm_waveform(trial_Range{ii}(1,1):...
            trial_Range{ii}(1,2)-1);
    end
    Vm_trialAve = Vm_trialAve./length(trial_nums);
    Vm_array{jj} = Vm_trialAve;
    dRange_array{jj} = disp_Range;
    f2 = figure;
    plot(disp_Range,Vm_trialAve,'k');
    hold on;
    xlim([-1 disp_Range(1,end)+1]);
    ylim([-(0.1*max(abs(Vm_trialAve))+max(abs(Vm_trialAve))) 0]);
    xlabel('Time');
    ylabel('Membrane potential (mV)');
    hold off;
    if save_it == 1,
        saveas(gcf,[pwd filesep 'Vm_trialAve_angle',...
            num2str(stimvalues(jj)), '.fig']);
    else
    end
    close(f2);
    %clear Vm_trialAve;  %may not be desirable under some circumstances,
                        %see: cycle phase response block
    %else
    %end
end

%%
stim_rasters = cell(length(stimvalues),1);
%part_rasters = cell(stimmatch_reps,length(stimvalues));
for b = 1:length(stimvalues),
    current_epochs = find(stimorder == b);
    align_on_ind = nStim_ON(current_epochs).*sampling_rate;
    align_off_ind = nStim_OFF(current_epochs).*sampling_rate;
    for a=1:length(current_epochs),
        rasters{a,1} = spike_locations(find(spike_locations>align_on_ind(a,1)&spike_locations<align_off_ind(a,1)))-(nStim_ON(current_epochs(1,a))*sampling_rate);
        %part_rasters{a,b} = rasters{a,1};
    end
    new_rasters = [];
    for a=1:length(current_epochs),
        new_rasters = [new_rasters;rasters{a,1}];
    end
    stim_rasters{b,1} = new_rasters;
end


%%
    %*******************************************************************
    %*******************************************************************
    %Run FFT on median filtered subthreshold responses.  Assess power in
    %low freq. and gamma bands
    %*******************************************************************
    %*******************************************************************
    %x_downsample = 10;     %makes sense for spiking analysis, not
    %necessary for spike-removed Vm
%trial_amp_coeff = cell(length(stimvalues)-1,stimmatch_reps);   ***
trial_amp_coeff = cell(length(stimvalues)-1,length(reps));
%trial_power_coeff = cell(length(stimvalues)-1,stimmatch_reps); ***
trial_power_coeff = cell(length(stimvalues)-1,length(reps));
series_amp_coeff = cell(length(stimvalues)-1,1);
series_power_coeff = cell(length(stimvalues)-1,1);
Vm_ext_array = cell2mat(Vm_trial_array');
for j2 = 1:length(stimvalues)-1,
    %for j3 = 1:stimmatch_reps,                 ***
    for j3 = 1:length(reps),
        %f_Vm = [Vm_trial_array{j2,j3}(1,1);...
        %    Vm_trial_array{j2,j3}(x_downsample:(end-(mod(length(Vm_trial_array{j2,j3}),x_downsample))),1)];
        series_f_Vm = Vm_ext_array(1:end,j2);
        f_Vm = Vm_trial_array{j2,j3}(1:end,1);
        ft_times = (1:length(f_Vm)).*sample_interval;
        ext_ft_times = (1:length(series_f_Vm)).*sample_interval;
        if use_windowing == 1,
            w_series = f_Vm.*hanning(length(f_Vm));
            %alt.: w_series = f_Vm.*blackman(length(f_Vm));
            ext_w_series = series_f_Vm.*hanning(length(series_f_Vm));
        else
            w_series = f_Vm;
            ext_w_series = series_f_Vm
        end
        fc = fft(w_series,length(ft_times));
        series_fc = fft(ext_w_series,length(ext_ft_times));
            
        %*********method below for use with spikes only***********
        %spike_counts = spiketimes2bins(w_series,ft_times);
        %[fc,freqs] = fouriercoeffs(spike_counts,sample_interval);
        %stimtrain_f{j2,j3} = freqs;
            
        trial_amp_coeff{j2,j3} = ...
            reshape(sqrt(fc.*conj(fc)),length(fc),1);
        trial_power_coeff{j2,j3} = ...
            reshape((fc.*conj(fc)),length(fc),1);
        freqs = sampling_rate*(0:length(ft_times)-1)/length(ft_times);
        series_amp_coeff{j2,1} = ...
            reshape((series_fc.*conj(series_fc)),length(series_fc),1);
        series_power_coeff{j2,1} = ...
            reshape(series_fc.*conj(series_fc),length(series_fc),1);
        ext_freq = ...
            sampling_rate*(0:length(ext_ft_times)-1)/length(ext_ft_times);
    end
end
ave_amp_coeff = cell(length(stimvalues)-1,1);
ave_power_coeff = cell(length(stimvalues)-1,1);
new_amp_coeff = cell2mat(reshape(trial_amp_coeff,1,size(trial_amp_coeff,2)*...
    size(trial_amp_coeff,1)));
for j4 = 1:length(stimvalues)-1,
    sum_amp_coeff = [];
    sum_power_coeff = [];
    %for j5 = 1:stimmatch_reps,             ***
    for j5 = 1:length(reps),
        sum_amp_coeff = [sum_amp_coeff,trial_amp_coeff{j4,j5}];
        sum_power_coeff = [sum_power_coeff,trial_power_coeff{j4,j5}];
    end
    %ave_amp_coeff{j4,1} = sum(sum_amp_coeff,2)./stimmatch_reps;    ***
    ave_amp_coeff{j4,1} = sum(sum_amp_coeff,2)./length(reps);
    %ave_power_coeff{j4,1} = sum(sum_power_coeff,2)./stimmatch_reps; ***
    ave_power_coeff{j4,1} = sum(sum_power_coeff,2)./length(reps);
end
    
    %plot results
    for j4 = 1:length(stimvalues)-1,
        fas = figure;
        h = plot(freqs,ave_amp_coeff{j4,1},'k');
        h.LineWidth = 2.5;
        xlabel('Frequency (Hz)');
        ylabel('Spectral power');
        xlim([-10 80]);
        title(['Amplitude spectrum of Vm traces,'...
            'Direction: ',num2str(stimvalues(1,j4))]);
        if save_it == 1,
            saveas(gcf,[pwd filesep 'FFTVm_angle',num2str(stimvalues(1,j4))]);
            close(fas);
        else
        end
    end
    
    %%
    %*********************************************************************
    %**********Compute autocorrelation of spike reponse,******************
    %**********then FFT of result -- for comparison with Vm spectrum******
    %*********************************************************************
autocorr_output = {};
for j7 = 1:length(stimvalues)-1,
    autocorr_output{j7,1} = struct('trial_corr',[],'trial_lags',[]);
end
modulation_band = [0 100];
window_size = stim_duration(1,1); %note- this does not accommodate blocks with variable cycle durations
bin_size = 0.001;
for j6 = 1:length(stimvalues)-1,
    if isempty(stim_rasters{j6,1}),
        trial_spike_times = NaN;
        continue;
    else
        trial_spike_times = stim_rasters{j6,1}/sampling_rate;
    end
    ach_bin_edges = [-(window_size):bin_size:window_size];
    N_ach = histc(trial_spike_times,ach_bin_edges);
    ach_bin_centers = (ach_bin_edges(1:end-1)+ach_bin_edges(2:end))/2;
    N_ach = N_ach(1:end-1);
    [corr,lags] = xcorr(N_ach,N_ach,round(window_size/bin_size));
    lags = lags*bin_size;
    autocorr_output{j6,1}.trial_corr = corr;
    autocorr_output{j6,1}.trial_lags = lags;
end
    
for j8 = 1:length(stimvalues)-1,
    f = figure;
    ach_h = bar(autocorr_output{j8,1}.trial_lags,...
        autocorr_output{j8,1}.trial_corr,'k');
    hold on;
    xlabel('Lags');
    ylabel('Count');
    xlim([-(stim_duration(1,1)) stim_duration(1,1)]);
    title(['Trial Spike-time autocorrelation, Direction angle ',...
        num2str(stimvalues(1,j8))]);
    if save_it == 1,
        saveas(gcf,[pwd filesep 'trial_spiketrain_autocorr_angle',...
            num2str(stimvalues(1,j8)),'.fig']);
        close(f);
    else
    end
end

%*****consider adding partial autocorrelation here*****
  
s_trial_amp_coeff = cell(length(stimvalues)-1,1);
s_trial_power_coeff = cell(length(stimvalues)-1,1);
    %now FFT 
for j9 = 1:length(stimvalues)-1,
    if isempty(autocorr_output{j9,1}.trial_corr)
        continue;
    else
        aft_times = (1:length(autocorr_output{j9,1}.trial_corr)).*...
            sample_interval;
    end
    if use_windowing == 1,
        series_ = reshape(autocorr_output{j9,1}.trial_corr,...
            length(autocorr_output{j9,1}.trial_corr),1);
        aw_series = ...
            series_.*hanning(length(autocorr_output{j9,1}.trial_corr));
    else
        aw_series = autocorr_output{j9,1}.trial_corr;
    end
    afc = fft(aw_series,length(aft_times));
    s_trial_amp_coeff{j9,1} = reshape(sqrt(afc.*conj(afc)),length(afc),1);
    s_trial_power_coeff{j9,1} = reshape((afc.*conj(afc)),length(afc),1);
    s_freqs = sampling_rate*(0:length(aft_times)-1)/length(aft_times);
end
    %sp_ave_amp_coeff = cell(length(stimvalues)-1,1);
    %sp_ave_power_coeff = cell(length(stimvalues)-1,1);
    %for j4=1:length(stimvalues)-1,
    %    sp_sum_amp_coeff = [];
    %    sp_sum_power_coeff = [];
    %    for j5=1:stimmatch_reps,
    %        sp_sum_amp_coeff = [sp_sum_amp_coeff;s_trial_amp_coeff{j4,j5}];
    %        sp_sum_power_coeff = [sp_sum_power_coeff;s_trial_power_coeff{j4,j5}];
    %    end
    %    sp_ave_amp_coeff{j4,1} = sum(sp_sum_amp_coeff,2)./stimmatch_reps;
    %    sp_ave_power_coeff{j4,1} = sum(sp_sum_power_coeff,2)./stimmatch_reps;
    %end
    %not necessary; averaged over autocorr step
        
    
    %plot spiking power spectra
for j9 = 1:length(stimvalues)-1,
    if isempty(s_trial_power_coeff{j9,1}),
        continue
    else
    end
    fps = figure;
    h = plot(s_freqs,s_trial_power_coeff{j9,1},'k');
    h.LineWidth = 2.5;
    xlabel('Frequency (Hz)');
    ylabel('Spectral power');
    xlim([-10 100]);
    title(['Power spectrum of spike response autocorrelation,'...
        'Direction: ',num2str(stimvalues(1,j9))]);
    if save_it == 1,
        saveas(gcf,[pwd filesep 'FFTspike_angle',...
            num2str(stimvalues(1,j9))]);
        close(fps);
    else
    end
end 


%%
        %*******Assess biophysical spike threshold using Azouz & Gray, 1999
        %*******methods****************************************************
        %******************************************************************
        
        %********SECTION RELOCATED ABOVE*********
        
       
        %***Create cell summary barplots of Spontaneous versus all 
        %***Stimulus-ON spike thresholds
        stim_th = [];
        spont_th = [];
        opening_set = find(threshold<(nStim_ON(1,1)/sample_interval));
        spont_locs = [];
        for qq = 1:length(nStim_ON),
            if ~isnan(nStim_OFF(qq,1)),
                stim_set = find((threshold>...
                    (nStim_ON(qq,1)/sample_interval))&...
                    (threshold<(nStim_OFF(qq,1)/sample_interval)));
                stim_th = [stim_th;Vm_th(stim_set,1)];
                if qq < length(nStim_ON),
                    spont_set = find((threshold>...
                        (nStim_OFF(qq,1)/sample_interval))&...
                        (threshold<(nStim_ON(qq+1,1)/sample_interval)));
                else
                    spont_set = find((threshold>...
                        (nStim_OFF(qq,1)/sample_interval)));
                end
                spont_th = [spont_th;Vm_th(spont_set,1)];
                spont_locs = [spont_locs;threshold(spont_set,1)];
            else
            end
        end
        
        %(NOTE - this is more of a back-up plot, hence the lack of save
        %command - see figure 6 for relevant data.
        f5 = figure;
        spont_x = 1.0+(rand(length(spont_th),1)-0.5);
        stim_x = 3.0+(rand(length(stim_th),1)-0.5);
        scatter(spont_x,spont_th,'filled');
        hold on;
        %sel = randi(length(spont_th),length(spont_th),1);
        scatter(stim_x,stim_th,'filled');
        %scatter(stim_x(sel,1),stim_th(sel,1),'filled');
        xlim([0 4]);
        ylim([-100 0]);
        m_spont = mean(spont_th);
        m_stim = mean(stim_th);
        sem_spont = std(spont_th)/sqrt(length(spont_th));
        sem_stim = std(stim_th)/sqrt(length(stim_th));
        errorbar(1.0,m_spont,sem_spont);
        bar(1.0,m_spont);
        errorbar(3.0,m_stim,sem_stim);
        bar(3.0,m_stim);
        ax = gca;
        ax.XTick = ([1 3]);
        ax.XTickLabel = {'Spont.','Stimulus'};
        ylabel('Membrane potential (mV)');
        close(f5);
        
        
        %follwing two blocks are transplanted from the VF-plot section; 
        %note pref/null combination indices are needed for the following section
        if ~exist('pref_stim','var')>0,
            try
                pref_stim = otpref_total;
            catch
                pref_stim = otpref_1;
            end
             null_stim = pref_stim+180;
             if null_stim>=360,
                 null_stim = null_stim-360;
             else
             end
         else
         end
         %***NOTE: POTENTIAL ERROR: some cells have different 
         %***variable names saved for pref_stim and null_stim.  
         %***This is an attempt to correct that, but may not catch all
         %cases.
         pref_ind = find(stimvalues==pref_stim);
         null_ind = find(stimvalues==null_stim);
         all_dirind = 1:length(stimvalues)-1;
         down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
         up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
         down_i = circshift(all_dirind,1,2);
         up_i = circshift(all_dirind,-1,2);
         pref_set = ...
             [down_(pref_ind);pref_stim;up_(pref_ind)];
         null_set = ...
             [down_(null_ind);null_stim;up_(null_ind)];
         pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
         null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
         
        
        %Partition thresholds by 3 categories:  Spont., Pref, Null
        %*****IMPORTANT: RESULTING ARRAY ORDER ALWAYS PROCEEDS IN THE FOLLOWING
        %*****ROW SEQUENCE - (THETA-1)*NUM REPS , THETA*NUM REPS , (THETA+1)*NUM REPS 
        %*****(WHERE THETA = PEAK PREF. DIRECTION) **********************
        pref_th = [];
        null_th = [];
        pref_set = {};
        null_set = {};
        for ii = 1:length(pref_set_index),
            pstim_locs = find(stimorder==pref_set_index(ii,1));
            for iii = 1:length(pstim_locs),
                pref_set{end+1,1} = find(threshold>...
                    (nStim_ON(pstim_locs(1,iii),1)/sample_interval)&...
                    (threshold<(nStim_OFF(pstim_locs(1,iii),1)/...
                    sample_interval)));
            end
            pref_th = [pref_th;Vm_th(cell2mat(pref_set),1)];
        end
        for ii = 1:length(null_set_index),
            nstim_locs = find(stimorder==null_set_index(ii,1));
            for iii = 1:length(nstim_locs),
                null_set{end+1,1} = find(threshold>...
                    (nStim_ON(nstim_locs(1,iii),1)/sample_interval)&...
                    (threshold<(nStim_OFF(nstim_locs(1,iii),1)/...
                    sample_interval)));
            end
            null_th = [null_th;Vm_th(cell2mat(null_set),1)];
        end
        
        f6 = figure;
        pref_x = 3.0+(rand(length(pref_th),1)-0.5);
        null_x = 5.0+(rand(length(null_th),1)-0.5);
        s = scatter(spont_x,spont_th,'filled');
        hold on;
        sp = scatter(pref_x,pref_th,'filled');
        sn = scatter(null_x,null_th,'filled');
        s.MarkerFaceColor = [0.2,0.2,0.2];
        sp.MarkerFaceColor = [0.8, 0, 0.5];
        sn.MarkerFaceColor = [0, 0, 0.9];
        xlim([0 6]);
        ylim([-100 0]);
        m_pref = mean(pref_th);
        m_null = mean(null_th);
        sd_pref = std(pref_th);
        sd_null = std(null_th);
        sem_pref = std(pref_th)/sqrt(length(pref_th));
        sem_null = std(null_th)/sqrt(length(null_th));
        errorbar(1.0,m_spont,sem_spont);
        bar(1.0,m_spont);
        errorbar(3.0,m_pref,sem_pref);
        bar(3.0,m_pref);
        errorbar(5.0,m_null,sem_null);
        bar(5.0,m_null);
        ax = gca;
        ax.XTick = ([1 3 5]);
        ax.XTickLabel = {'Spont.','Pref','Null'};
        ylabel('Membrane potential (mV)');
        title('Threshold potential by Input/Stim. type');
        if save_it == 1,
            saveas(gcf,[pwd filesep 'threshold_Comp_',char(directories),'.fig']);
            close(f6);
        else
        end
        
        %******** Now analyze relationship between Vrest and Vth - for every
        %interval between two spikes choose the subthreshold region that is more
        %than 4 ms after the previous and before the current threshold and take
        %the mean of that epoch.
        %
        %PREF SET
        %
        p_subth_mean = [];
        p_matched_th = [];
        p_el_record = [];
        p_subth_std = []; 
        p_matched_isi = {};
        matched_slope = [];
        for a = 1:length(pref_set),
            p_matched_isi{end+1,1} = diff(threshold(pref_set{a,1})).*...
                sample_interval;
            for b = 2:length(pref_set{a,1}),
                p_epoch_length = (threshold(pref_set{a,1}(b,1))-...
                    round(0.004/sample_interval))-...
                    (threshold(pref_set{a,1}(b-1,1))+round(0.004/sample_interval));
                slope_range = (threshold(pref_set{a,1}(b,1))-...
                    round(0.005/sample_interval)):threshold(pref_set{a,1}(b,1));
                %pf = polyfit(sample_interval*(1:length(slope_range))',Vm_waveform(slope_range)./1000,1);
                d_vec = first_deriv((Vm_waveform(slope_range)./...
                    1000),sample_interval);
                %matched_slope(end+1) = pf(1,1);
                matched_slope(end+1) = mean(d_vec);
                if p_epoch_length <= 0,
                    p_subth_mean(end+1,1) = NaN;
                    p_subth_std(end+1,1) = NaN;
                    p_matched_th(end+1,1) = NaN;
                else
                    p_subth_epoch = [(threshold(pref_set{a,1}(b-1,1))+...
                        round(0.004/sample_interval)):1:...
                        (threshold(pref_set{a,1}(b,1))-...
                        round(0.004/sample_interval))];
                    p_subth_mean(end+1,1) = mean(Vm_waveform(p_subth_epoch));
                    p_subth_std(end+1,1) = std(Vm_waveform(p_subth_epoch));
                    p_matched_th(end+1,1) = Vm_th(pref_set{a,1}(b,1),1);
                end
                p_el_record(end+1) = p_epoch_length;
            end
        end
        p_allspike_AHP_amp = [];
        p_trialspike_AHP_amp = cell(length(pref_set),1);
        for a1 = 1:length(pref_set),
            p_trial_AHP_amps = [];
            for b1 = 1:length(pref_set{a1,1})-1,  %***leaving out the last spike in trial
                p_th_epoch = [threshold(pref_set{a1,1}(b1,1)):1:...
                    threshold(pref_set{a1,1}(b1+1,1))];
                p_trialspike_AHP_amp{a1,1}(end+1,1) = ...
                    Vm_waveform(p_th_epoch(1,1))-min(Vm_waveform(p_th_epoch));
                p_trial_AHP_amps(end+1,1) = ...
                    Vm_waveform(p_th_epoch(1,1))-min(Vm_waveform(p_th_epoch));
                p_allspike_AHP_amp(end+1,1) = ...
                    Vm_waveform(p_th_epoch(1,1))-min(Vm_waveform(p_th_epoch));
            end
            p_trialmean_AHP_amp(a1,1) = mean(p_trial_AHP_amps);
            p_trialsd_AHP_amp(a1,1) = std(p_trial_AHP_amps);
        end
        mean_p_AHP = mean(p_allspike_AHP_amp);
        %trial_mean_p_AHP = 
        mean_p_rest = nanmean(p_subth_mean);
        p_jump = abs(p_matched_th-p_subth_mean);
        fpr = figure;
        scatter(p_subth_mean,p_matched_th,'k');
        hold on;
        xlabel('Pre-spike resting potential (mV)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN VREST AND VTH (PREF. stim. responses only)');
        %ylim([-100 -50]);
        %xlim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'Vrest_Vth_PREF.fig']);
            close(fpr);
        else
        end
        fpi = figure;
        scatter(p_subth_mean,p_el_record.*sample_interval,'k');
        hold on;
        xlabel('Pre-spike resting potential (mV)');
        ylabel('Subthreshold epoch duration (seconds)');
        title('RELATIONSHIP BETWEEN VREST AND EFFECTIVE ISI (PREF. STIMS.)');
        %xlim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vrest.fig']);
            close(fpi);
        else
        end
        fpa = figure;
        scatter(p_el_record.*sample_interval,p_matched_th,'k');
        hold on;
        xlabel('Subthreshold epoch duration (seconds)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN EFFECTIVE ISI AND FOLLOWING VTH (PREF. STIMS.)');
        %ylim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vth_PREF.fig']);
            close(fpa);
        else
        end
        fp3 = figure;
        scatter3(p_matched_th,p_subth_mean,p_el_record.*sample_interval,'k');
        hold on;
        xlabel('Following spike threshold (mV)');
        ylabel('Pre-spike resting potential (mV)');
        zlabel('Effective ISI (seconds)');
        title('PREF');
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vrest_Vth_PREFcombined.fig']);
            close(fp3);
        else
        end
        fpf = figure;
        scatter(p_subth_std,p_matched_th,'k');
        hold on;
        xlabel('Subthreshold standard deviation (mV)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN PRE-SPIKE FLUCTUATION AND THRESHOLD (PREF. STIMS.)');
        %ylim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'VmSD_Vth_PREF.fig']);
            close(fpf);
        else
        end
        fsl = figure;
        scatter(matched_slope,p_matched_th,'k');
        hold on;
        xlabel('Vm slope (mV/ms)');
        ylabel('Threshold (mV)');
        title('Vm slope 5ms before spike (PREF. STIMS.)');
        %ylim([-90 -20]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'Vmslope_Vth_PREF.fig']);
            close(fsl);
        else
        end
        %
        %NULL SET
        %
        n_subth_mean = [];
        n_matched_th = [];
        n_el_record = [];
        n_subth_std = [];
        n_matched_isi = [];
        nmatched_slope = [];
        for aa = 1:length(null_set),
            n_matched_isi{end+1,1} = ...
                diff(threshold(null_set{aa,1})).*sample_interval;
            for bb = 2:length(null_set{aa,1}),
                n_epoch_length = ...
                    (threshold(null_set{aa,1}(bb,1))-round(0.004/sample_interval))-...
                    (threshold(null_set{aa,1}(bb-1,1))+round(0.004/sample_interval));
                nslope_range = ...
                    (threshold(null_set{aa,1}(bb,1))-...
                    round(0.005/sample_interval)):threshold(null_set{aa,1}(bb,1));
                %pfn = polyfit(sample_interval*(1:length(nslope_range))',Vm_waveform(nslope_range)./1000,1);
                d_vec = ...
                    first_deriv((Vm_waveform(nslope_range)./1000),sample_interval);
                %nmatched_slope(end+1) = pfn(1,1);
                nmatched_slope(end+1) = mean(d_vec);
                if n_epoch_length <= 0,
                    n_subth_mean(end+1,1) = NaN;
                    n_matched_th(end+1,1) = NaN;
                    n_subth_std(end+1,1) = NaN;
                else
                    n_subth_epoch = [(threshold(null_set{aa,1}(bb-1,1))+...
                        round(0.004/sample_interval)):1:...
                        (threshold(null_set{aa,1}(bb,1))-...
                        round(0.004/sample_interval))];
                    n_subth_mean(end+1,1) = mean(Vm_waveform(n_subth_epoch));
                    n_subth_std(end+1,1) = std(Vm_waveform(n_subth_epoch));
                    n_matched_th(end+1,1) = Vm_th(null_set{aa,1}(bb,1),1);
                end
                n_el_record(end+1) = n_epoch_length;
            end
        end
        n_allspike_AHP_amp = [];
        n_trialspike_AHP_amp = cell(length(null_set),1);
        for a2 = 1:length(null_set),
            n_trial_AHP_amps = [];
            for b2 = 1:length(null_set{a2,1})-1,  %***leaving out the last spike in trial
                n_th_epoch = ...
                    [threshold(null_set{a2,1}(b2,1)):1:...
                    threshold(null_set{a2,1}(b2+1,1))];
                n_trialspike_AHP_amp{a2,1}(end+1,1) = ...
                    Vm_waveform(n_th_epoch(1,1))-min(Vm_waveform(n_th_epoch));
                n_trial_AHP_amps(end+1,1) = ...
                    Vm_waveform(n_th_epoch(1,1))-min(Vm_waveform(n_th_epoch));
                n_allspike_AHP_amp(end+1,1) = ...
                    Vm_waveform(n_th_epoch(1,1))-min(Vm_waveform(n_th_epoch));
            end
            n_trialmean_AHP_amp(a2,1) = mean(n_trial_AHP_amps);
            n_trialsd_AHP_amp(a2,1) = std(n_trial_AHP_amps);
        end
        mean_n_AHP = mean(n_allspike_AHP_amp);
        mean_n_rest = nanmean(n_subth_mean);
        n_jump = abs(n_matched_th-n_subth_mean);
        fnr = figure;
        scatter(n_subth_mean,n_matched_th,'k');
        hold on;
        xlabel('Pre-spike resting potential (mV)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN VREST AND VTH (NULL. stim. responses only)');
        %ylim([-100 -50]);
        %xlim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'Vrest_Vth_NULL.fig']);
            close(fnr);
        else
        end
        fni = figure;
        scatter(n_subth_mean,n_el_record.*sample_interval,'k');
        hold on;
        xlabel('Pre-spike resting potential (mV)');
        ylabel('Subthreshold epoch duration (seconds)');
        title('RELATIONSHIP BETWEEN VREST AND EFFECTIVE ISI (NULL STIMS.)');
        %xlim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vrest_NULL.fig']);
            close(fni);
        else
        end
        fna = figure;
        scatter(n_el_record.*sample_interval,n_matched_th,'k');
        hold on;
        xlabel('Subthreshold epoch duration (seconds)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN EFFECTIVE ISI AND FOLLOWING VTH (NULL STIMS.)');
        %ylim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vth_NULL.fig']);
            close(fna);
        else
        end
        fn3 = figure;
        scatter3(n_matched_th,n_subth_mean,n_el_record.*sample_interval,'k');
        hold on;
        xlabel('Following spike threshold (mV)');
        ylabel('Pre-spike resting potential (mV)');
        zlabel('Effective ISI (seconds)');
        title('NULL');
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vrest_Vth_NULLcombined.fig']);
            close(fn3);
        else
        end
        fnf = figure;
        scatter(n_subth_std,n_matched_th,'k');
        hold on;
        xlabel('Subthreshold standard deviation (mV)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN PRE-SPIKE FLUCTUATION AND THRESHOLD (NULL STIMS.)');
        %ylim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'VmSD_Vth_NULL.fig']);
            close(fnf);
        else
        end
        fsl = figure;
        scatter(nmatched_slope,n_matched_th,'k');
        hold on;
        xlabel('Vm slope (mV/ms)');
        ylabel('Threshold (mV)');
        title('Vm slope 5ms before spike (NULL STIMS.)');
        %ylim([-90 -20]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'Vmslope_Vth_NULL.fig']);
            close(fsl);
        else
        end
        %
        %SPONT. Set
        %
        sp_subth_mean = [];
        sp_matched_th = [];
        sel_record = diff(spont_locs);
        sp_subth_std = [];
        sp_allspike_AHP_amp = [];
        for c=2:length(spont_locs),
            if sel_record(c-1,1) > round(0.2/sample_interval),
                sp_subth_epoch = [(spont_locs(c,1)-...
                    round(0.1/sample_interval)):1:...
                    (spont_locs(c,1)-round(0.004/sample_interval))];
                sp_subth_mean(end+1,1) = mean(Vm_waveform(sp_subth_epoch));
                sp_subth_std(end+1,1) = std(Vm_waveform(sp_subth_epoch));
                sp_matched_th(end+1,1) = Vm_waveform(spont_locs(c,1),1);
                sp_allspike_AHP_amp(end+1,1) = ...
                    Vm_waveform(spont_locs(c-1,1),1)-...
                    min(Vm_waveform(sp_subth_epoch));
            else
                sp_epoch_length = ...
                    (spont_locs(c,1)-round(0.004/sample_interval))-...
                    (spont_locs(c-1,1)+round(0.004/sample_interval));
                if sp_epoch_length <= 0,
                    sp_subth_mean(end+1,1) = NaN;
                    sp_matched_th(end+1,1) = NaN;
                    sp_subth_std(end+1,1) = NaN;
                    sp_allspike_AHP_amp(end+1,1) = NaN;
                else
                    sp_subth_epoch = ...
                        [(spont_locs(c-1,1)+...
                        round(0.004/sample_interval)):1:...
                        (spont_locs(c,1)-round(0.004/sample_interval))];
                    sp_subth_mean(end+1,1) = mean(Vm_waveform(sp_subth_epoch));
                    sp_subth_std(end+1,1) = std(Vm_waveform(sp_subth_epoch));
                    sp_matched_th(end+1,1) = Vm_waveform(spont_locs(c,1),1);
                    if c <= length(spont_locs),
                        sp_allspike_AHP_amp(end+1,1) = ...
                            Vm_waveform(spont_locs(c-1,1),1)-...
                            min(Vm_waveform(spont_locs(c-1,1):...
                            (spont_locs(c-1,1)+round(0.2/sample_interval)),1));
                    else
                        sp_allspike_AHP_amp(end+1,1) = ...
                            Vm_waveform(spont_locs(c-1,1),1)-...
                            min(Vm_waveform(spont_locs(c-1,1):length(Vm_waveform)));
                    end
                end
            end 
        end
        mean_sp_rest = nanmean(sp_subth_mean);
        sp_jump = abs(sp_matched_th-sp_subth_mean);
        fr = figure;
        scatter(sp_subth_mean,sp_matched_th,'k');
        hold on;
        xlabel('Pre-spike resting potential (mV)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN VREST AND VTH (SPONT. SPIKING only)');
        %ylim([-100 -50]);
        %xlim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'Vrest_Vth_SPONT.fig']);
            close(fr);
        else
        end
        fi = figure;
        scatter(sp_subth_mean,sel_record.*sample_interval,'k');
        hold on;
        xlabel('Pre-spike resting potential (mV)');
        ylabel('Subthreshold epoch duration (seconds)');
        title('RELATIONSHIP BETWEEN VREST AND EFFECTIVE ISI (SPONT.)');
        %xlim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vrest_SPONT.fig']);
            close(fi);
        else
        end
        fa = figure;
        scatter(sel_record.*sample_interval,sp_matched_th,'k');
        hold on;
        xlabel('Subthreshold epoch duration (seconds)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN EFFECTIVE ISI AND FOLLOWING VTH (SPONT. SPIKES)');
        %ylim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vth_SPONT.fig']);
            close(fa);
        else
        end
        fsp3 = figure;
        scatter3(sp_matched_th,sp_subth_mean,sel_record.*sample_interval,'k');
        hold on;
        xlabel('Following spike threshold (mV)');
        ylabel('Pre-spike resting potential (mV)');
        zlabel('Effective ISI (seconds)');
        title('SPONTANEOUS');
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'ISI_Vrest_Vth_SPONTcombined.fig']);
            close(fsp3);
        else
        end
        fnf = figure;
        scatter(sp_subth_std,sp_matched_th,'k');
        hold on;
        xlabel('Subthreshold standard deviation (mV)');
        ylabel('Following spike threshold (mV)');
        title('RELATIONSHIP BETWEEN PRE-SPIKE FLUCTUATION AND THRESHOLD (SPONT. SPIKES)');
        %ylim([-100 -50]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'VmSD_Vth_SPONT.fig']);
            close(fnf);
        else
        end
        
        
        
%%
    %***************************************************************
    %*****Calculate all spike shape parameters
    %***************************************************************
    
    %*****************FOR PREF. STIM.:     *************************
    pref_spike_params = struct('new_window',{},'pref_dVdt',{},...
        'int_pref_dVdt',{},'max_dVdt',{},...
        'peak_Vm',{},'spike_amplitude',{},...
        'th_spike_width',{},'spike_halfwidth',{});
    for i2 = 1:length(pref_set),
        for j = 1:length(pref_set{i2,1}),
            re_half_win = round(0.010/sample_interval);
            re_pref_triggered_trace = ...
                waveform_dnlp(threshold(pref_set{i2,1}(j,1))-...
                re_half_win:threshold(pref_set{i2,1}(j,1),1)+re_half_win,1);
            re_windowed_t = (-(re_half_win):re_half_win).*sample_interval;
            pref_spike_params{i2,j}.pref_dVdt = ...
                gradient(re_pref_triggered_trace,re_windowed_t);
            pref_spike_params{i2,j}.int_pref_dVdt = ...
                interp(gradient(re_pref_triggered_trace,re_windowed_t),10);
            pref_spike_params{i2,j}.max_dVdt = ...
                max(gradient(re_pref_triggered_trace,re_windowed_t));
            pref_spike_params{i2,j}.peak_Vm = ...
                max(re_pref_triggered_trace);
            pref_spike_params{i2,j}.spike_amplitude = ...
                max(re_pref_triggered_trace)-...
                (waveform_dnlp(threshold(pref_set{i2,1}(j,1))));
            half_height = abs(max(re_pref_triggered_trace)-...
                (waveform_dnlp(threshold(pref_set{i2,1}(j,1)))))/2;
            %resample spike trace for improved accuracy (at 10x sampling
            %rate)
            resampled_trace = interp(re_pref_triggered_trace,10);
            resampled_windowed_t = interp(re_windowed_t,10);
            r_up_draft = ...
                resampled_trace(find(gradient(resampled_trace,...
                resampled_windowed_t)>0),1);
            r_down_draft = ...
                resampled_trace(find(gradient(resampled_trace,...
                resampled_windowed_t)<=0),1);
            [~,upward_ttime] = ...
                min(abs(r_up_draft-...
                (waveform_dnlp(threshold(pref_set{i2,1}(j,1))))));
            [~,down_ttime] = ...
                min(abs(r_down_draft-...
                (waveform_dnlp(threshold(pref_set{i2,1}(j,1))))));
            pref_spike_params{i2,j}.th_spike_width = ...
                (abs(length(r_up_draft)-...
                upward_ttime)+abs(down_ttime))*(sample_interval/10);
            [~,upward_htime] = ...
                min(abs(r_up_draft-...
                (waveform_dnlp(threshold(pref_set{i2,1}(j,1)))+half_height)));
            [~,down_htime] = ...
                min(abs(r_down_draft-...
                (waveform_dnlp(threshold(pref_set{i2,1}(j,1)))+half_height)));
            pref_spike_params{i2,j}.spike_halfwidth = ...
                (abs(length(r_up_draft)-...
                upward_htime)+abs(down_htime))*(sample_interval/10);
            pref_spike_params{i2,j}.new_window = ...
                re_pref_triggered_trace;
            pref_spike_params{i2,j}.resampled_window = ...
                resampled_trace;
        end
    end
    %make shape parameter phase plots
    %1. from raw traces
    ccode_stims = 1;
    color_set = {[0,0,1],[0,0,0],[1,0,0]};
    c = 0;
    for i3 = 1:length(pref_set),
        %if mod(i3,stimmatch_reps) == 1,            ***
        if mod(i3,length(reps)) == 1,
            c = c+1; 
            fff = figure; 
        else
        end
        for j = 1:length(pref_set{i3,1}),
            phase_port = ...
                plot(pref_spike_params{i3,j}.new_window,...
                pref_spike_params{i3,j}.pref_dVdt);
            if ccode_stims == 1,
                set(phase_port,'Color',color_set{c});
            else
            end
            hold on;
            xlabel('Voltage (mV)');
            ylabel('dV / dt');
            %ylim([-1.5e5 4e5]);
            %xlim([-100 30]);
        end
        if ccode_stims == 1,
            title(['PREF MODE:  orientation ',num2str(stimvalues(pref_set_index(c,1)))]);
        else
        end
        if save_it == 1,
            saveas(gcf,[pwd filesep 'PhasePort_PREF_ori',num2str(stimvalues(pref_set_index(c,1))),'.fig']);
            %close(fff);
        else
        end
        hold off;
    end
    %2. from interpolated traces (resampled at 10x)
    c = 0;
    for i4 = 1:length(pref_set),
        %if mod(i4,stimmatch_reps) == 1,        ***
        if mod(i4,length(reps)) == 1,
            c = c+1; 
            fff = figure; 
        else
        end
        for j = 1:length(pref_set{i4,1}),
            phase_port = ...
                plot(pref_spike_params{i4,j}.resampled_window,...
                pref_spike_params{i4,j}.int_pref_dVdt);
            if ccode_stims == 1,
                set(phase_port,'Color',color_set{c});
            else
            end
            hold on;
            xlabel('Voltage (mV)');
            ylabel('dV / dt');
            %ylim([-2.0e5 4e5]);
            %xlim([-100 50]);
        end
        if ccode_stims == 1,
            title(['PREF MODE (interpolated): orientation ',num2str(stimvalues(pref_set_index(c,1)))]);
        else
        end
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'InterpPhasePort_PREF_ori',num2str(stimvalues(pref_set_index(c,1))),'.fig']);
            %close(fff);
        else
        end
    end
    
    %******************FOR NULL STIM.:   ************************
    null_spike_params = struct('new_window',{},'null_dVdt',...
        {},'int_null_dVdt',{},'max_dVdt',{},...
        'peak_Vm',{},'spike_amplitude',{},...
        'th_spike_width',{},'spike_halfwidth',{});
    for i5 = 1:length(null_set),
        for j = 1:length(null_set{i5,1}),
            re_half_win = round(0.010/sample_interval);
            re_null_triggered_trace = ...
                waveform_dnlp(threshold(null_set{i5,1}(j,1))-...
                re_half_win:threshold(null_set{i5,1}(j,1),1)+re_half_win,1);
            re_windowed_t = (-(re_half_win):re_half_win).*sample_interval;
            null_spike_params{i5,j}.null_dVdt = ...
                gradient(re_null_triggered_trace,re_windowed_t);
            null_spike_params{i5,j}.int_null_dVdt = ...
                interp(gradient(re_null_triggered_trace,re_windowed_t),10);
            null_spike_params{i5,j}.max_dVdt = ...
                max(gradient(re_null_triggered_trace,re_windowed_t));
            null_spike_params{i5,j}.peak_Vm = ...
                max(re_null_triggered_trace);
            null_spike_params{i5,j}.spike_amplitude = ...
                max(re_null_triggered_trace)-...
                (waveform_dnlp(threshold(null_set{i5,1}(j,1))));
            half_height = abs(max(re_null_triggered_trace)-...
                (waveform_dnlp(threshold(null_set{i5,1}(j,1)))))/2;
            %resample spike trace for improved accuracy (at 10x sampling
            %rate)
            resampled_trace = interp(re_null_triggered_trace,10);
            resampled_windowed_t = interp(re_windowed_t,10);
            r_up_draft = resampled_trace(find(gradient(resampled_trace,...
                resampled_windowed_t)>0),1);
            r_down_draft = resampled_trace(find(gradient(resampled_trace,...
                resampled_windowed_t)<=0),1);
            [~,upward_ttime] = ...
                min(abs(r_up_draft-...
                (waveform_dnlp(threshold(null_set{i5,1}(j,1))))));
            [~,down_ttime] = ...
                min(abs(r_down_draft-...
                (waveform_dnlp(threshold(null_set{i5,1}(j,1))))));
            null_spike_params{i5,j}.th_spike_width = ...
                (abs(length(r_up_draft)-upward_ttime)+...
                abs(down_ttime))*(sample_interval/10);
            [~,upward_htime] = min(abs(r_up_draft-...
                (waveform_dnlp(threshold(null_set{i5,1}(j,1)))+half_height)));
            [~,down_htime] = min(abs(r_down_draft-...
                (waveform_dnlp(threshold(null_set{i5,1}(j,1)))+half_height)));
            null_spike_params{i5,j}.spike_halfwidth = ...
                (abs(length(r_up_draft)-upward_htime)+...
                abs(down_htime))*(sample_interval/10);
            null_spike_params{i5,j}.new_window = re_null_triggered_trace;
            null_spike_params{i5,j}.resampled_window = resampled_trace;
        end
    end
    %make shape parameter phase plots
    %1. from raw traces
    ccode_stims = 1;
    color_set = {[0,0,1],[0,0,0],[1,0,0]};
    c = 0;
    for i6 = 1:length(null_set),
        %if mod(i6,stimmatch_reps) == 1,        ***
        if mod(i6,length(reps)) == 1,
            c = c+1; 
            fff = figure; 
        else
        end
        for j = 1:length(null_set{i6,1}),
            phase_port = plot(null_spike_params{i6,j}.new_window,...
                null_spike_params{i6,j}.null_dVdt);
            if ccode_stims == 1,
                set(phase_port,'Color',color_set{c});
            else
            end
            hold on;
            xlabel('Voltage (mV)');
            ylabel('dV / dt');
            %ylim([-1.5e5 4e5]);
            %xlim([-100 30]);
        end
        if ccode_stims == 1,
            title(['NULL MODE: orientation ',num2str(stimvalues(null_set_index(c,1)))]);
        else
        end
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'PhasePort_NULL_ori',num2str(stimvalues(null_set_index(c,1))),'.fig']);
            %close(fff);
        else
        end
    end
    %2. from interpolated traces (resampled at 10x)
    c = 0;
    for i7 = 1:length(null_set),
        %if mod(i7,stimmatch_reps) == 1, c = c+1; fhs = figure; else end
        if mod(i7,length(reps)) == 1, c = c+1; fhs = figure; else end
        for j = 1:length(null_set{i7,1}),
            phase_port = plot(null_spike_params{i7,j}.resampled_window,...
                null_spike_params{i7,j}.int_null_dVdt);
            if ccode_stims == 1,
                set(phase_port,'Color',color_set{c});
            else
            end
            hold on;
            xlabel('Voltage (mV)');
            ylabel('dV / dt');
            %ylim([-2.0e5 4e5]);
            %xlim([-100 50]);
        end
        if ccode_stims == 1,
            title(['NULL MODE (interpolated): orientation ',num2str(stimvalues(null_set_index(c,1)))]);
        else
        end
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep 'InterpPhasePort_NULL_ori',num2str(stimvalues(null_set_index(c,1))),'.fig']);
            %close(fhs);
        else
        end
    end
    
    %******************FOR SPONT. ACTIVITY SPIKES**********************
    spont_spike_params = struct('new_window',{},'spont_dVdt',...
        {},'int_spont_dVdt',{},'max_dVdt',{},...
        'peak_Vm',{},'spike_amplitude',{},...
        'th_spike_width',{},'spike_halfwidth',{});
    for i8 = 1:length(spont_locs),
        re_half_win = round(0.010/sample_interval);
        re_spont_triggered_trace = waveform_dnlp(spont_locs(i8,1)-...
                re_half_win:spont_locs(i8,1)+re_half_win,1);
        re_windowed_t = (-(re_half_win):re_half_win).*sample_interval;
        spont_spike_params{i8,1}.spont_dVdt = ...
            gradient(re_spont_triggered_trace,re_windowed_t);
        spont_spike_params{i8,1}.int_spont_dVdt = ...
            interp(gradient(re_spont_triggered_trace,re_windowed_t),10);
        spont_spike_params{i8,1}.max_dVdt = ...
            max(gradient(re_spont_triggered_trace,re_windowed_t));
        spont_spike_params{i8,1}.peak_Vm = max(re_spont_triggered_trace);
        spont_spike_params{i8,1}.spike_amplitude = ...
            max(re_spont_triggered_trace)-(waveform_dnlp(spont_locs(i8,1)));
        half_height = abs(max(re_spont_triggered_trace)-...
            (waveform_dnlp(spont_locs(i8,1))))/2;
        %resample spike trace for improved accuracy (at 10x sampling
        %rate)
        resampled_trace = interp(re_spont_triggered_trace,10);
        resampled_windowed_t = interp(re_windowed_t,10);
        r_up_draft = resampled_trace(find(gradient(resampled_trace,...
            resampled_windowed_t)>0),1);
        r_down_draft = resampled_trace(find(gradient(resampled_trace,...
            resampled_windowed_t)<=0),1);
        [~,upward_ttime] = ...
            min(abs(r_up_draft-(waveform_dnlp(spont_locs(i8,1)))));
        [~,down_ttime] = ...
            min(abs(r_down_draft-(waveform_dnlp(spont_locs(i8,1)))));
        spont_spike_params{i8,1}.th_spike_width = ...
            (abs(length(r_up_draft)-upward_ttime)+...
            abs(down_ttime))*(sample_interval/10);
        [~,upward_htime] = ...
            min(abs(r_up_draft-...
            (waveform_dnlp(spont_locs(i8,1)))+half_height));
        [~,down_htime] = min(abs(r_down_draft-...
            (waveform_dnlp(spont_locs(i8,1)))+half_height));
        spont_spike_params{i8,1}.spike_halfwidth = ...
            (abs(length(r_up_draft)-upward_htime)+...
            abs(down_htime))*(sample_interval/10);
        spont_spike_params{i8,1}.new_window = re_spont_triggered_trace;
        spont_spike_params{i8,1}.resampled_window = resampled_trace;
    end
    %make shape parameter phase plots
    %1. from raw traces    
    fshs = figure; 
    for i9 = 1:length(spont_locs),
        phase_port = plot(spont_spike_params{i9,1}.new_window,...
            spont_spike_params{i9,1}.spont_dVdt,'k');
        hold on;
        xlabel('Voltage (mV)');
        ylabel('dV / dt');
        %ylim([-1.5e5 4e5]);
        %xlim([-100 30]);
        title('Spontaneous spiking');
    end
    hold off;
    if save_it == 1,
        saveas(gcf,[pwd filesep 'PhasePort_SPONT.fig']);
        close(fshs);
    else
    end
    %2. from interpolated traces (resampled at 10x)
    fssh = figure;
    for ii2 = 1:length(spont_locs),
        phase_port = plot(spont_spike_params{ii2,1}.resampled_window,...
            spont_spike_params{ii2,1}.int_spont_dVdt,'k');
        hold on;
        xlabel('Voltage (mV)');
        ylabel('dV / dt');
        %ylim([-2.0e5 4e5]);
        %xlim([-100 50]);
        title('Spontaneous spiking (interpolated)');
    end
    hold off;
    if save_it == 1,
        saveas(gcf,[pwd filesep 'InterpPhasePort_SPONT.fig']);
        close(fssh);
    else
    end
    
    %Summary collected subplots for all sets
    allpref_maxd = {};
    allpref_peakVm = {};
    allpref_amp = {};
    allpref_thwidth = {};
    allpref_halfwidth = {};
    for ii3 = 1:length(pref_set),
        for j = 1:length(pref_set{ii3,1}),
            allpref_maxd{end+1,1} = pref_spike_params{ii3,j}.max_dVdt;
            allpref_peakVm{end+1,1} = pref_spike_params{ii3,j}.peak_Vm;
            allpref_amp{end+1,1} = pref_spike_params{ii3,j}.spike_amplitude;
            allpref_thwidth{end+1,1} = ...
                pref_spike_params{ii3,j}.th_spike_width;
            allpref_halfwidth{end+1,1} = ...
                pref_spike_params{ii3,j}.spike_halfwidth;
        end
    end
    allnull_maxd = {};
    allnull_peakVm = {};
    allnull_amp = {};
    allnull_thwidth = {};
    allnull_halfwidth = {};
    for ii4 = 1:length(null_set),
        for j = 1:length(null_set{ii4,1}),
            allnull_maxd{end+1,1} = null_spike_params{ii4,j}.max_dVdt;
            allnull_peakVm{end+1,1} = null_spike_params{ii4,j}.peak_Vm;
            allnull_amp{end+1,1} = null_spike_params{ii4,j}.spike_amplitude;
            allnull_thwidth{end+1,1} = ...
                null_spike_params{ii4,j}.th_spike_width;
            allnull_halfwidth{end+1,1} = ...
                null_spike_params{ii4,j}.spike_halfwidth;
        end
    end
    allspont_maxd = {};
    allspont_peakVm = {};
    allspont_amp = {};
    allspont_thwidth = {};
    allspont_halfwidth = {};
    for ii5 = 1:length(spont_spike_params),
        allspont_maxd{end+1,1} = spont_spike_params{ii5,1}.max_dVdt;
        allspont_peakVm{end+1,1} = spont_spike_params{ii5,1}.peak_Vm;
        allspont_amp{end+1,1} = spont_spike_params{ii5,1}.spike_amplitude;
        allspont_thwidth{end+1,1} = ...
            spont_spike_params{ii5,1}.th_spike_width;
        allspont_halfwidth{end+1,1} = ...
            spont_spike_params{ii5,1}.spike_halfwidth;
    end
    spont_x = 1.0+(rand(length(spont_spike_params),1)-0.5);
    pref_x = 3.0+(rand(length(allpref_amp),1)-0.5);
    null_x = 5.0+(rand(length(allnull_amp),1)-0.5);
    pref_rex = 3.0+(rand(length(allpref_halfwidth),1)-0.5);
    null_rex = 5.0+(rand(length(allnull_halfwidth),1)-0.5);
    pref_aex = 3.0+(rand(length(p_allspike_AHP_amp),1)-0.5);
    null_aex = 5.0+(rand(length(n_allspike_AHP_amp),1)-0.5);
    ff = figure;
    %generate max-dV/dt subfigure
    subplot(2,3,1);
    scatter(spont_x,cell2mat(allspont_maxd),25,'k','filled');
    hold on;
    scatter(pref_x,cell2mat(allpref_maxd),25,'r','filled');
    scatter(null_x,cell2mat(allnull_maxd),25,'b','filled');
    ax = gca;
    ax.XTick = ([1 3 5]);
    ax.XTickLabel = {'Spont.','Pref.','Null'};
    ylabel('Maximum dV / dt');
    hold off;
    %generate peak Vm subfigure
    subplot(2,3,2);
    scatter(spont_x,cell2mat(allspont_peakVm),25,'k','filled');
    hold on;
    scatter(pref_x,cell2mat(allpref_peakVm),25,'r','filled');
    scatter(null_x,cell2mat(allnull_peakVm),25,'b','filled');
    ax = gca;
    ax.XTick = ([1 3 5]);
    ax.XTickLabel = {'Spont.','Pref.','Null'};
    ylabel('Peak potential (mV)');
    hold off;
    %generate spike amplitude subfigure
    subplot(2,3,3);
    scatter(spont_x,cell2mat(allspont_amp),25,'k','filled');
    hold on;
    scatter(pref_x,cell2mat(allpref_amp),25,'r','filled');
    scatter(null_x,cell2mat(allnull_amp),25,'b','filled');
    ax = gca;
    ax.XTick = ([1 3 5]);
    ax.XTickLabel = {'Spont.','Pref.','Null'};
    ylabel('Spike amplitude (mV)');
    hold off;
    %generate threshold-to-threshold spike width subfigure
    subplot(2,3,4);
    scatter(spont_x,cell2mat(allspont_thwidth),25,'k','filled');
    hold on;
    scatter(pref_rex,cell2mat(allpref_thwidth),25,'r','filled');
    scatter(null_rex,cell2mat(allnull_thwidth),25,'b','filled');
    ax = gca;
    ax.XTick = ([1 3 5]);
    ax.XTickLabel = {'Spont.','Pref.','Null'};
    ylabel('Threshold spike width (seconds)');
    hold off;
    %generate spike half-width subfigure
    subplot(2,3,5);
    scatter(spont_x,cell2mat(allspont_halfwidth),25,'k','filled');
    hold on;
    scatter(pref_rex,cell2mat(allpref_halfwidth),25,'r','filled');
    scatter(null_rex,cell2mat(allnull_halfwidth),25,'b','filled');
    ax = gca;
    ax.XTick = ([1 3 5]);
    ax.XTickLabel = {'Spont.','Pref.','Null'};
    ylabel('Spike half-width (seconds)');
    hold off;
    %generate AHP trough amplitude subfigure
    subplot(2,3,6);
    scatter(spont_x(1:end-1,1),sp_allspike_AHP_amp,25,'k','filled');
    hold on;
    scatter(pref_aex,p_allspike_AHP_amp,25,'r','filled');
    scatter(null_aex,n_allspike_AHP_amp,25,'b','filled');
    ax = gca;
    ax.XTick = ([1 3 5]);
    ax.XTickLabel = {'Spont.','Pref.','Null'};
    ylabel('AHP trough amplitude (mV)');
    hold off;
    if save_it == 1,
        saveas(gcf,[pwd filesep 'Spike_Params_',char(directories),'.fig']);
        close(ff);
    else
    end
    
    
    
    
%%
%*******************Generate Vm-to-FR plot for cell******************
%********************************************************************
if default > 0,
    VF_type = default;
else
    VF_type = input(['Enter ''1'' for conventional V-F plot, ''2''',...
    ' for maximally sampled/instantaneous FR V-F plot: ']);
end
pref_ind = find(stimvalues==pref_stim);
null_ind = find(stimvalues==null_stim);

switch VF_type
    
    case 1
            
        [f,f_temp] = temporal_conv_VF(Vm_array,...
            be_array,stimvalues,sample_interval,psth_array,save_it);
            
    case 2
        
        [FR_inst,opt_w] = rasters2instFR(stimvalues,stimorder,...
            stim_duration,spike_locations,t_vec,sampling_rate,nStim_ON,...
            nStim_OFF,reps_stimmatch,save_it,use_global_filter,filter_type,...
            pref_ind);          %***
        
         %%%(OPTIONAL) Select only adjacent directions at peak and null...
         %for V-F analysis (starting with flat criterion of adjacent...
         %directions only; SEE NOTE BELOW)**************
         
         %***||||||||||||||||||||||||||||||||||||||||||||||||||||||
         %***INCORPORATE CODE SELECTING FLANKING DIRECTIONS HERE***
         %***Will require moving up Fourier transform analysis, OR 
         %***moving down VF plot block.
         %***||||||||||||||||||||||||||||||||||||||||||||||||||||||
         
         try
             stim_tf = p.tFrequency;
         catch
             stim_tf = tempFrequency;
         end
            
         criterion = 6;
         
         if use_flanktest == 1,
            [Vm_ori_th,Vm_ori_accept,sp_ori_th,sp_ori_accept] = ...
                run_flankertest(pref_stim,null_stim,stimvalues,...
                s_trial_power_coeff,series_power_coeff,freq,s_freq,...
                stim_tf,stim_duration,reps_stimmatch,ave_trials,...
                criterion);             %***
         
         
         else
            all_dirind = 1:length(stimvalues)-1;
            down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
            up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
            down_i = circshift(all_dirind,1,2);
            up_i = circshift(all_dirind,-1,2);
            pref_set = [down_(pref_ind);pref_stim;up_(pref_ind)];
            null_set = [down_(null_ind);null_stim;up_(null_ind)];
            pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
            null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
         end
         
         %modify existing Vm traces
         for ii = 1:length(Vm_array),
             new_Vm_array{ii,1} = Vm_array{ii,1}(1:...
                 (length(Vm_array{ii,1})-(round(h_margins*sampling_rate))));
         end
         
         %********************************************************
         %METHOD 1 (incl. subthreshold), UNWEIGHTED
         %***If plan to use (i.e. ever debug) use_flanktest "if"-block above,
         %this block needs to be moved to the "else" sub-block
         array_disp = abs(length(FR_inst{1,1})-length(new_Vm_array{1,1}));
         if length(FR_inst{1,1})==length(new_Vm_array{1,1}),
             coll_pref_Vm = [new_Vm_array{pref_set_index(1,1),1}(1:end-1,1);...
                 new_Vm_array{pref_set_index(2,1),1}(1:end-1,1);...
                 new_Vm_array{pref_set_index(3,1),1}(1:end-1,1)];
             coll_pref_FR = [FR_inst{pref_set_index(1,1),1}(1:end-1,1);...
                 FR_inst{pref_set_index(2,1),1}(1:end-1,1);...
                 FR_inst{pref_set_index(3,1),1}(1:end-1,1)];
             coll_null_Vm = [new_Vm_array{null_set_index(1,1),1}(1:end-1,1);...
                 new_Vm_array{null_set_index(2,1),1}(1:end-1,1);...
                 new_Vm_array{null_set_index(3,1),1}(1:end-1,1)];
             coll_null_FR = [FR_inst{null_set_index(1,1),1}(1:end-1,1);...
                 FR_inst{null_set_index(2,1),1}(1:end-1,1);...
                 FR_inst{null_set_index(3,1),1}(1:end-1,1)];
         else
             for n = 1:length(FR_inst),
                FR_inst{n,1} = reshape(FR_inst{n,1},length(FR_inst{n,1}),1);
             end
             coll_pref_Vm = [new_Vm_array{pref_set_index(1,1),1}(1:end-array_disp,1);...
                 new_Vm_array{pref_set_index(2,1),1}(1:end-array_disp,1);...
                 new_Vm_array{pref_set_index(3,1),1}(1:end-array_disp,1)];
             coll_pref_FR = [FR_inst{pref_set_index(1,1),1};...
                 FR_inst{pref_set_index(2,1),1};...
                 FR_inst{pref_set_index(3,1),1}];
             coll_null_Vm = [new_Vm_array{null_set_index(1,1),1}(1:end-array_disp,1);...
                 new_Vm_array{null_set_index(2,1),1}(1:end-array_disp,1);...
                 new_Vm_array{null_set_index(3,1),1}(1:end-array_disp,1)];
             coll_null_FR = [FR_inst{null_set_index(1,1),1};...
                 FR_inst{null_set_index(2,1),1};...
                 FR_inst{null_set_index(3,1),1}];
         end
         
          %*****************PREF unbinned*******************************
         stimset = 1;
         z_scale = 1;
         match_PN = 1;  %indicates binning for pref and null plots should be matched
         [ f9,ph,bin_edges ] = VF_densitymap(coll_pref_Vm,coll_null_Vm,...
             coll_pref_FR,adapt_bins,stimset,z_scale,match_PN,...
             pref_set_index,null_set_index);
         cb = colorbar;
         ylabel(cb,'log count of time-step occurances');
         if save_it == 1,
             title(['Maximally-sampled V-F plot'...
                 'for collated PREF direction']);
             saveas(gcf,[pwd filesep 'subTH_rawVF_density_PREF.fig']);
             %close(f9);
         else
         end
         if fit_it == 1,
             %if save_it == 1,
             %    f10 = openfig('subTH_rawVF_density_PREF.fig');
             %    hold on;
             %else
             %end
             sub_elem_pcoll = find(coll_pref_FR==0);
             %Vth_est_pcoll = mode(coll_pref_Vm(sub_elem_pcoll));
             Vth_est_pcoll = mean(coll_pref_Vm(sub_elem_pcoll));
             weight = 1;
             [np_x,pm_out,pcl,pgof,pfitinfo] = vf_powerfit(coll_pref_Vm,...
                 coll_pref_FR,anchor_Vth,Vth_est_pcoll,weight);
             p_z_max = max(max(get(ph,'ZData')));
             np_h = line(np_x,pm_out,p_z_max*ones(length(np_x),1));
             np_h.LineWidth = 3.0;
             np_h.Color = [0 0 0];
             title({'PREF Max-sampled V-F plot ',...
                'a*rectify(Vm-Vth)^b fit parameters:  Vth = ',...
                num2str(Vth_est_pcoll),...
                ' a = ', num2str(pcl.a), ' b = ' num2str(pcl.b)});
             mn_line = line([m_pref m_pref],[0 max(pm_out)]);
             %make local regression smoothing comparison
             %[x,inds] = sort(coll_pref_Vm);
             %yy = smooth(x,coll_pref_FR(inds),5000,'lowess');
             %plot(x,yy,'b-');
             %make gaussian fits across y, smooth means in x bins
             x_bins = [];
             bin_ymean = [];
             bin_ysd = [];
             bin_yse = [];
             bin_count = [];
             for i3 = 2:length(bin_edges{1,1}),
                 if i3 < length(bin_edges{1,1}),
                     bin_inds = find(coll_pref_Vm>=bin_edges{1,1}(1,i3-1)&...
                         coll_pref_Vm<bin_edges{1,1}(1,i3));
                 else
                     bin_inds = find(coll_pref_Vm>=bin_edges{1,1}(1,i3-1)&...
                         coll_pref_Vm<=bin_edges{1,1}(1,i3));
                 end
                 bin_count(end+1,1) = length(bin_inds);
                 if ~isempty(bin_inds),
                     bin_ymean(end+1,1) = nanmean(coll_pref_FR(bin_inds,1));
                     bin_ysd(end+1,1) = std(coll_pref_FR(bin_inds,1));
                     bin_yse(end+1,1) = std(coll_pref_FR(bin_inds,1))./sqrt(length(bin_inds));
                     x_bins = [x_bins;i3-1];
                 else
                 end
             end
             bin_centers = (bin_edges{1,1}(1,2:end)+...
                 bin_edges{1,1}(1,1:end-1))/2;
             scatter3(bin_centers(x_bins),bin_ymean,...
                 (abs(p_z_max+1).*ones(length(bin_ymean),1)));
             for i5 = 1:length(bin_ymean),
                 eb = line([bin_centers(x_bins(i5,1)) bin_centers(x_bins(i5,1))],...
                     [bin_ymean(i5,1)-bin_ysd(i5,1) bin_ymean(i5,1)+bin_ysd(i5,1)]);
                 eb.LineWidth = 2.0;
                 eb.Color = 'r';
             end
             legend('density','fit','measured mean Vth',...
                 'X-binned FR means','Location','NorthEastOutside');
             hold off;
          else
          end
          if save_it == 1,
              if filter_type == 2,
                  saveas(gcf,[pwd filesep 'subTH_maxsampled_VF_density_PREF.fig']);
              elseif filter_type == 3,
                  saveas(gcf,[pwd filesep 'boxcar_VF_density_PREF.fig']);
              else
              end
             close(f9);
          else
          end
          %create secondary plot to overlay X-binned pref and null curves
          fel = figure;
          bin_centers = (bin_edges{1,1}(1,2:end)+...
              bin_edges{1,1}(1,1:end-1))/2;
          scatter(bin_centers(x_bins),bin_ymean,'filled','k');
          for i5 = 1:length(bin_ymean),
              eb = line([bin_centers(x_bins(i5,1)) bin_centers(x_bins(i5,1))],...
                  [bin_ymean(i5,1)-bin_yse(i5,1) bin_ymean(i5,1)+bin_yse(i5,1)]);
              eb.LineWidth = 2.0;
              eb.Color = 'r';
          end
          xlabel('Membrane potential (mV)');
          ylabel('Firing rate');
          hold on;
          
          
          
          %plot comparison of FR = 0 and FR > 0 Vm-distributions
          sub_pref_ind = find(coll_pref_FR == 0);
          sub_pref_Vm = coll_pref_Vm(sub_pref_ind);
          prefh = histc(sub_pref_Vm,bin_edges{1,1});
          prefh = prefh(1:end-1);
          fc = figure;
          bar(bin_centers,bin_count);
          grid on;
          grid minor;
          title('Pref response Vm distribution');
          ylabel('Counts');
          xlabel('Potential (mV)');
          hold on;
          bar(bin_centers,prefh,'m');
          legend('All Vm','Subthreshold Vm only',...
              'Location','EastOutside');
          if save_it == 1,
              saveas(gcf,[pwd filesep 'PREF_subVm_vs_VmDist','.fig']);
              close(fc);
          else
          end
          
            %***********************
            %***now NULL unbinned***
            %***********************
          z_scale = 1;
          stimset = null_set_index;
          [ f21,nh,bin_edges_ ] = VF_densitymap(coll_null_Vm,coll_pref_Vm,...
              coll_null_FR,adapt_bins,stimset,z_scale,match_PN,...
              pref_set_index,null_set_index);
          ncb = colorbar;
          ylabel(ncb,'log count of time-step occurances');
          if save_it == 1,
              title(['Maximally-sampled V-F plot'...
                  'for collated NULL direction']);
              saveas(gcf,[pwd filesep 'subTH_rawVF_density_NULL.fig']);
              %close(f21);
          else
          end
          if fit_it == 1,
              %if save_it == 1,
              %    f22 = openfig('subTH_rawVF_density_NULL.fig');
              %    hold on;
              %else 
              %end
              sub_elem_ncoll = find(coll_null_FR==0);
              %Vth_est_ncoll = mode(coll_null_Vm(sub_elem_ncoll));
              Vth_est_ncoll = mean(coll_null_Vm(sub_elem_ncoll));
              weight = 1;
              [nn_x,nm_out,ncl,ngof,nfitinfo] = vf_powerfit(coll_null_Vm,...
                  coll_null_FR,anchor_Vth,Vth_est_ncoll,weight);
              %alt: [nn_x,nm_out,ncl,ngof,nfitinfo] = vf_linfit(coll_null_Vm,coll_null_FR,anchor_Vth,Vth_est_ncoll,weight);
              n_z_max = max(max(get(nh,'ZData')));
              nn_h = line(nn_x,nm_out,n_z_max*ones(length(nn_x),1));
              nn_h.LineWidth = 3.0;
              nn_h.Color = [0 0 0];
              title({'NULL Max-sampled V-F plot ',...
                  'a*rectify(Vm-Vth)^b fit parameters:  Vth ='...
                  num2str(Vth_est_ncoll)...
                  ' a =' num2str(ncl.a) ' b =' num2str(ncl.b)});
              mn_line = line([m_null m_null],[0 200]);
              %make local regression smoothing comparison
              %[xx,inds2] = sort(coll_null_Vm);
              %yyy = smooth(xx,coll_null_FR(inds2),15000,'lowess');
              %plot(xx,yyy,'b-');
              %make gaussian fits across y, smooth means in x bins
              x_bins = [];
              bin_ymean = [];
              bin_ysd = [];
              bin_yse = [];
              bin_count = [];
              for i4 = 2:length(bin_edges{1,1}),
                  if i4 < length(bin_edges{1,1}),
                      bin_inds = find(coll_null_Vm>=bin_edges{1,1}(1,i4-1)&...
                      coll_null_Vm<bin_edges{1,1}(1,i4));
                  else
                      bin_inds = find(coll_null_Vm>=bin_edges{1,1}(1,i4-1)&...
                      coll_null_Vm<=bin_edges{1,1}(1,i4));
                  end
                  bin_count(end+1,1) = length(bin_inds);
                  if ~isempty(bin_inds),
                      bin_ymean(end+1,1) = nanmean(coll_null_FR(bin_inds,1));
                      bin_ysd(end+1,1) = std(coll_null_FR(bin_inds,1));
                      bin_yse(end+1,1) = std(coll_null_FR(bin_inds,1))./sqrt(length(bin_inds));
                      x_bins = [x_bins;i4-1];
                  else
                  end
              end
              bin_centers = (bin_edges{1,1}(1,2:end)+...
              bin_edges{1,1}(1,1:end-1))/2;
              scatter3(bin_centers(x_bins),bin_ymean,...
              (abs(n_z_max+1).*ones(length(bin_ymean),1)));
              for i6 = 1:length(bin_ymean),
                  eb = line([bin_centers(x_bins(i6,1)) bin_centers(x_bins(i6,1))],...
                  [bin_ymean(i6,1)-bin_ysd(i6,1) bin_ymean(i6,1)+bin_ysd(i6,1)]);
                  eb.LineWidth = 2.0;
                  eb.Color = 'b';
              end
            legend('density','fit','measured mean Vth',...
                'X-binned FR means','Location','NorthEastOutside');
            hold off;
        else
        end
        if save_it == 1,
            if filter_type == 2,
                saveas(gcf,[pwd filesep...
                'subTH_maxsampled_VF_density_NULL.fig']);      
            elseif filter_type == 3,
                  saveas(gcf,[pwd filesep 'boxcar_VF_density_PREF.fig']);
            else
            end
            close(f21);
        else
        end
        %complete secondary plot
        figure(fel);
        bin_centers = (bin_edges{1,1}(1,2:end)+...
            bin_edges{1,1}(1,1:end-1))/2;
        scatter(bin_centers(x_bins),bin_ymean,'filled','k');
        for i6 = 1:length(bin_ymean),
            eb = line([bin_centers(x_bins(i6,1)) bin_centers(x_bins(i6,1))],...
                [bin_ymean(i6,1)-bin_yse(i6,1) bin_ymean(i6,1)+bin_yse(i6,1)]);
            eb.LineWidth = 2.0;
            eb.Color = 'b';
        end
        title('X-binned Firing rate means and errors, PREF and NULL');
        hold off;
        if save_it == 1,
            if filter_type == 2,
                saveas(gcf,[pwd filesep 'Xbin_VF_overlay.fig']);
            elseif filter_type == 3,
                saveas(gcf,[pwd filesep 'boxcar_Xbin_VFoverlay.fig']);
            else
            end
            close(fel);
        else
        end
        
        
         %plot comparison of FR = 0 and FR > 0 Vm-distributions
        sub_null_ind = find(coll_null_FR == 0);
        sub_null_Vm = coll_null_Vm(sub_null_ind);
        nullh = histc(sub_null_Vm,bin_edges{1,1});
        nullh = nullh(1:end-1);
        ffn = figure;
        bar(bin_centers,bin_count);
        grid on;
        grid minor;
        title('Null response Vm distribution');
        ylabel('Counts');
        xlabel('Potential (mV)');
        hold on;
        bar(bin_centers,nullh,'c');
        legend('All Vm','Subthreshold Vm only','Location','NorthEast');
        if save_it == 1,
            saveas(gcf,[pwd filesep 'PREF_subVm_vs_VmDist','.fig']);
            close(ffn);
        else
        end
        
end




%%
%********************
%assess F vs Vm noise
%Method 1: FR at last timestep in 5-ms bins matched to Vm SD over same time bin
%****
%pref
for iii = 1:round(length(coll_pref_FR)/floor(0.005*sampling_rate)),
    p_bin_edges(iii,1) = (floor(0.005*sampling_rate)*iii)+1;
end
p_bin_edges = [1;p_bin_edges];
for n = 1:length(p_bin_edges),
    if n<length(p_bin_edges)-1,
        p_bin_inds = p_bin_edges(n):(p_bin_edges(n+1)-1);
    else
        p_bin_inds = p_bin_edges(n):length(coll_pref_FR);
    end
    p_bin_FRmed(n,1) = median(coll_pref_FR(p_bin_inds));
    p_bin_FRmode(n,1) = mode(coll_pref_FR(p_bin_inds));
    p_bin_noise(n,1) = std(coll_pref_Vm(p_bin_inds));
    %supplemental (see third plot below)
    early_p = round(length(p_bin_inds)/2);
    late_p = early_p+1;
    early_p_bin_noise(n,1) = std(coll_pref_FR(p_bin_inds(1:early_p)));
    late_p_bin_FRmode(n,1) = mode(coll_pref_FR(p_bin_inds(late_p:end)));
end
%****
%null
for ii = 1:round(length(coll_null_FR)/floor(0.005*sampling_rate)),
    n_bin_edges(ii,1) = (floor(0.005*sampling_rate)*ii)+1;
end
n_bin_edges = [1;n_bin_edges];
for nn = 1:length(n_bin_edges),
    if nn<length(n_bin_edges)-1,
        n_bin_inds = n_bin_edges(nn):(n_bin_edges(nn+1)-1);
    else
        n_bin_inds = n_bin_edges(nn):length(coll_null_FR);
    end
    n_bin_FRmed(nn,1) = median(coll_null_FR(n_bin_inds));
    n_bin_FRmode(nn,1) = mode(coll_null_FR(n_bin_inds));
    n_bin_noise(nn,1) = std(coll_null_Vm(n_bin_inds));
    %supplemental (see third plot below)
    early_n = round(length(n_bin_inds)/2);
    late_n = early_n+1;
    early_n_bin_noise(nn,1) = std(coll_null_FR(n_bin_inds(1:early_n)));
    late_n_bin_FRmode(nn,1) = mode(coll_null_FR(n_bin_inds(late_n:end)));
end
%generate corresponding mode scatterplot
ffno = figure;
ps = scatter(p_bin_noise,p_bin_FRmode,'filled','r');
hold on;
pn = scatter(n_bin_noise,n_bin_FRmode,'filled','b');
xlabel('Vm noise, SD (mV)');
ylabel('Firing rate, bin mode');
legend('Pref','Null','Location','SouthEast');
title('5-ms binned FR mode vs. Vm noise');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'FRmode_v_noise_binned.fig']);
    close(ffno);
else
end
%generate corresponding median scatterplot
ffne = figure;
ps = scatter(p_bin_noise,p_bin_FRmed,'filled','r');
hold on;
pn = scatter(n_bin_noise,n_bin_FRmed,'filled','b');
xlabel('Vm noise, SD (mV)');
ylabel('Firing rate, bin median');
legend('Pref','Null','Location','SouthEast');
title('5-ms binned FR median vs. Vm noise');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'FRmed_v_noise_binned.fig']);
    close(ffne);
else
end
%generate supplemental plot dividing early bin noise from late bin FR
ffns = figure;
ps = scatter(early_p_bin_noise,late_p_bin_FRmode,'filled','r');
hold on;
pn = scatter(early_n_bin_noise,late_n_bin_FRmode,'filled','b');
xlabel('Vm noise, SD (mV) - early bin only');
ylabel('Firing rate, LATE BIN MODE');
legend('Pref','Null','Location','SouthEast');
title('5-ms bins: late-bin FR mode vs. early-bin Vm noise');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'lateFRmode_v_earlyNoise_binned.fig']);
    close(ffns);
else
end

%generate x-binned error bar plot
if mod(ceil(max([p_bin_noise;n_bin_noise])*10),2),
    max_xbin = ((ceil(max([p_bin_noise;n_bin_noise])*10))+1)/10;
else
    max_xbin = (ceil(max([p_bin_noise;n_bin_noise])*10))/10;
end
xbin_range = 0:0.2:max_xbin;
p_xbin_locs = cell(length(xbin_range)-1,1);
n_xbin_locs = cell(length(xbin_range)-1,1);
for b = 1:length(xbin_range)-1,
    p_xbin_locs{b,1} = ...
        find(p_bin_noise>xbin_range(1,b)&p_bin_noise<xbin_range(1,b+1));
    n_xbin_locs{b,1} = ...
        find(n_bin_noise>xbin_range(1,b)&n_bin_noise<xbin_range(1,b+1));
end
for bb = 1:length(xbin_range)-1,
    p_new_FR_mean(bb,1) = mean(p_bin_FRmode(p_xbin_locs{bb,1}));
    p_new_FR_sem(bb,1) = ...
        std(p_bin_FRmode(p_xbin_locs{bb,1}))./...
        sqrt(length(p_xbin_locs{bb,1}));
    n_new_FR_mean(bb,1) = mean(n_bin_FRmode(n_xbin_locs{bb,1}));
    n_new_FR_sem(bb,1) = ...
        std(n_bin_FRmode(n_xbin_locs{bb,1}))./...
        sqrt(length(n_xbin_locs{bb,1}));
end
xbin_centers = (xbin_range(2:end)+xbin_range(1:end-1))/2;
ffnx = figure;
scatter(xbin_centers,p_new_FR_mean,'filled','r');
hold on;
for d = 1:length(xbin_centers),
    p_xeb = line([xbin_centers(1,d) xbin_centers(1,d)],...
        [(p_new_FR_mean(d,1)-...
        p_new_FR_sem(d,1)) (p_new_FR_mean(d,1)+...
        p_new_FR_sem(d,1))]);
    p_xeb.LineWidth = 2.0;
    p_xeb.Color = 'r';
end
scatter(xbin_centers,n_new_FR_mean,'filled','b');
for dd = 1:length(xbin_centers),
    n_xeb = line([xbin_centers(1,dd) xbin_centers(1,dd)],...
        [(n_new_FR_mean(dd,1)-...
        n_new_FR_sem(dd,1)) (n_new_FR_mean(dd,1)+...
        n_new_FR_sem(dd,1))]);
    n_xeb.LineWidth = 2.0;
    n_xeb.Color = 'b';
end
xlabel('Vm noise, SD (mV)');
ylabel('Firing rate');
%legend('Pref','Null','Location','SouthEast');
title('X-binned (0.2 mV) mean FR and SEM');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'FR_v_noise_XbinSEM.fig']);
    close(ffnx);
else
end
%generate additional supplemental plot using demeaned firing rate (from
%means binned in x according to above; mean FR for each bin is subtracted
%from each FR datapoint in each bin, and then the binned-mean of resulting
%datapoints is calculated (basically residual plot from "quick-and-dirty"
%type of local regression)
p_subt_FR = cell(length(xbin_range)-1,1);
n_subt_FR = cell(length(xbin_range)-1,1);
for k = 1:length(xbin_range)-1,
    p_subt_FR{k,1} = p_bin_FRmode(p_xbin_locs{k,1})-p_new_FR_mean(k,1);
    n_subt_FR{k,1} = n_bin_FRmode(n_xbin_locs{k,1})-n_new_FR_mean(k,1);
end
for k = 1:length(xbin_range)-1,
    mp_subt_FR(k,1) = mean(p_subt_FR{k,1});
    mn_subt_FR(k,1) = mean(n_subt_FR{k,1});
end
ffnss = figure;
ps = scatter(xbin_centers,mp_subt_FR,'filled','r');
hold on;
pn = scatter(xbin_centers,mn_subt_FR,'filled','b');
xlabel('Vm noise, SD (mV)');
ylabel('Mean (Firing rate - mean(FR))');
legend('Pref','Null','Location','SouthEast');
title('X-binned (0.2 mV) demeaned Firing rate vs Noise');
hold off
if save_it == 1,
    saveas(gcf,[pwd filesep 'demeanedFR_v_noise_Xbin.fig']);
    close(ffnss);
else
end


%Method 2: Sliding 5-ms bin captures FR at every timestep; matched to Vm SD
%from same 5-ms sliding bin
l_bound = ceil(0.005*sampling_rate);
win_len = floor(0.005*sampling_rate);
%pref
p_new_FR = coll_pref_FR(l_bound:length(coll_pref_FR));
for m = 1:length(p_new_FR),
    p_win_sd(m,1) = std(coll_pref_Vm(m:m+win_len));
end
%null
n_new_FR = coll_null_FR(l_bound:length(coll_null_FR));
for mm = 1:length(n_new_FR),
    n_win_sd(mm,1) = std(coll_null_Vm(mm:mm+win_len));
end
%generate corresponding sliding window scatterplot
ffnw = figure;
ps = scatter(p_win_sd,p_new_FR,'filled','r');
hold on;
pn = scatter(n_win_sd,n_new_FR,'filled','b');
xlabel('Vm noise, SD (mV)');
ylabel('Firing rate');
legend('Pref','Null','Location','SouthEast');
title('Sliding window: 5-ms VmSD noise matched to inst. FR');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep 'slidingWin_instFR_v_noise.fig']);
    close(ffnw);
else
end






%%
%********************Firing rate autocorrelation**************************
%*************************************************************************
%ultimately compares actual pref and null FR autocorr to FR autocorr
%predicted by data fits
autocorr_FR = {};
bin_FR_ACH = 0;    % ***NEED TO ENTER THIS INTO FUNCION INPUTS***
for m2 = 1:length(stimvalues),
    autocorr_FR{m2,1} = struct('stim_FR_corr',[],'stim_FR_lags',[]);
end
fa_window_size = stim_duration(1,1); 
if bin_FR_ACH == 1,
    bin_size = 0.005;
    for m3 = 1:length(stimvalues),
        if isempty(FR_inst{m3,1}),
            continue;
        else
        end
        win_ind = ...
            (1:floor(bin_size/sample_interval):floor(fa_window_size/sample_interval))';
        for m4 = 1:length(win_ind),
            if m4 == length(win_ind),
                f_mb(m4,1) = mean(FR_inst{m3,1}(win_ind(m4,1):end));
            else
                f_mb(m4,1) = mean(FR_inst{m3,1}(win_ind(m4,1):win_ind(m4+1,1)));
            end
        end
        [corr,lags] = xcorr(f_mb,f_mb,length(win_ind),'coeff');
        lags = lags*bin_size;
        autocorr_FR{m3,1}.stim_FR_corr = corr;
        autocorr_FR{m3,1}.stim_FR_lags = lags;
        clear f_mb corr lags;
    end
else
    bin_size = 0.001;
    for m3 = 1:length(stimvalues),
        [corr,lags] = ...
            xcorr(FR_inst{m3,1},FR_inst{m3,1},...
            round(fa_window_size/bin_size),'coeff');
        lags = lags*bin_size;
        autocorr_FR{m3,1}.stim_FR_corr = corr;
        autocorr_FR{m3,1}.stim_FR_lags = lags;
    end
end
    
for m5 = 1:length(stimvalues),
    f = figure;
    ach_f = ...
        bar(autocorr_FR{m5,1}.stim_FR_lags,...
        autocorr_FR{m5,1}.stim_FR_corr,'k');
    hold on;
    xlabel('Lags (seconds)');
    ylabel('Correlation');
    %xlim([-(stim_duration(1,1)) stim_duration(1,1)]);
    title(['Firing rate autocorrelation, Direction angle ',...
        num2str(stimvalues(1,m5))]);
    if save_it == 1,
        if bin_FR_ACH == 1,
            saveas(gcf,[pwd filesep 'binned_FR_autocorr_angle',...
                num2str(stimvalues(1,m5)),'.fig']);
            close(f);
        else
            saveas(gcf,[pwd filesep 'FR_autocorr_angle',...
                num2str(stimvalues(1,m5)),'.fig']);
            close(f);
        end
    else
    end
end

%%
 %
    %
    %
    %
    %Generate raw plot overlays of Vm and spiking traces to show divergence
    %in pref and null transfer functions
    f_o = figure;
    if detrend == 1,
        plot(t_vec,sub_wv,'r');
    else
        plot(t_vec,waveform_dnlp,'r');
    end
    hold on;
    plot(t_vec,Vm_waveform,'b');
    pref_markers = nStim_ON(find(stimvalues(stimorder)==pref_stim));
    if (pref_stim+180)>360
        null_markers = ...
            nStim_ON(find(stimvalues(stimorder)==(pref_stim-180)));
    else
        null_markers = ...
            nStim_ON(find(stimvalues(stimorder)==(pref_stim+180)));
    end
    sz = 50;
    mh = scatter(pref_markers,repmat([-30],...
        length(pref_markers),1),sz,'filled','d');
    mh.MarkerFaceColor = [0 0 0];
    mh.MarkerEdgeColor = 'k';
    hold on;
    mmh = scatter(null_markers,repmat([-30],...
        length(null_markers),1),sz,'filled','v');
    mmh.MarkerFaceColor = [0 1 0];
    mmh.MarkerEdgeColor = 'g';
    legend('Spike trace','Vm trace','Pref on trials',...
        'Null on trials','Location','NorthEast');
    xlabel('Time vector (seconds)');
    ylabel('Membrane potential (mV)');
    if save_it == 1,
        saveas(gcf,[pwd filesep 'FullTrace_SpikeVm_Overlay.fig']);
        close(f_o);
    else
    end
    
    
    
    
    %%
     %**********************CREATE CYCLE TRACES****************************
    %*********************************************************************
    format long
    frameRate = 1/60;
    for i2 = 1:length(stimvalues)-1,
        N_trialFrames = length(mti2_{1,i2}.frameTimes);
        N_cycleFrames = N_trialFrames/nCycles(1,1);
        pad_nFrames = 3;
        c_rep_locs = find(stimorder == i2);
        for k2 = 1:length(c_rep_locs),
            for ii2 = 1:nCycles(1,1),
                c_light = ...
                    round((nStim_ON(c_rep_locs(1,k2),1)+...
                    (N_cycleFrames*frameRate*(ii2-1)))/sample_interval);
                c_a = ...
                    c_light-round((pad_nFrames*frameRate)/sample_interval);
                c_light_off = round((nStim_ON(c_rep_locs(1,k2),1)+...
                    ((N_cycleFrames/2)*frameRate*(ii2)))/sample_interval);
                c_b = round((nStim_ON(c_rep_locs(1,k2),1)+...
                    (N_cycleFrames*frameRate*ii2))/sample_interval)+...
                    round((pad_nFrames*frameRate)/sample_interval);
                cycle_array_tvec{i2,ii2,k2} = ...
                    [0:sample_interval:sample_interval*(c_b-c_a)];
                if detrend == 1,
                    cycle_array_Sp{i2,ii2,k2} = sub_wv(c_a:c_b,1);
                else
                    cycle_array_Sp{i2,ii2,k2} = waveform_dnlp(c_a:c_b,1);
                end
                cycle_array_Vm{i2,ii2,k2} = Vm_waveform(c_a:c_b,1);
                c_a_array(i2,ii2,k2) = c_a*sample_interval;
                c_b_array(i2,ii2,k2) = c_b*sample_interval;
                rel_light_array(i2,ii2,k2) = (c_light-c_a)*sample_interval;
                abs_light_array(i2,ii2,k2) = c_light;
                rel_light_off_array(i2,ii2,k2) = ...
                    (c_light_off-c_a)*sample_interval;
                abs_light_off_array(i2,ii2,k2) = c_light_off;
            end
        end
    end
    %use full trace to visually inspect light on locations
    figure;
    if detrend == 1,
        plot(t_vec,sub_wv,'k');
    else
        plot(t_vec,waveform_dnlp,'k');
    end
    hold on;
    for k3 = 1:size(rel_light_array,1),
        c_rep_locs = find(stimorder == k3);
        for k4 = 1:length(c_rep_locs),
            for k5 = 1:nCycles(1,1),
                line([(abs_light_array(k3,k5,k4)*...
                    sample_interval) (abs_light_array(k3,k5,k4)*...
                    sample_interval)],[-150 50]);
            end
        end
    end
    %find cycle CV of ISIs   (consider adding CV2 here to account for
    %'burstiness'; see Holt et al. 1996)
    stim_collected_cycles = cell(length(stimvalues)-1,1);
    %within_cycle_CV2 = NaN(length(stimvalues)-1,nCycles(1,1),stimmatch_reps);  ***
    within_cycle_CV2 = NaN(length(stimvalues)-1,nCycles(1,1),length(reps));
    for k3 = 1:length(stimvalues)-1,
        stacked_stimcycle = [];
        c_rep_locs = find(stimorder == k3);
        collected_rep_ISIs = [];
        for k4 = 1:length(c_rep_locs),
            collected_trial_ISIs = [];
            for k5 = 1:nCycles(1,1),
                if k5 < nCycles(1,1),
                    s_i = ...
                        spike_locations(find((abs_light_array(k3,k5,k4)<...
                        spike_locations)&(spike_locations<...
                        abs_light_array(k3,k5+1,k4))));
                    stacked_stimcycle = ...
                        [stacked_stimcycle;(s_i-abs_light_array(k3,k5,k4))];
                    N_cycle_spikes(k3,k5,k4) = length(s_i);
                else
                    s_i = ...
                        spike_locations(find((spike_locations>...
                        abs_light_array(k3,k5,k4))&(spike_locations<...
                        abs_light_array(k3,k5,k4)+...
                        ((N_cycleFrames*frameRate)/sample_interval))));
                    stacked_stimcycle = ...
                        [stacked_stimcycle;(s_i-abs_light_array(k3,k5,k4))];
                    N_cycle_spikes(k3,k5,k4) = length(s_i);
                end
                cycle_ISIs = diff(s_i).*sample_interval;
                if ~isempty(cycle_ISIs),
                    within_cycle_CV(k3,k5,k4) = ...
                        std(cycle_ISIs)/mean(cycle_ISIs);
                    within_cycle_CV2(k3,k5,k4) = ...
                        nanmean((2.*abs(cycle_ISIs(2:end)-cycle_ISIs(1:end-1)))./...
                        (cycle_ISIs(1:end-1)+cycle_ISIs(2:end)));
                    collected_trial_ISIs = [collected_trial_ISIs;cycle_ISIs];
                else
                end
            end
            collected_rep_ISIs = [collected_rep_ISIs;collected_trial_ISIs];
        end
        collected_ISIs{k3,1} = collected_rep_ISIs;
        bystim_cycle_CV{k3,1} = ...
            std(collected_rep_ISIs)/mean(collected_rep_ISIs);
        stim_collected_cycles{k3,1} = sort(stacked_stimcycle);
    end
    f = figure;
    for i8 = 1:length(collected_ISIs),
        scatter(repmat(stimvalues(1,i8),...
            length(collected_ISIs{i8,1}),1),collected_ISIs{i8,1}(:),'k');
        hold on;
    end
    xlim([0 360]);
    xlabel('Direction angle');
    ylabel('Cycle-framed ISIs (sec.)');
    if save_it == 1,
        saveas(gcf,[pwd filesep 'CycleISIDists_byStim.fig']);
        close(f);
    else
    end
    %REGARDING CV2 MEASUREMENT - begin with time-dep. CV2 (a value that
    %varies in real-time and is a comparison of adjacent ISIs), get a
    %single CV2 measurement for each cycle by averaging time-dep. CV2s over
    %the cycle; then split analysis into median CV2 across reps. (i.e.
    %matched-cycles) or through trial (i.e. matched-trials)
    %plot CV2 by stim
    f = figure;
    matched_trialMed = nanmedian(nanmedian(within_cycle_CV2,2),3);
    matched_cycleMed = nanmedian(nanmedian(within_cycle_CV2,3),2);
    unmatched_Med = nanmedian(reshape(within_cycle_CV2,...
        size(within_cycle_CV2,1),size(within_cycle_CV2,2)*size(within_cycle_CV2,3)),2);
    for n = 1:length(collected_ISIs),
        for n1 = 1:length(c_rep_locs),
            for n2 = 1:nCycles(1,1),
                scatter(stimvalues(1,n),...
                    within_cycle_CV2(n,n2,n1),'k');
                hold on;
            end
        end
        scatter(stimvalues(1,n),matched_cycleMed(n,1),'filled','b');
        hold on;
        scatter(stimvalues(1,n),matched_trialMed(n,1),'filled','m');
        scatter(stimvalues(1,n),unmatched_Med(n,1),'filled','r');
    end
    xlim([0 360]);
    ylim([0 2]);
    xlabel('Direction angle');
    ylabel('Cycle CV2');
    title('CV2 by stim.; all cycles');
    if save_it == 1,
        saveas(gcf,[pwd filesep 'CV2_bystim_medians.fig']);
        close(f);
    else
    end
    %plot CV2 by cycle
    f = figure;
    partial_cycleMed = nanmedian(within_cycle_CV2,3);
    cycle_vals = 1:nCycles(1,1);
    for n = 1:length(collected_ISIs),
        plot(cycle_vals,partial_cycleMed(n,:));
        hold on;
    end
    xlabel('Cycle number');
    ylabel('Cycle-matched median CV2');
    legend(cellstr(num2str(stimvalues(1,1:length(stimvalues)-1)')),...
        'Location','EastOutside');
    if save_it == 1,
        saveas(gcf,[pwd filesep 'CV2_bycycle_medians.fig']);
        close(f);
    else
    end
    
    
   
    
    
    %*****calculate cycle-based spike-train autocorrelation
    autocorr_output = {};
    for i9 = 1:length(stimvalues)-1,
        autocorr_output{i9,1} = struct('cycle_corr',[],'cycle_lags',[]);
    end
    window_size = stim_duration(1,1)/10; %note- this does not accommodate blocks with variable cycle durations
    bin_size = 0.001;
    for i11 = 1:length(stimvalues)-1,
        if isempty(stim_collected_cycles{i11,1}),
            cycle_spike_times = NaN;
            continue;
        else
            cycle_spike_times = stim_collected_cycles{i11,1}/sampling_rate;
        end
        ach_bin_edges = [-(window_size):bin_size:window_size];
        N_ach = histc(cycle_spike_times,ach_bin_edges);
        ach_bin_centers = (ach_bin_edges(1:end-1)+ach_bin_edges(2:end))/2;
        N_ach = N_ach(1:end-1);
        [corr,lags] = xcorr(N_ach,N_ach,round(window_size/bin_size));
        lags = lags*bin_size;
        autocorr_output{i11,1}.cycle_corr = corr;
        autocorr_output{i11,1}.cycle_lags = lags;
    end
    %generate plots
     for i12 = 1:length(stimvalues)-1,
        f = figure;
        ach_h = ...
            bar(autocorr_output{i12,1}.cycle_lags,...
            autocorr_output{i12,1}.cycle_corr,'k');
        hold on;
        xlabel('Lags');
        ylabel('Count');
        xlim([-(stim_duration(1,1)/10) stim_duration(1,1)/10]);
        title(['Cycle-averaged Spike-train autocorrelation, Direction angle ',...
            num2str(stimvalues(1,i12))]);
        if save_it == 1,
            saveas(gcf,[pwd filesep 'cycle_spiketrain_autocorr_angle',...
                num2str(stimvalues(1,i12)),'.fig']);
            close(f);
        else
        end
     end
     
     %******Cycle-based Vm plots******************************
     %********************************************************
     allCycle_meanVm = cellfun(@mean, cycle_array_Vm);
     allCycle_SDVm = cellfun(@std, cycle_array_Vm);
     stim_meanVm = mean(mean(allCycle_meanVm,3),2);
     matchedCycle_meanVm = mean(allCycle_meanVm,3);
     %matchedCycle_SDVm{v,v1} = ...
     %    cellfun(@(x,y) cycle_array_Vm{x,y,:},...
     %    1:size(cycle_array_Vm),'UniformOutput',0);
     for v = 1:length(stimvalues)-1,
         coll_Vmcycles = [];
         for v1=1:nCycles(1,1),
             coll_Vmreps = [];
             for v2=1:length(c_rep_locs),
                coll_Vmreps = [coll_Vmreps;cycle_array_Vm{v,v1,v2}(:)];
             end
             cycle_SDVm(v,v1) = std(coll_Vmreps);
             coll_Vmcycles = [coll_Vmcycles;coll_Vmreps];
         end
         stim_SDVm(v,1) = std(coll_Vmcycles);
     end
     %plot Vm SD by cycle
     fn = figure;
     cycle_num = 1:size(cycle_SDVm,2);
     if ~exist('pref_stim','var'),
         pref_stim = otpref_total;
     else
     end
     if (pref_stim+180)>=360
         nullval = (pref_stim+180)-360;
     else
         nullval = pref_stim+180;
     end
     for a = 1:size(cycle_SDVm,1),
         nl = plot(cycle_num,cycle_SDVm(a,:));
         if a == find(stimvalues==pref_stim)||a == find(stimvalues==nullval),
             nl.LineWidth = 3.0;
         else
         end
         hold on;
     end
     legend(cellstr(num2str(stimvalues')),'Location','EastOutside');
     xlabel('Cycle number');
     ylabel('Cycle Vm SD (mV)');
     if save_it == 1,
         saveas(gcf,[pwd filesep 'VmSD_bycycle.fig']);
         close(fn);
     else
     end
     %plot Vm SD as function of Vm mean
     fs = figure;
     for aa = 1:size(matchedCycle_meanVm,1),
         sm = scatter(matchedCycle_meanVm(aa,:),cycle_SDVm(aa,:));
         if aa == find(stimvalues==pref_stim)||aa == find(stimvalues==nullval),
             sm.MarkerFaceColor = 'flat';
         else
         end
         hold on;
     end
     legend(cellstr(num2str(stimvalues')),'Location','EastOutside');
     xlabel('Cycle mean Vm (mV)');
     ylabel('Cycle Vm SD (mV)');
     if save_it == 1,
         saveas(gcf,[pwd filesep 'VmMean_v_VmSD.fig']);
         close(fs);
     else
     end
     
     
     
     
     
    %%
    clear a a1 a2 aa ach_f ach_h b b1 b2 bb c c_a c_b dd
    clear fas fc fel ff fff ffn ffne ffno ffns ffnss ffnw ffnx fhs fi
    clear fn fn3 fna fnf fni fnr fp3 fpa fpf fpi fpr fps fr fs 
    clear i11 i12 i2 i3 i4 i5 i6 i7 i8 i9 ii ii2 ii3 ii4 ii5 iii
    clear j j2 j3 j4 j5 j6 j7 j8 j9 jj k k2 k3 k4 k5
    clear m m2 m3 m5 mm n n1 n2 nn q qq sub_loc
    %*********************************************
    %******SAVE WORKSPACE FOR THIS CELL
    save([char(directories),'_allVars_',char(date)]);
    %*********************************************
    
    if overwrite == 0,
        cd ..
    else
    end
    
end