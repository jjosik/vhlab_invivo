function [ freq_vec,ampsp_vec,psd_vec,autocorr_output ] = ...
    process_spiketrace_FFT( spike_trace,t_vec,sampling_rate,auto_detect,...
    detrend,use_windowing,display_spikes )
%**********************
%PROCESS_SPIKETRACE_FFT - 
%   INPUTS -    SPIKE_TRACE: 1D raw vector containing values for voltage at
%                   every timestep in spiketrace (units should be in mV)
%               T_VEC: 1D raw vector containing timestamps in seconds (i.e.
%                   1 second/sampling_rate)
%               SAMPLING_RATE: rate of signal acquisition (=1 / TIMESTEP)
%               AUTO_DETECT: Set to '1' if coarse threshold criterion
%                   should be set to default -20 mV, or '0' if user wishes to
%                   select criterion using GUI.
%               DETREND: Set to '1' if spike_trace should be detrended
%                   before further processing (see function
%                   DETREND_SPIKETRACE.m)
%               USE_WINDOWING: Set to '1' if Fourier transform evaluation 
%                   should be averaged over weighted windows; current
%                   default set to Hanning window-type.
%               DISPLAY_SPIKES: when set to '1', will diplay plot of
%                   complete spike trace with spike raster locations;
%                   allows for visual inspection of spike identifications.
%
%   OUTPUTS -   FREQ_VEC: "x-axis" output of Fourier analysis
%               AMPSP_VEC: Amplitude spectrum ("y-axis") output of Fourier
%                   analysis (NOT recommended for spike data)   
%               PSD_VEC: Power spectral density ("y-axis") output of spike Fourier analysis
%                   (computed by FFT of autocorrelogram)
%               AUTOCORR_OUTPUT: Corresponding autocorrelation traces used
%                   in the FFT computation described above.  The variable is a
%                   structure that contains 1.) a TRIAL_CORR vector variable (y) and
%                   2.) a TRIAL_LAGS vector variable (x) for autocorrelation plotting.  

sample_interval = 1/sampling_rate;
reset2phys = 1;     %default - should be set up as input to PROCESS_SPIKETRACE_FFT function

if auto_detect == 1,
    th_ = -20.0;   %in mV
else
    f = figure;
    plot(t_vec,spike_trace,'k');
    hold on;
    display(['Use cursor to select coarse detection threshold for the experiment trace.  Press ''Return'' when finished']);
    [~,th_] = ginput;
end
    
spike_locations = [];
spike_stepLocations = find(spike_trace>th_);
if isempty(spike_stepLocations),
    msgID = 'MATLAB:th_test';
    msgtext = 'No spikes were detected in this dataset.';
    ME = MException(msgID,msgtext);
    throw(ME)
else
end
onsets = diff(spike_stepLocations);
sp_edge = find(onsets>1);
ind_edge = [spike_stepLocations(1,1);spike_stepLocations(sp_edge(:,1)+1,1)];

for ii = 1:length(ind_edge),
    if ii<length(ind_edge),
        subset = find((spike_stepLocations>ind_edge(ii,1))&...
            (spike_stepLocations<ind_edge(ii+1,1)));
    else
        subset = find(spike_stepLocations>ind_edge(ii,1));
    end
    total_subset = [ind_edge(ii,1);spike_stepLocations(subset,1)];
    [~,sub_loc] = max(spike_trace(total_subset(1:end),1));
    spike_locations(ii,1) = ind_edge(ii,1)+(sub_loc-1);
end


if display_spikes == 1,
    if ishandle('f'),
        gcf;
        scatter(spike_locations(:,1)./sampling_rate,...
            zeros(length(spike_locations),1),'r');
        title('Scatterpoints show spike peak raster locations');
        if max(spike_trace) < 0,
            ylim([-100 10]);
        else
            ylim([-100 max(spike_trace)+0.1*abs(max(spike_trace))]);
        end
    else
        plot(t_vec,spike_trace,'k');
        hold on;
        scatter(spike_locations(:,1)./sampling_rate,...
            zeros(length(spike_locations),1),'r');
        title('Scatterpoints show spike peak raster locations');
        if max(spike_trace) < 0,
            ylim([-100 10]);
        else
            ylim([-100 (max(spike_trace)+0.1*abs(max(spike_trace)))]);
        end
    end
else
end

%*******Assess biophysical spike threshold using Azouz & Gray, 1999
%*******methods****************************************************
%******************************************************************
        
vt_slope = gradient(spike_trace,sample_interval); 
slope_criterion = 0.033; %fraction of peak dV/dt calibrated by visual inspection to match threshold
ex_flag = [];   %Spikes that fall on the edge of acquisition or selected trace window, or otherwise have
                %incomplete waveforms are excluded.
peak_slope = zeros(length(spike_locations),1);
slope_Vm = zeros(length(spike_locations),1);
th_loc = zeros(length(spike_locations),1);
Vm_th = zeros(length(spike_locations),1);
for q = 1:length(spike_locations),
    search_size = ceil(0.004/sample_interval);
    if search_size > spike_locations(q,1),
        ex_flag = [ex_flag;q];
        continue;
    else
    end
    search_pad = (spike_locations(q,1)-search_size):...
        spike_locations(q,1)-1;
    [max_slope,peak_ind] = max(vt_slope(reshape(search_pad,...
        length(search_pad),1),1));
    peak_slope(q,1) = search_pad(1,1)+(peak_ind-1);
    new_search = search_pad(1,1):peak_slope(q,1);
    slope_Vm(q,1) = spike_trace((search_pad(1,1)+(peak_ind-1)),1);
    [~,th_ind] = min(abs((slope_criterion*max_slope)-...
        vt_slope(new_search,1)));
    th_loc(q,1) = search_pad(1,1)+(th_ind-1);
    Vm_th(q,1) = spike_trace((search_pad(1,1)+(th_ind-1)),1);
end
if any(ex_flag),
    th_loc(ex_flag) = [];
    Vm_th(ex_flag) = [];
else
end

if detrend == 1,   
    t_vec = reshape(t_vec,1,length(t_vec));
    [new_trace] = detrend_spiketrace(spike_trace,t_vec,th_loc,Vm_th,reset2phys);
    
else
    new_trace = spike_trace;
end

%**********************************************************************
%**********Compute autocorrelation of spike response,******************
%**********then FFT to get power spectrum *****************************
autocorr_output = struct('trial_corr',[],'trial_lags',[]);
half_window = length(t_vec);
bin_size = 0.002;
ach_bin_edges = [-(half_window)/sampling_rate:bin_size:half_window/sampling_rate];
N_ach = histc(spike_locations./sampling_rate,ach_bin_edges);
N_ach = N_ach(1:end-1);
[corr,lags] = xcorr(N_ach,N_ach,round((half_window/sampling_rate)/bin_size));
lags = lags.*bin_size;
autocorr_output.trial_corr = corr;
autocorr_output.trial_lags = lags;

aft_times = (1:length(autocorr_output.trial_corr)).*sample_interval;
if use_windowing == 1,  %recommended if spectra will be averaged over cycles or trials
    series_ = reshape(autocorr_output.trial_corr,...
        length(autocorr_output.trial_corr),1);
    aw_series = series_.*hanning(length(autocorr_output.trial_corr));
else
    aw_series = autocorr_output.trial_corr;
end
afc = fft(aw_series,length(aft_times));
amp_coeff = reshape(sqrt(afc.*conj(afc)),length(afc),1);
power_coeff = reshape((afc.*conj(afc)),length(afc),1);
s_freqs = sampling_rate*(0:length(aft_times)-1)/length(aft_times);

if mod(length(s_freqs),2) ~= 0
    ampsp_vec =  amp_coeff(1:(length(amp_coeff)-1)/2,1);
    psd_vec = power_coeff(1:(length(power_coeff)-1)/2,1);
    freq_vec = s_freqs(1,1:(length(s_freqs)-1)/2);
else
    ampsp_vec = amp_coeff(1:(length(amp_coeff)/2),1);
    psd_vec = power_coeff(1:(length(power_coeff)/2),1);
    freq_vec = s_freqs(1,1:(length(s_freqs)/2));
end


end



