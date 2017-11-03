function [ new_trace,sub_tr,f,f_c ] = detrend_spiketrace( spike_trace,t_vec,th_time,th_value,reset2phys )
%DETREND_SPIKETRACE -  Summary of this function goes here
%
%   INPUTS -    SPIKE_TRACE: 1D raw vector containing values for voltage at
%                   every timestep in spiketrace (units should be in mV)
%               T_VEC: 1D raw vector containing timestamps in seconds (i.e.
%                   1 second/sampling_rate)
%               TH_TIME: absolute time indices indicating locations of
%                   all individual spike thresholds in the T_VEC vector
%               TH_VALUE: voltage value (in units of mV) of all individual 
%                   thresholds matched by index to onsets in TH_TIME
%               RESET2PHYS: when set to '1' makes an additional correction
%                   adding back to the simple detrended trace the median of the
%                   first 5-msecs in the original trace, i.e. returns the
%                   baseline to "physiological range".
%
%   OUTPUTS -   NEW_TRACE: 1D vector of *final* corrected voltage values
%                   (i.e, y-axis)
%               F:  figure handle for the plot of the original uncorrected trace
%                   (made using the the original T_VEC on time axis)
%               F_C: figure handle for the plot of the new corrected trace
%                   (also made using the original T_VEC)

default_span = 51;
check_interval = diff(t_vec(1,th_time(:,1)));
if max(check_interval) < 5000,    %NOTE: this condition may not be necessary; need to examine over many traces
    span = default_span;
else
    span = 501;         
end

ytr_line = smooth(t_vec(1,th_time(:,1)),th_value(:,1),span,'lowess');
Vq = interp1(t_vec(1,th_time(:,1)),ytr_line,t_vec(1,:),'linear','extrap');
sub_tr = spike_trace-Vq';
if reset2phys == 1,
    initial_Vm = median(spike_trace(1:default_span));
    new_trace = sub_tr+initial_Vm;
else
    new_trace = sub_tr;
end

%display 'before' trace with trendline
f = figure;
plot(t_vec,spike_trace,'k');
hold on;
scatter(t_vec(1,th_time(:,1)),th_value(:,1),'x','r');
plot(t_vec,Vq,'k--');
xlabel('Time (seconds)');
ylabel('Potential (mV)');
title(' BEFORE trace with trendline');
hold off;

%display 'after' trace 
f_c = figure;
plot(t_vec,new_trace,'k');
hold on;
scatter(t_vec(1,th_time(:,1)),th_value(:,1)-Vq(1,th_time(:,1))'+initial_Vm,'x','r');
xlabel('Time (seconds)');
ylabel('Potential (mV)');
title('AFTER trace');


end

