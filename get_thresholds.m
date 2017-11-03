function [ th_loc,Vm_th,new_trace,sub_wv,max_dvdt ] = ...
    get_thresholds( spike_trace,t_vec,spike_locations,sample_interval,use_detrend,reset2phys )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
    max_dvdt = [];
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
        max_dvdt = [max_dvdt;max_slope(1,1)];
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

    if use_detrend == 1,   
        t_vec = reshape(t_vec,1,length(t_vec));
        [new_trace,sub_wv] = detrend_spiketrace(spike_trace,t_vec,th_loc,Vm_th,reset2phys);
    
    else
        new_trace = spike_trace;
        sub_wv = spike_trace;
    end
    
end

