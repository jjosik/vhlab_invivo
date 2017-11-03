function [ output_args ] = drift_test( data_file,pre_process )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load(data_file);
sample_interval = spike2data_Ch1.interval;
sampling_rate = round(1/sample_interval);
waveform = spike2data_Ch1.values.*100;  %in mV
t_vec = [sample_interval:sample_interval:(length(waveform)*sample_interval)];

switch pre_process

    case 0  %with no pre-processing
        f = figure;
        plot(t_vec,waveform,'k');
        
        display_spikes = 1;
        auto_detect = 0;
        [spike_locations,h_margins] = find_spikes(waveform,t_vec,...
        sampling_rate,display_spikes,auto_detect);
        
        %biophysical threshold
        %first with detrending and reset to physiological range
        reset2phys = 1;
        use_detrend = 1;
        [th_loc,Vm_th,new_trace,sub_wv] = get_thresholds(waveform,t_vec,...
            spike_locations,sample_interval,use_detrend,reset2phys);
        %now again without detrending or reset
        reset2phys = 0;
        use_detrend = 0;
        [th_loc_,Vm_th_,new_trace_,sub_wv_] = get_threshold(waveform,t_vec,...
            spike_locations,sample_interval,use_detrend,reset2phys);











    case 1

        f = figure;
        plot(t_vec,waveform,'k');


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
            waveform = waveform(1:round(x_break/sample_interval));
            close(f);
        else
        end

        

end

