function [ Vm_waveform,Vm_array,Vm_trial_array,Vm_array_smooth ] = ...
    extract_Vm( waveform_dnlp,sub_wv,sample_interval,h_margins,stimvalues,...
    nStim_ON,nStim_OFF,stimorder,reps,use_detrend,smooth_Vm,save_it)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%*****Remove spikes to get raw Vm traces and trial-averaged Vm trace**
    %*********************************************************************
    win_d = ceil(0.01/sample_interval); 
    %alt: win_d = 41;  %corresponds to the 4-ms used in Priebe et al. 2004
    if mod(win_d,2) == 0,
        win_d = win_d+1;
    else
    end
    if use_detrend == 1,
        Vm_waveform = medfilt1(sub_wv,win_d);
    else
        Vm_waveform = medfilt1(waveform_dnlp,win_d);
    end
    trial_Range = cell(3,1);
    Vm_array = cell(length(stimvalues),1);
    Vm_trial_array = cell(length(stimvalues),length(reps));
    dRange_array = cell(length(stimvalues),1);
    for jj = 1:length(stimvalues),
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
        %smooth traces
        %win = round(0.05/sample_interval);
        %if mod(win,2)==0,
        %    win = win + 1;
        %else
        %end
        %Vm_trialAve_smooth = sgolayfilt(Vm_trialAve,5,win);
        %Vm_array_smooth{jj} = Vm_trialAve_smooth;
    end
        
        %smooth Vm traces
        %method 2:
        
    Vm_array_smoothed = cell(length(stimvalues),1);
    win = 0.05;
    t_vec = disp_Range;
    Xn = [];
    for qq = 2:length(t_vec),
        Xn(end+1) = mean([t_vec(1,qq) t_vec(1,qq-1)]);
    end
    t_index = t_vec/sample_interval;
    for jj = 1:length(stimvalues),
        if ~isempty(Vm_array{jj}),
            pad_len = round(win/sample_interval/2);
            padded_Vm = [fliplr(Vm_array{jj}(1:pad_len));...
                Vm_array{jj};fliplr(Vm_array{jj}(length(Vm_array{jj})-pad_len:end))];
            Yn = padded_Vm;
            g_kernel = gausswin(round(win/sample_interval));
            norm_g_kernel = g_kernel./sum(g_kernel);
            temp_array_smooth = conv(Yn,norm_g_kernel,'same');
            Vm_array_smooth{jj,1} = temp_array_smooth(pad_len+1:length(temp_array_smooth)-pad_len-1); 
        else
            Vm_array_smooth{jj,1} = zeros(length(t_vec)-1,1);
        end
        
    end
    
    smooth_Vm = 1;
    t_vec = disp_Range;
    if smooth_Vm == 1,
        for jj = 1:length(stimvalues)-1,
            f = figure;
            plot(t_vec,Vm_array_smooth{jj},'k');
            hold on;
            xlabel('Time (seconds)');
            ylabel('Membrane potential (mV)');
            hold off;
            if save_it == 1,
                saveas(gcf,[pwd filesep 'Gauss_smoothed_trialAve_Vm_',...
                    num2str(stimvalues(1,jj)),'.fig']);
                close(f);
            else
            end
        end
    else
    end
        
end



