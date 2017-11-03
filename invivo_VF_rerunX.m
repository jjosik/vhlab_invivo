function [ Vth_est_coll,sta_vm_pref,sta_vm_null,rebin_sub_trial_array ] = ...
    invivo_VF_rerunX( save_it,rebin_trials,overwrite,model,collate_stims,...
    use_pre_analyzed,display_spikes,reset2phys,auto_detect,use_detrend,use_trial_baseline,...
    sub_folder,cur_dir,default,use_global_filter,filter_type,adapt_bins,...
    anchor_Vth,fit_it,code,get_VmSTA,stack_rasters,smooth_Vm,supp_singletrial_FR)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%   varargout is rebin_part_rasters variable
try
    load([sub_folder filesep cur_dir filesep 'analyzed_spike.mat'],'reps');
catch
    load([pwd filesep 'analyzed_spike.mat'],'reps');
end
try 
    load([sub_folder filesep cur_dir filesep 'analyzed_spike.mat'],'tempFrequency');
catch
    load([pwd filesep 'analyzed_spike.mat'],'tempFrequency');
end
if use_pre_analyzed == 1,
    saved_analysis = uigetfile('*.mat');
    if saved_analysis == 0,
        m_pref = NaN;
        m_null = NaN;
        Vth_est_pcoll = NaN;
        Vth_est_ncoll = NaN;
        %p_subth_CV = NaN;save_it = 1;
        %n_subth_CV = NaN;
        p_subth_std = NaN;
        n_subth_std = NaN;
        return;
    else
    end
    load(saved_analysis);
    %if overwrite == 0,
    %    d = dir(pwd);
    %    isub = [d(:).isdir];
    %    sub_d = {d(isub).name}';
    %    target = ['rerun_',char(datetime('today'))];
    %    if ~any(strcmp(target,sub_d)),
    %        mkdir(['rerun_',char(datetime('today'))]);
    %        cd(['rerun_',char(datetime('today'))]);
    %    else
    %        cd(['rerun_',char(datetime('today'))]);
    %    end
    %else
    %end
else
    smooth_method = 0;
    spike_prefile = 'analyzed_spike.mat';
    stims = 'stims.mat';
    stimsv = 'stimvalues.mat';
    stimso = 'stimorder.mat';
    data = 'data.mat';
        
    [waveform_dnlp,sample_interval,stimvalues,stimorder,nStim_ON,nStim_OFF,...
        otpref_total,stim_duration,reps_stimmatch,p] = prep_rawspikes(reps,spike_prefile,...
        stims,stimsv,stimso,data);
    
     %if overwrite == 0,
    %    d = dir(pwd);
    %    isub = [d(:).isdir];
    %    sub_d = {d(isub).name}';
    %    target = ['rerun_',char(datetime('today'))];
    %    if ~any(strcmp(target,sub_d)),
    %        mkdir(['rerun_',char(datetime('today'))]);
    %        cd(['rerun_',char(datetime('today'))]);
    %    else
    %        cd(['rerun_',char(datetime('today'))]);
    %    end
    %else
    %end
    %%
    sampling_rate = round(1/sample_interval);
    t_vec = sample_interval:sample_interval:...
        (length(waveform_dnlp)*sample_interval);
    if ~exist('spike_trace'),
        spike_trace = waveform_dnlp;
    else
    end
    
    [spike_locations,h_margins] = find_spikes(spike_trace,t_vec,...
        sampling_rate,display_spikes,auto_detect);
    
    %*******Assess biophysical spike threshold using Azouz & Gray, 1999
    %*******methods****************************************************
    %******************************************************************
        
   [th_loc,Vm_th,new_trace,sub_wv,max_dvdt] = get_thresholds(spike_trace,t_vec,...
       spike_locations,sample_interval,use_detrend,reset2phys);
   
   %%
    %*****Remove spikes to get raw Vm traces and trial-averaged Vm trace**
    %*********************************************************************
    rebin_size = 0.05;                                  %* added 5/29 for bootstrap predictions
    rebin_length = round(rebin_size/sample_interval);   %* added 5/29 for bootstrap predictions
    
    [Vm_waveform,Vm_array,Vm_trial_array,Vm_array_smooth] = extract_Vm(waveform_dnlp,sub_wv,sample_interval,...
        h_margins,stimvalues,nStim_ON,nStim_OFF,stimorder,reps,use_detrend,smooth_Vm,save_it);
    
    [trial_baseline,trial_baseline_array,corr_thloc_array,corr_Vmth_array] =...
        trial_baseline_correction(Vm_waveform,t_vec,sample_interval,...
        stimvalues,stimorder,nStim_ON,nStim_OFF,th_loc,Vm_th);
    
    process_Vm_array = cell(size(Vm_trial_array,1),1);
    sub_trial_array = cell(size(Vm_trial_array,1),size(Vm_trial_array,2));
    if use_trial_baseline == 1,
        for i = 1:size(Vm_trial_array,1)-1,
            for j = 1:size(Vm_trial_array,2),
                range = length(Vm_trial_array{i,j});
                sub_trial_array{i,j} = Vm_trial_array{i,j}-trial_baseline_array(i,j);
                temp_trial_rebin = [];                  %*
                y = 1;                                  %*
                for z = rebin_length:rebin_length:length(sub_trial_array{i,j}),                 %*
                    temp_trial_rebin = [temp_trial_rebin;nanmean(sub_trial_array{i,j}(y:z))];   %*
                    y = z;                                                                      %*
                end
                rebin_sub_trial_array{i,j} = temp_trial_rebin;                                  %*
            end
            process_Vm_trialAve = zeros(range,1);
            for jj = 1:size(Vm_trial_array,2),
                process_Vm_trialAve = process_Vm_trialAve + sub_trial_array{i,jj};
            end
            process_Vm_trialAve = process_Vm_trialAve./size(Vm_trial_array,2);
            if smooth_Vm == 1,
                switch smooth_method
                    case 0
                        [Vm_trialAve_smooth] = Vm_GaussianSmooth(process_Vm_trialAve,...
                            sample_interval);
                    case 1
                        [Vm_trialAve_smooth] = Vm_SGsmooth(process_Vm_trialAve,...
                            sample_interval);
                end
                process_Vm_array{i} = Vm_trialAve_smooth;
            else
                process_Vm_array{i} = process_Vm_trialAve;
            end
        end
    else
        process_Vm_array = Vm_array;
    end
    %get corrected threshold means, by trial and by cell
    repmean_corr_th = cellfun(@mean, corr_Vmth_array);
    trialmean_corr_th = nanmean(repmean_corr_th,2);
    cellmean_corr_th = nanmean(trialmean_corr_th);
    
    
    
    
    %%
   
    stim_rasters = align_rasters(spike_locations,stimvalues,stimorder,...
        nStim_ON,nStim_OFF,sampling_rate);
    
    %%
    %if collate_stims==0,
        %ID 'Pref' and 'Null' stimulus directions
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
    pref_ind = find(stimvalues==pref_stim);
    null_ind = find(stimvalues==null_stim);
    %else
    %end
    
    %%
    if default > 0,
        VF_type = default;
    else
        VF_type = input(['Enter ''1'' for conventional V-F plot, ''2''',...
            ' for maximally sampled/instantaneous FR V-F plot: ']);
    end
    
    switch VF_type
        
        case 1
            
           [f,f_temp] = temporal_conv_VF(process_Vm_array,...
                be_array,stimvalues,sample_interval,psth_array,save_it);
            
        case 2
           if supp_singletrial_FR == 1,
               stack_rasters = 0;
               [FR_inst,opt_w,part_rasters] = rasters2instFR_global3(stimvalues,stimorder,...
                    stim_duration,spike_locations,t_vec,sampling_rate,nStim_ON,...
                    nStim_OFF,reps_stimmatch,save_it,use_global_filter,filter_type,...
                    pref_ind,null_ind,stack_rasters);
                trial_FR_inst = FR_inst;
           else
               [FR_inst,opt_w,part_rasters] = rasters2instFR_global3(stimvalues,stimorder,...
                    stim_duration,spike_locations,t_vec,sampling_rate,nStim_ON,...
                    nStim_OFF,reps_stimmatch,save_it,use_global_filter,filter_type,...
                    pref_ind,null_ind,stack_rasters);
           end
            if ~exist('tempFrequency','var'),
                tempFrequency = p.tFrequency;
            else
            end
            [Vth_est_coll,params,gof_obj,allfits_info] = fitVF_maxsamp_Xbins(FR_inst,process_Vm_array,...
                stimvalues,model,h_margins,sampling_rate,tempFrequency,...
                adapt_bins,anchor_Vth,fit_it,save_it,code,use_trial_baseline,...
                cellmean_corr_th,trialmean_corr_th,pref_ind,null_ind);
            % NOTE - below code added temporarily  5/29/17 only to run
            % binning: not needed for regular VF analysis
            % analysis for bootstrapped predictions
            %*METHOD 1 - from rasters
            %for n = 1:size(part_rasters,2),
            %    for m = 1:size(part_rasters,1),
            %        s_trial_rebin = [];
            %        y = 0;
            %        for z = rebin_length:rebin_length:round(stim_duration(1,1)/sample_interval),
            %            s_trial_rebin = [s_trial_rebin;(length(find(part_rasters{m,n}>y & part_rasters{m,n}<z))/rebin_size)];
            %            y = z;
            %        end
            %        rebin_part_rasters{m,n} = s_trial_rebin;
            %    end
            %end
            %METHOD 2 - from rasters, similar to METHOD 1, but using mean
            %ISI per bin
            %for n = 1:size(part_rasters,2),
            %    for m = 1:size(part_rasters,1),
            %        mi_trial_rebin = [];
            %        for z = rebin_length:rebin_length:round(stim_duration(1,1)/sample_interval),
            %            bin_set = find(part_rasters{m,n}>y & part_rasters{m,n}<z);
            %            if length(bin_set) == 0,
            %                mi_trial_rebin = [mi_trial_rebin;0];
            %            elseif length(bin_set) == 1,
            %                mi_trial_rebin = [mi_trial_rebin;(1/rebin_size)];
            %            elseif length(bin_set) > 1,
            %                mi_trial_rebin = [mi_trial_rebin;1/(nanmean(diff(bin_set.*sample_interval)))];
            %            else
            %            end
            %        end
            %        rebin_part_mISI{m,n} = mi_trial_rebin;
            %    end
            %end
            %METHOD 3 - from single trial FR estimates
            %for nn = 1:size(trial_FR_inst,1),
            %    for mm = 1:size(trial_FR_inst,2),
            %        fr_trial_rebin = [];
            %        yy = 1;
            %        for zz = rebin_length:rebin_length:round(stim_duration(1,1)/sample_interval),
            %            fr_trial_rebin = [fr_trial_rebin;nanmean(trial_FR_inst{nn,mm}(yy:zz,1))];
            %            yy = zz;
            %        end
            %        rebin_trial_FR{nn,mm} = fr_trial_rebin;
            %    end
            %end
            
            
            
            if save_it == 1,
                rbfilename = strcat(code,'_preVF_variables.mat');
                save(rbfilename);
            else
            end
            
            if rebin_trials == 1,
            %    sta_vm_pref = [];
            %    sta_vm_null = [];
                Vth_est_coll = [];
            else
            end
                    
            
    end
    
    if get_VmSTA == 1,
       
        %create collated stimulus indexing that excludes spontaneous trials
        all_dirind = 1:length(stimvalues)-1;
        down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
        up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
        down_i = circshift(all_dirind,1,2);
        up_i = circshift(all_dirind,-1,2);
        pref_set = [down_(pref_ind);pref_stim;up_(pref_ind)];
        null_set = [down_(null_ind);null_stim;up_(null_ind)];
        pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
        null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
        
        cell_spikewindows = struct('pref',{},'null',{});
        use_flankers = 1; %might take this out; if not add to higher level params.
        if use_flankers == 1,
            %compute spike_triggered average Vm (PREF)
            sta_win = -round(0.05/sample_interval):1:round(0.05/sample_interval);
            sta_vm_pre = zeros(length(sta_win),1);
            spike_count = 0;
            for ii = 1:length(pref_set_index),
                pstim_loc = pref_set_index(ii,1);
                abstim_loc = find(stimorder==pstim_loc);
                for iii = 1:reps_stimmatch,
                    for s = 1:length(corr_thloc_array{pstim_loc,iii}),
                        spike_count = spike_count + 1;
                        st_Vmind = corr_thloc_array{pstim_loc,iii}(s,1);
                        sub_ind = st_Vmind-round(nStim_ON(abstim_loc(1,iii),1)/sample_interval);
                        st_Vmrange = [(sub_ind-round(0.05/sample_interval)):1:...
                            (sub_ind+round(0.05/sample_interval))];
                        if any(st_Vmrange<0),
                            spike_count = spike_count-1;
                            continue;
                        else
                        end
                        st_Vmsamp = sub_trial_array{pstim_loc,iii}(st_Vmrange);
                        cell_spikewindows.pref{end+1,1} = st_Vmsamp;
                        sta_vm_pre = sta_vm_pre + st_Vmsamp;
                    end
                end
            end
            sta_vm_pref = sta_vm_pre./spike_count;
            t_win = -0.05:(0.05/round(0.05/sample_interval)):0.05;
            f30 = figure;
            plot(t_win,sta_vm_pref,'r');
            xlabel('Time (sec)');
            ylabel('Membrane potential (mV)');
            hold on;
            
            %compute spike-triggered average Vm (NULL)
            sta_vm_pre = zeros(length(sta_win),1);
            spike_count_ = 0;
            for jj = 1:length(null_set_index),
                nstim_loc = null_set_index(jj,1);
                abstim_loc = find(stimorder==nstim_loc);
                for jjj = 1:reps_stimmatch,
                    for s = 1:length(corr_thloc_array{nstim_loc,jjj}),
                        spike_count_ = spike_count_ + 1;
                        st_Vmind = corr_thloc_array{nstim_loc,jjj}(s,1);
                        sub_ind_ = st_Vmind-round(nStim_ON(abstim_loc(1,jjj),1)/sample_interval);
                        st_Vmrange = [(sub_ind_-round(0.05/sample_interval)):1:...
                            (sub_ind_+round(0.05/sample_interval))];
                        if any(st_Vmrange<0),
                            spike_count = spike_count-1;
                            continue;
                        else
                        end
                        st_Vmsamp = sub_trial_array{nstim_loc,jjj}(st_Vmrange);
                        cell_spikewindows.null{end+1,1} = st_Vmsamp;
                        sta_vm_pre = sta_vm_pre + st_Vmsamp;
                    end
                end
            end
            sta_vm_null = sta_vm_pre./spike_count_;
            plot(t_win,sta_vm_null,'b');
            legend('Pref','Null','Location','NorthWest');
            if save_it == 1,
                saveas(gcf,[pwd filesep 'baselinecorr_STAVm.fig']);
                close(f30);
                save('bc10_STA_variables.mat','sta_vm_pref','sta_vm_null',...
                    't_win','sample_interval');
                save('bc10_cell_spikewindows.mat','cell_spikewindows');
            else
            end
        else
            %...under construction...
            
        end
        
    else
        sta_vm_pref = [];
        sta_vm_null = [];
    end
    
end
        
    cd ..
end



