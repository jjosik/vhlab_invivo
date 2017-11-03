function [ FR_inst,opt_w,part_rasters ] = ...
    rasters2instFR_global3( stimvalues,stimorder,stim_duration,spike_locations,...
    t_vec,sampling_rate,nStim_ON,nStim_OFF,stimmatch_reps,save_it,...
    use_global_filter,filter_type,pref_ind,null_ind,stack_rasters )

%DESCRIPTION - this version reduces the number of options down to 2 (also removes 
%some additional fluff):
%   1. a basic gaussian kernel smoother with window of pre-determined size
%   2. the computed optimal gaussian filter (still uses global option based
%   on collating side-by-side rather than stacking)
% ********NOTE rasters2instFR_global4.m now has the exact functionality that
% global3 had before stack_rasters option was added (used in all analysis up to ~4/12/17).
%INPUTS - 
 
%OUTPUTS -

sample_interval = 1/sampling_rate;

%*******Stack spiketrains by stim (needed for both methods following)**
    stim_rasters = cell(length(stimvalues),1);
    part_raster = cell(stimmatch_reps,length(stimvalues));
    for b = 1:length(stimvalues),
        current_epochs = find(stimorder == b);
        align_on_ind = nStim_ON(current_epochs).*sampling_rate;
        align_off_ind = nStim_OFF(current_epochs).*sampling_rate;
        for a=1:length(current_epochs),
            rasters{a,1} = spike_locations(find(spike_locations>align_on_ind(a,1)&spike_locations<align_off_ind(a,1)))-(nStim_ON(current_epochs(1,a))*sampling_rate);
            part_rasters{a,b} = rasters{a,1};
        end
        new_rasters = [];
        for a=1:length(current_epochs),
            new_rasters = [new_rasters;rasters{a,1}];
        end
        stim_rasters{b,1} = new_rasters;
    end
    if use_global_filter == 1,
        p_epochs = find(stimorder == pref_ind);
        align_off_ind = nStim_OFF(p_epochs).*sampling_rate;
        align_on_ind = nStim_ON(p_epochs).*sampling_rate;
        n_offset = mean(align_off_ind-align_on_ind);
        pn_rasters = [stim_rasters{pref_ind,1}(:);(n_offset+stim_rasters{null_ind,1}(:))];
    else
    end
    
    switch filter_type
    
        case 1      %flat window size for gaussian smoothing
            
            switch stack_rasters
                
                case 1
                    
                    opt_w = 0.05        
                    stacked_t_vec = sample_interval:sample_interval:stim_duration(1,1);
                    Xn = [];
                        for qq = 2:length(stacked_t_vec),
                            Xn(end+1) = mean([stacked_t_vec(1,qq) stacked_t_vec(1,qq-1)]);
                        end
                        stacked_t_index = stacked_t_vec/sample_interval;
                        FR_inst = cell(length(stim_rasters),1);
                        for bb=1:length(stim_rasters),
                            display(['Computing firing rate estimates for stim. orientation ',num2str(stimvalues(1,bb)),'...']);
                            Yn = zeros(length(stacked_t_vec)-1,1);
                            if ~isempty(stim_rasters{bb,1}),
                                for q = 1:length(stim_rasters{bb,1}),
                                    range = length(find(stacked_t_index<stim_rasters{bb,1}(q,1)));
                                    if range <= 1,
                                        try
                                            Yn(length(find(stacked_t_index<stim_rasters{bb,1}(q,1))),1) = 1;
                                        catch
                                            Yn(length(find(stacked_t_index<stim_rasters{bb,1}(q,1))+1),1) = 1;
                                        end
                                    else
                                        Yn((length(find(stacked_t_index<stim_rasters{bb,1}(q,1)))-1),1) = 1;
                                    end
                                end
                                g_kernel = gausswin(round(opt_w*sampling_rate));
                                norm_g_kernel = g_kernel./sum(g_kernel);
                                FR_inst{bb,1} = (1/(sample_interval*stimmatch_reps)).*conv(Yn,norm_g_kernel,'same');
                        %ALT. METHOD - DO NOT USE/NOT FULLY FUNCTIONAL - Correcting
                        %for shift instead
                        %***NOTE: methods for use with convolution function yield
                        %reasonable output but it is shifted in time with respect
                        %to input (a fatal error in this analysis); use equivalence
                        %between multiplication in frequency domain and convolution
                        %to avoid the shift as follows
                        %Yn = resample(Yn,1,10);
                        %len = length(Yn);
                        %pad_len = max(1:len+3*(opt_w(bb,1)/sample_interval));
                        %int = 2^nextpow2(pad_len);
                        %nY = fft(Yn,int);
                        %f = (0:int-1)/int;
                        %f = [-f(1:int/2+1) f(int/2:-1:2)];
                        %fd_g_kernel = exp(-(1/2)*(2*pi*f*opt_w(bb,1)).^2);
                        %FR_inst{bb,1} = ifft(nY*fd_g_kernel,int);
                        %END OF ALT. METHOD
                            else
                                FR_inst{bb,1} = zeros(length(stacked_t_vec)-1,1);
                            end
                    %plot stacked rasters and inst. FR trace
                    %FOR ALT METHOD ONLY: down_t = resample(stacked_t_vec,1,10);
                            tick_ht = 1.0;
                            f2 = figure;
                            ax_upper = subplot(2,1,1);
                            for c = 1:stimmatch_reps,
                                for cc = 1:length(part_rasters{c,bb}),
                                    disp_raster = line([part_rasters{c,bb}(cc,1) part_rasters{c,bb}(cc,1)],[(tick_ht*c)-tick_ht,(tick_ht*c)]);
                                    hold on;
                                    set(disp_raster,'Color',[0 0 0],'LineWidth',1.5);
                                end
                                ylim([-0.1 (tick_ht*(c+1))]);
                                xlim([0 length(stacked_t_vec)-1]);
                                xlabel('Timestep');
                                title(['Raster and inst. firing rate for angle ',num2str(stimvalues(1,c))]);
                                hold off;
                            end
                            ax_lower = subplot(2,1,2);
                    %FOR ALT METHOD ONLY: plot(down_t,FR_inst{bb,1},'k');
                            plot(stacked_t_vec(1,1:end-1),FR_inst{bb,1},'k');
                            hold on;
                            xlabel('Time (seconds)');
                            ylabel('Inst. firing rate');
                            hold off;
                            if save_it == 1,
                                saveas(gcf,[pwd filesep 'smoothed_inst_FR_rasters_',...
                                    num2str(stimvalues(1,bb)),'.fig']);
                                close(f2);
                            else
                            end
                        end
                        
                case 0
                    
                    opt_w = 0.05        
                    train_t_vec = sample_interval:sample_interval:stim_duration(1,1);
                    Xn = [];
                    for qq = 2:length(train_t_vec),
                        Xn(end+1) = mean([train_t_vec(1,qq) train_t_vec(1,qq-1)]);
                    end
                    train_t_index = train_t_vec/sample_interval;
                    FR_inst = cell(size(part_rasters,2),size(part_rasters,1));
                    for bb=1:size(part_rasters,2),
                        for cc = 1:size(part_rasters,1),
                            display(['Computing firing rate estimates for stim. orientation ',num2str(stimvalues(1,bb)),', rep. #',num2str(cc),'...']);
                            Yn = zeros(length(train_t_vec)-1,1);
                            if ~isempty(part_rasters{cc,bb}),
                                for q = 1:length(part_rasters{cc,bb}),
                                    range = length(find(train_t_index<part_rasters{cc,bb}(q,1)));
                                    if range <= 1,
                                        try
                                            Yn(length(find(train_t_index<part_rasters{cc,bb}(q,1))),1) = 1;
                                        catch
                                            Yn(length(find(train_t_index<part_rasters{cc,bb}(q,1))+1),1) = 1;
                                        end
                                    else
                                        Yn((length(find(train_t_index<part_rasters{cc,bb}(q,1)))-1),1) = 1;
                                    end
                                end
                                g_kernel = gausswin(round(opt_w*sampling_rate));
                                norm_g_kernel = g_kernel./sum(g_kernel);
                                FR_inst{bb,cc} = (1/sample_interval).*conv(Yn,norm_g_kernel,'same');
                            else
                                FR_inst{bb,cc} = zeros(length(train_t_vec)-1,1);
                            end
                            %tick_ht = 1.0;
                            %f2 = figure;
                            %for c = 1:stimmatch_reps,
                            %    ax_upper = subplot(2,3,c);
                            %    for cc = 1:length(part_rasters{c,bb}),
                            %        disp_raster = line([part_rasters{c,bb}(cc,1) part_rasters{c,bb}(cc,1)],[0 tick_ht]);
                            %        hold on;
                            %        set(disp_raster,'Color',[0 0 0],'LineWidth',1.5);
                            %        ylim([-0.1 (tick_ht*2)]);
                            %        xlim([0 length(train_t_vec)-1]);
                            %        xlabel('Timestep');
                            %        title(['Raster and inst. firing rate for angle ',num2str(stimvalues(1,bb)),', rep. ',num2str(c)]);
                            %        hold off;
                            %    end
                            %    ax_lower = subplot(2,3,c+stimmatch_reps);
                            %    plot(train_t_vec(1,1:end-1),FR_inst{bb,c},'k');
                            %    hold on;
                            %    xlabel('Time (seconds)');
                            %    ylabel('Inst. firing rate');
                            %    hold off;
                            %end
                            %if save_it == 1,
                            %    saveas(gcf,[pwd filesep 'unstacked_flatwin_inst_FR_rasters_',...
                            %        num2str(stimvalues(1,bb)),'.fig']);
                            %    close(f2);
                            %else
                            %end
                            %block above removed 
                        end
                    end
                    
            end
                            
                    
            
               
        case 2      %Bandwidth-optimized Gaussian kernel with moving average; can be switched to less 
                %optimal global kernel

        %**********************
        %*Type 2: Band-optimized kernel method (fixed bandwidth) (algorithm from Shimazaki & Shinomoto, J. Comp. Neurosci, 2010) 
        %*Here we use the closed-form analytic result for gaussian kernel
        %windows together with brute-force parameter search.  
        %***A global filter size is extracted by collating end-to-end in time the stacked
        %trial rasters for 'pref' and 'null' conditions. (in effect when
        %use_global_filter = 1).
        w_vec = 0.005:0.001:0.5;
        opt_w = zeros(length(stim_rasters),1);
        switch use_global_filter
            
            case 0
                
                for aa=1:length(stim_rasters),
                    display(['Computing cost function for stim. orientation ',num2str(stimvalues(1,aa)),'...']);
                    %if isempty(stim_rasters{aa,1},
                    %    continue;
                    %else
                    for n=1:length(w_vec),
                        gauss_width = round(w_vec(1,n)*sampling_rate);
                        inner_sum = 0;
                        for j = 2:length(stim_rasters{aa,1}),
                            for iii = 1:j-1
                                inner_sum = inner_sum + (exp(-(((1/sampling_rate)*stim_rasters{aa,1}(iii,1)-(1/sampling_rate)*stim_rasters{aa,1}(j,1))^2)./(4*w_vec(1,n)^2)))-...
                                    ((2*sqrt(2)).*(exp(-(((1/sampling_rate)*stim_rasters{aa,1}(iii,1)-(1/sampling_rate)*stim_rasters{aa,1}(j,1))^2)./(2*w_vec(1,n)^2))));
                            end
                        end
                        C_n(n,1) = (stimmatch_reps/w_vec(1,n))+(2/w_vec(1,n)).*(inner_sum);  %cost function
                    end
                    [l_cost,n_ind] = min(C_n);      %find global minimum of cost function
                    opt_w(aa,1) = w_vec(1,n_ind);
                    ff = figure;
                    plot(w_vec(1,:),C_n(:,1),'k');
                    hold on;
                    xlim([0 w_vec(1,end)]);
                    ylim([min(C_n)-0.15*abs(min(C_n)) max(C_n)+0.15*abs(max(C_n))]);
                    h_xm = line([opt_w(aa,1) opt_w(aa,1)],[min(C_n)-0.15*abs(min(C_n)) max(C_n)+0.15*abs(max(C_n))]);
                    h_ym = line([0 w_vec(1,end)],[l_cost l_cost]);
                    h_xm.LineStyle = '--';
                    h_ym.LineStyle = '--';
                    title(['Cost function and optimal gaussian window width for ',num2str(stimvalues(1,aa)),' degree stim.']);
                    if save_it == 1,
                        saveas(gcf,[pwd filesep 'gausswin_costfunction_',num2str(stimvalues(1,aa)),'.fig']);
                        close(ff);
                    else
                    end
                end
                
                
            case 1
                
                opt_w = [];
                display('Computing global cost function for combined pref/null stim. orientations ...');
                for n=1:length(w_vec),
                        gauss_width = round(w_vec(1,n)*sampling_rate);
                        inner_sum = 0;
                        for j = 2:length(pn_rasters),
                            for iii = 1:j-1
                                inner_sum = inner_sum + (exp(-(((1/sampling_rate)*pn_rasters(iii,1)-(1/sampling_rate)*pn_rasters(j,1))^2)./(4*w_vec(1,n)^2)))-...
                                    ((2*sqrt(2)).*(exp(-(((1/sampling_rate)*pn_rasters(iii,1)-(1/sampling_rate)*pn_rasters(j,1))^2)./(2*w_vec(1,n)^2))));
                            end
                        end
                        C_n(n,1) = (stimmatch_reps/w_vec(1,n))+(2/w_vec(1,n)).*(inner_sum);  %cost function
                end
                [l_cost,n_ind] = min(C_n);      %find global minimum of cost function
                opt_w = w_vec(1,n_ind);
                %opt_w = 0.05;       %***IMPORTANT -- TEMP CHANGE FOR TESTING
                g_opt_w = opt_w;    %for passing back into general function
                ff = figure;
                plot(w_vec(1,:),C_n(:,1),'k');
                hold on;
                xlim([0 w_vec(1,end)]);
                ylim([min(C_n)-0.15*abs(min(C_n)) max(C_n)+0.15*abs(max(C_n))]);
                h_xm = line([opt_w opt_w],[min(C_n)-0.15*abs(min(C_n)) max(C_n)+0.15*abs(max(C_n))]);
                h_ym = line([0 w_vec(1,end)],[l_cost l_cost]);
                h_xm.LineStyle = '--';
                h_ym.LineStyle = '--';
                title(['Cost function and optimal gaussian window width for combined pref/null stims.']);
                if save_it == 1,
                    saveas(gcf,[pwd filesep '50_gausswin_costfunction_global.fig']);
                    close(ff);
                else
                end
                
        end
                
                
        stacked_t_vec = sample_interval:sample_interval:stim_duration(1,1);
        Xn = [];
        for qq = 2:length(stacked_t_vec),
            Xn(end+1) = mean([stacked_t_vec(1,qq) stacked_t_vec(1,qq-1)]);
        end
        stacked_t_index = stacked_t_vec/sample_interval;
        FR_inst = cell(length(stim_rasters),1);
        for bb=1:length(stim_rasters),
            display(['Computing firing rate estimates for stim. orientation ',num2str(stimvalues(1,bb)),'...']);
            Yn = zeros(length(stacked_t_vec)-1,1);
            if ~isempty(stim_rasters{bb,1}),
                for q = 1:length(stim_rasters{bb,1}),
                    range = length(find(stacked_t_index<stim_rasters{bb,1}(q,1)));
                    if range <= 1,
                        try
                            Yn(length(find(stacked_t_index<stim_rasters{bb,1}(q,1))),1) = 1;
                        catch
                            Yn(length(find(stacked_t_index<stim_rasters{bb,1}(q,1))+1),1) = 1;
                        end
                    else
                        Yn((length(find(stacked_t_index<stim_rasters{bb,1}(q,1)))-1),1) = 1;
                    end
                end
                g_kernel = gausswin(round(opt_w(bb,1)*sampling_rate));
                norm_g_kernel = g_kernel./sum(g_kernel);
                FR_inst{bb,1} = (1/(sample_interval*stimmatch_reps)).*conv(Yn,norm_g_kernel,'same');
                        %ALT. METHOD - DO NOT USE/NOT FULLY FUNCTIONAL - Correcting
                        %for shift instead
                        %***NOTE: methods for use with convolution function yield
                        %reasonable output but it is shifted in time with respect
                        %to input (a fatal error in this analysis); use equivalence
                        %between multiplication in frequency domain and convolution
                        %to avoid the shift as follows
                        %Yn = resample(Yn,1,10);
                        %len = length(Yn);
                        %pad_len = max(1:len+3*(opt_w(bb,1)/sample_interval));
                        %int = 2^nextpow2(pad_len);
                        %nY = fft(Yn,int);
                        %f = (0:int-1)/int;
                        %f = [-f(1:int/2+1) f(int/2:-1:2)];
                        %fd_g_kernel = exp(-(1/2)*(2*pi*f*opt_w(bb,1)).^2);
                        %FR_inst{bb,1} = ifft(nY*fd_g_kernel,int);
                        %END OF ALT. METHOD
            else
                FR_inst{bb,1} = zeros(length(stacked_t_vec)-1,1);
            end
            %plot stacked rasters and inst. FR trace
            %FOR ALT METHOD ONLY: down_t = resample(stacked_t_vec,1,10);
            tick_ht = 1.0;
            f6 = figure;
            ax_upper = subplot(2,1,1);
            for c = 1:stimmatch_reps,
                for cc = 1:length(part_rasters{c,bb}),
                    disp_raster = line([part_rasters{c,bb}(cc,1) part_rasters{c,bb}(cc,1)],[(tick_ht*c)-tick_ht,(tick_ht*c)]);
                    hold on;
                    set(disp_raster,'Color',[0 0 0],'LineWidth',1.5);
                end
                ylim([-0.1 (tick_ht*(c+1))]);
                xlim([0 length(stacked_t_vec)-1]);
                xlabel('Timestep');
                title(['Raster and inst. firing rate for angle ',num2str(stimvalues(1,c))]);
                hold off;
            end
            ax_lower = subplot(2,1,2);
            %FOR ALT METHOD ONLY: plot(down_t,FR_inst{bb,1},'k');
            plot(stacked_t_vec(1,1:end-1),FR_inst{bb,1},'k');
            hold on;
            xlabel('Time (seconds)');
            ylabel('Inst. firing rate');
            hold off;
            if save_it == 1,
                if use_global_filter == 1,
                    saveas(gcf,[pwd filesep 'globalFilt_instFR_rasters_',...
                        num2str(stimvalues(1,bb)),'.fig']);
                else
                    saveas(gcf,[pwd filesep 'inst_FR_rasters_',...
                        num2str(stimvalues(1,bb)),'.fig']);
                end
                close(f6);
            else
            end
        end
        
    end
    
end

