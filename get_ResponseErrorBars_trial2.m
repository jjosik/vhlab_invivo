function [ output_args ] = get_ResponseErrorBars_trial2( allow_flip )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
% ALLOW_FLIP should be set to '1' if simulation results are permitted to
% switch response categories, i.e. if PREF/NULL is an absolute designator for
% Dir1/Dir2, then allow_flip = 0; if PREF/NULL is relative in the usual sense
% then allow_flip = 1.
load('collected_params.mat');

[folder,directories] = walk_directories();
cd(folder);

o_frdist_p_meas = {};
o_frdist_p_pred = {};
y_frdist_p_meas = {};
y_frdist_p_meas = {};
o_frmean_p_meas = [];
o_frmean_p_pred = [];
o_frse_p_meas = [];
o_frse_p_pred = [];
y_frmean_p_meas = [];
y_frmean_p_pred = [];
y_frse_p_meas = [];
y_frse_p_pred = [];
o_frdist_n_meas = {};
o_frdist_n_pred = {};
y_frdist_n_meas = {};
y_frdist_n_pred = {};
y_frdist_n_meas = {};
y_frdist_p_pred = {};
o_frmean_n_meas = [];
o_frmean_n_pred = [];
o_frse_n_meas = [];
o_frse_n_pred = [];
y_frmean_n_meas = [];
y_frmean_n_pred = [];
y_frse_n_meas = [];
y_frse_n_pred = [];
o_mean_actFR_dsi = [];
y_mean_actFR_dsi = [];
o_se_actFR_dsi = [];
y_se_actFR_dsi = [];
o_samp_actFR_dsi = [];
y_samp_actFR_dsi = [];
o_samp_modelFR_dsi = [];
y_samp_modelFR_dsi = [];



for i = 1:size(directories,1),
    cd(strtrim(directories(i,:)));
    sub_folder = uigetdir();
    d = dir(sub_folder);
    isub = [d(:).isdir];
    sel_Snames = {d(isub).name}';
    sel_Snames(ismember(sel_Snames,{'.','..','Codes'})) = [];
    [s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',sel_Snames);
    sub_dir = char(sel_Snames(s));
    [~,subfile,~] = fileparts(sub_folder);
    for j = 1:length(sub_dir),
        cur_dir = sub_dir(j,:);
        cd(cur_dir);
        load('analyzed_Vm.mat','total_vm_byAngle','stimvalues');
        load('analyzed_spike.mat','total_spike_byAngle');
        %load('na10M1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        load('GsmoothM1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        load('noflankers_alltrial_variables.mat','opref_cell_FR_alltrials','onull_cell_FR_alltrials',...
            'ypref_cell_FR_alltrials','ynull_cell_FR_alltrials','opref_cell_Vm_alltrials',...
            'onull_cell_Vm_alltrials','ypref_cell_Vm_alltrials','ynull_cell_Vm_alltrials',...
            'old_cellnames','young_cellnames');
        if strcmp(subfile,'old'),
            pref_Vm = opref_cell_Vm_alltrials(j,1);
            pref_spike = opref_cell_FR_alltrials(j,1);
            null_Vm = onull_cell_Vm_alltrials(j,1);
            null_spike = onull_cell_FR_alltrials(j,1);
        else
            pref_Vm = ypref_cell_Vm_alltrials(j,1);
            pref_spike = ypref_cell_FR_alltrials(j,1);
            null_Vm = ynull_cell_Vm_alltrials(j,1);
            null_spike = ynull_cell_FR_alltrials(j,1);
            
        end
        nsims = 1000;
        
        samp_dist_n_meas = [];
        samp_dist_n_pred = [];
        o_samp_err_mps = [];
        o_samp_err_mns = [];
        y_samp_err_mps = [];
        y_samp_err_mns = [];
        samp_spike_dsi = [];
        % vm
        p_running_ave_vm = zeros(length(pref_Vm{1,1}(1,1)),1);
        n_running_ave_vm = zeros(length(null_Vm{1,1}(1,1)),1);
        for n = 1:nsims,
            %if allow_flip == 1,
                samp_pref_vm = datasample(pref_Vm,size(pref_Vm,1)*size(pref_Vm,2));
                linearize_psamp = cell2mat(reshape(samp_pref_vm{1,1},1,size(samp_pref_vm{1,1},1)*size(samp_pref_vm{1,1},2)));
                p_running_ave_vm = p_running_ave_vm + sum(linearize_psamp,2);
                samp_null_vm = datasample(null_Vm,size(null_Vm,1)*size(null_Vm,2));
                linearize_nsamp = cell2mat(reshape(samp_null_vm{1,1},1,size(samp_null_vm{1,1},1)*size(samp_null_vm{1,1},2)));
                n_running_ave_vm = n_running_ave_vm + sum(linearize_nsamp,2);
                %apply to model
                
            %else
            %end
        end
        p_sim_ave_vm = n_running_ave_vm./(nsims*size(linearize_psamp,2));
        n_sim_ave_vm = p_running_ave_vm./(nsims*size(linearize_nsamp,2));
        %apply to model
        if strcmp(subfile,'old'),
            o_samp_pref_modelFR = rectify(old_params{j,1}.tanh_params.a +...
                old_params{j,1}.tanh_params.b.*tanh((p_sim_ave_vm-old_params{j,1}.tanh_params.c)./...
                old_params{j,1}.tanh_params.d));
            o_samp_null_modelFR = rectify(old_params{j,1}.tanh_params.a +...
                old_params{j,1}.tanh_params.b.*tanh((n_sim_ave_vm-old_params{j,1}.tanh_params.c)./...
                old_params{j,1}.tanh_params.d));
            o_samp_modelFR_dsi = [o_samp_modelFR_dsi;abs(mean(o_samp_pref_modelFR)-...
                        mean(o_samp_null_modelFR))/max([mean(o_samp_pref_modelFR);mean(o_samp_null_modelFR)])];
        else
            y_samp_pref_modelFR = rectify(young_params{j,1}.tanh_params.a +...
                young_params{j,1}.tanh_params.b.*tanh((p_sim_ave_vm-young_params{j,1}.tanh_params.c)./...
                young_params{j,1}.tanh_params.d));
            y_samp_null_modelFR = rectify(young_params{j,1}.tanh_params.a +...
                young_params{j,1}.tanh_params.b.*tanh((n_sim_ave_vm-young_params{j,1}.tanh_params.c)./...
                young_params{j,1}.tanh_params.d));
            y_samp_modelFR_dsi = [y_samp_modelFR_dsi;abs(mean(y_samp_pref_modelFR)-...
                mean(y_samp_null_modelFR))/max([mean(y_samp_pref_modelFR);mean(y_samp_null_modelFR)])];
        end
        % measured fr
        p_running_ave_fr = zeros(length(pref_spike{1,1}(1,1)),1);
        n_running_ave_fr = zeros(length(null_spike{1,1}(1,1)),1);
        
        for n = 1:nsims,
            if allow_flip == 1,
                samp_pref_FR = datasample(pref_spike,size(pref_spike,1)*size(pref_spike,2));
                linFR_psamp = cell2mat(reshape(samp_pref_FR{1,1},1,size(samp_pref_FR{1,1},1)*size(samp_pref_FR{1,1},2)));
                p_running_ave_fr = p_running_ave_fr + sum(linFR_psamp,2);
                samp_null_FR = datasample(null_spike,size(null_spike,1)*size(null_spike,2));
                linFR_nsamp = cell2mat(reshape(samp_null_FR{1,1},1,size(samp_null_FR{1,1},1)*size(samp_null_FR{1,1},2)));
                n_running_ave_fr = n_running_ave_fr + sum(linFR_nsamp,2);
                %get trial actual DSI
                if strcmp(subfile,'old'),
                    int_sim_pref_FR = sum(linFR_psamp,2)/size(linFR_psamp,2);
                    int_sim_null_FR = sum(linFR_nsamp,2)/size(linFR_nsamp,2);
                    o_samp_actFR_dsi = [o_samp_actFR_dsi;abs(mean(int_sim_pref_FR)-...
                        mean(int_sim_null_FR))/max([mean(int_sim_pref_FR);mean(int_sim_null_FR)])];
                else
                    int_sim_pref_FR = sum(linFR_psamp,2)/size(linFR_psamp,2);
                    int_sim_null_FR = sum(linFR_nsamp,2)/size(linFR_nsamp,2);
                    y_samp_actFR_dsi = [y_samp_actFR_dsi;abs(mean(int_sim_pref_FR)-...
                        mean(int_sim_null_FR))/max([mean(int_sim_pref_FR);mean(int_sim_null_FR)])];
                end
            else
            end
        end
               
        p_sim_ave_fr = p_running_ave_fr/(nsims*size(linFR_psamp,2));
        n_sim_ave_fr = n_running_ave_fr/(nsims*size(linFR_nsamp,2));
        if strcmp(subfile,'old'),
            o_frdist_p_meas{end+1,1} = p_sim_ave_fr;
            o_frmean_p_meas(end+1,1) = mean(p_sim_ave_fr);
            o_frse_p_meas(end+1,1) = std(p_sim_ave_fr)/(sqrt(length(p_sim_ave_fr)));
            o_frdist_p_pred{end+1,1} = o_samp_pref_modelFR;
            o_frmean_p_pred(end+1,1) = mean(o_samp_pref_modelFR);
            o_frse_p_pred(end+1,1) = std(o_samp_pref_modelFR)/(sqrt(length(o_samp_pref_modelFR)));
            o_frdist_n_meas{end+1,1} = n_sim_ave_fr;
            o_frmean_n_meas(end+1,1) = mean(n_sim_ave_fr);
            o_frse_n_meas(end+1,1) = std(n_sim_ave_fr)/(sqrt(length(n_sim_ave_fr)));
            o_frdist_n_pred{end+1,1} = o_samp_null_modelFR;
            o_frmean_n_pred(end+1,1) = mean(o_samp_null_modelFR);
            o_frse_n_pred(end+1,1) = std(o_samp_null_modelFR)/(sqrt(length(o_samp_null_modelFR)));
            o_mean_actFR_dsi(end+1,1) = mean(o_samp_actFR_dsi);
            o_se_actFR_dsi(end+1,1) = std(o_samp_actFR_dsi)/(sqrt(length(o_samp_actFR_dsi)));
        else
            y_frdist_p_meas{end+1,1} = p_sim_ave_fr;
            y_frmean_p_meas(end+1,1) = mean(p_sim_ave_fr);
            y_frse_p_meas(end+1,1) = std(p_sim_ave_fr)/(sqrt(length(p_sim_ave_fr)));
            y_frdist_p_pred{end+1,1} = y_samp_pref_modelFR;
            y_frmean_p_pred(end+1,1) = mean(y_samp_pref_modelFR);
            y_frse_p_pred(end+1,1) = std(y_samp_pref_modelFR)/(sqrt(length(y_samp_pref_modelFR)));
            y_frdist_n_meas{end+1,1} = n_sim_ave_fr;
            y_frmean_n_meas(end+1,1) = mean(n_sim_ave_fr);
            y_frse_n_meas(end+1,1) = std(n_sim_ave_fr)/(sqrt(length(n_sim_ave_fr)));
            y_frdist_n_pred{end+1,1} = y_samp_null_modelFR;
            y_frmean_n_pred(end+1,1) = mean(y_samp_null_modelFR);
            y_frse_n_pred(end+1,1) = std(y_samp_null_modelFR)/(sqrt(length(y_samp_null_modelFR)));
            y_mean_actFR_dsi(end+1,1) = mean(y_samp_actFR_dsi);
            y_se_actFR_dsi(end+1,1) = std(y_samp_actFR_dsi)/(sqrt(length(y_samp_actFR_dsi)));
        end
        
        clear pref_ind null_ind;
        cd ../
    end
    cd ../
end

%alt. actual spike DSI
load('variables.mat','spike_dsi_old','spike_dsi_young');
spike_dsi_old = spike_dsi_old';
spike_dsi_young = spike_dsi_young';
                
  
%Predicted vs. Measured FR Response
    f = figure;
    subplot(1,2,1);
    scatter(y_frmean_n_meas,y_frmean_n_pred,'b','filled');
    hold on;
    scatter(y_frmean_p_meas,y_frmean_p_pred,'r','filled');
    for m = 1:length(y_frmean_n_meas),
        a = line([y_frmean_n_meas(m)-y_frse_n_meas(m) y_frmean_n_meas(m)+y_frse_n_meas(m)],...
            [y_frmean_n_pred(m) y_frmean_n_pred(m)]);
        a.Color = [0 0 1];
        a.LineWidth = 2.0;
        e = line([y_frmean_n_meas(m) y_frmean_n_meas(m)],[rectify(y_frmean_n_pred(m)-y_frse_n_pred(m)) y_frmean_n_pred(m)+y_frse_n_pred(m)]);
        e.Color = [0 0 1];
        e.LineWidth = 2.0;
    end
    for m = 1:length(y_frmean_p_meas),
        a = line([y_frmean_p_meas(m)-y_frse_p_meas(m) y_frmean_p_meas(m)+y_frse_p_meas(m)],...
            [y_frmean_p_pred(m) y_frmean_p_pred(m)]);
        a.Color = [1 0 0];
        a.LineWidth = 2.0;
        e = line([y_frmean_p_meas(m) y_frmean_p_meas(m)],[rectify(y_frmean_p_pred(m)-y_frse_p_pred(m)) y_frmean_p_pred(m)+y_frse_p_pred(m)]);
        e.Color = [1 0 0];
        e.LineWidth = 2.0;
    end
    fullset = [y_frmean_n_meas+y_frse_n_meas;y_frmean_p_meas+y_frse_p_meas;y_frmean_n_pred+y_frse_n_pred;y_frmean_p_pred+y_frse_p_pred];
    unl = line([0 max(fullset)],[0 max(fullset)]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Measured Response (Hz)');
    ylabel('Predicted Response (Hz)');
    try
        [nb,nb_ci,nresid,nrint,nstats] = ...
            regress(y_frmean_n_pred,horzcat(ones(length(y_frmean_n_meas),1),y_frmean_n_meas));
    catch
        [nslope,noffset,nslope_ci,nresid,nrint,nstats] = ...
            quickregression(y_frmean_n_meas,y_frmean_n_pred,0.05);
    end
    [nslope,noffset,nslope_ci,nresid,nrint,qnstats] = ...
            quickregression(reshape(y_frmean_n_meas,size(y_frmean_n_pred)),y_frmean_n_pred,0.05);
    try
        [pb,pb_ci,presid,print,pstats] = ...
            regress(y_frmean_p_pred,horzcat(ones(length(y_frmean_p_meas),1),y_frmean_p_meas));
    catch
        [pslope,poffset,pslope_ci,presid,print,pstats] = ...
            quickregression(y_frmean_p_meas,y_frmean_p_pred,0.05);
    end
    [pslope,poffset,pslope_ci,presid,print,qpstats] = ...
            quickregression(reshape(y_frmean_p_meas,size(y_frmean_n_pred)),y_frmean_n_pred,0.05);
    legend(['Null, R^2 = ',num2str(nstats(1,1)),', p = ',num2str(nstats(1,3)),', Slope CI: [',num2str(nslope_ci(1,1)),' ',num2str(nslope_ci(1,2)),']'],...
    ['Pref, R^2 = ',num2str(pstats(1,1)),', p = ',num2str(pstats(1,3)),', Slope CI: [',num2str(pslope_ci(1,1)),' ',num2str(pslope_ci(1,2)),']'],...
    'Location','NorthWest');
    title('YOUNG');
    hold off;
    subplot(1,2,2);
    scatter(o_frmean_n_meas,o_frmean_n_pred,'b','filled');
    hold on;
    scatter(o_frmean_p_meas,o_frmean_p_pred,'r','filled');
    for m = 1:length(o_frmean_n_meas),
        a = line([o_frmean_n_meas(m)-o_frse_n_meas(m) o_frmean_n_meas(m)+o_frse_n_meas(m)],...
            [o_frmean_n_pred(m) o_frmean_n_pred(m)]);
        a.Color = [0 0 1];
        a.LineWidth = 2.0;
        e = line([o_frmean_n_meas(m) o_frmean_n_meas(m)],[rectify(o_frmean_n_pred(m)-o_frse_n_pred(m)) o_frmean_n_pred(m)+o_frse_n_pred(m)]);
        e.Color = [0 0 1];
        e.LineWidth = 2.0;
    end
    for m = 1:length(o_frmean_p_meas),
        a = line([o_frmean_p_meas(m)-o_frse_p_meas(m) o_frmean_p_meas(m)+o_frse_p_meas(m)],...
            [o_frmean_p_pred(m) o_frmean_p_pred(m)]);
        a.Color = [1 0 0];
        a.LineWidth = 2.0;
        e = line([o_frmean_p_meas(m) o_frmean_p_meas(m)],[rectify(o_frmean_p_pred(m)-o_frse_p_pred(m)) o_frmean_p_pred(m)+o_frse_p_pred(m)]);
        e.Color = [1 0 0];
        e.LineWidth = 2.0;
    end
    fullset = [o_frmean_n_meas+o_frse_n_meas;o_frmean_p_meas+o_frse_p_meas;o_frmean_n_pred+o_frse_n_pred;o_frmean_p_pred+o_frse_p_pred];
    unl = line([0 max(fullset)],[0 max(fullset)]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Measured Response (Hz)');
    ylabel('Predicted Response (Hz)');
    try
        [nb,nb_ci,nresid,nrint,nstats_] = ...
            regress(o_frmean_n_pred,horzcat(ones(length(o_frmean_n_meas),1),o_frmean_n_meas));
    catch
        [nslope,noffset,nslope_ci,nresid,nrint,nstats_] = ...
            quickregression(o_frmean_n_meas,o_frmean_n_pred,0.05);
    end
     [nslope,noffset,nslope_ci,nresid,nrint,qnstats_] = ...
            quickregression(reshape(o_frmean_n_meas,size(o_frmean_n_pred)),o_frmean_n_pred,0.05);
    try
        [pb,pb_ci,presid,print,pstats_] = ...
            regress(o_frmean_p_pred,horzcat(ones(length(o_frmean_p_meas),1),o_frmean_p_meas));
    catch
        [pslope,poffset,pslope_ci,presid,print,qpstats_] = ...
            quickregression(reshape(o_frmean_p_meas,size(o_frmean_p_pred)),o_frmean_p_pred,0.05);
    end
    legend(['Null, R^2 = ',num2str(nstats_(1,1)),', p = ',num2str(nstats_(1,3)),', Slope CI: [',num2str(nslope_ci(1,1)),' ',num2str(nslope_ci(1,2)),']'],...
    ['Pref, R^2 = ',num2str(pstats_(1,1)),', p = ',num2str(pstats_(1,3)),', Slope CI: [',num2str(pslope_ci(1,1)),' ',num2str(pslope_ci(1,2)),']'],...
    'Location','NorthWest');
    title('OLD');
    hold off;
    %*******
    %Predicted vs. Measured DSI
    %*******
    f30 = figure;
    %scatter(y_mean_actFR_dsi,y_samp_modelFR_dsi,'g','filled');
    scatter(spike_dsi_young,y_samp_modelFR_dsi,'g','filled');
    hold on;
    %scatter(o_mean_actFR_dsi,o_samp_modelFR_dsi,'m','filled');
    scatter(spike_dsi_old,o_samp_modelFR_dsi,'m','filled');
    %for m = 1:length(y_mean_actFR_dsi),
    %    e = line([y_mean_actFR_dsi(1,m) y_mean_actFR_dsi(1,m)],...
    %        [m_model_ydsi(m)-se_model_ydsi(m) m_model_ydsi(m)+se_model_ydsi(m)]);
    %    e.Color = [0 1 0];
    %    e.LineWidth = 2.0;
    %end
    %for m = 1:length(spike_dsi_old),
    %    e = line([spike_dsi_old(1,m) spike_dsi_old(1,m)],...
    %        [m_model_odsi(m)-se_model_odsi(m) m_model_odsi(m)+se_model_odsi(m)]);
    %    e.Color = [0.6 0 0.6];
    %    e.LineWidth = 2.0;
    %end
    unl = line([0 1.2],[0 1.2]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Experimental DSI');
    ylabel('Predicted DSI');
    try
        [yb,yb_ci,yresid,yrint,ystats] = ...
            regress(y_samp_modelFR_dsi,horzcat(ones(length(y_mean_actFR_dsi),1),y_mean_actFR_dsi));
    catch
        [yslope,yoffset,yslope_ci,yresid,yrint,ystats] = ...
            quickregression(y_mean_actFR_dsi',y_samp_modelFR_dsi,0.05);
    end
    [yslope,yoffset,yslope_ci,yresid,yrint,qystats] = ...
            quickregression(reshape(y_mean_actFR_dsi,size(y_samp_modelFR_dsi)),y_samp_modelFR_dsi,0.05);
    try
        [ob,ob_ci,oresid,orint,ostats] = ...
            regress(o_samp_modelFR_dsi,horzcat(ones(length(o_mean_actFR_dsi),1),o_mean_actFR_dsi));
    catch
        [oslope,ooffset,oslope_ci,oresid,orint,ostats] = ...
            quickregression(reshape(o_mean_actFR_dsi,size(o_samp_modelFR_dsi)),o_samp_modelFR_dsi,0.05);
    end
    [oslope,ooffset,oslope_ci,oresid,orint,qostats] = ...
            quickregression(reshape(o_mean_actFR_dsi,size(o_samp_modelFR_dsi)),o_samp_modelFR_dsi,0.05);
    legend(['YOUNG, R^2 = ',num2str(ystats(1,1)),', p = ',num2str(ystats(1,3)),', Young slope CI: [',num2str(yslope_ci(1,1)),' ',num2str(yslope_ci(1,2)),']'],...
    ['OLD, R^2 = ',num2str(ostats(1,1)),', p = ',num2str(ostats(1,3)),', Old slope CI: [',num2str(oslope_ci(1,1)),' ',num2str(oslope_ci(1,2)),']'],...
    'Location','NorthWest');
    title('Observed vs. Predicted DSI');
    hold off;




end

