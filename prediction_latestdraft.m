

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
o_actual_FRpref = [];
o_actual_FRnull = [];
y_actual_FRpref = [];
y_actual_FRnull = [];

se_opv = [];
se_ops = [];
se_onv = [];
se_ons = [];
se_ypv = [];
se_yps = [];
se_ynv = [];
se_yns = [];
se_odsi = [];
se_ydsi = [];
o_m_mps = [];
o_m_mns = [];
o_se_mps = [];
o_se_mns = [];
y_m_mps = [];
y_m_mns = [];
y_se_mps = [];
y_se_mns = [];
m_model_odsi = [];
m_model_ydsi = [];
se_model_odsi = [];
se_model_ydsi = [];
spike_dsi_old = [];
spike_dsi_young = [];
spike_sedsi_old = [];
spike_sedsi_young = [];


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
        load('analyzed_spike.mat','total_spike_byAngle','sampleRate');
        %load('na10M1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        load('GsmoothM1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        load('alltrial_variables.mat','opref_cell_FR_alltrials','onull_cell_FR_alltrials',...
            'ypref_cell_FR_alltrials','ynull_cell_FR_alltrials','opref_cell_Vm_alltrials',...
            'onull_cell_Vm_alltrials','ypref_cell_Vm_alltrials','ynull_cell_Vm_alltrials',...
            'old_cellnames','young_cellnames');
        sampling_rate = sampleRate;
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
        ms_block = 0.03;
        block_size = ceil(ms_block*sampling_rate);
        N_blocks = 1000;
        samp_err_pv = [];
        samp_err_ps = [];
        samp_err_nv = [];
        samp_err_ns = [];
        o_samp_err_mps = [];
        o_samp_err_mns = [];
        y_samp_err_mps = [];
        y_samp_err_mns = [];
        samp_spike_dsi = [];
        o_samp_modelFR_dsi = [];
        y_samp_modelFR_dsi = [];
        mean_pv = sum(cell2mat(reshape(pref_Vm{1,1},1,...
            size(pref_Vm{1,1},1)*size(pref_Vm{1,1},2))),2);
        mean_pv = mean_pv./length(pref_Vm{1,1});
        mean_ps = sum(cell2mat(reshape(pref_spike{1,1},1,...
            size(pref_spike{1,1},1)*size(pref_spike{1,1},2))),2);
        mean_ps = mean_ps./length(pref_spike{1,1});
        mean_nv = sum(cell2mat(reshape(null_Vm{1,1},1,...
            size(null_Vm{1,1},1)*size(null_Vm{1,1},2))),2);
        mean_nv = mean_nv./length(null_Vm{1,1});
        mean_ns = sum(cell2mat(reshape(null_spike{1,1},1,...
        size(null_spike{1,1},1)*size(null_spike{1,1},2))),2);
        mean_ns = mean_ns./length(null_spike{1,1});
        block_range = 1:(length(mean_pv)-block_size);
        
        for n = 1:length(block_range),
            block_edges{n,1} = [n:n+block_size];
        end
        block_samp = cell(N_blocks,1);
        
         block_indices = 1:N_blocks;
        sample_indices = datasample(block_indices,length(block_indices));
        for b = 1:N_blocks,
            pv_block_samp{b,1} = mean_pv(block_edges{sample_indices(1,b)},1);
            ps_block_samp{b,1} = mean_ps(block_edges{sample_indices(1,b)},1);
            nv_block_samp{b,1} = mean_nv(block_edges{sample_indices(1,b)},1);
            ns_block_samp{b,1} = mean_ns(block_edges{sample_indices(1,b)},1);
        end
        samp_pref_blockmv = [];
        samp_pref_blockms = [];
        samp_null_blockmv = [];
        samp_null_blockms = [];
        %get randomized block samples to apply to measured data evaluation
        for s = 1:length(sample_indices),
            if mean(pv_block_samp{s,1}) > mean(nv_block_samp{s,1}),
                samp_pref_blockmv = [samp_pref_blockmv;mean(pv_block_samp{s,1})];
                samp_null_blockmv = [samp_null_blockmv;mean(nv_block_samp{s,1})];
            else
                samp_pref_blockmv = [samp_pref_blockmv;mean(nv_block_samp{s,1})];
                samp_null_blockmv = [samp_null_blockmv;mean(pv_block_samp{s,1})];
            end
            if mean(ps_block_samp{s,1}) > mean(ns_block_samp{s,1}),
                samp_pref_blockms = [samp_pref_blockms;mean(ps_block_samp{s,1})];
                samp_null_blockms = [samp_null_blockms;mean(ns_block_samp{s,1})];
            else
                samp_pref_blockms = [samp_pref_blockms;mean(ns_block_samp{s,1})];
                samp_null_blockms = [samp_null_blockms;mean(ps_block_samp{s,1})];
            end
        end
        if strcmp(subfile,'old'),
            spike_dsi_old(end+1,1) = mean(samp_pref_blockms./samp_null_blockms);
            spike_sedsi_old(end+1,1) = std(samp_pref_blockms./samp_null_blockms)/sqrt(length(samp_pref_blockms));
        else
            spike_dsi_young(end+1,1) = mean(samp_pref_blockms./samp_null_blockms);
            spike_sedsi_young(end+1,1) = std(samp_pref_blockms./samp_null_blockms)/sqrt(length(samp_null_blockms));
        end
        
        %samp_err_pv = mean(samp_pref_blockmv);
        %samp_err_ps = mean(samp_pref_blockms);
        %samp_err_nv = mean(samp_null_blockmv);
        %samp_err_ns = mean(samp_null_blockms);
        %if strcmp(subfile,'old'),
        %    se_opv = sqrt((1/(N_blocks-1))*sum((samp_pref_blockmv-samp_err_pv).^2));
        %    se_ops = sqrt((1/(N_blocks-1))*sum((samp_pref_blockms-samp_err_ps).^2));
        %    se_onv = sqrt((1/(N_blocks-1))*sum((samp_null_blockmv-samp_err_nv).^2));
        %    se_ons = sqrt((1/(N_blocks-1))*sum((samp_null_blockms-samp_err_ns).^2));
        %else
        %    se_ypv = sqrt((1/(N_blocks-1))*sum((samp_pref_blockmv-samp_err_pv).^2));
        %    se_yps = sqrt((1/(N_blocks-1))*sum((samp_pref_blockms-samp_err_ps).^2));
        %    se_ynv = sqrt((1/(N_blocks-1))*sum((samp_null_blockmv-samp_err_nv).^2));
        %    se_yns = sqrt((1/(N_blocks-1))*sum((samp_null_blockms-samp_err_ns).^2));
        %end
        
        %get consecutive bins to apply to model predictions
        ord_blocks = block_size:block_size:length(pref_Vm{1,1}{1,1});
        ord_blocks = [0,ord_blocks];
        for n = 1:length(ord_blocks)-1,
            pv_block_ord{n,1} = mean_pv(ord_blocks(1,n)+1:ord_blocks(1,n+1));
            ps_block_ord{n,1} = mean_ps(ord_blocks(1,n)+1:ord_blocks(1,n+1));
            nv_block_ord{n,1} = mean_nv(ord_blocks(1,n)+1:ord_blocks(1,n+1));
            ns_block_ord{n,1} = mean_ns(ord_blocks(1,n)+1:ord_blocks(1,n+1));
        end
        
        %get data sample for model
        ord_indices = 1:length(ord_blocks)-1;
        %ord_loc = datasample(ord_indices,length(ord_indices));
        samp_pref_modelmv = [];
        samp_null_modelmv = [];
        samp_pref_modelms = [];
        samp_null_modelms = [];
        for nn = 1:nsims,
            ord_loc = datasample(ord_indices,length(ord_indices));
            new_pv_block = cell(1,length(ord_loc));
            new_ps_block = cell(1,length(ord_loc));
            new_nv_block = cell(1,length(ord_loc));
            new_ns_block = cell(1,length(ord_loc));
            for ns = 1:length(ord_loc)
                new_pv_block{1,ns} = pv_block_ord{ord_loc(1,ns),1};
                new_ps_block{1,ns} = ps_block_ord{ord_loc(1,ns),1};
                new_nv_block{1,ns} = nv_block_ord{ord_loc(1,ns),1};
                new_ns_block{1,ns} = ns_block_ord{ord_loc(1,ns),1};
                
            end
            samp_pref_modelmv = [samp_pref_modelmv;mean(mean(cell2mat(new_pv_block),2))];
            samp_pref_modelms = [samp_pref_modelms;mean(mean(cell2mat(new_ps_block),2))];
            samp_null_modelmv = [samp_null_modelmv;mean(mean(cell2mat(new_nv_block),2))];
            samp_null_modelms = [samp_null_modelms;mean(mean(cell2mat(new_ns_block),2))];
        end
        
        if strcmp(subfile,'old'),
            o_samp_pref_modelFR = rectify(old_params{j,1}.tanh_params.a +...
                old_params{j,1}.tanh_params.b.*tanh((samp_pref_modelmv-old_params{j,1}.tanh_params.c)./...
                old_params{j,1}.tanh_params.d));
            o_samp_null_modelFR = rectify(old_params{j,1}.tanh_params.a +...
                old_params{j,1}.tanh_params.b.*tanh((samp_null_modelmv-old_params{j,1}.tanh_params.c)./...
                old_params{j,1}.tanh_params.d));
            if mean(o_samp_pref_modelFR) > mean(o_samp_null_modelFR),
                o_samp_err_mps = [o_samp_err_mps;mean(o_samp_pref_modelFR)];
                o_samp_err_mns = [o_samp_err_mns;mean(o_samp_null_modelFR)];
            else
                o_samp_err_mps = [o_samp_err_mps;mean(o_samp_null_modelFR)];
                o_samp_err_mns = [o_samp_err_mns;mean(o_samp_pref_modelFR)];
            end
            o_samp_modelFR_dsi = [o_samp_modelFR_dsi;abs(mean(o_samp_pref_modelFR)-...
                mean(o_samp_null_modelFR))/max([mean(o_samp_pref_modelFR);mean(o_samp_null_modelFR)])];
        else
            y_samp_pref_modelFR = rectify(young_params{j,1}.tanh_params.a +...
                    young_params{j,1}.tanh_params.b.*tanh((samp_pref_modelmv-young_params{j,1}.tanh_params.c)./...
                    young_params{j,1}.tanh_params.d));
            y_samp_null_modelFR = rectify(young_params{j,1}.tanh_params.a +...
                    young_params{j,1}.tanh_params.b.*tanh((samp_null_modelmv-young_params{j,1}.tanh_params.c)./...
                    young_params{j,1}.tanh_params.d));
            if mean(y_samp_pref_modelFR) > mean(y_samp_null_modelFR),
                    y_samp_err_mps = [y_samp_err_mps;mean(y_samp_pref_modelFR)];
                    y_samp_err_mns = [y_samp_err_mns;mean(y_samp_null_modelFR)];
            else
                    y_samp_err_mps = [y_samp_err_mps;mean(y_samp_null_modelFR)];
                    y_samp_err_mns = [y_samp_err_mns;mean(y_samp_pref_modelFR)];
            end
            y_samp_modelFR_dsi = [y_samp_modelFR_dsi;abs(mean(y_samp_pref_modelFR)-...
                mean(y_samp_null_modelFR))/max([mean(y_samp_pref_modelFR);mean(y_samp_null_modelFR)])];
        end
        %actual DSI
        samp_spike_dsi = [samp_spike_dsi;abs(mean(samp_pref_blockms)-mean(samp_null_blockms))...
            /max([mean(samp_pref_blockms);mean(samp_null_blockms)])];
        %predicted DSI
        
        
        samp_err_pv = mean(samp_pref_blockmv);
        samp_err_ps = mean(samp_pref_blockms);
        samp_err_nv = mean(samp_null_blockmv);
        samp_err_ns = mean(samp_null_blockms);
        if strcmp(subfile,'old'),
            se_opv(end+1,1) = sqrt((1/(nsims-1))*sum((samp_pref_blockmv-samp_err_pv).^2));
            se_ops(end+1,1) = sqrt((1/(nsims-1))*sum((samp_pref_blockms-samp_err_ps).^2));
            se_onv(end+1,1) = sqrt((1/(nsims-1))*sum((samp_null_blockmv-samp_err_nv).^2));
            se_ons(end+1,1) = sqrt((1/(nsims-1))*sum((samp_null_blockms-samp_err_ns).^2));
            se_odsi(end+1,1) = sqrt((1/(nsims-1))*sum((samp_spike_dsi-mean(samp_spike_dsi)).^2));
            o_m_mps(end+1,1) = mean(o_samp_err_mps);
            o_m_mns(end+1,1) = mean(o_samp_err_mns);
            o_se_mps(end+1,1) = sqrt((1/(length(ord_blocks)-1))*sum((o_samp_err_mps-mean(o_samp_err_mps)).^2));
            o_se_mns(end+1,1) = sqrt((1/(length(ord_blocks)-1))*sum((o_samp_err_mns-mean(o_samp_err_mns)).^2));
            m_model_odsi(end+1,1) = mean(o_samp_modelFR_dsi);
            se_model_odsi(end+1,1) = sqrt((1/(length(ord_blocks)-1))*sum((o_samp_modelFR_dsi-mean(o_samp_modelFR_dsi)).^2));
            o_actual_FRpref(end+1,1) = mean(samp_pref_blockms);
            o_actual_FRnull(end+1,1) = mean(samp_null_blockms);
        else
            se_ypv(end+1,1) = sqrt((1/(N_blocks-1))*sum((samp_pref_blockmv-samp_err_pv).^2));
            se_yps(end+1,1) = sqrt((1/(N_blocks-1))*sum((samp_pref_blockms-samp_err_ps).^2));
            se_ynv(end+1,1) = sqrt((1/(N_blocks-1))*sum((samp_null_blockmv-samp_err_nv).^2));
            se_yns(end+1,1) = sqrt((1/(N_blocks-1))*sum((samp_null_blockms-samp_err_ns).^2));
            se_ydsi(end+1,1) = sqrt((1/(nsims-1))*sum((samp_spike_dsi-mean(samp_spike_dsi)).^2));
            y_m_mps(end+1,1) = mean(y_samp_err_mps);
            y_m_mns(end+1,1) = mean(y_samp_err_mns);
            y_se_mps(end+1,1) = sqrt((1/(length(ord_blocks)-1))*sum((y_samp_err_mps-mean(y_samp_err_mps)).^2));
            y_se_mns(end+1,1) = sqrt((1/(length(ord_blocks)-1))*sum((y_samp_err_mns-mean(y_samp_err_mns)).^2));
            m_model_ydsi(end+1,1) = mean(y_samp_modelFR_dsi);
            se_model_ydsi(end+1,1) = sqrt((1/(length(ord_blocks)-1))*sum((y_samp_modelFR_dsi-mean(y_samp_modelFR_dsi)).^2));
            y_actual_FRpref(end+1,1) = mean(samp_pref_blockms);
            y_actual_FRnull(end+1,1) = mean(samp_null_blockms);
        end
        cd ../
    end
    cd ../
end
       
        
        
        
        %%%**************
        %%%PLOTTING
        %%%**************
if allow_flip == 1,
    %Predicted vs. Measured FR Response
    f = figure;
    subplot(1,2,1);
    scatter(y_actual_FRnull,y_m_mns,'b','filled');
    hold on;
    scatter(y_actual_FRpref,y_m_mps,'r','filled');
    for m = 1:length(y_actual_FRnull),
        a = line([y_actual_FRnull(m)-se_yns(m) y_actual_FRnull(m)+se_yns(m)],...
            [y_m_mns(m) y_m_mns(m)]);
        a.Color = [0 0 1];
        a.LineWidth = 2.0;
        e = line([y_actual_FRnull(m) y_actual_FRnull(m)],[rectify(y_m_mns(m)-y_se_mns(m)) y_m_mns(m)+y_se_mns(m)]);
        e.Color = [0 0 1];
        e.LineWidth = 2.0;
    end
    for m = 1:length(y_actual_FRpref),
        a = line([y_actual_FRpref(m)-se_yps(m) y_actual_FRpref(m)+se_yps(m)],...
            [y_m_mps(m) y_m_mps(m)]);
        a.Color = [1 0 0];
        a.LineWidth = 2.0;
        e = line([y_actual_FRpref(m) y_actual_FRpref(m)],[rectify(y_m_mps(m)-y_se_mps(m)) y_m_mps(m)+y_se_mps(m)]);
        e.Color = [1 0 0];
        e.LineWidth = 2.0;
    end
    fullset = [y_actual_FRnull+se_yns;y_actual_FRpref+se_yps;y_m_mns+y_se_mns;y_m_mps+y_se_mps];
    unl = line([0 max(fullset)],[0 max(fullset)]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Measured Response (Hz)');
    ylabel('Predicted Response (Hz)');
    try
        [nb,nb_ci,nresid,nrint,nstats] = ...
            regress(y_m_mns,horzcat(ones(length(y_actual_FRnull),1),y_actual_FRnull));
    catch
        [nslope,noffset,nslope_ci,nresid,nrint,nstats] = ...
            quickregression(y_actual_FRnull,y_m_mns,0.05);
    end
    [nslope,noffset,nslope_ci,nresid,nrint,qnstats] = ...
            quickregression(reshape(y_actual_FRnull,size(y_m_mns)),y_m_mns,0.05);
    try
        [pb,pb_ci,presid,print,pstats] = ...
            regress(y_m_mps,horzcat(ones(length(y_actual_FRpref),1),y_actual_FRpref));
    catch
        [pslope,poffset,pslope_ci,presid,print,pstats] = ...
            quickregression(y_actual_FRpref,y_m_mps,0.05);
    end
    [pslope,poffset,pslope_ci,presid,print,qpstats] = ...
            quickregression(reshape(y_actual_FRpref,size(y_m_mps)),y_m_mps,0.05);
    legend(['Null, R^2 = ',num2str(nstats(1,1)),', p = ',num2str(nstats(1,3)),', Slope CI: [',num2str(nslope_ci(1,1)),' ',num2str(nslope_ci(1,2)),']'],...
    ['Pref, R^2 = ',num2str(pstats(1,1)),', p = ',num2str(pstats(1,3)),', Slope CI: [',num2str(pslope_ci(1,1)),' ',num2str(pslope_ci(1,2)),']'],...
    'Location','NorthWest');
    title('YOUNG');
    hold off;
    subplot(1,2,2);
    scatter(o_actual_FRnull,o_m_mns,'b','filled');
    hold on;
    scatter(o_actual_FRpref,o_m_mps,'r','filled');
    for m = 1:length(o_actual_FRnull),
        a = line([o_actual_FRnull(m)-se_ons(m) o_actual_FRnull(m)+se_ons(m)],...
            [o_m_mns(m) o_m_mns(m)]);
        a.Color = [0 0 1];
        a.LineWidth = 2.0;
        e = line([o_actual_FRnull(m) o_actual_FRnull(m)],[rectify(o_m_mns(m)-o_se_mns(m)) o_m_mns(m)+o_se_mns(m)]);
        e.Color = [0 0 1];
        e.LineWidth = 2.0;
    end
    for m = 1:length(o_actual_FRpref),
        a = line([o_actual_FRpref(m)-se_ops(m) o_actual_FRpref(m)+se_ops(m)],...
            [o_m_mps(m) o_m_mps(m)]);
        a.Color = [1 0 0];
        a.LineWidth = 2.0;
        e = line([o_actual_FRpref(m) o_actual_FRpref(m)],[rectify(o_m_mps(m)-o_se_mps(m)) o_m_mps(m)+o_se_mps(m)]);
        e.Color = [1 0 0];
        e.LineWidth = 2.0;
    end
    fullset = [o_actual_FRnull+se_ons;o_actual_FRpref+se_ops;o_m_mns+o_se_mns;o_m_mps+o_se_mps];
    unl = line([0 max(fullset)],[0 max(fullset)]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Measured Response (Hz)');
    ylabel('Predicted Response (Hz)');
    try
        [nb,nb_ci,nresid,nrint,nstats_] = ...
            regress(o_m_mns,horzcat(ones(length(o_actual_FRnull),1),o_actual_FRnull));
    catch
        [nslope,noffset,nslope_ci,nresid,nrint,nstats_] = ...
            quickregression(o_actual_FRnull,o_m_mns,0.05);
    end
     [nslope,noffset,nslope_ci,nresid,nrint,qnstats_] = ...
            quickregression(reshape(o_actual_FRnull,size(o_m_mns)),o_m_mns,0.05);
    try
        [pb,pb_ci,presid,print,pstats_] = ...
            regress(o_m_mps,horzcat(ones(length(o_actual_FRpref),1),o_actual_FRpref));
    catch
        [pslope,poffset,pslope_ci,presid,print,qpstats_] = ...
            quickregression(reshape(o_actual_FRpref,size(o_m_mps)),o_m_mps,0.05);
    end
    legend(['Null, R^2 = ',num2str(nstats_(1,1)),', p = ',num2str(nstats_(1,3)),', Slope CI: [',num2str(nslope_ci(1,1)),' ',num2str(nslope_ci(1,2)),']'],...
    ['Pref, R^2 = ',num2str(pstats_(1,1)),', p = ',num2str(pstats_(1,3)),', Slope CI: [',num2str(pslope_ci(1,1)),' ',num2str(pslope_ci(1,2)),']'],...
    'Location','NorthWest');
    title('OLD');
    hold off;
    %Predicted vs. Measured DSI
    f30 = figure;
    scatter(spike_dsi_young,m_model_ydsi,'g','filled');
    hold on;
    scatter(spike_dsi_old,m_model_odsi,'m','filled');
    for m = 1:length(spike_dsi_young),
        e = line([spike_dsi_young(1,m) spike_dsi_young(1,m)],...
            [m_model_ydsi(m)-se_model_ydsi(m) m_model_ydsi(m)+se_model_ydsi(m)]);
        e.Color = [0 1 0];
        e.LineWidth = 2.0;
        em = line([spike_dsi_young(m,1)-spike_sedsi_young(m,1) spike_dsi_young(m,1)+...
            spike_sedsi_young(m,1)],[m_model_ydsi(m) m_model_ydsi(m)]);
        em.Color = [0 1 0];
        em.LineWidth = 2.0;
    end
    for m = 1:length(spike_dsi_old),
        e = line([spike_dsi_old(1,m) spike_dsi_old(1,m)],...
            [m_model_odsi(m)-se_model_odsi(m) m_model_odsi(m)+se_model_odsi(m)]);
        e.Color = [0.6 0 0.6];
        e.LineWidth = 2.0;
        em = line([spike_dsi_old(m,1)-spike_sedsi_old(m,1) spike_dsi_old(m,1)+...
            spike_sedsi_old(m,1)],[m_model_odsi(m) m_model_odsi(m)]);
        em.Color = [0.6 0 0.6];
        em.LineWidth = 2.0;
    end
    unl = line([0 1.2],[0 1.2]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Experimental DSI');
    ylabel('Predicted DSI');
    try
        [yb,yb_ci,yresid,yrint,ystats] = ...
            regress(m_model_ydsi,horzcat(ones(length(spike_dsi_young),1),spike_dsi_young));
    catch
        [yslope,yoffset,yslope_ci,yresid,yrint,ystats] = ...
            quickregression(spike_dsi_young',m_model_ydsi,0.05);
    end
    [yslope,yoffset,yslope_ci,yresid,yrint,qystats] = ...
            quickregression(reshape(spike_dsi_young,size(m_model_ydsi)),m_model_ydsi,0.05);
    try
        [ob,ob_ci,oresid,orint,ostats] = ...
            regress(m_model_odsi,horzcat(ones(length(spike_dsi_old),1),spike_dsi_old));
    catch
        [oslope,ooffset,oslope_ci,oresid,orint,ostats] = ...
            quickregression(reshape(spike_dsi_old,size(m_model_odsi)),m_model_odsi,0.05);
    end
    [oslope,ooffset,oslope_ci,oresid,orint,qostats] = ...
            quickregression(reshape(spike_dsi_old,size(m_model_odsi)),m_model_odsi,0.05);
    legend(['YOUNG, R^2 = ',num2str(ystats(1,1)),', p = ',num2str(ystats(1,3)),', Young slope CI: [',num2str(yslope_ci(1,1)),' ',num2str(yslope_ci(1,2)),']'],...
    ['OLD, R^2 = ',num2str(ostats(1,1)),', p = ',num2str(ostats(1,3)),', Old slope CI: [',num2str(oslope_ci(1,1)),' ',num2str(oslope_ci(1,2)),']'],...
    'Location','NorthWest');
    title('Observed vs. Predicted DSI');
    hold off;
    
else
    
    f = figure;
    subplot(1,2,1);
    scatter(y_actual_FRnull,y_predicted_FRnull,'b','filled');
    hold on;
    scatter(y_actual_FRpref,y_predicted_FRpref,'r','filled');
    for m = 1:length(y_actual_FRnull),
        a = line([y_actual_FRnull(m)-se_yns(m) y_actual_FRnull(m)+se_yns(m)],[y_predicted_FRnull(m) y_predicted_FRnull(m)]);
        a.Color = [0 0 1];
        e = line([y_actual_FRnull(m) y_actual_FRnull(m)],[y_pre_FRnull_down(m) y_pre_FRnull_up(m)]);
        e.Color = [0 0 1];
    end
    for m = 1:length(y_actual_FRpref),
        a = line([y_actual_FRpref(m)-se_yps(m) y_actual_FRpref(m)+se_yps(m)],[y_predicted_FRpref(m) y_predicted_FRpref(m)]);
        a.Color = [1 0 0];
        e = line([y_actual_FRpref(m) y_actual_FRpref(m)],[y_pre_FRpref_down(m) y_pre_FRpref_up(m)]);
        e.Color = [1 0 0];
    end
    fullset = [y_actual_FRnull;y_predicted_FRnull;y_actual_FRpref;y_predicted_FRpref];
    unl = line([0 max(fullset)],[0 max(fullset)]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Measured Response (Hz)');
    ylabel('Predicted Response (Hz)');
    try
        [nb,nb_ci,nresid,nrint,nstats] = ...
            regress(y_predicted_FRnull,horzcat(ones(length(y_actual_FRnull),1),y_actual_FRnull));
    catch
        [nslope,noffset,nslope_ci,nresid,nrint,nstats] = ...
            quickregression(y_actual_FRnull,y_predicted_FRnull,0.05);
    end
    try
        [pb,pb_ci,presid,print,pstats] = ...
            regress(y_predicted_FRpref,horzcat(ones(length(y_actual_FRpref),1),y_actual_FRpref));
    catch
        [pslope,poffset,pslope_ci,presid,print,pstats] = ...
        quickregression(y_actual_FRpref,y_predicted_FRpref,0.05);
    end
    legend(['Null, R^2 = ',num2str(nstats(1,1)),', p = ',num2str(nstats(1,3))],...
        ['Pref, R^2 = ',num2str(pstats(1,1)),', p = ',num2str(pstats(1,3))],'Location','NorthWest');
    title('YOUNG');
    hold off;
    subplot(1,2,2);
    scatter(o_actual_FRnull,o_predicted_FRnull,'b','filled');
    hold on;
    scatter(o_actual_FRpref,o_predicted_FRpref,'r','filled');
    for m = 1:length(o_actual_FRnull),
        a = line([o_actual_FRnull(m)-se_ons(m) o_actual_FRnull(m)+se_ons(m)],[o_predicted_FRnull(m) o_predicted_FRnull(m)]);
        a.Color = [0 0 1];
        e = line([o_actual_FRnull(m) o_actual_FRnull(m)],[o_pre_FRnull_down(m) o_pre_FRnull_up(m)]);
        e.Color = [0 0 1];
    end
    for m = 1:length(o_actual_FRpref),
        a = line([o_actual_FRpref(m)-se_ops(m) o_actual_FRpref(m)+se_ops(m)],[o_predicted_FRpref(m) o_predicted_FRpref(m)]);
        a.Color = [1 0 0];
        e = line([o_actual_FRpref(m) o_actual_FRpref(m)],[o_pre_FRpref_down(m) o_pre_FRpref_up(m)]);
        e.Color = [1 0 0];
    end
    o_fullset = [o_actual_FRnull;o_predicted_FRnull;o_actual_FRpref;o_predicted_FRpref];
    unl = line([0 max(o_fullset)],[0 max(o_fullset)]);
    unl.LineStyle = '--';
    unl.Color = [0 0 0];
    xlabel('Measured Response (Hz)');
    ylabel('Predicted Response (Hz)');
    try
        [nb,nb_ci,nresid,nrint,nstats_] = ...
            regress(o_predicted_FRnull,horzcat(ones(length(o_actual_FRnull),1),o_actual_FRnull));
    catch
        [nslope,noffset,nslope_ci,nresid,nrint,nstats_] = ...
            quickregression(o_actual_FRnull,o_predicted_FRnull,0.05);
    end
    try
        [pb,pb_ci,presid,print,pstats_] = ...
            regress(o_predicted_FRpref,horzcat(ones(length(o_actual_FRpref),1),o_actual_FRpref));
    catch
        [pslope,poffset,pslope_ci,presid,print,pstats_] = ...
            quickregression(o_actual_FRpref,o_predicted_FRpref,0.05);
    end
    legend(['Null, R^2 =',num2str(nstats_(1,1)),', p = ',num2str(nstats_(1,3))],...
        ['Pref, R^2 =',num2str(pstats_(1,1)),', p = ',num2str(pstats_(1,3))],'Location','NorthWest');
    title('OLD');
    hold off;

end

