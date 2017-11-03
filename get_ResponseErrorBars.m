function [ output_args ] = get_ResponseErrorBars( allow_flip )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
% ALLOW_FLIP should be set to '1' if simulation results are permitted to
% switch response categories, i.e. if PREF/NULL is an absolute designator for
% Dir1/Dir2, then allow_flip = 0; if PREF/NULL is relative in the usual sense
% then allow_flip = 1.
load('collected_params.mat');

[folder,directories] = walk_directories();
cd(folder);

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
    sub_dir = char(sel_Snames(s))
    [~,subfile,~] = fileparts(sub_folder);
    for j = 1:length(sub_dir),
        cur_dir = sub_dir(j,:);
        cd(cur_dir);
        load('analyzed_Vm.mat','total_vm_byAngle','stimvalues');
        load('analyzed_spike.mat','total_spike_byAngle');
        %load('na10M1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        load('GsmoothM1_F_BL_X_collstim_Xbinfitdata.mat','pref_ind','null_ind');
        all_dirind = 1:length(stimvalues)-1;
        down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
        up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
        down_i = circshift(all_dirind,1,2);
        up_i = circshift(all_dirind,-1,2);
        %pref_set = [down_(pref_ind);pref_stim;up_(pref_ind)];
        %null_set = [down_(null_ind);null_stim;up_(null_ind)];
        pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
        null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
        pref_Vm = reshape(total_vm_byAngle(pref_set_index(:,1),:),...
            (length(pref_set_index)*size(total_vm_byAngle,2)),1);
        pref_spike = reshape(total_spike_byAngle(pref_set_index(:,1),:),...
            (length(pref_set_index)*size(total_vm_byAngle,2)),1);
        null_Vm = reshape(total_vm_byAngle(null_set_index(:,1),:),...
            (length(null_set_index)*size(total_vm_byAngle,2)),1);
        null_spike = reshape(total_spike_byAngle(null_set_index(:,1),:),...
            (length(null_set_index)*size(total_spike_byAngle,2)),1);
        nsims = 10000;
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
        for n = 1:nsims,
            if allow_flip == 1,
                samp_pref_vm = datasample(pref_Vm,length(pref_Vm));
                samp_null_vm = datasample(null_Vm,length(null_Vm));
                if mean(samp_pref_vm) > mean(samp_null_vm),
                    samp_err_pv = [samp_err_pv;mean(samp_pref_vm)];
                    samp_err_nv = [samp_err_nv;mean(samp_null_vm)];
                else
                    samp_err_pv = [samp_err_pv;mean(samp_null_vm)];
                    samp_err_nv = [samp_err_nv;mean(samp_pref_vm)];
                end
                if strcmp(subfile,'old'),
                    o_samp_pref_modelFR = rectify(old_params{j,1}.tanh_params.a +...
                        old_params{j,1}.tanh_params.b.*tanh((samp_pref_vm-old_params{j,1}.tanh_params.c)./...
                        old_params{j,1}.tanh_params.d));
                    o_samp_null_modelFR = rectify(old_params{j,1}.tanh_params.a +...
                        old_params{j,1}.tanh_params.b.*tanh((samp_null_vm-old_params{j,1}.tanh_params.c)./...
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
                        young_params{j,1}.tanh_params.b.*tanh((samp_pref_vm-young_params{j,1}.tanh_params.c)./...
                        young_params{j,1}.tanh_params.d));
                    y_samp_null_modelFR = rectify(young_params{j,1}.tanh_params.a +...
                        young_params{j,1}.tanh_params.b.*tanh((samp_null_vm-young_params{j,1}.tanh_params.c)./...
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
                
                %Alt.
                samp_pref_spike = datasample(pref_spike,length(pref_spike));
                samp_null_spike = datasample(null_spike,length(null_spike));
                if mean(samp_pref_spike) > mean(samp_null_spike),
                    samp_err_ps = [samp_err_ps;mean(samp_pref_spike)];
                    samp_err_ns = [samp_err_ns;mean(samp_null_spike)];
                else
                    samp_err_ps = [samp_err_ps;mean(samp_null_spike)];
                    samp_err_ns = [samp_err_ns;mean(samp_pref_spike)];
                end
                samp_spike_dsi = [samp_spike_dsi;abs(mean(samp_pref_spike)-mean(samp_null_spike))...
                        /max([mean(samp_pref_spike);mean(samp_null_spike)])];
            else
                samp_pref_vm = datasample(pref_Vm,length(pref_Vm));
                samp_err_pv = [samp_err_pv;mean(samp_pref_vm)];
                samp_pref_spike = datasample(pref_spike,length(pref_spike));
                samp_err_ps = [samp_err_ps;mean(samp_pref_spike)];
                samp_null_vm = datasample(null_Vm,length(null_Vm));
                samp_err_nv = [samp_err_nv;mean(samp_null_vm)];
                samp_null_spike = datasample(null_spike,length(null_spike));
                samp_err_ns = [samp_err_ns;mean(samp_null_spike)];
            end
        end
        if strcmp(subfile,'old'),
            se_opv(end+1,1) = (1/(nsims-1))*sum((samp_err_pv-mean(samp_err_pv)).^2);
            se_ops(end+1,1) = (1/(nsims-1))*sum((samp_err_ps-mean(samp_err_ps)).^2);
            se_onv(end+1,1) = (1/(nsims-1))*sum((samp_err_nv-mean(samp_err_nv)).^2);
            se_ons(end+1,1) = (1/(nsims-1))*sum((samp_err_ns-mean(samp_err_ns)).^2);
            se_odsi(end+1,1) = (1/(nsims-1))*sum((samp_spike_dsi-mean(samp_spike_dsi)).^2);
            o_m_mps(end+1,1) = mean(o_samp_err_mps);
            o_m_mns(end+1,1) = mean(o_samp_err_mns);
            o_se_mps(end+1,1) = (1/(nsims-1))*sum((o_samp_err_mps-mean(o_samp_err_mps)).^2);
            o_se_mns(end+1,1) = (1/(nsims-1))*sum((o_samp_err_mns-mean(o_samp_err_mns)).^2);
            m_model_odsi(end+1,1) = mean(o_samp_modelFR_dsi);
            se_model_odsi(end+1,1) = (1/(nsims-1))*sum((o_samp_modelFR_dsi-mean(o_samp_modelFR_dsi)).^2);
        else
            se_ypv(end+1,1) = (1/(nsims-1))*sum((samp_err_pv-mean(samp_err_pv)).^2);
            se_yps(end+1,1) = (1/(nsims-1))*sum((samp_err_ps-mean(samp_err_ps)).^2);
            se_ynv(end+1,1) = (1/(nsims-1))*sum((samp_err_nv-mean(samp_err_nv)).^2);
            se_yns(end+1,1) = (1/(nsims-1))*sum((samp_err_ns-mean(samp_err_ns)).^2);
            se_ydsi(end+1,1) = (1/(nsims-1))*sum((samp_spike_dsi-mean(samp_spike_dsi)).^2);
            y_m_mps(end+1,1) = mean(y_samp_err_mps);
            y_m_mns(end+1,1) = mean(y_samp_err_mns);
            y_se_mps(end+1,1) = (1/(nsims-1))*sum((y_samp_err_mps-mean(y_samp_err_mps)).^2);
            y_se_mns(end+1,1) = (1/(nsims-1))*sum((y_samp_err_mns-mean(y_samp_err_mns)).^2);
            m_model_ydsi(end+1,1) = mean(y_samp_modelFR_dsi);
            se_model_ydsi(end+1,1) = (1/(nsims-1))*sum((y_samp_modelFR_dsi-mean(y_samp_modelFR_dsi)).^2);
        end
            
        
        %following section incorrect since error bars are on single-cell
        %points
        %if strcmp(subfile,'old'),
        %    for sp = 1:length(pref_set_index),
        %        o_pref_Vm = [o_pref_Vm;reshape(total_vm_byAngle(pref_set_index(sp,1),:),...
        %            size(total_vm_byAngle,2),1)];
        %        o_pref_spike = [o_pref_spike;reshape(total_spike_byAngle(pref_set_index(sp,1),:),...
        %            size(total_spike_byAngle,2),1)];
        %    end
        %    for sn = 1:length(null_set_index),
        %        o_null_Vm = [o_null_Vm;reshape(total_vm_byAngle(null_set_index(sn,1),:),...
        %            size(total_vm_byAngle,2),1)];
        %        o_null_spike = [o_null_spike;reshape(total_spike_byAngle(null_set_index(sn,1),:),...
        %            size(total_spike_byAngle,2),1)];
        %    end
        %else
        %    for sp = 1:length(pref_set_index),
        %        y_pref_Vm = [y_pref_Vm;reshape(total_vm_byAngle(pref_set_index(sp,1),:),...
        %            size(total_vm_byAngle,2),1)];
        %        y_pref_spike = [y_pref_spike;reshape(total_spike_byAngle(pref_set_index(sp,1),:),...
        %            size(total_spike_byAngle,2),1)];
        %    end
        %    for sn = 1:length(null_set_index),
        %        y_null_Vm = [y_null_Vm;reshape(total_vm_byAngle(null_set_index(sn,1),:),...
        %            size(total_vm_byAngle,2),1)];
        %        y_null_spike = [y_null_spike;reshape(total_spike_byAngle(null_set_index(sn,1),:),...
        %            size(total_spike_byAngle,2),1)];
        %    end
        %end
        clear pref_ind null_ind;
        cd ../
    end
    cd ../
end

%nsims = 10000;
%samp_err_opv = [];
%samp_err_ops = [];
%samp_err_onv = [];
%samp_err_ons = [];
%samp_err_ypv = [];
%samp_err_yps = [];
%samp_err_ynv = [];
%samp_err_yns = [];
%
%for n = 1:nsims,
%    samp_opref_vm = datasample(o_pref_Vm,length(o_pref_Vm));
%    samp_err_opv = [samp_err_opv;mean(samp_opref_vm)];
%    samp_opref_spike = datasample(o_pref_spike,length(o_pref_spike));
%    samp_err_ops = [samp_err_ops;mean(samp_opref_spike)];
%    samp_onull_vm = datasample(o_null_Vm,length(o_null_Vm));
%    samp_err_onv = [samp_err_onv;mean(samp_onull_vm)];
%    samp_onull_spike = datasample(o_null_spike,length(o_null_spike));
%    samp_err_ons = [samp_err_ons;mean(samp_onull_spike)];
%    samp_ypref_vm = datasample(y_pref_Vm,length(y_pref_Vm));
%    samp_err_ypv = [samp_err_ypv;mean(samp_ypref_vm)];
%    samp_ypref_spike = datasample(y_pref_spike,length(y_pref_spike));
%    samp_err_yps = [samp_err_yps;mean(samp_ypref_spike)];
%    samp_ynull_vm = datasample(y_null_Vm,length(y_null_Vm));
%    samp_err_ynv = [samp_err_ynv;mean(samp_ynull_vm)];
%    samp_ynull_spike = datasample(y_null_spike,length(y_null_spike));
%    samp_err_yns = [samp_err_yns;mean(samp_ynull_spike)];
%end
%
%se_opv = (1/(nsims-1))*sum((samp_err_opv-mean(samp_err_opv)).^2);
%se_ops = (1/(nsims-1))*sum((samp_err_ops-mean(samp_err_ops)).^2);
%se_onv = (1/(nsims-1))*sum((samp_err_onv-mean(samp_err_onv)).^2);
%se_ons = (1/(nsims-1))*sum((samp_err_ons-mean(samp_err_ons)).^2);
%se_ypv = (1/(nsims-1))*sum((samp_err_ypv-mean(samp_err_ypv)).^2);
%se_yps = (1/(nsims-1))*sum((samp_err_yps-mean(samp_err_yps)).^2);
%se_ynv = (1/(nsims-1))*sum((samp_err_ynv-mean(samp_err_ynv)).^2);
%se_yns = (1/(nsims-1))*sum((samp_err_yns-mean(samp_err_yns)).^2);

load('variables.mat');


y_predicted_FRnull = [];
y_predicted_FRpref = [];
y_pre_FRnull_up = [];
y_pre_FRnull_down = [];
y_pre_FRpref_up = [];
y_pre_FRpref_down = [];
y_predicted_DSI = [];
y_actual_FRnull = [];
y_actual_FRpref = [];
yAR_order = [];

o_predicted_FRnull = [];
o_predicted_FRpref = [];
o_pre_FRnull_up = [];
o_pre_FRnull_down = [];
o_pre_FRpref_up = [];
o_pre_FRpref_down = [];
o_predicted_DSI = [];
o_actual_FRnull = [];
o_actual_FRpref = [];
oAR_order = [];

if allow_flip == 1,
    for j = 1:length(young_cells),
        translate_name = num2str(str2num(regexprep(young_cells{j,1},'cell_','')));
        yAR_arrayloc = find(ismember(cell_young,translate_name));
        yAR_order(end+1,1) = yAR_arrayloc;
        present_Vmnull = vm_null_young(1,yAR_arrayloc);
        present_Vmpref = vm_pref_young(1,yAR_arrayloc);
        y_actual_FRnull(end+1,1) = spike_null_young(1,yAR_arrayloc);
        y_actual_FRpref(end+1,1) = spike_pref_young(1,yAR_arrayloc);
    end
    
    for k = 1:length(old_cells),
        translate_name = num2str(str2num(regexprep(old_cells{k,1},'cell_','')));
        oAR_arrayloc = find(ismember(cell_old,translate_name));
        oAR_order(end+1,1) = oAR_arrayloc;
        o_actual_FRnull(end+1,1) = spike_null_old(1,oAR_arrayloc);
        o_actual_FRpref(end+1,1) = spike_pref_old(1,oAR_arrayloc);
    end
else
    for j = 1:length(young_cells),
        translate_name = num2str(str2num(regexprep(young_cells{j,1},'cell_','')));
        yAR_arrayloc = find(ismember(cell_young,translate_name));
        yAR_order(end+1,1) = yAR_arrayloc;
        present_Vmnull = vm_null_young(1,yAR_arrayloc);
        present_Vmpref = vm_pref_young(1,yAR_arrayloc);
        y_predicted_FRnull(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
            young_params{j,1}.tanh_params.b.*tanh((present_Vmnull-young_params{j,1}.tanh_params.c)./...
            young_params{j,1}.tanh_params.d));
        y_pre_FRnull_up(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
            young_params{j,1}.tanh_params.b.*tanh(((present_Vmnull+se_ynv(j,1))-young_params{j,1}.tanh_params.c)./...
            young_params{j,1}.tanh_params.d));
        y_pre_FRnull_down(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
            young_params{j,1}.tanh_params.b.*tanh(((present_Vmnull-se_ynv(j,1))-young_params{j,1}.tanh_params.c)./...
            young_params{j,1}.tanh_params.d));
        y_predicted_FRpref(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
            young_params{j,1}.tanh_params.b.*tanh((present_Vmpref-young_params{j,1}.tanh_params.c)./...
            young_params{j,1}.tanh_params.d));
        y_pre_FRpref_up(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
            young_params{j,1}.tanh_params.b.*tanh(((present_Vmpref+se_ypv(j,1))-young_params{j,1}.tanh_params.c)./...
            young_params{j,1}.tanh_params.d));
        y_pre_FRpref_down(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
            young_params{j,1}.tanh_params.b.*tanh(((present_Vmpref-se_ypv(j,1))-young_params{j,1}.tanh_params.c)./...
            young_params{j,1}.tanh_params.d));
        y_predicted_DSI(end+1,1) = (y_predicted_FRpref(j,1)-y_predicted_FRnull(j,1))/y_predicted_FRpref(j,1);
        y_actual_FRnull(end+1,1) = spike_null_young(1,yAR_arrayloc);
        y_actual_FRpref(end+1,1) = spike_pref_young(1,yAR_arrayloc);
    end


    for k = 1:length(old_cells),
        translate_name = num2str(str2num(regexprep(old_cells{k,1},'cell_','')));
        oAR_arrayloc = find(ismember(cell_old,translate_name));
        oAR_order(end+1,1) = oAR_arrayloc;
        present_Vmnull_ = vm_null_old(1,oAR_arrayloc);
        present_Vmpref_ = vm_pref_old(1,oAR_arrayloc);
        o_predicted_FRnull(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
            old_params{k,1}.tanh_params.b.*tanh((present_Vmnull_-old_params{k,1}.tanh_params.c)./...
            old_params{k,1}.tanh_params.d));
        o_pre_FRnull_up(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
            old_params{k,1}.tanh_params.b.*tanh(((present_Vmnull_+se_onv(k,1))-old_params{k,1}.tanh_params.c)./...
            old_params{k,1}.tanh_params.d));
        o_pre_FRnull_down(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
            old_params{k,1}.tanh_params.b.*tanh(((present_Vmnull_-se_onv(k,1))-old_params{k,1}.tanh_params.c)./...
            old_params{k,1}.tanh_params.d));
        o_predicted_FRpref(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
            old_params{k,1}.tanh_params.b.*tanh((present_Vmpref_-old_params{k,1}.tanh_params.c)./...
            old_params{k,1}.tanh_params.d));
        o_pre_FRpref_up(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
            old_params{k,1}.tanh_params.b.*tanh(((present_Vmpref_+se_opv(k,1))-old_params{k,1}.tanh_params.c)./...
            old_params{k,1}.tanh_params.d));
        o_pre_FRpref_down(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
            old_params{k,1}.tanh_params.b.*tanh(((present_Vmpref_-se_opv(k,1))-old_params{k,1}.tanh_params.c)./...
            old_params{k,1}.tanh_params.d));
        o_predicted_DSI(end+1,1) = (o_predicted_FRpref(k,1)-o_predicted_FRnull(k,1))/o_predicted_FRpref(k,1);
        o_actual_FRnull(end+1,1) = spike_null_old(1,oAR_arrayloc);
        o_actual_FRpref(end+1,1) = spike_pref_old(1,oAR_arrayloc);
    end
end

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
    end
    for m = 1:length(spike_dsi_old),
        e = line([spike_dsi_old(1,m) spike_dsi_old(1,m)],...
            [m_model_odsi(m)-se_model_odsi(m) m_model_odsi(m)+se_model_odsi(m)]);
        e.Color = [0.6 0 0.6];
        e.LineWidth = 2.0;
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

end

