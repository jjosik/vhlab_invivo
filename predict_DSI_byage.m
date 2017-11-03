function [ output_args ] = predict_DSI_byage( plot_cells,bin_cellX,use_trial_baseline )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

sub_folder = uigetdir();
d = dir(sub_folder);
isub = [d(:).isdir];
sel_Snames = {d(isub).name}';
sel_Snames(ismember(sel_Snames,{'.','..','Codes'})) = [];
[s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',sel_Snames);

directories = char(sel_Snames(s));
cd(sub_folder);

old_means = [];
old_bins = [];
young_means = [];
young_bins = [];
for i = 1:size(directories,1),
    [collated_bin_ymean,collated_bin_centers,...
        case_simpleslope,case_maxderiv,case_xinter,...
        case_cellTH,case_prefTH,case_nullTH,params_array,cell_names] =...
        get_VFdata(plot_cells,bin_cellX);
    if strcmp(directories(i,:),'old'),
        old_params = params_array;
        old_cells = cell_names;
    else
        young_params = params_array;
        young_cells = cell_names;
    end
end
cd ../
save('collected_params.mat','old_params','old_cells','young_params','young_cells');
load('variables.mat');
y_predicted_FRnull = [];
y_predicted_FRpref = [];
y_predicted_DSI = [];
y_actual_FRnull = [];
y_actual_FRpref = [];
yAR_order = [];
for j = 1:length(young_cells),
    translate_name = num2str(str2num(regexprep(young_cells{j,1},'cell_','')));
    yAR_arrayloc = find(ismember(cell_young,translate_name));
    yAR_order(end+1,1) = yAR_arrayloc;
    present_Vmnull = vm_null_young(1,yAR_arrayloc);
    present_Vmpref = vm_pref_young(1,yAR_arrayloc);
    y_predicted_FRnull(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
        young_params{j,1}.tanh_params.b.*tanh((present_Vmnull-young_params{j,1}.tanh_params.c)./...
        young_params{j,1}.tanh_params.d));
    y_predicted_FRpref(end+1,1) = rectify(young_params{j,1}.tanh_params.a +...
        young_params{j,1}.tanh_params.b.*tanh((present_Vmpref-young_params{j,1}.tanh_params.c)./...
        young_params{j,1}.tanh_params.d));
    y_predicted_DSI(end+1,1) = (y_predicted_FRpref(j,1)-y_predicted_FRnull(j,1))/y_predicted_FRpref(j,1);
    y_actual_FRnull(end+1,1) = spike_null_young(1,yAR_arrayloc);
    y_actual_FRpref(end+1,1) = spike_pref_young(1,yAR_arrayloc);
end

o_predicted_FRnull = [];
o_predicted_FRpref = [];
o_predicted_DSI = [];
o_actual_FRnull = [];
o_actual_FRpref = [];
oAR_order = [];
for k = 1:length(old_cells),
    translate_name = num2str(str2num(regexprep(old_cells{k,1},'cell_','')));
    oAR_arrayloc = find(ismember(cell_old,translate_name));
    oAR_order(end+1,1) = oAR_arrayloc;
    present_Vmnull_ = vm_null_old(1,oAR_arrayloc);
    present_Vmpref_ = vm_pref_old(1,oAR_arrayloc);
    o_predicted_FRnull(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
        old_params{k,1}.tanh_params.b.*tanh((present_Vmnull_-old_params{k,1}.tanh_params.c)./...
        old_params{k,1}.tanh_params.d));
    o_predicted_FRpref(end+1,1) = rectify(old_params{k,1}.tanh_params.a +...
        old_params{k,1}.tanh_params.b.*tanh((present_Vmpref_-old_params{k,1}.tanh_params.c)./...
        old_params{k,1}.tanh_params.d));
    o_predicted_DSI(end+1,1) = (o_predicted_FRpref(k,1)-o_predicted_FRnull(k,1))/o_predicted_FRpref(k,1);
    o_actual_FRnull(end+1,1) = spike_null_old(1,oAR_arrayloc);
    o_actual_FRpref(end+1,1) = spike_pref_old(1,oAR_arrayloc);
end

%plot DSI actual vs. predicted 
f = figure;
subplot(2,2,1);
scatter(spike_dsi_young(reshape(yAR_order,1,length(yAR_order))),y_predicted_DSI,'g','filled');
hold on;
unl = line([0 1.5],[0 1.5]);
unl.LineStyle = '--';
unl.Color = [0 0 0];
xlabel('Measured DSI');
ylabel('Predicted DSI');
try
    sdy = spike_dsi_young(reshape(yAR_order,length(yAR_order),1));
    [yb,yb_ci,yresid,yrint,ystats] = ...
        regress(y_predicted_DSI,horzcat(ones(length(sdy),1),sdy'));
catch
    [yslope,yoffset,yslope_ci,yresid,yrint,ystats] = ...
        quickregression(spike_dsi_young(reshape(yAR_order,length(yAR_order),1)),...
        y_predicted_DSI,0.05);
end
title(['YOUNG, R^2 = ',num2str(ystats(1,1)),', p = ',num2str(ystats(1,3))]);
hold off;
subplot(2,2,2);
scatter(spike_dsi_old(reshape(oAR_order,1,length(oAR_order))),o_predicted_DSI,'m','filled');
hold on;
unl = line([0 1.5],[0 1.5]);
unl.LineStyle = '--';
unl.Color = [0 0 0];
xlabel('Measured DSI');
ylabel('Predicted DSI');
try
    sdo = spike_dsi_old(reshape(oAR_order,length(oAR_order),1));
    [ob,ob_ci,oresid,orint,ostats] = ...
        regress(o_predicted_DSI,horzcat(ones(length(sdo),1),sdo'));
catch
    [ob,ob_ci,oresid,orint,ostats] = ...
        quickregression(spike_dsi_old(reshape(oAR_order,length(oAR_order),1)),...
        o_predicted_DSI,0.05);
end
title(['OLD, R^2 = ',num2str(ostats(1,1)),', p = ',num2str(ostats(1,3))]);
hold off;
subplot(2,2,3);
scatter(y_actual_FRnull,y_predicted_FRnull,'b','filled');
hold on;
scatter(y_actual_FRpref,y_predicted_FRpref,'r','filled');
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
subplot(2,2,4);
scatter(o_actual_FRnull,o_predicted_FRnull,'b','filled');
hold on;
scatter(o_actual_FRpref,o_predicted_FRpref,'r','filled');
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
