function [ Vth_est_coll,params,gof_obj,allfits_info ] = fitVF_maxsamp_Xbins( FR_inst,Vm_array,stimvalues,...
    model,h_margins,sampling_rate,tempFrequency,adapt_bins,anchor_Vth,...
    fit_it,save_it,code,use_trial_baseline,cellmean_corr_th,...
    trialmean_corr_th,pref_ind,null_ind)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

params = struct('power_params',[],'linear_params',[],'tanh_params',[]);
gof_obj = struct('power_gof',[],'linear_gof',[],'tanh_gof',[]);
allfits_info = struct('power_info',[],'linear_info',[],'tanh_info',[]);
try
    stim_tf = p.tFrequency;
catch
    stim_tf = tempFrequency;
end

all_dirind = 1:length(stimvalues)-1;
 %modify existing Vm traces
for ii = 1:length(Vm_array),
    new_Vm_array{ii,1} = Vm_array{ii,1}(1:...
        (length(Vm_array{ii,1})-(round(h_margins*sampling_rate))));
end

array_disp = abs(length(FR_inst{1,1})-length(new_Vm_array{1,1}));
if ~(length(FR_inst{1,1})==length(new_Vm_array{1,1})),
    coll_FR = [];
    for n = 1:length(FR_inst)-1,
        coll_FR = [coll_FR;reshape(FR_inst{n,1},length(FR_inst{n,1}),1)];
    end
    coll_Vm = [];
    for m = 1:length(new_Vm_array)-1,
        coll_Vm = [coll_Vm;new_Vm_array{m,1}(1:end-array_disp,1)];
    end
else
    coll_FR = [];
    for n = 1:length(FR_inst)-1,
        coll_FR = [coll_FR;reshape(FR_inst{n,1},length(FR_inst{n,1}),1)];
    end
    coll_Vm = [];
    for m = 1:length(new_Vm_array)-1,
        coll_Vm = [coll_Vm;new_Vm_array{m,1}];
    end
end

f = figure;
stimset = 1;
scale_z = 1;
match_bins = 1;
if adapt_bins == 1,
    loc_min = min(coll_Vm);
    loc_max = max(coll_Vm);
    low = round(loc_min/10)*10;
    high = low+50;
    if loc_max>high,
        high = ceil(loc_max/10)*10;
    else
    end
    v_bin_edges = low:1:high;
    if length(v_bin_edges)==51,
        f_bin_edges = 0:5:250;
    else
        if mod(250,length(v_bin_edges)-1)==0,
            f_bin_edges = 0:250/(length(v_bin_edges)-1):250;
        else
            link_div = ceil(250/length(v_bin_edges)-1);
            f_bin_edges = 0:link_div:(link_div*(length(v_bin_edges)-1));
        end
    end
    bin_edges = {v_bin_edges f_bin_edges};
    subplot(2,2,1);
    if iscell(coll_Vm),
        hc = hist3([coll_Vm{stimset,1}(1:end-1,1),coll_FR{stimset,1}(1:end,1)],'Edges',bin_edges);
    else
        hc = hist3([coll_Vm(1:end-1,1),coll_FR(1:end-1)],'Edges',bin_edges);
    end
else
    N_bins = [50,50];
    subplot(2,2,1);
    if iscell(coll_Vm),
        hc = hist3([coll_Vm{stimset,1}(1:end-1,1),coll_FR{stimset,1}(1:end,1)],N_bins);
    else
        hc = hist3([coll_Vm(1:end-1,1),coll_FR(1:end-1,1)],N_bins);
    end
end
if scale_z == 1,
    hc1 = 1+log10(hc');
else
    hc1 = hc';
end
hc1(size(hc,1)+1,size(hc,2)+1) = 0;
if iscell(coll_Vm),
    xc = linspace(min(coll_Vm{stimset,1}(1:end-1,1)),max(coll_Vm{stimset,1}(1:end-1,1)),size(hc,1)+1);
    yc = linspace(min(coll_FR{stimset,1}(1:end,1)),max(coll_FR{stimset,1}(1:end,1)),size(hc,1)+1);
else
    coll_FR = reshape(coll_FR,1,length(coll_FR));
    xc = linspace(min(coll_Vm(1:end-1,1)),max(coll_Vm(1:end-1,1)),size(hc,1)+1);
    yc = linspace(min(coll_FR(1,1:end-1)),max(coll_FR(1,1:end-1)),size(hc,1)+1);
end
h = pcolor(xc,yc,hc1);
xlabel('Membrane potential (mV)');
ylabel('Inst. Firing rate');
axis square;
hold on;
    
cb = colorbar;
ylabel(cb,'log count of time-steps per joint event bin');
if fit_it == 1,
    sub_elem = find(coll_FR==0);
    Vth_est_coll = mean(coll_Vm(sub_elem));
    weight = 1;
    %run the three fits consecutively:
    %1. rectified power law
    [n_x,m_out,cl,gof,fitinfo] = vf_powerfit_ALT(coll_Vm,coll_FR,anchor_Vth,...
        Vth_est_coll,weight);
    %2. rectilinear
    [ln_x,lm_out,lcl,lgof,lfitinfo] = vf_linfit_ALT(coll_Vm,coll_FR,anchor_Vth,...
        Vth_est_coll,weight);
    %3. hyperbolic tangent
    [s_x,tm_out,tcl,tgof,tfitinfo] = alt_tanhfit(coll_Vm,coll_FR,use_trial_baseline);
    z_max = max(max(get(h,'ZData')));
    n_h = line(n_x,m_out,z_max*ones(length(n_x),1));
    n_h.LineWidth = 3.0;
    n_h.Color = [0 0 0];
    n_h.LineStyle = '-';
    ln_h = line(ln_x,lm_out,z_max*ones(length(ln_x),1));
    ln_h.LineWidth = 3.0;
    ln_h.Color = [0 0 0];
    ln_h.LineStyle = '--';
    s_h = line(s_x,tm_out,z_max*ones(length(s_x),1));
    s_h.LineWidth = 3.0;
    s_h.Color = [0 0 0];
    s_h.LineStyle = ':';
    %legend('Density','Rectified power law','Rectilinear','Tanh','Location','NorthEastOutside');
    %th_line = line([cellmean_corr_th cellmean_corr_th],[0 max(m_out)]);
    %th_line.LineWidth = 2.0;
    %th_line.Color = [0 0 1];
    hold off;
    
    %X-binned plots
    x_bins = [];
    bin_ymean = [];
    bin_ysd = [];
    bin_yse = [];
    bin_count = [];
    for i = 2:length(bin_edges{1,1}),
        if i < length(bin_edges{1,1}),
            bin_inds = find(coll_Vm>=bin_edges{1,1}(1,i-1)&...
                coll_Vm<bin_edges{1,1}(1,i));
        else
            bin_inds = find(coll_Vm>=bin_edges{1,1}(1,i-1)&...
                coll_Vm<=bin_edges{1,1}(1,i));
        end
        bin_count(end+1,1) = length(bin_inds);
        if ~isempty(bin_inds),
            bin_ymean(end+1,1) = nanmean(coll_FR(1,bin_inds));
            bin_ysd(end+1,1) = std(coll_FR(1,bin_inds));
            bin_yse(end+1,1) = std(coll_FR(1,bin_inds))./sqrt(length(bin_inds));
            x_bins = [x_bins;i-1];
        else
        end
    end
    bin_centers = (bin_edges{1,1}(1,2:end)+bin_edges{1,1}(1,1:end-1))/2;
    subplot(2,2,2);
    scatter(bin_centers(x_bins),bin_ymean,'r');
    for j = 1:length(bin_ymean),
        eb = line([bin_centers(x_bins(j,1)) bin_centers(x_bins(j,1))],...
            [bin_ymean(j,1)-bin_ysd(j,1) bin_ymean(j,1)+bin_ysd(j,1)]);
        eb.LineWidth = 2.0;
        eb.Color = [0.4 0.4 0.4];
    end
    hold on;
    %fit rectified power law to means
    m_anchor_Vth = 0;
    m_Vth_est = NaN;
    m_weight = 1;
    [mn_x,mm_out,mcl,mgof,mfitinfo] = vf_powerfit_ALT(bin_centers(x_bins)',bin_ymean,m_anchor_Vth,...
        m_Vth_est,m_weight);
    mn_h = line(mn_x,mm_out);
    mn_h.LineWidth = 3.0;
    mn_h.Color = [0 0 0];
    %cell_th_line = line([cellmean_corr_th cellmean_corr_th],[-10 70]);
    %cell_th_line.LineWidth = 2.0;
    %cell_th_line.Color = [0 0 0];
    %pc_th_line = line([trialmean_corr_th(pref_ind,1) trialmean_corr_th(pref_ind,1)],...
    %    [-10 70]);
    %pc_th_line.LineWidth = 2.0;
    %pc_th_line.Color = [1 0 0];
    %nc_th_line = line([trialmean_corr_th(null_ind,1) trialmean_corr_th(null_ind,1)],...
    %    [-10 70]);
    %nc_th_line.LineWidth = 2.0;
    %nc_th_line.Color = [0 0 1];
    xlabel('Membrane potential (mV)');
    ylabel('Mean Firing rate');
    title(['Power law fit: ',num2str(mcl.a),'.*rectify(x - ',num2str(mcl.b),').^^',num2str(mcl.c)]);
    hold off;
    subplot(2,2,3);
    scatter(bin_centers(x_bins),bin_ymean,'r');
    for j = 1:length(bin_ymean),
        eb = line([bin_centers(x_bins(j,1)) bin_centers(x_bins(j,1))],...
            [bin_ymean(j,1)-bin_ysd(j,1) bin_ymean(j,1)+bin_ysd(j,1)]);
        eb.LineWidth = 2.0;
        eb.Color = [0.4 0.4 0.4];
    end
    hold on;
    %rectilinear fit to means
    [ml_x,ml_out,mlcl,mlgof,mlfitinfo] = vf_linfit_ALT(bin_centers(x_bins)',bin_ymean,m_anchor_Vth,...
        m_Vth_est,m_weight);
    ml_h = line(ml_x,ml_out);
    ml_h.LineWidth = 3.0;
    ml_h.Color = [0 0 0];
    %cell_th_line = line([cellmean_corr_th cellmean_corr_th],[-10 70]);
    %cell_th_line.LineWidth = 2.0;
    %cell_th_line.Color = [0 0 0];
    %pc_th_line = line([trialmean_corr_th(pref_ind,1) trialmean_corr_th(pref_ind,1)],...
    %    [-10 70]);
    %pc_th_line.LineWidth = 2.0;
    %pc_th_line.Color = [1 0 0];
    %nc_th_line = line([trialmean_corr_th(null_ind,1) trialmean_corr_th(null_ind,1)],...
    %    [-10 70]);
    %nc_th_line.LineWidth = 2.0;
    %nc_th_line.Color = [0 0 1];
    xlabel('Membrane potential (mV)');
    ylabel('Mean Firing rate');
    title(['Rectilinear fit: ',num2str(mlcl.a),'.*rectify(x - ',num2str(mlcl.b),')']);
    hold off;
    subplot(2,2,4);
    scatter(bin_centers(x_bins),bin_ymean,'r');
    for j = 1:length(bin_ymean),
        eb = line([bin_centers(x_bins(j,1)) bin_centers(x_bins(j,1))],...
            [bin_ymean(j,1)-bin_ysd(j,1) bin_ymean(j,1)+bin_ysd(j,1)]);
        eb.LineWidth = 2.0;
        eb.Color = [0.4 0.4 0.4];
    end
    hold on;
    %fit hyperbolic tangent to means
    [ms_x,mt_out,mtcl,mtgof,mtfitinfo] = alt_tanhfit(bin_centers(x_bins)',bin_ymean,use_trial_baseline);
    ms_h = line(ms_x,mt_out);
    ms_h.LineWidth = 3.0;
    ms_h.Color = [0 0 0];
    %cell_th_line = line([cellmean_corr_th cellmean_corr_th],[-10 70]);
    %cell_th_line.LineWidth = 2.0;
    %cell_th_line.Color = [0 0 0];
    %pc_th_line = line([trialmean_corr_th(pref_ind,1) trialmean_corr_th(pref_ind,1)],...
    %    [-10 70]);
    %pc_th_line.LineWidth = 2.0;
    %pc_th_line.Color = [1 0 0];
    %nc_th_line = line([trialmean_corr_th(null_ind,1) trialmean_corr_th(null_ind,1)],...
    %    [-10 70]);
    %nc_th_line.LineWidth = 2.0;
    %nc_th_line.Color = [0 0 1];
    xlabel('Membrane potential (mV)');
    ylabel('Mean Firing rate');
    title(['Saturating fit: ',num2str(mtcl.a),' + ',num2str(mtcl.b),'.*tanh((x - ',num2str(mtcl.c),')./',num2str(mtcl.d),')']);
    hold off;
    if save_it == 1,
        saveas(gcf,[pwd filesep strcat(code,'_VF_subplots_Vm10.fig')]);
        close(f);
    else
    end
    
    prefmean_corr_th = trialmean_corr_th(pref_ind,1);
    nullmean_corr_th = trialmean_corr_th(null_ind,1);
    
    params.power_params = mcl;
    params.linear_params = mlcl;
    params.tanh_params = mtcl;
    gof_obj.power_gof = mgof;
    gof_obj.linear_gof = mlgof;
    gof_obj.tanh_gof = mtgof;
    allfits_info.power_info = mfitinfo;
    allfits_info.linear_info = mlfitinfo;
    allfits_info.tanh_info = mtfitinfo;
else
end

%plot comparison of FR = 0 and FR > 0 Vm-distributions
sub_ind = find(coll_FR == 0);
sub_Vm = coll_Vm(sub_ind);
sh = histc(sub_Vm,bin_edges{1,1});
sh = sh(1:end-1);
fc = figure;
bar(bin_centers,bin_count);
grid on;
grid minor;
title('Response Vm distribution');
ylabel('Counts');
xlabel('Potential (mV)');
hold on;
bar(bin_centers,sh,'c');
legend('All Vm','Subthreshold Vm only','Location','EastOutside');
if save_it == 1,
    saveas(gcf,[pwd filesep strcat(code,'_COLL_subVm_vs_VmDist','.fig')]);
    close(fc);
else
end

if save_it == 1,
    vffilename = strcat(code,'_collstim_Xbinfitdata.mat');
    save(vffilename);
else
end
    
    
    
    


end

