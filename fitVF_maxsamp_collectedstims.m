function [ Vth_est_coll ] = fitVF_maxsamp_collectedstims( FR_inst,Vm_array,stimvalues,...
    model,h_margins,sampling_rate,tempFrequency,adapt_bins,anchor_Vth,fit_it,save_it,code)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
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

stimset = 1;
z_scale = 1;
match_bins = 1;
[ f80,ph,bin_edges ] = VF_densitymap_C(coll_Vm,coll_FR,adapt_bins,...
    stimset,z_scale,match_bins,h_margins);
cb = colorbar;
ylabel(cb,'log count of time-step occurances');
if save_it == 1,
    title(['Maximally-sampled V-F plot'...
        'for all collated stim. directions']);
    saveas(gcf,[pwd filesep strcat(code,'subTH_rawVF_density_COLL.fig')]);
    %close(f9);
else
end
if fit_it == 1,
    sub_elem = find(coll_FR==0);
    Vth_est_coll = mean(coll_Vm(sub_elem));
    weight = 1;
    switch model
        case 0
            [n_x,m_out,cl,gof,fitinfo] = vf_powerfit(coll_Vm,coll_FR,...
                anchor_Vth,Vth_est_coll,weight);
        case 1
            [n_x,m_out,cl,gof,fitinfo] = vf_linfit(coll_Vm,coll_FR,...
                anchor_Vth,Vth_est_coll,weight);
    end
    z_max = max(max(get(ph,'ZData')));
    n_h = line(n_x,m_out,z_max*ones(length(n_x),1));
    n_h.LineWidth = 3.0;
    n_h.Color = [0 0 0];
    switch model
        case 0
            title({'All stims. max-sampled V-F plot ',...
                'a*rectify(Vm-Vth)^b fit parameters:  Vth = ',...
                num2str(Vth_est_coll),...
                ' a = ', num2str(cl.a), ' b = ' num2str(cl.b)});
        case 1
            title({'All stims. max-sampled V-F plot ',...
                'a*rectify(Vm-Vth) fit parameters: Vth = ',...
                num2str(Vth_est_coll),' a = ', num2str(cl.a)});
    end
    %mn_line = line([m_pref m_pref],[0 max(pm_out)]);  blocked 3/17
    %make local regression smoothing comparison
    %[x,inds] = sort(coll_pref_Vm);
    %yy = smooth(x,coll_pref_FR(inds),5000,'lowess');
    %plot(x,yy,'b-');
    %make gaussian fits across y, smooth means in x bins
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
            bin_ymean(end+1,1) = nanmean(coll_FR(bin_inds,1));
            bin_ysd(end+1,1) = std(coll_FR(bin_inds,1));
            bin_yse(end+1,1) = std(coll_FR(bin_inds,1))./sqrt(length(bin_inds));
            x_bins = [x_bins;i-1];
        else
        end
    end
    bin_centers = (bin_edges{1,1}(1,2:end)+bin_edges{1,1}(1,1:end-1))/2;
    scatter3(bin_centers(x_bins),bin_ymean,...
        (abs(z_max+1).*ones(length(bin_ymean),1)));
    for j = 1:length(bin_ymean),
        eb = line([bin_centers(x_bins(j,1)) bin_centers(x_bins(j,1))],...
            [bin_ymean(j,1)-bin_ysd(j,1) bin_ymean(j,1)+bin_ysd(j,1)]);
        eb.LineWidth = 2.0;
        eb.Color = 'r';
    end
    legend('density','fit','1 mV-binned mean Vm','1 mV-binned Vm SD',...
        'Location','NorthEastOutside');
    hold off;
else
end
if save_it == 1,
    saveas(gcf,[pwd filesep strcat(code,'subTH_maxsampled_VF_density_COLL.fig')]);
    close(f80);
else
end
%create secondary plot of X-binned means
fel = figure;
bin_centers = (bin_edges{1,1}(1,2:end)+...
    bin_edges{1,1}(1,1:end-1))/2;
scatter(bin_centers(x_bins),bin_ymean,'filled','k');
for i5 = 1:length(bin_ymean),
    eb = line([bin_centers(x_bins(i5,1)) bin_centers(x_bins(i5,1))],...
        [bin_ymean(i5,1)-bin_yse(i5,1) bin_ymean(i5,1)+bin_yse(i5,1)]);
    eb.LineWidth = 2.0;
    eb.Color = 'k';
end
xlabel('Membrane potential (mV)');
ylabel('Firing rate');
if save_it == 1,
    saveas(gcf,[pwd filesep strcat(code,'_Xbinned_VF_COLL.fig')]);
    close(fel);
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
    saveas(gcf,[pwd filesep strcat(code,'COLL_subVm_vs_VmDist','.fig')]);
    close(fc);
else
end

if save_it == 1,
    vffilename = strcat(code,'_collstim_fitdata.mat');
    save(vffilename);
else
end

end

