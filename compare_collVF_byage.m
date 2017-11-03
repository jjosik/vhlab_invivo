function [ ogof,ygof ] = compare_collVF_byage( plot_cells,bin_cellX,use_trial_baseline,use_gridsearch )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
%NOTE - need to fix this directory/subdirectory walk; currently 'old' and
%'young' directories have to run separately (i.e. twice manually), with the second
%run skipping over the array initialization (***the proper structure is
%already in use in the STA codes)

[folder,directories] = walk_directories();
cd(folder);

old_means = [];
old_bins = [];
young_means = [];
young_bins = [];
o_coll_Vm = [];
o_coll_FR = [];
y_coll_Vm = [];
y_coll_FR = [];
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
    collated_bin_ymean = [];
    collated_bin_centers = [];
    if bin_cellX,
        case_simpleslope = [];
        case_maxderiv = [];
        case_xinter = [];
        case_cellTH = [];
        case_prefTH = [];
        case_nullTH = [];
        params_array = {};
        cell_names = {};
    else
    end
    for j = 1:size(sub_dir,1),
        cur_dir = sub_dir(j,:);
        cd(cur_dir);
      
        %[collated_bin_ymean,collated_bin_centers,...
        %    case_simpleslope,case_maxderiv,case_xinter,...
        %    case_cellTH,case_prefTH,case_nullTH] = get_VFdata(plot_cells,bin_cellX);
        dd = dir(fullfile(cur_dir,'*.mat'));
        cfilename = uigetfile('*.mat');
        if cfilename == 0,
            cd ..
            continue
        else
        end
        if bin_cellX == 1,
        load([sub_folder filesep cur_dir filesep cfilename],...
            'coll_FR','coll_Vm','bin_edges','cellmean_corr_th',...
            'prefmean_corr_th','nullmean_corr_th','params',...
            'gof_obj','allfits_info','ms_x','mt_out');
        else
            load([sub_folder filesep cur_dir filesep cfilename],...
            'coll_FR','coll_Vm','bin_edges','cl');
        end
        temp_ymeans = [];
        temp_xbins = [];
        bin_count = [];
        for j = 2:length(bin_edges{1,1}),
            if j < length(bin_edges{1,1}),
                bin_inds = find(coll_Vm>=bin_edges{1,1}(1,j-1)&...
                    coll_Vm<bin_edges{1,1}(1,j));
            else
                bin_inds = find(coll_Vm>=bin_edges{1,1}(1,j-1)&...
                    coll_Vm<=bin_edges{1,1}(1,j));
            end
            bin_count(end+1,1) = length(bin_inds);
            if ~isempty(bin_inds),
                temp_ymeans(end+1,1) = nanmean(coll_FR(bin_inds));
                temp_xbins = [temp_xbins;j-1];
            else
            end
        end
        gen_bin_centers = (bin_edges{1,1}(1,2:end)+bin_edges{1,1}(1,1:end-1))/2;
        if bin_cellX == 1,
            bin_centers = gen_bin_centers(temp_xbins);
        else
            bin_centers = gen_bin_centers(temp_xbins)-cl.Vth_est;
        end
        collated_bin_ymean = [collated_bin_ymean;reshape(temp_ymeans,length(temp_ymeans),1)];
        collated_bin_centers = [collated_bin_centers;reshape(bin_centers,length(bin_centers),1)];
        
        if bin_cellX,
        
        %parameter arrays
            case_simpleslope(end+1,1) = params.tanh_params.b/params.tanh_params.d;
            case_maxderiv(end+1,1) = max(gradient(mt_out,mean(diff(ms_x))));
            case_xinter(end+1,1) = params.tanh_params.c+...
                params.tanh_params.d*(atanh(-(params.tanh_params.a/params.tanh_params.b)));
            case_cellTH(end+1,1) = cellmean_corr_th;
            case_prefTH(end+1,1) = prefmean_corr_th;
            case_nullTH(end+1,1) = nullmean_corr_th;
            params_array{end+1,1} = params;
            cell_names{end+1,1} = directories(i,:);
        else
        end
        if strcmp(cellstr(directories(i,:)),'old'),
            o_coll_Vm = [o_coll_Vm;coll_Vm];
            o_coll_FR = [o_coll_FR,coll_FR];
        else
            y_coll_Vm = [y_coll_Vm;coll_Vm];
            y_coll_FR = [y_coll_FR,coll_FR];
        end
        
        cd ..
    end
    
    if strcmp(cellstr(directories(i,:)),'old'),
        old_means = collated_bin_ymean;
        old_bins = floor(collated_bin_centers)+0.5;
        old_simpleslope = case_simpleslope;
        old_maxderiv = case_maxderiv;
        old_xinter = case_xinter;
        old_cellTH = case_cellTH;
        old_prefTH = case_prefTH;
        old_nullTH = case_nullTH;
        %o_coll_Vm = [o_coll_Vm;coll_Vm];
        %o_coll_FR = [o_coll_FR,coll_FR];
    else
        young_means = collated_bin_ymean;
        young_bins = floor(collated_bin_centers)+0.5;
        young_simpleslope = case_simpleslope;
        young_maxderiv = case_maxderiv;
        young_xinter = case_xinter;
        young_cellTH = case_cellTH;
        young_prefTH = case_prefTH;
        young_nullTH = case_nullTH;
        %y_coll_Vm = [y_coll_Vm;coll_Vm];
        %y_coll_FR = [y_coll_FR,coll_FR];
    end
    clear coll_Vm coll_FR;
    cd ..
end


common_oldbins = unique(old_bins);
common_youngbins = unique(young_bins);
total_ymeans = [];
total_yse = [];
y_set = cell(length(common_youngbins),1);
for j = 1:length(common_youngbins),
    inst_loc = find(young_bins==common_youngbins(j,1));
    if isempty(inst_loc),
        total_ymeans = [total_ymeans;NaN];
        total_yse = [total_yse;NaN];
    else
        total_ymeans = [total_ymeans;mean(young_means(inst_loc))];
        total_yse = [total_yse;std(young_means(inst_loc))./sqrt(length(young_means(inst_loc)))];
        y_set{j,1} = young_means(inst_loc);
    end
end

total_omeans = [];
total_ose = [];
o_set = cell(length(common_oldbins),1);
for k = 1:length(common_oldbins),
    inst_loc_ = find(old_bins==common_oldbins(k,1));
    if isempty(inst_loc_),
        total_omeans = [total_omeans;NaN];
        total_ose = [total_ose;NaN];
    else
        total_omeans = [total_omeans;mean(old_means(inst_loc_))];
        total_ose = [total_ose;std(old_means(inst_loc_))./sqrt(length(old_means(inst_loc_)))];
        o_set{k,1} = old_means(inst_loc_);
    end
end

f1 = figure;
scatter(common_youngbins,total_ymeans,'r','filled');
hold on;
scatter(common_oldbins,total_omeans,'b','filled');
errorbar(common_youngbins,total_ymeans,total_yse,'.');
errorbar(common_oldbins,total_omeans,total_ose,'.');
xlabel('Relative Vm (mV)');
ylabel('Firing rate');
legend('young','old','Location','NorthWest');
common_bins = unique([common_youngbins;common_oldbins]);
y_multcomp = [];
y_multcomp_p = [];
bonf_alpha = 0.05/length(common_bins);
for n = 1:length(common_bins),
    if ~ismember(common_bins(n,1),common_youngbins)||~ismember(common_bins(n,1),common_oldbins),
        y_multcomp(end+1,1) = NaN;
        continue;
    else
        xy = find(common_youngbins==common_bins(n,1));
        xo = find(common_oldbins==common_bins(n,1));
        g1 = y_set{xy,1};
        g2 = o_set{xo,1};
        [h,p] = ttest2(g1,g2,'alpha',bonf_alpha);
        y_multcomp(end+1,1) = h;
        y_multcomp_p(end+1,1) = p;
    end
    clear g1 g2
end
y_multcomp(isnan(y_multcomp)) = 0;
start_bin = min(common_bins);
aster = find(y_multcomp);
ledger = 0.0;
if ~isempty(aster),
    scatter(start_bin+(aster-1),ledger.*ones(length(aster),1),'*');
else
end
%tanh fits
if use_gridsearch,
    reps = 500;
    use_random = 1;
    reinit_N = 5;
    [o_x,o_out,ocl,rmse,si] = tanfit_gridsearch(common_oldbins,total_omeans,...
        reinit_N,use_random,reps,use_trial_baseline);
    o_h = line(o_x,o_out);
    o_h.LineWidth = 3.0;
    o_h.Color = [0 0 1];
    [y_x,y_out,ycl,rmse,si] = tanfit_gridsearch(common_youngbins,total_ymeans,...
        reinit_N,use_random,reps,use_trial_baseline);
    y_h = line(y_x,y_out);
    y_h.LineWidth = 3.0;
    y_h.Color = [1 0 0];
    %title(['Old fit: ',num2str(ocl.a),' + ',num2str(ocl.b),'.*((Vm - ',num2str(ocl.c),')./',num2str(ocl.d),') ; Young fit: ',...
    %    num2str(ycl.a),' + ',num2str(ycl.b),'.*((Vm - ',num2str(ycl.c),')./',num2str(ycl.d),')']);
    hold off;
else
    [o_x,o_out,ocl,ogof,ofitinfo] = alt_tanhfit(common_oldbins,total_omeans,use_trial_baseline);
    o_h = line(o_x,o_out);
    o_h.LineWidth = 3.0;
    o_h.Color = [0 0 1];
    [y_x,y_out,ycl,ygof,yfitinfo] = alt_tanhfit(common_youngbins,total_ymeans,use_trial_baseline);
    y_h = line(y_x,y_out);
    y_h.LineWidth = 3.0;
    y_h.Color = [1 0 0];
    title(['Old fit: ',num2str(ocl.a),' + ',num2str(ocl.b),'.*((Vm - ',num2str(ocl.c),')./',num2str(ocl.d),') ; Young fit: ',...
        num2str(ycl.a),' + ',num2str(ycl.b),'.*((Vm - ',num2str(ycl.c),')./',num2str(ycl.d),')']);
    hold off;
end

f3 = figure;
ho = scatterhist(o_coll_Vm,o_coll_FR);
f4 = figure;
hy = scatterhist(y_coll_Vm,y_coll_FR);
cond_o_Vm = o_coll_Vm(find(o_coll_FR>0));
cond_o_FR = o_coll_FR(find(o_coll_FR>0));
cond_y_Vm = y_coll_Vm(find(y_coll_FR));
cond_y_FR = y_coll_FR(find(y_coll_FR));
f5 = figure;
hoc = scatterhist(cond_o_Vm,cond_o_FR);
f6 = figure;
hyc = scatterhist(cond_y_Vm,cond_y_FR);

keyboard

%collected parameter plots and stats
%old_nullTH(isnan(old_nullTH)) = 0; %*can't do this; it effectively means threshold is
%equal to baseline; must use nanmean/nanstd in null cases.  
%young_nullTH(isnan(young_nullTH)) = 0;
f2 = figure;
subplot(2,2,1);
scatter((1.5+(0.15-(rand(size(young_simpleslope))*0.3))).*ones(size(young_simpleslope)),young_simpleslope,'c','filled');
hold on;
yb = errorbar(1.5,mean(young_simpleslope),std(young_simpleslope)./sqrt(length(young_simpleslope)));
yb.Color = [0 0 0];
yb.LineWidth = 3.0;
scatter((3.0+(0.15-(rand(size(old_simpleslope))*0.3))).*ones(size(old_simpleslope)),old_simpleslope,'b','filled');
ob = errorbar(3.0,mean(old_simpleslope),std(old_simpleslope)./sqrt(length(old_simpleslope)));
ob.Color = [0 0 0];
ob.LineWidth = 3.0;
boxwhisker(young_simpleslope,1.5);
boxwhisker(old_simpleslope,3.0);
set(gca,'xtick',[1.5 3],'XTickLabel',{'Young','Old'});
xlim([0 4.5]);
[h,~] = kstest(young_simpleslope);
[h1,~] = kstest(old_simpleslope);
if h||h1,
    %[~,pk] = kstest2(young_simpleslope,old_simpleslope);
    [~,pk] = ranksum(young_simpleslope,old_simpleslope);
    p = pk;
else
    %[~,p_] = ttest2(young_simpleslope,old_simpleslope);
    [~,p_] = ranksum(young_simpleslope,old_simpleslope);
    p = p_;
end
title(['Tanh fit gain from parameters, WRS test: p = ',num2str(p)]);
ylabel('Slope ( Hz / mV )');
xlabel('Age group');
hold off;
subplot(2,2,2);
scatter((1.5+(0.15-(rand(size(young_maxderiv))*0.3))).*ones(size(young_maxderiv)),young_maxderiv,'c','filled');
hold on;
yb = errorbar(1.5,mean(young_maxderiv),std(young_maxderiv)./sqrt(length(young_maxderiv)));
yb.Color = [0 0 0];
yb.LineWidth = 3.0;
scatter((3.0+(0.15-(rand(size(old_maxderiv))*0.3))).*ones(size(old_maxderiv)),old_maxderiv,'b','filled');
ob = errorbar(3.0,mean(old_maxderiv),std(old_maxderiv)./sqrt(length(old_maxderiv)));
ob.Color = [0 0 0];
ob.LineWidth = 3.0;
boxwhisker(young_maxderiv,1.5);
boxwhisker(old_maxderiv,3.0);
set(gca,'xtick',[1.5 3],'XTickLabel',{'Young','Old'});
xlim([0 4.5]);
[h,~] = kstest(young_maxderiv);
[h1,~] = kstest(old_maxderiv);
if h||h1,
    %[~,pka] = kstest2(young_maxderiv,old_maxderiv);
    [~,pka] = ranksum(young_maxderiv,old_maxderiv);
    pa = pka;
else
    %[~,pa_] = ttest2(young_maxderiv,old_maxderiv);
    [~,pa_] = ranksum(young_maxderiv,old_maxderiv);
    pa = pa_;
end
title(['Tanh fit gain from max. derivative, WRS test: p =',num2str(pa)]);
ylabel('Slope ( Hz / mV )');
xlabel('Age group');
hold off;
subplot(2,2,3);
scatter((1.5+(0.15-(rand(size(young_xinter))*0.3))).*ones(size(young_xinter)),real(young_xinter),'c','filled');
hold on;
yb = errorbar(1.5,mean(real(young_xinter)),std(real(young_xinter))./sqrt(length(young_xinter)));
yb.Color = [0 0 0];
yb.LineWidth = 3.0;
scatter((3.0+(0.15-(rand(size(old_xinter))*0.3))).*ones(size(old_xinter)),real(old_xinter),'b','filled');
ob = errorbar(3.0,mean(real(old_xinter)),std(real(old_xinter))./sqrt(length(old_xinter)));
ob.Color = [0 0 0];
ob.LineWidth = 3.0;
boxwhisker(real(young_xinter),1.5);
boxwhisker(real(old_xinter),3.0);
set(gca,'xtick',[1.5 3],'XTickLabel',{'Young','Old'});
xlim([0 4.5]);
[h,~] = kstest(real(young_xinter));
[h1,~] = kstest(real(old_xinter));
if h||h1,
    %[~,pkb] = kstest2(real(young_xinter),real(old_xinter));
    [~,pkb] = ranksum(real(young_xinter),real(old_xinter));
    pb = pkb;
else
    %[~,pb_] = kstest2(real(young_xinter),real(old_xinter));
    [~,pb_] = ranksum(real(young_xinter),real(old_xinter));
    pb = pb_;
end
title(['Fit-based Vm-intercept, WRS test: p = ',num2str(pb)]);
ylabel('Relative Vm (mV)');
xlabel('Age group');
hold off;
subplot(2,2,4);
scatter((0.5+(0.1-(rand(size(young_nullTH))*0.2))).*ones(size(young_nullTH)),young_nullTH,'b','filled');
hold on;
ynb = errorbar(0.5,nanmean(young_nullTH),nanstd(young_nullTH)./sqrt(length(young_nullTH)));
ynb.Color = [0 0 0];
ynb.LineWidth = 3.0;
scatter((1.0+(0.1-(rand(size(young_prefTH))*0.2))).*ones(size(young_prefTH)),young_prefTH,'r','filled');
ypb = errorbar(1.0,mean(young_prefTH),std(young_prefTH)./sqrt(length(young_prefTH)));
ypb.Color = [0 0 0];
ypb.LineWidth = 3.0;
scatter((1.5+(0.1-(rand(size(young_cellTH))*0.2))).*ones(size(young_cellTH)),young_cellTH,'k','filled');
yb = errorbar(1.5,mean(young_cellTH),std(young_cellTH)./sqrt(length(young_cellTH)));
yb.Color = [0 0 0];
yb.LineWidth = 3.0;
scatter((2.5+(0.1-(rand(size(old_nullTH))*0.2))).*ones(size(old_nullTH)),old_nullTH,'b','filled');
onb = errorbar(2.5,nanmean(old_nullTH),nanstd(old_nullTH)./sqrt(length(old_nullTH)));
onb.Color = [0 0 0];
onb.LineWidth = 3.0;
scatter((3.0+(0.1-(rand(size(old_prefTH))*0.2))).*ones(size(old_prefTH)),old_prefTH,'r','filled');
opb = errorbar(3.0,mean(old_prefTH),std(old_prefTH)./sqrt(length(old_prefTH)));
opb.Color = [0 0 0];
opb.LineWidth = 3.0;
scatter((3.5+(0.1-(rand(size(old_cellTH))*0.2))).*ones(size(old_cellTH)),old_cellTH,'k','filled');
ob = errorbar(3.5,mean(old_cellTH),std(old_cellTH)./sqrt(length(old_cellTH)));
ob.Color = [0 0 0];
ob.LineWidth = 3.0;
boxwhisker(young_nullTH,0.5);
boxwhisker(young_prefTH,1.0);
boxwhisker(young_cellTH,1.5);
boxwhisker(old_nullTH,2.5);
boxwhisker(old_prefTH,3.0);
boxwhisker(old_cellTH,3.5);
xlim([0 4]);
set(gca,'xtick',[0.5 1.0 1.5 2.5 3.0 3.5],'XTickLabel',{'YN','YP','Young','ON','OP','Old'});
[h2,~] = kstest(young_nullTH); [h3,~] = kstest(young_prefTH); [h4,~] = kstest(old_nullTH);
[h5,~] = kstest(old_prefTH); [h6,~] = kstest(old_prefTH); [h7,~] = kstest(young_cellTH);
[h8,~] = kstest(old_cellTH);
if h7||h8,
    [ah,ap] = kstest2(young_cellTH,old_cellTH);
else
    %[ah,ap] = ttest2(young_cellTH,old_cellTH);
    [ah,ap] = ranksum(young_cellTH,old_cellTH);
end
h_results = [h2;h3;h4;h5;h6];
data = [young_nullTH;young_prefTH;old_nullTH;old_prefTH];
if any(h_results),
    groups = [ones(length(young_nullTH),1);2.*ones(length(young_prefTH),1);...
        3.*ones(length(old_nullTH),1);4.*ones(length(old_prefTH),1)];
    [pg,anovatab,stats] = kruskalwallis(data,groups);
else
    groups = [[ones(length(young_nullTH),1),ones(length(young_nullTH),1)];...
        [ones(length(young_prefTH),1),2.*ones(length(young_prefTH),1)];...
        [2.*ones(length(old_nullTH),1),ones(length(old_nullTH),1)];...
        [2.*ones(length(old_prefTH),1),2.*ones(length(old_prefTH),1)]];
    [pg,anovatab,stats] = anovan(data,groups);
end
%[yh,yp] = ttest2(young_nullTH,young_prefTH);
%[oh,op] = ttest2(old_nullTH,old_prefTH);
%[nh,np] = ttest2(young_nullTH,old_nullTH);
%[ph,pp] = ttest2(young_prefTH,old_prefTH);
%title({['dV/dt method spike thresholds, YN-YP WRS test: p = ',num2str(yp)],...
%    ['ON-OP t-test: p = ',num2str(op)],['Young-old WRS test: p = ',num2str(ap)],...
%    ['YP-OP t-test: p = ',num2str(pp),'; YN-ON WRS test: p = ',num2str(np)]});
if any(h_results),
    figure(f2);
    title(['dV/dt method spike thresholds, Kruskal-Wallis: p = ',num2str(pg)]);
else
    figure(f2);
    title(['dV/dt method spike thresholds, Unbalanced Design 2-way ANOVA: p = ',num2str(pg)]);
end
ylabel('Relative Vm (mV)');
xlabel('Age/Stim. category');
hold off;

save('collected_params_vars.mat');


        



end

