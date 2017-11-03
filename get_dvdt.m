function [ output_args ] = get_dvdt(  )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


[folder,directories] = walk_directories();
cd(folder);

%old_null_maxslope = {};
%young_pref_maxslope = {};
%young_null_maxslope = {};
old_maxslope = {};
young_maxslope = {};
cellmean_old = [];
cellmean_young = [];

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
        try
            load('DVDT_GsmoothM1_F_BL_X_rebinned_trial_responses.mat','pref_ind','null_ind',...
                 'max_dvdt');
        catch
            try
                load('GsmoothM1_F_BL_X_rebinned_trial_responses.mat','pref_ind','null_ind',...
                     'max_dvdt');
            catch
                load('GsmoothM1_F_BL_X_collstim_Xbinfitdata','pref_ind','null_ind',...
                     'max_dvdt');
            end
        end
        all_dirind = 1:length(stimvalues)-1;
        down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
        up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
        down_i = circshift(all_dirind,1,2);
        up_i = circshift(all_dirind,-1,2);
        %pref_set = [down_(pref_ind);pref_stim;up_(pref_ind)];
        %null_set = [down_(null_ind);null_stim;up_(null_ind)];
        pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
        null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
        if strcmp(subfile,'old'),
            old_maxslope{end+1,1} = max_dvdt;
            cellmean_old = [cellmean_old;mean(max_dvdt)];
        else
            young_maxslope{end+1,1} = max_dvdt;
            cellmean_young = [cellmean_young;mean(max_dvdt)];
        end
        clear pref_ind null_ind
        cd ../
    end
    cd ../
end
    
old_dvdt = cell2mat(old_maxslope)./10^3;
young_dvdt = cell2mat(young_maxslope)./10^3;
cellmean_old = cellmean_old./10^3;
cellmean_young = cellmean_young./10^3;

%pooled
%histograms
f1 = figure;
bin_edges = [50:5:400];
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
oh = histc(old_dvdt,bin_edges);
yh = histc(young_dvdt,bin_edges);
oh = oh(1:end-1);
yh = yh(1:end-1);
bar(bin_centers,oh,'m');
hold on;
bar(bin_centers,yh,'g');
alpha(0.5);
xlabel('Max dV / dt (mV / Sec.)');
ylabel('Frequency');
legend('Experienced','Naive');
%cumulative
f2 = figure;
[ox,oy] = cumhist(old_dvdt,[0 400],1.0);
plot(ox,oy,'m-');
hold on;
[yx,yy] = cumhist(young_dvdt,[0 400],0.1);
plot(yx,yy,'g-');
xlabel('Max dV/dt (mV/sec.)');
ylabel('Proportion');
legend('Experienced','Naive','Location','NorthWest');
%stats
[p,h] = ranksum(old_dvdt,young_dvdt);
%bar
f3 = figure;

%on cell mean dvdt
%bars
f4 = figure;
scatter(1.5+(0.25.*(rand(length(cellmean_young),1)-0.5)),cellmean_young,'g');
hold on;
scatter(3.0+(0.25.*(rand(length(cellmean_old),1)-0.5)),cellmean_old,'m');
xlim([0 4.5]);
boxwhisker(cellmean_young,1.5);
boxwhisker(cellmean_old,3.0);
ax = gca;
set(ax,'XTick',[1.5,3],'XTickLabels',{'Naive','Experienced'});
xlabel('Age/Experience');
ylabel('Max dV/dt (mV/sec.)');
legend('Naive','Experienced','Location','NorthWest');
%stats
[pp,hh] = ranksum(cellmean_young,cellmean_old);




