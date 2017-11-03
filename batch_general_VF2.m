function [ p_collated_bin_ymean,p_collated_bin_centers,...
    n_collated_bin_ymean,n_collated_bin_centers ] = batch_general_VF2( plot_cells,bin_x )
%BATCH_GENERAL_VF2 - intermediate function under COMPARE_COLLVF_2.m
%that accesses PREF/NULL VF response data from user selected .mat file and 
%ensures that they are aligned to matching bins. 

%select cell folder
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

p_collated_bin_ymean = [];
p_collated_bin_centers = [];
n_collated_bin_ymean = [];
n_collated_bin_centers = [];
if plot_cells == 1,
    fp = figure;
    fn = figure;
else
end
for i = 1:size(directories,1),
    cd(directories(i,:));
    dd = dir(fullfile(directories(i,:),'*.mat'));
    cfilename = uigetfile('*.mat');
    if cfilename == 0,
        cd ..
        continue
    else
    end
    if bin_x == 1,
        load([sub_folder filesep directories(i,:) filesep cfilename],...
            'coll_FR','coll_Vm','bin_edges',
    load([sub_folder filesep directories(i,:) filesep cfilename],...
        'coll_pref_FR','coll_null_FR','coll_pref_Vm','coll_null_Vm','bin_edges','pcl','ncl');
    %pref
    temp_ymeans = [];
    temp_xbins = [];
    bin_count = [];
    for j = 2:length(bin_edges{1,1}),
        if j < length(bin_edges{1,1}),
            p_bin_inds = find(coll_pref_Vm>=bin_edges{1,1}(1,j-1)&...
                coll_pref_Vm<bin_edges{1,1}(1,j));
        else
            p_bin_inds = find(coll_pref_Vm>=bin_edges{1,1}(1,j-1)&...
                coll_pref_Vm<=bin_edges{1,1}(1,j));
        end
        bin_count(end+1,1) = length(p_bin_inds);
        if ~isempty(p_bin_inds),
            temp_ymeans(end+1,1) = nanmean(coll_pref_FR(p_bin_inds,1));
            temp_xbins = [temp_xbins;j-1];
        else
        end
    end
    gen_bin_centers = (bin_edges{1,1}(1,2:end)+bin_edges{1,1}(1,1:end-1))/2;
    bin_centers = gen_bin_centers(temp_xbins)-pcl.Vth_est;
    p_collated_bin_ymean = [p_collated_bin_ymean;reshape(temp_ymeans,length(temp_ymeans),1)];
    p_collated_bin_centers = [p_collated_bin_centers;reshape(bin_centers,length(bin_centers),1)];

    %null
    ntemp_ymeans = [];
    ntemp_xbins = [];
    nbin_count = [];
    for k = 2:length(bin_edges{1,1}),
        if k < length(bin_edges{1,1}),
            n_bin_inds = find(coll_null_Vm>=bin_edges{1,1}(1,k-1)&...
                coll_null_Vm<bin_edges{1,1}(1,k));
        else
            n_bin_inds = find(coll_null_Vm>=bin_edges{1,1}(1,k-1)&...
                coll_null_Vm<=bin_edges{1,1}(1,k));
        end
        nbin_count(end+1,1) = length(n_bin_inds);
        if ~isempty(n_bin_inds),
            ntemp_ymeans(end+1,1) = nanmean(coll_null_FR(n_bin_inds,1));
            ntemp_xbins = [ntemp_xbins;k-1];
        else
        end
    end
    nbin_centers = gen_bin_centers(ntemp_xbins)-ncl.Vth_est;
    n_collated_bin_ymean = [n_collated_bin_ymean;reshape(ntemp_ymeans,length(ntemp_ymeans),1)];
    n_collated_bin_centers = [n_collated_bin_centers;reshape(nbin_centers,length(nbin_centers),1)];
    cd ..
end
end

