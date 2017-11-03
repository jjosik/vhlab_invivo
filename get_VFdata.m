function [ collated_bin_ymean,collated_bin_centers,...
    case_simpleslope,case_maxderiv,case_xinter,case_cellTH,case_prefTH,...
    case_nullTH,params_array,cell_names] = get_VFdata( plot_cells,bin_cellX )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

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

collated_bin_ymean = [];
collated_bin_centers = [];
if plot_cells == 1,
    ff = figure;
    %...under construction
else
end
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
for i = 1:size(directories,1),
    cd(directories(i,:));
    dd = dir(fullfile(directories(i,:),'*.mat'));
    cfilename = uigetfile('*.mat');
    if cfilename == 0,
        cd ..
        continue
    else
    end
    if bin_cellX == 1,
        load([sub_folder filesep directories(i,:) filesep cfilename],...
            'coll_FR','coll_Vm','bin_edges','cellmean_corr_th',...
            'prefmean_corr_th','nullmean_corr_th','params',...
            'gof_obj','allfits_info','ms_x','mt_out');
    else
        load([sub_folder filesep directories(i,:) filesep cfilename],...
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
    
    cd ..
    
end
end



