function [ th_selected_kernel,flag_peak ] = ...
    alt_singlecluster_bypeaks( unclustered_kernel,...
    mean_sx,mean_ty,single_by_size,roi_threshold,input_binary)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

interim = unclustered_kernel-min(min(unclustered_kernel));
[bb,L] = bwboundaries(interim,4);
si = cellfun(@size, bb, 'UniformOutput',false);
ssi = cellfun(@max, si, 'UniformOutput',false);
if single_by_size == 1,
    [~,ci] = max(cell2mat(ssi));
else
    re_ind = sub2ind(size(unclustered_kernel),mean_ty,mean_sx);
    ci = L(re_ind);
end

if ci == 0,
    flag_peak = 1;
    selected_kernel = [];
else
    flag_peak = 0;
    rois = regionprops(L,'PixelIdxList');
    a = rois(ci).PixelIdxList;
    sub_k = zeros(size(unclustered_kernel));
    sub_k(a(1:end)) = 1;
    new_msk = unclustered_kernel.*sub_k;
    selected_kernel = new_msk;
end

if roi_threshold > 0 && roi_threshold < 1,
    if ~isempty(selected_kernel),
        [tempset,~] = sort(selected_kernel,'descend');
        tempset(find(tempset==0)) = [];
        sub_size = ceil(roi_threshold*length(tempset));
        sub_set = find(abs(selected_kernel)>tempset(sub_size));
        temp_selected_kernel = zeros(size(selected_kernel));
        temp_selected_kernel(sub_set) = 1;
        th_selected_kernel = temp_selected_kernel.*selected_kernel;
        baseline_correct = min(abs(th_selected_kernel(find(th_selected_kernel))));
        th_selected_kernel = ...
            heaviside(th_selected_kernel).*(th_selected_kernel-baseline_correct);
    else
        th_selected_kernel = unclustered_kernel;
    end
else
    if input_binary == 1,
        th_selected_kernel = selected_kernel;
    else
        baseline_correct = min(abs(selected_kernel(find(selected_kernel))));
        th_selected_kernel = ...
        heaviside(selected_kernel).*(selected_kernel-baseline_correct);
    end
end
    

        
        


end

