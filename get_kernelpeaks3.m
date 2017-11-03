function [ N_sx,N_ty,N_peaks,M_sx,M_ty,M_peaks ] = ...
    get_kernelpeaks3( kernel,o_kernel,o_kernel_norm,N,M,stim_signal,varargin )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

global scale;
%set up threshold to ID bootstrap ROIs
if nargin < 6,
    alpha = 0.05;
else
    alpha = varargin{1};
end
[bs,bsimg] = rc_bootstrap_set(stim_signal,500,alpha);
%partition into positive and negative ROI subfields
p_kernel = find(bs==max(max(bs)));
n_kernel = find(bs==min(min(bs)));
sub_p = zeros(size(bs));
sub_n = zeros(size(bs));
sub_p(p_kernel) = 1;
sub_n(n_kernel) = 1;
fp_kernel = o_kernel_norm.*sub_p;
temp_fn_kernel = -1*(o_kernel.*sub_n);
fn_kernel = rescale(temp_fn_kernel,scale.*[0 1],[0 abs(min(min(o_kernel_norm)))]);
%fn_kernel = -1*(o_kernel_norm.*sub_n);

%SEARCH one ON subfield
if N == 1,
    N_peaks = max(max(kernel));
    N_ind = find(N_peaks==kernel);
    [N_ty,N_sx] = ind2sub(size(kernel),N_ind);
elseif N > 1,     % or, multiple ON subfields
    nv_ = [];
    nloc_p_ = [];
    p_ = [];
    for p = 1:size(kernel,2),
        [nv,nloc_p] = findpeaks(kernel(:,p));
        nv_ = [nv_;nv];
        nloc_p_ = [nloc_p_;nloc_p];
        p_ = [p_;p.*ones(length(nv),1)];
    end
    [new_nv,nvi] = sort(nv_,1,'descend');
    test_N_sx = p_(nvi);
    test_N_ty = nloc_p_(nvi);
    test_N_peaks = new_nv(:);
    abs_ind_ = sub2ind(size(kernel),test_N_ty,test_N_sx);
    
    %resample kernel for first peak subunit ROI
    %for j = 1:size(sub_p,1),
    %    temp_(j,:) = resample(sub_p(j,:),100,size(sub_p,2),1);
   %end
    %for j = 1:100,
    %    resampled_sub_p(:,j) = resample(temp_(:,j),100,size(sub_p,1),1);
    %end
    %binary_resub_p = (resampled_sub_p ~= 0);
    b = imresize(sub_p,[100 100],'nearest');
    binary_resub_p = b;
    
    sub_kernel_ = cell(N+10,1);
    single_by_size = 0;
    roi_threshold = -1;
    input_binary = 1;
    [nsub_kernel_,flag_peak_] = alt_singlecluster_bypeaks(binary_resub_p,...
        test_N_sx(1,1),test_N_ty(1,1),single_by_size,roi_threshold,input_binary);
    sub_kernel_{1,1} = nsub_kernel_;
    alt_N_sx(1,1) = test_N_sx(1,1);
    alt_N_ty(1,1) = test_N_ty(1,1);
    alt_N_peaks(1,1) = test_N_peaks(1,1);
    
    %now check subsequent putative peaks against ROIS matched to previously
    %accepted peaks
    for n = 2:N+10,
        if ~ismember(abs_ind_(n),find(sub_kernel_{n-1,1})),
            roi_threshold = -1;
            input_binary = 1;
            [psub_kernel_,flag_peak_] = ...
                alt_singlecluster_bypeaks(binary_resub_p,test_N_sx(n,1),...
                test_N_ty(n,1),single_by_size,roi_threshold,input_binary);
            sub_kernel_{n,1} = psub_kernel_;
            alt_N_sx(n,1) = test_N_sx(n,1);
            alt_N_ty(n,1) = test_N_ty(n,1);
            alt_N_peaks(n,1) = test_N_peaks(n,1);
        else
            %if peak is rejected keep the same kernel for next search
            sub_kernel_{n,1} = sub_kernel_{n-1,1}; 
            alt_N_sx(n,1) = NaN;
            alt_N_ty(n,1) = NaN;
            alt_N_peaks(n,1) = NaN;
        end
    end
    alt_N_sx = alt_N_sx(~isnan(alt_N_sx));
    alt_N_ty = alt_N_ty(~isnan(alt_N_ty));
    alt_N_peaks = alt_N_peaks(~isnan(alt_N_peaks));
    N_sx = alt_N_sx(1:N);
    N_ty = alt_N_ty(1:N);
    N_peaks = alt_N_peaks(1:N);
end

%SEARCH one OFF subfield
if M == 1,
    M_peaks = min(min(kernel));
    M_ind = find(M_peaks==kernel);
    [M_ty,M_sx] = ind2sub(size(kernel),M_ind);
elseif M > 1,     % or, multiple OFF subfields
    inv_kernel = (-1.*kernel)+abs(min(min(kernel)))+255;
    mv_ = [];
    mloc_p_ = [];
    s_ = [];
    for s = 1:size(inv_kernel,2),
        [mv,mloc_p] = findpeaks(inv_kernel(:,s));
        mv_ = [mv_;mv];
        mloc_p_ = [mloc_p_;mloc_p];
        s_ = [s_;s.*ones(length(mv),1)];
    end
    [new_mv,mvi] = sort(mv_,1,'descend');
    test_M_sx = s_(mvi);
    test_M_ty = mloc_p_(mvi);
    test_M_peaks = new_mv(:);
    abs_ind = sub2ind(size(inv_kernel),test_M_ty,test_M_sx);
    
    %for i=1:size(sub_n,1),
    %    temp(i,:) = resample(sub_n(i,:),100,size(sub_n,2),1);
    %end
    %for i = 1:100,
    %    resampled_sub_n(:,i) = resample(temp(:,i),100,size(sub_n,1),1);
    %end
    %binary_resub_n = (resampled_sub_n ~= 0);
    b_ = imresize(sub_n,[100 100],'nearest');
    binary_resub_n = b_;
    
    sub_kernel = cell(M+10,1);
    single_by_size = 0;
    roi_threshold = -1;
    input_binary = 1;
    [nsub_kernel,flag_peak] = alt_singlecluster_bypeaks(binary_resub_n,...
        test_M_sx(1,1),test_M_ty(1,1),single_by_size,roi_threshold,input_binary);
    sub_kernel{1,1} = nsub_kernel;
    alt_M_sx(1,1) = test_M_sx(1,1);
    alt_M_ty(1,1) = test_M_ty(1,1);
    alt_M_peaks(1,1) = test_M_peaks(1,1);
    
    for m = 2:M+10,
        if ~ismember(abs_ind(m),find(sub_kernel{m-1,1})),
            roi_threshold = -1;
            input_binary = 1;
            [psub_kernel,flag_peak] = ...
                alt_singlecluster_bypeaks(binary_resub_n,test_M_sx(m,1),...
                test_M_ty(m,1),single_by_size,roi_threshold,input_binary);
            sub_kernel{m,1} = psub_kernel;
            alt_M_sx(m,1) = test_M_sx(m,1);
            alt_M_ty(m,1) = test_M_ty(m,1);
            alt_M_peaks(m,1) = test_M_peaks(m,1);
        else
            sub_kernel{m,1} = sub_kernel{m-1,1};
            alt_M_sx(m,1) = NaN;
            alt_M_ty(m,1) = NaN;
            alt_M_peaks(m,1) = NaN;
        end
    end
    alt_M_sx = alt_M_sx(~isnan(alt_M_sx));
    alt_M_ty = alt_M_ty(~isnan(alt_M_ty));
    alt_M_peaks = alt_M_peaks(~isnan(alt_M_peaks));
    M_sx = alt_M_sx(1:M);
    M_ty = alt_M_ty(1:M);
    M_peaks = alt_M_peaks(1:M);
end
    
    
    
end

