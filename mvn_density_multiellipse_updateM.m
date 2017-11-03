function [ ON_ellipse_array,OFF_ellipse_array,...
    singleunit_params,interunit_params] = ...
    mvn_density_multiellipse_updateM( resampled_kernel,o_kernel,pixels_square,watch_fit,...
    save_it,theta_per_bar,sampleRate,show_contours,roi_method,output_vm,varargin )

%******'_updateM' suffix indicates latest March 2017 update


%******INPUTS - 
%******     RESAMPLED_KERNEL - mean kernel at final resolution to be
%               analyzed; in general input STRFs are typically resampled to 50x50 or
%               100x100 to smooth out mismatched sampling between temporal
%               and spatial dimensions.
%           O_KERNEL - "original" mean kernel at original resolution 
%               (typically saved under analysis files as 'Kernel_mean')
%           PIXELS_SQUARE - simple scalar value of square pixel
%               resolution used in the corresponding resampled kernel (i.e.
%               if kernel is 50x50, pixels_square should be set to 50).
%               NOTE: this is initially set in the batch analyze file and
%               is used by that part of the code to correctly ID the
%               resampled kernel .mat file.  So don't change it here.  
%           WATCH_FIT - should be set to '1' to display window showing
%               initial fit conditions and fit parameters/results in progress.
%           SAVE_IT - set to '1' if workspace and all plots should be saved
%               to disk.
%           THETA_PER_BAR - gives the degrees of visual angle to columns 
%               conversion factor for the spatial dimension.  This is
%               passed to the function from the .mat files opened in batch
%               analysis.
%           SAMPLERATE - simply sampling rate, the rate at which discrete observations are
%               obtained during the recording.  Obtained from the 
%               'analyzed_xcorr_Vm.mat' file.  Typically passed from the
%               batch analysis function.
%           SHOW_CONTOURS - shows output STRF images with fits superimposed
%               on isocontours of the responses defined in z.
%           ROI_METHOD - determines which of two region-of-interest (ROI) ID
%               methods should be utilized.  Set to '1', 
%******THRESHOLDING METHODS -- 
%         ON_OFF_PIXELSORT function (to use set variable pixelsort_ON = 1)
%           Thresholds all data in image to include upper 80% for each
%           ON/OFF valence BEFORE ROIs are ID'd.  
%         ROI_THRESHOLD (set to a value between 0 and 1) sets threshold
%         level of data within each ROI
%******ROI ID METHODS -- 
%         ROI_METHOD = 1, uses k-means clustering to repartition kernel by
%         z-values.
%         ROI_METHOD = 2, uses stimulus bootstrap to ID "significant
%         pixels" first.  Creates ROI masks for both response valences.  
%       IN GENERAL, primary ROIs are ID'd through overlap with global
%       peaks by valence.  These primary ROIs are then used to exclude
%       spurious local extrema within the boundaries of the region.
%       Secondary ROIs are only retained for analysis IFF their total area
%       > 50% of the area of the primary ROI.  The number of subunits to be
%       fit are also constrained by user input at the beginning of the
%       program.  

sampleInterval = 1/sampleRate;

figure(10);
imagesc(resampled_kernel);
N = input(['Enter number of ON subunits: ']);
M = input(['Enter number of OFF subunits: ']);
% 
global scalemax scalemin scale;
scalemax = max(resampled_kernel(:));
scalemin = abs(min(resampled_kernel(:)));
%scale = max(scalemax,scalemin);
scale = scalemax+scalemin;

xgrid = 1:size(resampled_kernel,2);
tgrid = 1:size(resampled_kernel,1);
[X,Y] = meshgrid(xgrid,tgrid);

%if nargin < 11,
%    new_scalemax = 255; %default value in the absence of user input
%else
%    new_scalemax = varargin{1};
%end
kernel_norm = rescale(resampled_kernel,[-scalemin scalemax],[0 scale]);
o_kernel_norm = rescale(o_kernel,[-scalemin scalemax],[0 scale]);
half_scale = scalemin;
%if min(resampled_kernel(:)) == -(scalemin),
%    adj_range = [-(scalemin/scale) 1];
%else
%    adj_range = [-1 (scalemax/scale)];
%end
%kernel_norm = rescale(resampled_kernel,scale.*adj_range,[0 new_scalemax]);
%o_kernel_norm = rescale(o_kernel,scale.*adj_range,[0 new_scalemax]);
%***alt. method: following 3 lines set pivot point of ON/OFF subunit inversion 
%***to value equivalent to that in the original scaling; default sets 
%***"zero" to exactly half of the new scale range.***
%offset_prop = (scalemax+(-scalemin))/(scalemax+scalemin);
%zero_adjust = round(offset_prop*new_scalemax);
%half_scale = new_scalemax/2 + zero_adjust;


%get putative peaks as estimates for means (inspired by Leong et al, J.
%Neurosci, 2016)
%try
%    [N_sx,N_ty,N_peaks,M_sx,M_ty,M_peaks] = get_kernelpeaks(kernel_norm,N,M);
%catch
%    [N_sx,N_ty,N_peaks,M_sx,M_ty,M_peaks] = get_kernelpeaks2(kernel_norm,N,M,X,Y);
%end
alpha_pid = 0.05;
stim_signal = output_vm.xc_stimsignal;
[N_sx,N_ty,N_peaks,M_sx,M_ty,M_peaks] = get_kernelpeaks3(kernel_norm,...
    o_kernel,o_kernel_norm,N,M,stim_signal,alpha_pid);   %this method uses bootstrap ROIs to
                                            %exclude spurious local extrema

%default value:
pixelsort_ON = 0;
if pixelsort_ON == 1,
%ID upper 80 percentile of respective data to separate ON and OFF valences
    [p_bin_sub,n_bin_sub] = ON_OFF_pixelsort(resampled_kernel);
else
    p_bin_sub = ones(size(kernel_norm));
    n_bin_sub = ones(size(kernel_norm));
end
                                            
%default: roi_method = 1                                           
switch  roi_method
    
    case 1

        %use K-means clustering to repartition kernel by intensities; k=3 will 
        %implicitly partition raw data matched to either (1)ON, (2)OFF subfields,
        % or (3)noisy values near baseline (***NOTE: this result is NOT clustering by 
        %spatial location so data from different subfields with the same valence
        %will be classified together at this stage)
        options = statset('Display','final');
        nk = kmeans(kernel_norm(:),3,'Distance','sqeuclidean','Replicates',5,...
            'Options',options);
        zclustered_k = reshape(nk,size(kernel_norm));

        p_s_kernel = kernel_norm;
        n_s_kernel = kernel_norm;
        background_s_kernel = kernel_norm;
        querykernel_a = max(kernel_norm(find(zclustered_k==1)));
        querykernel_b = max(kernel_norm(find(zclustered_k==2)));
        querykernel_c = max(kernel_norm(find(zclustered_k==3)));
        qk = [querykernel_a;querykernel_b;querykernel_c];
        [~,qi] = sort(qk,'descend');
        base_perm = [1;2;3];
        n_s_kernel(find(zclustered_k~=base_perm(qi(3)))) = 0;
        background_s_kernel(find(zclustered_k~=base_perm(qi(2)))) = 0;
        p_s_kernel(find(zclustered_k~=base_perm(qi(1)))) = 0;
        p_s_kernel = p_s_kernel.*p_bin_sub;
        n_s_kernel = n_s_kernel.*n_bin_sub;
        
        push_method = 1;
        
    case 2
        
        alpha_roi = 0.05;
        [bs,bsimg] = rc_bootstrap_set(output_vm.xc_stimsignal,500,alpha_roi);
        p_kernel = find(bs==max(max(bs)));
        n_kernel = find(bs==min(min(bs)));
        sub_p = zeros(size(bs));
        sub_n = zeros(size(bs));
        sub_p(p_kernel) = 1;
        sub_n(n_kernel) = 1;
        binary_resub_p = imresize(sub_p,[pixels_square pixels_square],'nearest');
        binary_resub_n = imresize(sub_n,[pixels_square pixels_square],'nearest');
        fp_kernel = kernel_norm.*binary_resub_p.*p_bin_sub;
        fn_kernel = -1.*(resampled_kernel.*binary_resub_n.*n_bin_sub);
        %fn_kernel = rescale(fn_kernel,scale.*[0 1],[0 new_scalemax]);
        p_trace = heaviside(del2(binary_resub_p));
        n_trace = heaviside(del2(binary_resub_n));
        
        
        push_method = 2;
        
end


%fit positive subfields; first recluster spatially in x-y if N > 1

if push_method ~= 2,
    ms_kernel = p_s_kernel+sign(background_s_kernel).*half_scale+...
        (sign(n_s_kernel).*half_scale);
else
end
single_by_size = 0;
roi_threshold = 0.95;
input_binary = 0;
coll_sel_kernel = cell(N,1);
if N > 1,
    for ii = 1:N,
        if push_method ~= 2,
            [sel_kernel,~] = alt_singlecluster_bypeaks(ms_kernel,N_sx(ii,1),...
                N_ty(ii,1),single_by_size,roi_threshold);
            coll_sel_kernel{ii,1} = sel_kernel;
        else
            [sel_kernel,~] = alt_singlecluster_bypeaks(fp_kernel,N_sx(ii,1),...
                N_ty(ii,1),single_by_size,roi_threshold,input_binary);
            coll_sel_kernel{ii,1} = sel_kernel;
        end
        %refit = 0;
        [fit,p_array,opt_array,sse_array,x_track,y_track,mu,sigma,ellipse_array] = ...
            mvn_density_singleellipseALT(kernel_norm,sel_kernel,...
            X,Y,N_sx(ii,1),N_ty(ii,1),N_peaks(ii,1),watch_fit);
        [xi,yi] = find(~(opt_array-min(min(min(opt_array)))));
        [xi,yi] = find(~(sse_array-min(min(min(sse_array)))));
        ON_mu_array{ii,1} = mu;
        ON_sigma_array{ii,1} = sigma;
        ON_ellipse_array{ii,1} = ellipse_array{xi,yi};
    end
else
    if push_method ~= 2,
        [sel_kernel,~] = alt_singlecluster_bypeaks(ms_kernel,N_sx,...
                N_ty,single_by_size,roi_threshold,input_binary);
        coll_sel_kernel = sel_kernel;
    else
        [sel_kernel,~] = alt_singlecluster_bypeaks(fp_kernel,N_sx,...
                N_ty,single_by_size,roi_threshold,input_binary);
        coll_sel_kernel = sel_kernel;
    end
    %refit = 0;
    [fit,p_array,opt_array,sse_array,x_track,y_track,mu,sigma,ellipse_array] = ...
        mvn_density_singleellipseALT(kernel_norm,sel_kernel,...
        X,Y,N_sx,N_ty,N_peaks,watch_fit);
    [xi,yi] = find(~(opt_array-min(min(min(opt_array)))));
    [xi,yi] = find(~(sse_array-min(min(min(sse_array)))));
    ON_ellipse_array = ellipse_array{xi,yi};
end

%fit negative subfields; first reclustering spatially in x-y if M > 1
if push_method ~= 2,
    n_ms_kernel = half_scale+(((half_scale-n_s_kernel)./half_scale).*...
        (half_scale.*sign(n_s_kernel)));
    n_bg_kernel = (1+((half_scale-background_s_kernel)./half_scale)).*...
        (half_scale.*sign(background_s_kernel));
else
end
single_by_size = 0;
roi_threshold = 0.95;
input_binary = 0;
coll_nsel_kernel = cell(M,1);
if M > 1,
    for j=1:M,
        if push_method ~= 2,
            [nsel_kernel,~] = alt_singlecluster_bypeaks(n_ms_kernel,M_sx(j,1),...
                M_ty(j,1),single_by_size,roi_threshold,input_binary);
            coll_nsel_kernel{j,1} = nsel_kernel;
        else
            [nsel_kernel,~] = alt_singlecluster_bypeaks(fn_kernel,M_sx(j,1),...
                M_ty(j,1),single_by_size,roi_threshold,input_binary);
            coll_nsel_kernel{j,1} = nsel_kernel;
        end
        %refit = 0;
        [nfit,nparray,nopt_array,nsse_array,mx_track,my_track,m_mu,m_sigma,...
            n_ellipse_array] = mvn_density_singleellipseALT(kernel_norm,...
            nsel_kernel,X,Y,M_sx(j,1),M_ty(j,1),...
            M_peaks(j,1),watch_fit);
        [xj,yj] = find(~(nopt_array-min(min(min(nopt_array)))));
        [xj,yj] = find(~(nsse_array-min(min(min(nsse_array)))));
        OFF_mu_array{j,1} = m_mu;
        OFF_sigma_array{j,1} = m_sigma;
        OFF_ellipse_array{j,1} = n_ellipse_array{xj,yj};
        clear nsel_kernel;
    end
else
    if push_method ~= 2,
        [nsel_kernel,~] = alt_singlecluster_bypeaks(n_ms_kernel,M_sx,...
                M_ty,single_by_size,roi_threshold,input_binary);
        coll_nsel_kernel = nsel_kernel;
    else
        [nsel_kernel,~] = alt_singlecluster_bypeaks(fn_kernel,M_sx,...
                M_ty,single_by_size,roi_threshold,input_binary);
        coll_nsel_kernel = nsel_kernel;
    end
    %refit = 0;
    [nfit,np_array,nopt_array,nsse_array,mx_track,my_track,m_mu,m_sigma,...
        n_ellipse_array] = ...
        mvn_density_singleellipseALT(kernel_norm,nsel_kernel,X,Y,M_sx,M_ty,...
        M_peaks,watch_fit);
    [xj,yj] = find(~(nopt_array-min(min(min(nopt_array)))));
    [xj,yj] = find(~(nsse_array-min(min(min(nsse_array)))));
    OFF_ellipse_array = n_ellipse_array{xj,yj};
end

%generate elliptical error overlay
figure(2);
%imagesc(kernel_norm);
imagesc(resampled_kernel);
x_scale = 1;        %x scale corrected in resample file
y_scale = (size(o_kernel,1)*sampleInterval)/100; %y scale was not
colormap(redblue);
set(gca,'Ydir','normal');
hold on;
if N > 1,
    for p = 1:N,
        h = plot(ON_ellipse_array{p,1}.plot_ellipse(1,:)./x_scale,...
            flipud(ON_ellipse_array{p,1}.plot_ellipse(2,:)),'k');
        h.LineWidth = 2.0;
    end
else
    h = plot(ON_ellipse_array.plot_ellipse(1,:)./x_scale,...
        flipud(ON_ellipse_array.plot_ellipse(2,:)),'k');
    h.LineWidth = 2.0;
end
xlabel('Space (deg)');
ylabel('Time (sec)');
xlim([1 size(resampled_kernel,2)/x_scale]);
ylim([1 size(resampled_kernel,1)]);
hold on;
if M > 1,
    for np = 1:M,
        hh = plot(OFF_ellipse_array{np,1}.plot_ellipse(1,:)./x_scale,...
            flipud(OFF_ellipse_array{np,1}.plot_ellipse(2,:)),'k--');
        hh.LineWidth = 2.0;
    end
else
    hh = plot(OFF_ellipse_array.plot_ellipse(1,:)./x_scale,...
       flipud(OFF_ellipse_array.plot_ellipse(2,:)),'k--');
    hh.LineWidth = 2.0;
end
colorbar;
a = gca;
newx = round((a.XTick.*theta_per_bar*x_scale),1);
newy = round((a.YTick.*y_scale),3);
set(gca,'XTickLabel',newx);
set(gca,'YTickLabel',newy);
%legend('ON subunits','OFF subunits','Location','EastOutside');
%default:
show_contours = 0;
if show_contours == 1,
    contour(resampled_kernel);
else
end
if save_it == 1,
    cd figures
    if show_contours == 1,
        saveas(gcf,[pwd filesep 'TH80_fitted_resample50STRF_wC.fig']);
    else
        saveas(gcf,[pwd filesep 'TH80_fitted_resample50STRF.fig']);
    end
    close(gcf);
    cd ..
else
end

%generate STRF showing ROIs and selection criteria   
figure(3);
%kernel_renorm = rescale(kernel_norm,...
%    [min(min(kernel_norm)) max(max(kernel_norm))],[-1 1]);
resamp_renorm = rescale(resampled_kernel,...
    [min(min(resampled_kernel)) max(max(resampled_kernel))],[-1 1]);
trace_knorm = resamp_renorm;
loc_ptrace = find(p_trace);
loc_ntrace = find(n_trace);
trace_knorm(loc_ntrace) = min(min(resamp_renorm));
trace_knorm(loc_ptrace) = max(max(resamp_renorm));
% find final ROI regions
if iscell(coll_sel_kernel),
    for k = 1:length(coll_sel_kernel),
        subp_trace = heaviside(del2(coll_sel_kernel{k,1}));
        subloc_ptrace = find(subp_trace);
        %trace_sub(subloc_ptrace) = max(max(resamp_renorm));
    end
else
    subp_trace = heaviside(del2(coll_sel_kernel));
    subloc_ptrace = find(subp_trace);
end
if iscell(coll_nsel_kernel),
    for kk = 1:length(coll_nsel_kernel),
        subn_trace = heaviside(del2(coll_nsel_kernel{kk,1}));
        subloc_ntrace = find(subn_trace);
        %trace_sub(subloc_ntrace) = min(min(resamp_renorm));
    end
else
    subn_trace = heaviside(del2(coll_nsel_kernel));
    subloc_ntrace = find(subn_trace);
end
sub_trace = resamp_renorm;
sub_trace(subloc_ntrace) = min(min(resamp_renorm));
sub_trace(subloc_ptrace) = max(max(resamp_renorm));
imagesc(trace_knorm);
x_scale = 1;        %x scale corrected in resample file
y_scale = (size(o_kernel,1)*sampleInterval)/100; %y scale was not
colormap(redblue);
set(gca,'Ydir','normal');
hold on;
sz = 50;
scatter(N_sx,N_ty,sz,'y','filled');
scatter(M_sx,M_ty,sz,'k','filled');
xlabel('Space (deg)');
ylabel('Time (sec)');
xlim([1 size(resampled_kernel,2)/x_scale]);
ylim([1 size(resampled_kernel,1)]);
cb = colorbar;
cb.Label.String = 'Correlation (Vm * contrast)';
a = gca;
newx = round((a.XTick.*theta_per_bar*x_scale),1);
newy = round((a.YTick.*y_scale),3);
set(gca,'XTickLabel',newx);
set(gca,'YTickLabel',newy);
title('Bootstrap ROI selection with detected peaks');
if save_it == 1,
    cd figures
    saveas(gcf,[pwd filesep 'ROI_SelectionCriteria50.fig']);
    close(gcf);
    cd ..
else
end
%Generate Final selection ROI overlay
figure(4);  
imagesc(sub_trace);
x_scale = 1;        %x scale corrected in resample file
y_scale = (size(o_kernel,1)*sampleInterval)/100; %y scale was not
colormap(redblue);
set(gca,'Ydir','normal');
hold on;
sz = 50;
scatter(N_sx,N_ty,sz,'y','filled');
scatter(M_sx,M_ty,sz,'k','filled');
xlabel('Space (deg)');
ylabel('Time (sec)');
xlim([1 size(resampled_kernel,2)/x_scale]);
ylim([1 size(resampled_kernel,1)]);
cb = colorbar;
cb.Label.String = 'Correlation (Vm * contrast)';
a = gca;
newx = round((a.XTick.*theta_per_bar*x_scale),1);
newy = round((a.YTick.*y_scale),3);
set(gca,'XTickLabel',newx);
set(gca,'YTickLabel',newy);
title('Final thresholded ROI selection with detected peaks');
if save_it == 1,
    cd figures
    saveas(gcf,[pwd filesep 'ROI_FinalSelectionCrit50.fig']);
    close(gcf);
    cd ..
else
end



%****************************
%Calculate subunit parameters
%****************************
%***IMPORTANT***

singleunit_params = struct('subunit_temporal_onset',[],'temporal_onset',[],'time_to_peak',[],...
    'spatial_extent',[],'temporal_extent',[],'peak_amplitude',[],...
    'ST_velocity',[],'ST_tilt',[]);
interunit_params = struct('spatial_offset',[],'on2off_velocity',[]);
%first translate coordinates into time and visual space
cN_sx = N_sx.*theta_per_bar;
cN_ty = N_ty.*y_scale;
cM_sx = M_sx.*theta_per_bar;
cM_ty = M_ty.*y_scale;
%1. Spatial offset
for i = 1:N,
    for j = 1:M,
        spatial_offset{i,j} = abs(cN_sx(i,1)-cM_sx(j,1));
    end
end
interunit_params.spatial_offset = spatial_offset;
%2. Temporal onset
subunit_temporal_onset = cell(N+M,1);
singlearray_onsets = [];
for n = 1:N+M,
    if n<=N&&N==1,
        subunit_temporal_onset{n,1} = ...
            (min(min(ON_ellipse_array.plot_ellipse(2,:))).*y_scale)-...
            0.1;
        singlearray_onsets = [singlearray_onsets;subunit_temporal_onset{n,1}];
    else
    end
    if n<=N&&N>1,
        subunit_temporal_onset{n,1} = ...
            (min(min(ON_ellipse_array{n,1}.plot_ellipse(2,:))).*...
            y_scale)-0.1;
        singlearray_onsets = [singlearray_onsets;subunit_temporal_onset{n,1}];
    else
    end
    if n==N+1,
        try
            subunit_temporal_onset{n,1} = ...
                (min(min(OFF_ellipse_array.plot_ellipse(2,:))).*y_scale)-0.1;
        catch
             subunit_temporal_onset{n,1} = ...
                (min(min(OFF_ellipse_array{1,1}.plot_ellipse(2,:))).*y_scale)-0.1;
        end
        singlearray_onsets = [singlearray_onsets;subunit_temporal_onset{n,1}];
    else
    end
    if n>N+1,
        subunit_temporal_onset{n,1} = ...
            (min(min(OFF_ellipse_array{n-N,1}.plot_ellipse(2,:))).*y_scale)-...
            0.1;
        singlearray_onsets = [singlearray_onsets;subunit_temporal_onset{n,1}];
    else
    end
end
singleunit_params.temporal_onset = min(singlearray_onsets);
singleunit_params.subunit_temporal_onset = subunit_temporal_onset;
%3. Time to peak
time_to_peak = cell(N+M,1);
for n = 1:N+M,
    if n<=N,
        time_to_peak{n,1} = cN_ty(n)-0.1;
    else
        time_to_peak{n,1} = cM_ty(n-N)-0.1;
    end
end
singleunit_params.time_to_peak = time_to_peak;
%4. Spatial extent
spatial_extent = cell(N+M,1);
for n = 1:N+M,
    if n<=N&&N==1,
        spatial_extent{n,1} = ...
            (abs(max(max(ON_ellipse_array.plot_ellipse(1,:)))-...
            min(min(ON_ellipse_array.plot_ellipse(1,:))))).*theta_per_bar;
    else
    end
    if n<=N&&N>1,
        spatial_extent{n,1} = ...
            (abs(max(max(ON_ellipse_array{n,1}.plot_ellipse(1,:)))-...
            min(min(ON_ellipse_array{n,1}.plot_ellipse(1,:))))).*theta_per_bar;
    else
    end
    if n==N+1,
        try
            spatial_extent{n,1} = ...
                abs(max(max(OFF_ellipse_array.plot_ellipse(1,:)))-...
                min(min(OFF_ellipse_array.plot_ellipse(1,:)))).*theta_per_bar;
        catch
            spatial_extent{n,1} = ...
                (abs(max(max(OFF_ellipse_array{1,1}.plot_ellipse(1,:)))-...
                min(min(OFF_ellipse_array{1,1}.plot_ellipse(1,:))))).*theta_per_bar;
        end
    else
    end
    if n>N+1,
        spatial_extent{n,1} = ...
            (abs(max(max(OFF_ellipse_array{n-N,1}.plot_ellipse(1,:)))-...
            min(min(OFF_ellipse_array{n-N,1}.plot_ellipse(1,:))))).*theta_per_bar;
    else
    end
end
singleunit_params.spatial_extent = spatial_extent;
%5. Temporal extent
temporal_extent = cell(N+M,1);
for n = 1:N+M,
    if n<=N&&N==1,
        temporal_extent{n,1} = ...
            abs(max(max(ON_ellipse_array.plot_ellipse(2,:))).*y_scale-...
            min(min(ON_ellipse_array.plot_ellipse(2,:))).*y_scale);
    else
    end
    if n<=N&&N>1
        temporal_extent{n,1} = ...
            abs(max(max(ON_ellipse_array{n,1}.plot_ellipse(2,:))).*y_scale-...
            min(min(ON_ellipse_array{n,1}.plot_ellipse(2,:))).*y_scale);
    else
    end
    if n==N+1,
        try
            temporal_extent{n,1} = ...
                abs(max(max(OFF_ellipse_array.plot_ellipse(2,:))).*y_scale-...
                min(min(OFF_ellipse_array.plot_ellipse(2,:))).*y_scale);
        catch
            temporal_extent{n,1} = ...
                abs(max(max(OFF_ellipse_array{1,1}.plot_ellipse(2,:))).*y_scale-...
                min(min(OFF_ellipse_array{1,1}.plot_ellipse(2,:))).*y_scale);
        end
    else
    end
    if n>N+1,
        temporal_extent{n,1} = ...
            abs(max(max(OFF_ellipse_array{n-N,1}.plot_ellipse(2,:))).*y_scale-...
            min(min(OFF_ellipse_array{n-N,1}.plot_ellipse(2,:))).*y_scale);
    else
    end
end
singleunit_params.temporal_extent = temporal_extent;
%6. peak amplitude
peak_amplitude = cell(N+M,1);
for n = 1:N+M,
    if n<=N,
        peak_amplitude{n,1} = N_peaks(n);
    else
        peak_amplitude{n,1} = M_peaks(n-N);
    end
end
singleunit_params.peak_amplitude = peak_amplitude;
%7. S-T tilts and S_T velocity
ST_velocity = cell(N+M,1);
for n = 1:N+M,
    ST_velocity{n,1} = ...
        (singleunit_params.spatial_extent{n,1}/theta_per_bar)/(singleunit_params.temporal_extent{n,1}/y_scale);
end
singleunit_params.ST_velocity = ST_velocity;
ST_tilt = cell(N+M,1);
for n = 1:N+M,
    if n == 1,
        try
            ST_tilt{n,1} = 90 - acosd(((singleunit_params.temporal_extent{n,1}/y_scale)/2)/ON_ellipse_array.major);
        catch
            ST_tilt{n,1} = 90 - acosd(((singleunit_params.temporal_extent{n,1}/y_scale)/2)/ON_ellipse_array{n,1}.major);
        end
    else 
    end
    if n<=N&&N>1,
        ST_tilt{n,1} = 90 - acosd(((singleunit_params.temporal_extent{n,1}/y_scale)/2)/ON_ellipse_array{n,1}.major);
    else
    end
    if n==N+1,
        try
            ST_tilt{n,1} = 90 - acosd(((singleunit_params.temporal_extent{n,1}/y_scale)/2)/OFF_ellipse_array.major);
        catch
            ST_tilt{n,1} = 90 - acosd(((singleunit_params.temporal_extent{n,1}/y_scale)/2)/OFF_ellipse_array{n-N,1}.major);
        end
    else
    end
    if n>N+1,
        ST_tilt{n,1} = 90 - acosd(((singleunit_params.temporal_extent{n,1}/y_scale)/2)/OFF_ellipse_array{n-N,1}.major);
    else
    end
end
singleunit_params.ST_tilt = ST_tilt;
%8. ON-OFF transition velocity (peak-to-peak)
on2off_velocity = cell(N,M);
for n = 1:N,
    for m = 1:M,
        on2off_velocity{n,m} = ...
            abs(N_sx(n,1)-M_sx(m,1))/abs(N_ty(n,1)-M_ty(m,1));
    end
end
interunit_params.on2off_velocity = on2off_velocity;

if save_it == 1,
    save('STRF50_structure_TH80CI90.mat');
else
end


end

