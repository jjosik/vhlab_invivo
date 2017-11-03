function [  ] = strfshape_batchanalyze_RESAMPLE( pixels_square,varargin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%   Use: run_from_file = 1 to open saved mat-file analyses for comparison
%   -or- run_from_file = 0 to run all new fits
%   DEFAULT value is 0.

save_it = 1;
watch_fit = 1;
overwrite = 0;
show_contours = 0;
roi_method = 2;
if nargin < 2,
    run_from_file = 0;
else
    run_from_file = varargin{1};
end

%top = ('/Users/vhlab/Documents/intracellular_data/');
%p = addpath(genpath(top));

%select cells to analyze
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

switch run_from_file
    
    case 0
        
        [~,subfile,~] = fileparts(sub_folder);
        if strcmp('young',subfile),
            y_coll_subtemponset = [];
            y_coll_temponset = [];
            y_coll_time2peak = [];
            y_coll_spatialext = [];
            y_coll_tempext = [];
            y_coll_peakamp = [];
            y_coll_STtilt = [];
            y_coll_STvelocity = [];
            y_coll_spoffset = [];
            y_coll_onoffvel = [];
            y_cell_id = [];
            Iy_cell_id = [];
        else
        end
        if strcmp('old',subfile),
            o_coll_subtemponset = [];
            o_coll_temponset = [];
            o_coll_time2peak = [];
            o_coll_spatialext = [];
            o_coll_tempext = [];
            o_coll_peakamp = [];
            o_coll_STtilt = [];
            o_coll_STvelocity = [];
            o_coll_spoffset = [];
            o_coll_onoffvel = [];
            o_cell_id = [];
            Io_cell_id = [];
        else
        end

        for i = 1:size(directories,1),
            cd(directories(i,:));
            refile = ['resampled_kernel_',num2str(pixels_square) '.mat'];
            load([sub_folder filesep directories(i,:) filesep 'stims.mat']);
            load([sub_folder filesep directories(i,:) filesep 'data.mat']);
            load([sub_folder filesep directories(i,:) filesep 'stimorder.mat']);
            load([sub_folder filesep directories(i,:) filesep 'stimvalues.mat']);
            load([sub_folder filesep directories(i,:) filesep 'analyzed_xcorr_Vm.mat']);
            load([sub_folder filesep directories(i,:) filesep refile]);
            [ON_ellipse_array,OFF_ellipse_array,singleunit_params,interunit_params] = ...
                mvn_density_multiellipse_resampled(resampled_kernel_50,...
                Kernel_mean,pixels_square,watch_fit,save_it,resampled_theta_per_bar,...
                sampleRate,show_contours,roi_method,output_vm);
    %collate subunit params by age groups
            if strcmp('young',subfile),
                y_coll_subtemponset = [y_coll_subtemponset;cell2mat(singleunit_params.subunit_temporal_onset)];
                y_coll_temponset = [y_coll_temponset;singleunit_params.temporal_onset];
                y_coll_time2peak = [y_coll_time2peak;cell2mat(singleunit_params.time_to_peak)];
                y_coll_spatialext = [y_coll_spatialext;cell2mat(singleunit_params.spatial_extent)];
                y_coll_tempext = [y_coll_tempext;cell2mat(singleunit_params.temporal_extent)];
                y_coll_peakamp = [y_coll_peakamp;cell2mat(singleunit_params.peak_amplitude)];
                y_coll_STtilt = [y_coll_STtilt;cell2mat(singleunit_params.ST_tilt)];
                y_coll_STvelocity = [y_coll_STtilt;cell2mat(singleunit_params.ST_velocity)];
                y_cell_id = [y_cell_id;directories(i,:)];
            else
            end
            if strcmp('old',subfile),
             
                o_coll_subtemponset = [o_coll_subtemponset;cell2mat(singleunit_params.subunit_temporal_onset)];
                o_coll_temponset = [o_coll_temponset;singleunit_params.temporal_onset];
                o_coll_time2peak = [o_coll_time2peak;cell2mat(singleunit_params.time_to_peak)];
                o_coll_spatialext = [o_coll_spatialext;cell2mat(singleunit_params.spatial_extent)];
                o_coll_tempext = [o_coll_tempext;cell2mat(singleunit_params.temporal_extent)];
                o_coll_peakamp = [o_coll_peakamp;cell2mat(singleunit_params.peak_amplitude)];
                o_coll_STtilt = [o_coll_STtilt;cell2mat(singleunit_params.ST_tilt)];
                o_coll_STvelocity = [o_coll_STvelocity;cell2mat(singleunit_params.ST_velocity)];
                o_cell_id = [o_cell_id;directories(i,:)];
            else
            end
            if strcmp('young',subfile),
                y_coll_spoffset = [y_coll_spoffset;...
                    reshape(cell2mat(interunit_params.spatial_offset),...
                    size(interunit_params.spatial_offset,1)*size(interunit_params.spatial_offset,2),1)];
                y_coll_onoffvel = [y_coll_onoffvel;...
                    reshape(cell2mat(interunit_params.on2off_velocity),...
                    size(interunit_params.on2off_velocity,1)*size(interunit_params.on2off_velocity,2),1)];
                Iy_cell_id = [Iy_cell_id;directories(i,:)];
            else
            end
            if strcmp('old',subfile),
                o_coll_spoffset = [o_coll_spoffset;...
                    reshape(cell2mat(interunit_params.spatial_offset),...
                    size(interunit_params.spatial_offset,1)*size(interunit_params.spatial_offset,2),1)];
                o_coll_onoffvel = [o_coll_onoffvel;...
                    reshape(cell2mat(interunit_params.on2off_velocity),...
                    size(interunit_params.on2off_velocity,1)*size(interunit_params.on2off_velocity,2),1)];
                Io_cell_id = [Io_cell_id;directories(i,:)];
            else
            end
    
            cd ..
        end

    case 1
        
        [~,subfile,~] = fileparts(sub_folder);
        
        if strcmp('young',subfile),
            y_coll_subtemponset = [];
            y_coll_temponset = [];
            y_coll_time2peak = [];
            y_coll_spatialext = [];
            y_coll_tempext = [];
            y_coll_peakamp = [];
            y_coll_STtilt = [];
            y_coll_STvelocity = []
            y_coll_spoffset = [];
            y_coll_onoffvel = [];
            y_cell_id = [];
            Iy_cell_id = [];
        else
        end
        if strcmp('old',subfile),
            o_coll_subtemponset = [];
            o_coll_temponset = [];
            o_coll_time2peak = [];
            o_coll_spatialext = [];
            o_coll_tempext = [];
            o_coll_peakamp = [];
            o_coll_STtilt = [];
            o_coll_STvelocity = [];
            o_coll_spoffset = [];
            o_coll_onoffvel = [];
            o_cell_id = [];
            Io_cell_id = [];
        else
        end
        
        for i = 1:size(directories,1),
            cd(directories(i,:));
            if ~exist('STRF_structure_TH80CI90.mat','file');
                display(['File not found for ',directories(i,:),' ...']);
                continue;
            else
            end
            display(['Retrieving data for ',directories(i,:),' ...']);
            load([sub_folder filesep directories(i,:) filesep 'STRF_structure_TH80CI90.mat'],'singleunit_params');
            load([sub_folder filesep directories(i,:) filesep 'STRF_structure_TH80CI90.mat'],'interunit_params');
            if strcmp('young',subfile),
                y_coll_subtemponset = [y_coll_subtemponset;cell2mat(singleunit_params.subunit_temporal_onset)];
                y_coll_temponset = [y_coll_temponset;singleunit_params.temporal_onset];
                y_coll_time2peak = [y_coll_time2peak;cell2mat(singleunit_params.time_to_peak)];
                y_coll_spatialext = [y_coll_spatialext;cell2mat(singleunit_params.spatial_extent)];
                y_coll_tempext = [y_coll_tempext;cell2mat(singleunit_params.temporal_extent)];
                y_coll_peakamp = [y_coll_peakamp;cell2mat(singleunit_params.peak_amplitude)];
                y_coll_STtilt = [y_coll_STtilt;cell2mat(singleunit_params.ST_tilt)];
                y_coll_STvelocity = [y_coll_STvelocity;cell2mat(singleunit_params.ST_velocity)];
                y_cell_id = [y_cell_id;directories(i,:)];
            else
            end
            if strcmp('old',subfile),
                y_coll_subtemponset = [o_coll_subtemponset;cell2mat(singleunit_params.subunit_temporal_onset)];
                o_coll_temponset = [o_coll_temponset;singleunit_params.temporal_onset];
                o_coll_time2peak = [o_coll_time2peak;cell2mat(singleunit_params.time_to_peak)];
                o_coll_spatialext = [o_coll_spatialext;cell2mat(singleunit_params.spatial_extent)];
                o_coll_tempext = [o_coll_tempext;cell2mat(singleunit_params.temporal_extent)];
                o_coll_peakamp = [o_coll_peakamp;cell2mat(singleunit_params.peak_amplitude)];
                o_coll_STtilt = [o_coll_STtilt;cell2mat(singleunit_params.ST_tilt)];
                o_coll_STvelocity = [o_coll_STvelocity;cell2mat(singleunit_params.ST_velocity)];
                o_cell_id = [o_cell_id;directories(i,:)];
            else 
            end
            if strcmp('young',subfile),
                y_coll_spoffset = [y_coll_spoffset;...
                    reshape(cell2mat(interunit_params.spatial_offset),...
                    size(interunit_params.spatial_offset,1)*size(interunit_params.spatial_offset,2),1)];
                y_coll_onoffvel = [y_coll_onoffvel;...
                    reshape(cell2mat(interunit_params.on2off_velocity),...
                    size(interunit_params.on2off_velocity,1)*size(interunit_params.on2off_velocity,2),1)];
                Iy_cell_id = [Iy_cell_id;directories(i,:)];
            else
            end
            if strcmp('old',subfile),
                o_coll_spoffset = [o_coll_spoffset;...
                    reshape(cell2mat(interunit_params.spatial_offset),...
                    size(interunit_params.spatial_offset,1)*size(interunit_params.spatial_offset,2),1)];
                o_coll_onoffvel = [o_coll_onoffvel;...
                    reshape(cell2mat(interunit_params.on2off_velocity),...
                    size(interunit_params.on2off_velocity,1)*size(interunit_params.on2off_velocity,2),1)];
                Io_cell_id = [Io_cell_id;directories(i,:)];
            else
            end
            
            cd ..
        end
  
        
        if run_from_file==1&&strcmp('young',subfile),
            
        %histograms - young 
            figure;
            subplot(5,2,1);
            utrange = ...
                min(y_coll_subtemponset)-(0.1*max(y_coll_subtemponset)):...
                0.02:max(y_coll_subtemponset)+(0.1*max(y_coll_subtemponset));
            nh = histc(y_coll_subtemponset,utrange);
            nh = nh(1:end-1);
            bin_centers = utrange(2:end)+utrange(1:end-1)/2;
            h1 = bar(bin_centers,nh,'b');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([utrange(1) utrange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Subunit temporal onset');
            clear nh;
            subplot(5,2,2);
            trange = ...
                min(y_coll_temponset)-(0.1*max(y_coll_temponset)):...
                0.02:max(y_coll_temponset)+(0.1*max(y_coll_temponset));
            nh = histc(y_coll_temponset,trange);
            nh = nh(1:end-1);
            bin_centers = trange(2:end)+trange(1:end-1)/2;
            h1 = bar(bin_centers,nh,'b');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([trange(1) trange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Temporal onset');
            clear nh;
            subplot(5,2,3);
            y_coll_time2peak = unique(y_coll_time2peak);
            tprange = min(y_coll_time2peak)-(0.1*max(y_coll_time2peak)):...
                0.01:max(y_coll_time2peak)+(0.1*max(y_coll_time2peak));
            nh = histc(y_coll_time2peak,tprange);
            nh = nh(1:end-1);
            bin_centers = tprange(2:end)+tprange(1:end-1)/2;
            h2 = bar(bin_centers,nh,'b');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([tprange(1) tprange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Time to peak');
            clear nh;
            subplot(5,2,4);
            srange = min(y_coll_spatialext)-(0.1*max(y_coll_spatialext)):...
                1:max(y_coll_spatialext)+(0.1*max(y_coll_spatialext));
            nh = histc(y_coll_spatialext,srange);
            nh = nh(1:end-1);
            bin_centers = srange(2:end)+srange(1:end-1)/2;
            h3 = bar(bin_centers,nh,'b');
            xlabel('space (deg)');
            ylabel('counts');
            xlim([srange(1) srange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Spatial extent');
            clear nh;
            subplot(5,2,5);
            terange = min(y_coll_tempext)-(0.1*max(y_coll_tempext)):...
                0.01:max(y_coll_tempext)+(0.1*max(y_coll_tempext));
            nh = histc(y_coll_tempext,terange);
            nh = nh(1:end-1);
            bin_centers = terange(2:end)+terange(1:end-1)/2;
            h4 = bar(bin_centers,nh,'b');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([terange(1) terange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Temporal extent');
            clear nh;
            subplot(5,2,6);
            prange = min(y_coll_peakamp)-(0.1*max(y_coll_peakamp)):...
                10:max(y_coll_peakamp)+(0.1*max(y_coll_peakamp));
            nh = histc(y_coll_peakamp,prange);
            nh = nh(1:end-1);
            bin_centers = prange(2:end)+prange(1:end-1)/2;
            h5 = bar(bin_centers,nh,'b');
            xlabel('Amplitude (0-255)');
            ylabel('counts');
            xlim([prange(1) prange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Peak amplitude');
            clear nh;
            subplot(5,2,7);
            strange = min(y_coll_STtilt)-(0.1*max(y_coll_STtilt)):...
                10:max(y_coll_STtilt)+(0.1*max(y_coll_STtilt));
            nh = histc(y_coll_STtilt,strange);
            nh = nh(1:end-1);
            bin_centers = strange(2:end)+strange(1:end-1)/2;
            h6 = bar(bin_centers,nh,'b');
            xlabel('Abs angle (deg)');
            ylabel('counts');
            %xlim([strange(1) strange(end)]);
            %ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) spatiotemporal tilt');
            clear nh
            subplot(5,2,8);
            sorange = min(y_coll_spoffset)-(0.1*max(y_coll_spoffset)):...
                0.1:max(y_coll_spoffset)+(0.1*max(y_coll_spoffset));
            nh = histc(y_coll_spoffset,sorange);
            nh = nh(1:end-1);
            bin_centers = sorange(2:end)+sorange(1:end-1)/2;
            h7 = bar(bin_centers,nh,'b');
            xlabel('space (deg)');
            ylabel('counts');
            xlim([sorange(1) sorange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Interunit spatial offsets');
            clear nh;
            subplot(5,2,9);
            vrange =  min(y_coll_onoffvel)-(0.1*max(y_coll_onoffvel)):...
                0.01:max(y_coll_onoffvel)+(0.1*max(y_coll_onoffvel));
            nh = histc(y_coll_onoffvel,vrange);
            nh = nh(1:end-1);
            bin_centers = vrange(2:end)+vrange(1:end-1)/2;
            h8 = bar(bin_centers,nh,'b');
            xlabel('velocity (deg/sec)');
            ylabel('counts');
            xlim([vrange(1) vrange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Interunit ON-OFF velocity');
            clear nh;
            subplot(5,2,10);
            stvrange = min(y_coll_STvelocity)-(0.1*max(y_coll_STvelocity)):...
                0.01:max(y_coll_STvelocity)+(0.1*max(y_coll_STvelocity));
            nh = histc(y_coll_STvelocity,stvrange);
            nh = nh(1:end-1);
            bin_centers = stvrange(2:end)+stvrange(1:end-1)/2;
            h9 = bar(bin_centers,nh,'b');
            xlabel('velocity (deg/sec)');
            ylabel('counts');
            %xlim([stvrange(1) stvrange(end)]);
            %ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Young) Spatiotemporal velocity');
            clear nh;
            if save_it == 1,
                save('young_collated_TH80CI90.mat');
            else
            end
            

        else

        %histograms - old
            figure;
            subplot(5,2,1);
            utrange = ...
                min(y_coll_subtemponset)-(0.1*max(y_coll_subtemponset)):...
                0.02:max(y_coll_subtemponset)+(0.1*max(y_coll_subtemponset));
            nh = histc(y_coll_subtemponset,utrange);
            nh = nh(1:end-1);
            bin_centers = utrange(2:end)+utrange(1:end-1)/2;
            h1 = bar(bin_centers,nh,'m');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([utrange(1) utrange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Subunit temporal onset');
            clear nh;
            subplot(5,2,2);
            trange = ...
                min(o_coll_temponset)-(0.1*max(o_coll_temponset)):...
                0.02:max(o_coll_temponset)+(0.1*max(o_coll_temponset));
            nh = histc(o_coll_temponset,trange);
            nh = nh(1:end-1);
            bin_centers = trange(2:end)+trange(1:end-1)/2;
            h1 = bar(bin_centers,nh,'m');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([trange(1) trange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Temporal onset');
            clear nh;
            subplot(5,2,3);
            o_coll_time2peak = unique(o_coll_time2peak);
            tprange = min(o_coll_time2peak)-(0.1*max(o_coll_time2peak)):...
                0.01:max(o_coll_time2peak)+(0.1*max(o_coll_time2peak));
            nh = histc(o_coll_time2peak,tprange);
            nh = nh(1:end-1);
            bin_centers = tprange(2:end)+tprange(1:end-1)/2;
            h2 = bar(bin_centers,nh,'m');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([tprange(1) tprange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Time to peak');
            clear nh;
            subplot(5,2,4);
            srange = min(o_coll_spatialext)-(0.1*max(o_coll_spatialext)):...
                1:max(o_coll_spatialext)+(0.1*max(o_coll_spatialext));
            nh = histc(o_coll_spatialext,srange);
            nh = nh(1:end-1);
            bin_centers = srange(2:end)+srange(1:end-1)/2;
            h3 = bar(bin_centers,nh,'m');
            xlabel('space (deg)');
            ylabel('counts');
            xlim([srange(1) srange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Spatial extent');
            clear nh;
            subplot(5,2,5);
            terange = min(o_coll_tempext)-(0.1*max(o_coll_tempext)):...
                0.01:max(o_coll_tempext)+(0.1*max(o_coll_tempext));
            nh = histc(o_coll_tempext,terange);
            nh = nh(1:end-1);
            bin_centers = terange(2:end)+terange(1:end-1)/2;
            h4 = bar(bin_centers,nh,'m');
            xlabel('time (sec)');
            ylabel('counts');
            xlim([terange(1) terange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Temporal extent');
            clear nh;
            subplot(5,2,6);
            prange = min(abs(o_coll_peakamp-127.5))-(0.1*max(abs(o_coll_peakamp-127.5))):...
                10:max(abs(o_coll_peakamp-127.5))+(0.1*max(abs(o_coll_peakamp-127.5)));
            nh = histc(abs(o_coll_peakamp-127.5),prange);
            nh = nh(1:end-1);
            bin_centers = prange(2:end)+prange(1:end-1)/2;
            h5 = bar(bin_centers,nh,'m');
            xlabel('Amplitude (0-128)');
            ylabel('counts');
            xlim([prange(1) prange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Peak amplitude');
            clear nh;
            subplot(5,2,7);
            strange = min(o_coll_STtilt)-(0.1*max(o_coll_STtilt)):...
                2:max(o_coll_STtilt)+(0.1*max(o_coll_STtilt));
            nh = histc(o_coll_STtilt,strange);
            nh = nh(1:end-1);
            bin_centers = strange(2:end)+strange(1:end-1)/2;
            h6 = bar(bin_centers,nh,'m');
            xlabel('Abs angle (deg)');
            ylabel('counts');
            %xlim([strange(1) strange(end)]);
            %ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) spatiotemporal tilt');
            clear nh
            subplot(5,2,8);
            sorange = min(o_coll_spoffset)-(0.1*max(o_coll_spoffset)):...
                0.1:max(o_coll_spoffset)+(0.1*max(o_coll_spoffset));
            nh = histc(o_coll_spoffset,sorange);
            nh = nh(1:end-1);
            bin_centers = sorange(2:end)+sorange(1:end-1)/2;
            h7 = bar(bin_centers,nh,'m');
            xlabel('space (deg)');
            ylabel('counts');
            xlim([sorange(1) sorange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Interunit spatial offsets');
            clear nh;
            subplot(5,2,9);
            vrange =  min(o_coll_onoffvel)-(0.1*max(o_coll_onoffvel)):...
                0.01:max(o_coll_onoffvel)+(0.1*max(o_coll_onoffvel));
            nh = histc(o_coll_onoffvel,vrange);
            nh = nh(1:end-1);
            bin_centers = vrange(2:end)+vrange(1:end-1)/2;
            h8 = bar(bin_centers,nh,'m');
            xlabel('velocity (deg/sec)');
            ylabel('counts');
            xlim([vrange(1) vrange(end)]);
            ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Interunit ON-OFF velocity');
            clear nh;
            subplot(5,2,10);
            stvrange = min(o_coll_STvelocity)-(0.1*max(o_coll_STvelocity)):...
                0.0001:max(o_coll_STvelocity)+(0.1*max(o_coll_STvelocity));
            nh = histc(o_coll_STvelocity,stvrange);
            nh = nh(1:end-1);
            bin_centers = stvrange(2:end)+stvrange(1:end-1)/2;
            h9 = bar(bin_centers,nh,'m');
            xlabel('velocity (deg/sec)');
            ylabel('counts');
            %xlim([stvrange(1) stvrange(end)]);
            %ylim([0 max(nh)+(0.1*max(nh))]);
            title('(Old) Spatiotemporal velocity');
            clear nh;
            
            if save_it == 1,
                save('old_collated_TH80CI90.mat');
            else
            end
        end
        
end
        


%comparison plots
cd ../old;
load(['old' filesep 'old_collated_TH80CI90.mat']);
cd ../young;
load(['young' filesep 'young_collated_TH80CI90.mat']);
%Subunit temporal onset
figure;
subplot(5,2,1);
scatter(repmat(1.0,length(y_coll_subtemponset),1),y_coll_subtemponset,'b');
hold on;
b1 = bar(1.0,median(y_coll_subtemponset));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_subtemponset),1),o_coll_subtemponset,'m');
b2 = bar(3.0,median(o_coll_subtemponset));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Time (sec)');
[~,p] = kstest2(y_coll_temponset,o_coll_temponset);s
title(['Temporal onset, K-S test p = ',num2str(p)]);
hold off;
%Response temporal onset
subplot(5,2,2);
scatter(repmat(1.0,length(y_coll_temponset),1),y_coll_temponset,'b');
hold on;
b1 = bar(1.0,median(y_coll_temponset));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_temponset),1),o_coll_temponset,'m');
b2 = bar(3.0,median(o_coll_temponset));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Time (sec)');
[~,p] = kstest2(y_coll_temponset,o_coll_temponset);
title(['Temporal onset, K-S test p = ',num2str(p)]);
hold off;
%Time to peak (all subunits)
subplot(5,2,3);
scatter(repmat(1.0,length(y_coll_time2peak),1),y_coll_time2peak,'b');
hold on;
b1 = bar(1.0,median(y_coll_time2peak));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_time2peak),1),o_coll_time2peak,'m');
b2 = bar(3.0,median(o_coll_time2peak));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Time (sec)');
[~,p] = kstest2(y_coll_time2peak,o_coll_time2peak);
title(['Time to peak, K-S test p = ',num2str(p)]);
hold off;
%Spatial extent
subplot(5,2,4);
scatter(repmat(1.0,length(y_coll_spatialext),1),y_coll_spatialext,'b');
hold on;
b1 = bar(1.0,median(y_coll_spatialext));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_spatialext),1),o_coll_spatialext,'m');
b2 = bar(3.0,median(o_coll_spatialext));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Space (deg)');
[~,p] = kstest2(y_coll_spatialext,o_coll_spatialext);
title(['Spatial extent, K-S test p = ',num2str(p)]);
hold off;
%Temporal extent
subplot(5,2,5);
scatter(repmat(1.0,length(y_coll_tempext),1),y_coll_tempext,'b');
hold on;
b1 = bar(1.0,median(y_coll_tempext));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_tempext),1),o_coll_tempext,'m');
b2 = bar(3.0,median(o_coll_tempext));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Time (sec)');
[~,p] = kstest2(y_coll_tempext,o_coll_tempext);
title(['Temporal extent, K-S test p = ',num2str(p)]);
hold off;
%Peak amplitude
subplot(5,2,6);
scatter(repmat(1.0,length(y_coll_peakamp),1),y_coll_peakamp,'b');
hold on;
b1 = bar(1.0,median(y_coll_peakamp));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_peakamp),1),o_coll_peakamp,'m');
b2 = bar(3.0,median(o_coll_peakamp));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Amplitude ()');
[~,p] = kstest2(y_coll_peakamp,o_coll_peakamp);
title(['Peak amplitude, K-S test p = ',num2str(p)]);
hold off;
%ST tilt angle
subplot(5,2,7);
scatter(repmat(1.0,length(y_coll_STtilt),1),y_coll_STtilt,'b');
hold on;
b1 = bar(1.0,median(y_coll_STtilt));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_STtilt),1),o_coll_STtilt,'m');
b2 = bar(3.0,median(o_coll_STtilt));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Abs. angle (deg)');
[~,p] = kstest2(y_coll_STtilt,o_coll_STtilt);
title(['Spatiotemporal tilt, K-S test p = ',num2str(p)]);
hold off;
%Interunit spatial offset
subplot(5,2,8);
scatter(repmat(1.0,length(y_coll_spoffset),1),y_coll_spoffset,'b');
hold on;
b1 = bar(1.0,median(y_coll_spoffset));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_spoffset),1),o_coll_spoffset,'m');
b2 = bar(3.0,median(o_coll_spoffset));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Space (deg)');
[~,p] = kstest2(y_coll_spoffset,o_coll_spoffset);
title(['Interunit spatial offsets, K-S test p = ',num2str(p)]);
hold off;
%Interunit on-off velocity (peak-to-peak)
subplot(5,2,9);
scatter(repmat(1.0,length(y_coll_onoffvel),1),y_coll_onoffvel,'b');
hold on;
b1 = bar(1.0,median(y_coll_onoffvel));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_onoffvel),1),o_coll_onoffvel,'m');
b2 = bar(3.0,median(o_coll_onoffvel));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Velocity (deg/sec)');
[~,p] = kstest2(y_coll_onoffvel,o_coll_onoffvel);
title(['Interunit ON-OFF velocity, K-S test p = ',num2str(p)]);
hold off;
%ST velocity (deg/sec)
subplot(5,2,10);
scatter(repmat(1.0,length(y_coll_STvelocity),1),y_coll_STvelocity,'b');
hold on;
b1 = bar(1.0,median(y_coll_STvelocity));
b1.FaceColor = 'none';
scatter(repmat(3.0,length(o_coll_STvelocity),1),o_coll_STvelocity,'m');
b2 = bar(3.0,median(o_coll_STvelocity));
b2.FaceColor = 'none';
ax = gca;
ax.XTickLabel = {[],'Young',[],'Old',[]};
ylabel('Velocity (deg/sec)');
[~,p] = kstest2(y_coll_STvelocity,o_coll_STvelocity);
title(['Spatiotemporal velocity, K-S test p = ',num2str(p)]);
hold off;
if save_it == 1,
    cd ../
    saveas(gcf,[pwd filesep 'parameter_comp_TH80CI90.fig']);
else
end







end


