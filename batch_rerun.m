function [ output_args ] = batch_rerun( collate_ages,varargin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

save_it = 1;
overwrite = 0;
display_spikes = 0;
use_flankers = 1;
use_exclusive_AG = 0;  %will likely have complementary value to use_flankers (but not strictly necessary)
analyze_bycycle = 0;
auto_detect = 0; %set to 0 if user-defined coarse threshold should be queried for spike detection 
                 %(recommended in traces with electrical drift)
use_detrend = 1;
default = 2; %any value >0 will set VF plot type to a default, thereby preventing code from
             %being interrupted by query (which is generally desirable for
             %batch analyses)
use_global_filter = 0; %sub-option to be used in conjunction with filter_type = 2 (doesn't do anything otherwise)
filter_type = 1; 
%NOTE - filter_types:  if using rasters2instFR_global2, choices are:
%1=basic fixed gaussian smoother,2=optimal smoother,3=boxcar
%if using rasters2instFR_global3 (RECOMMENDED), choices are:
%1=fixed window gaussian smoother,2=optimal window smoother
use_flanktest = 0;
adapt_bins = 1;
anchor_Vth = 1;
fit_it = 1;
model = 0;  %set to '0' if rectified power law, '1' if rectilinear
code = '_M0_F_Bd';
if nargin < 2,
    use_pre_analyzed = 0;
else
    use_pre_analyzed = varargin{1};
end

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

o_cell_prefmean = [];
o_cell_nullmean = [];
o_cell_Vth_pref = [];
o_cell_Vth_null = [];
%o_cell_CV_pref = [];
%o_cell_CV_null = [];
o_cell_std_pref = [];
o_cell_std_null = [];
o_cell_AGprefmean = [];
o_cell_AGnullmean = [];
y_cell_prefmean = [];
y_cell_nullmean = [];
y_cell_Vth_pref = [];
y_cell_Vth_null = [];
%y_cell_CV_pref = [];
%y_cell_CV_null = [];
y_cell_std_pref = [];
y_cell_std_null = [];
y_cell_AGprefmean = [];
y_cell_AGnullmean = [];
cell_prefmean = [];
cell_nullmean = [];
cell_Vth_pref = [];
cell_Vth_null = [];
%cell_CV_pref = [];
%cell_CV_null = [];
cell_std_pref = [];
cell_std_null = [];
cell_AGprefmean = [];
cell_AGnullmean = [];

[~,subfile,~] = fileparts(sub_folder);
for i = 1:size(directories,1),
    cur_dir = directories(i,:);
    %try
        cd(directories(i,:));
    %catch
    %    cd ../
    %    cd(directories(i,:));
    %end
    load([sub_folder filesep directories(i,:) filesep 'stims.mat']);
    load([sub_folder filesep directories(i,:) filesep 'data.mat']);
    load([sub_folder filesep directories(i,:) filesep 'stimorder.mat']);
    load([sub_folder filesep directories(i,:) filesep 'stimvalues.mat']);
    load([sub_folder filesep directories(i,:) filesep 'analyzed_spike.mat']);
    [m_pref,m_null,AG_m_pref,AG_m_null,Vth_est_pcoll,Vth_est_ncoll,...
        p_subth_std,n_subth_std] = invivo_subth_rerun(save_it,overwrite,model,...
        use_pre_analyzed,display_spikes,use_flankers,use_exclusive_AG,...
        analyze_bycycle,auto_detect,use_detrend,use_trial_baseline,sub_folder,cur_dir,...
        default,use_global_filter,filter_type,use_flanktest,adapt_bins,...
        anchor_Vth,fit_it,code);
    if collate_ages == 1,
        cell_prefmean = [cell_prefmean;m_pref];
        cell_nullmean = [cell_nullmean;m_null];
        cell_Vth_pref = [cell_Vth_pref;Vth_est_pcoll];
        cell_Vth_null = [cell_Vth_null;Vth_est_ncoll];
        %cell_CV_pref = [cell_CV_pref;nanmean(p_subth_CV)];
        %cell_CV_null = [cell_CV_null;nanmean(n_subth_CV)];
        cell_AGprefmean = [cell_AGprefmean;AG_m_pref];
        cell_AGnullmean = [cell_AGnullmean;AG_m_null];
    else
        if strcmp('old',subfile),
            o_cell_prefmean = [o_cell_prefmean;m_pref];
            o_cell_nullmean = [o_cell_nullmean;m_null];
            o_cell_Vth_pref = [o_cell_Vth_pref;Vth_est_pcoll];
            o_cell_Vth_null = [o_cell_Vth_null;Vth_est_ncoll];
            %o_cell_CV_pref = [o_cell_CV_pref;nanmean(p_subth_CV)];
            %o_cell_CV_null = [o_cell_CV_null;nanmean(n_subth_CV)];
            o_cell_AGprefmean = [o_cell_AGprefmean;AG_m_pref];
            o_cell_AGnullmean = [o_cell_AGnullmean;AG_m_null];
        else
            y_cell_prefmean = [y_cell_prefmean;m_pref];
            y_cell_nullmean = [y_cell_nullmean;m_null];
            y_cell_Vth_pref = [y_cell_Vth_pref;Vth_est_pcoll];
            y_cell_Vth_null = [y_cell_Vth_null;Vth_est_ncoll];
            %y_cell_CV_pref = [y_cell_CV_pref;nanmean(p_subth_CV)];
            %y_cell_CV_null = [y_cell_CV_null;nanmean(n_subth_CV)];
            y_cell_AGprefmean = [y_cell_AGprefmean;AG_m_pref];
            y_cell_AGnullmean = [y_cell_AGnullmean;AG_m_null];
        end
    end
    cd ..

end

% young, pref/null biophysical threshold comparison
f1 = figure;
scatter(y_cell_prefmean,y_cell_nullmean,'k','filled');
hold on;
line([-80 -40],[-80 -40]);
xlabel('PREF biophys. threshold (mV)');
ylabel('NULL biophys. threshold (mV)');
title('Young PREF/NULL baseline-corrected biophys. threshold comparison by cell');
hold off;
% old, pref/null biophysical threshold comparison
f2 = figure;
scatter(o_cell_prefmean,o_cell_nullmean,'k','filled');
hold on;
line([-80 -40],[-80 -40]);
xlabel('PREF biophys. threshold (mV)');
ylabel('NULL biophys. threshold (mV)');
title('Old PREF/NULL baseline-corrected biophys. threshold comparison by cell');
hold off;
% combined pref/null biophysical threshold comparison
f3 = figure;
cell_prefmean = [y_cell_prefmean;o_cell_prefmean];
cell_nullmean = [y_cell_nullmean;o_cell_nullmean];
scatter(cell_prefmean,cell_nullmean,'k','filled');
hold on;
line([-80 -40],[-80 -40]);
xlabel('PREF biophys. threshold (mV)');
ylabel('NULL biophys. threshold (mV)');
title('Combined PREF/NULL biophys. threshold comparison by cell');
hold off;

%young, pref/null pre-stim. relative threshold comparison
f10 = figure;
scatter(y_cell_AGprefmean,y_cell_AGnullmean,'k','filled');
hold on;
line([-10 40],[-10 40]);
xlabel('PREF relative threshold (mV)');
ylabel('NULL relative threshold (mV)');
title('Young PREF/NULL pre-stim. relative threshold comparison across cells');
hold off;
%old, pref/null pre.-stim relative threshold comparison
f11 = figure;
scatter(o_cell_AGprefmean,o_cell_AGnullmean,'k','filled');
hold on;
line([-10 40],[-10 40]);
xlabel('PREF relative threshold (mV)');
ylabel('NULL relative threshold (mV)');
title('Old PREF/NULL pre-stim. relative threshold comparison across cells');
hold off;
%combined pref/null pre-stim. relative threshold comparison
f12 = figure;
scatter(cell_AGprefmean,cell_AGnullmean,'k','filled');
hold on;
line([-10 30],[-10 30]);
xlabel('PREF relative threshold (mV)');
ylabel('NULL relative threshold (mV)');
title('Combined PREF/NULL pre-stim. relative threshold comparison across cells');
hold off;

% young, pref/null fit threshold comparison
f4 = figure;
scatter(y_cell_Vth_pref,y_cell_Vth_null,'k','filled');
hold on;
line([-100 0],[-100 0]);
xlabel('PREF mean subthreshold Vm (mV)');
ylabel('NULL mean subthreshold Vm (mV)');
title('Young PREF/NULL fit subthreshold Vm comparison by cell');
hold off;
% old, pref/null fit subthreshold Vm comparison
f5 = figure;
scatter(o_cell_Vth_pref,o_cell_Vth_null,'k','filled');
hold on;
line([-100 0],[-100 0]);
xlabel('PREF mean subthreshold Vm (mV)');
ylabel('NULL mean subthreshold Vm (mV)');
title('Old PREF/NULL fit subthreshold Vm comparison by cell');
hold off;
% combined pref/null fit subthreshold Vm comparison
f6 = figure;
cell_Vth_pref = [y_cell_Vth_pref;o_cell_Vth_pref];
cell_Vth_null = [y_cell_Vth_null;o_cell_Vth_null];
scatter(cell_Vth_pref,cell_Vth_null,'k','filled');
hold on;
line([-100 0],[-100 0]);
xlabel('PREF mean subthreshold Vm (mV)');
ylabel('NULL mean subthreshold Vm (mV)');
title('Combined PREF/NULL fit subthreshold Vm comparison by cell');
hold off;

%consider adding pref-to-pref subth Vm to threshold difference (..also
%corr. null-to-null)

% young, pref/null CV comparison
%f7 = figure;
%scatter(y_cell_CV_pref,y_cell_CV_null,'k','filled');
%hold on;
%line([-15 5],[-15 5]);
%xlabel('PREF Vm CV');
%ylabel('NULL Vm CV');
%title('Young PREF/NULL Vm CV comparison by cell');
%hold off;
% old, pref/null CV comparison
%f8 = figure;
%scatter(o_cell_CV_pref,o_cell_CV_null,'k','filled');
%hold on;
%line([-15 5],[-15 5]);
%xlabel('PREF Vm CV');
%ylabel('NULL Vm CV');
%title('Old PREF/NULL Vm CV comparison by cell');
%hold off;
% combined pref/null CV comparison
%cell_CV_pref = [y_cell_CV_pref;o_cell_CV_pref];
%cell_CV_null = [y_cell_CV_null;o_cell_CV_null];
%f9 = figure;
%scatter(cell_CV_pref,cell_CV_null,'k','filled');
%hold on;
%line([-15 5],[-15 5]);
%xlabel('PREF Vm CV');
%ylabel('NULL Vm CV');
%title('Combined PREF/NULL Vm CV comparison by cell');

% young, pref/null std comparison
f101 = figure;
scatter(y_cell_std_pref,y_cell_std_null,'k','filled');
hold on;
line([-15 15],[-15 15]);
xlabel('PREF Vm STD');
ylabel('NULL Vm STD');
title('Young PREF/NULL Vm STD comparison by cell');
hold off;
% old, pref/null std comparison
f102 = figure;
scatter(o_cell_std_pref,o_cell_std_null,'k','filled');
hold on;
line([-15 15],[-15 15]);
xlabel('PREF Vm STD');
ylabel('NULL Vm STD');
title('Old PREF/NULL Vm STD comparison by cell');
hold off;
% combined pref/null std comparison
cell_CV_pref = [y_cell_std_pref;o_cell_std_pref];
cell_CV_null = [y_cell_std_null;o_cell_std_null];
f103 = figure;
scatter(cell_std_pref,cell_std_null,'k','filled');
hold on;
line([-15 15],[-15 15]);
xlabel('PREF Vm STD');
ylabel('NULL Vm STD');
title('Combined PREF/NULL Vm STD comparison by cell');



end

