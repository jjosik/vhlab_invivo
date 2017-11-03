function [ output_args ] = get_STAdata( plot_cells,target_FR )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


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

o_collated_sta_pref = {};
o_collated_sta_null = {};
o_cell_t_win = {};
o_cell_SR = [];
y_collated_sta_pref = {};
y_collated_sta_null = {};
y_cell_t_win = {};
y_cell_SR = [];
if plot_cells == 1,
    f = figure;
    %... under construction
else
end

[~,check_dir,~] = fileparts(sub_folder);
for i = 1:size(directories,1),
    cd(directories(i,:));
    dd = dir(fullfile(directories(i,:),'*.mat'));
    cfilename = uigetfile('*.mat');
    if cfilename == 0,
        cd ..
        continue
    else
    end
    load([sub_folder filesep directories(i,:) filesep cfilename],...
        'sta_vm_pref','sta_vm_null','t_win','sample_interval');
    if strcmp(check_dir,'old'),
        o_collated_sta_pref{end+1,1} = sta_vm_pref;
        o_collated_sta_null{end+1,1} = sta_vm_null;
        o_cell_t_win{end+1,1} = t_win;
        o_cell_SR(end+1,1) = 1/sample_interval;
    else
        y_collated_sta_pref{end+1,1} = sta_vm_pref;
        y_collated_sta_null{end+1,1} = sta_vm_null;
        y_cell_t_win{end+1,1} = t_win;
        y_cell_SR(end+1,1) = 1/sample_interval;
    end
    cd ../
end
cd ../
%first need to resample all traces
target_FR = 10000;  %target frame rate
o_collated_sta_pref_RES = cell(length(o_collated_sta_pref),1);
o_collated_sta_null_RES = cell(length(o_collated_sta_null),1);
for j = 1:length(o_collated_sta_pref),
    [p,q] = rat(target_FR / o_cell_SR(j,1));
    temp_res = resample(o_collated_sta_pref{j,1},p,q);
    o_collated_sta_pref_RES{j,1} = temp_res;
    temp_res_ = resample(o_collated_sta_null{j,1},p,q);
    o_collated_sta_null_RES{j,1} = temp_res_;
end
y_collated_sta_pref_RES = cell(length(y_collated_sta_pref),1);
y_collated_sta_null_RES = cell(length(y_collated_sta_null),1);
for k = 1:length(y_collated_sta_pref),
    [p,q] = rat(target_FR / y_cell_SR(k,1));
    temp_res = resample(y_collated_sta_pref{k,1},p,q);
    y_collated_sta_pref_RES{k,1} = temp_res;
    temp_res_ = resample(y_collated_sta_null{k,1},p,q);
    y_collated_sta_null_RES{k,1} = temp_res_;
end
%new_t_win = 1/target_FR:1/target_FR:(length(temp_res)/target_FR);
new_t_win = -(length(temp_res)/target_FR)/2:1/target_FR:((length(temp_res)/target_FR)/2-1/target_FR);
if save_it == 1,
    save('collected_VmSTA.mat','o_collated_sta_pref','o_collated_sta_null',...
        'o_cell_t_win','y_collated_sta_pref','y_collated_sta_null','y_cell_t_win',...
        'o_collated_sta_pref_RES','o_collated_sta_null_RES','y_collated_sta_pref_RES',...
        'y_collated_sta_null_RES','new_t_win');
else
end

%generate 2x2 plot of grand averages
o_AVEsta_pref = zeros(length(new_t_win),1);
o_AVEsta_null = zeros(length(new_t_win),1);
y_AVEsta_pref = zeros(length(new_t_win),1);
y_AVEsta_null = zeros(length(new_t_win),1);
for j = 1:length(o_collated_sta_pref_RES),
    o_AVEsta_pref = o_AVEsta_pref + o_collated_sta_pref_RES{j,1};
    o_AVEsta_null = o_AVEsta_null + o_collated_sta_null_RES{j,1};
end
o_AVEsta_pref = o_AVEsta_pref./length(o_collated_sta_pref_RES);
o_AVEsta_null = o_AVEsta_null./length(o_collated_sta_null_RES);
nan_flag = 0;
for k = 1:length(y_collated_sta_pref_RES),
    if ~isnan(y_collated_sta_null_RES{k,1}),
        y_AVEsta_pref = y_AVEsta_pref + y_collated_sta_pref_RES{k,1};
        y_AVEsta_null = y_AVEsta_null + y_collated_sta_null_RES{k,1};
    else
        nan_flag = nan_flag + 1;
        continue;
    end
end
y_AVEsta_pref = y_AVEsta_pref./length(y_collated_sta_pref_RES);
y_AVEsta_null = y_AVEsta_null./(length(y_collated_sta_null_RES)-nan_flag);

f40 = figure;
plot(new_t_win,o_AVEsta_pref,'r');
hold on;
plot(new_t_win,o_AVEsta_null,'b');
plot(new_t_win,y_AVEsta_pref,'r--');
plot(new_t_win,y_AVEsta_null,'b--');
xlabel('Time (seconds)');
ylabel('Membrane potential (mV)');
legend('Old PREF','Old NULL','Young PREF','Young NULL','Location','NorthWest');
hold off;
xlim([-0.05 0.05]);

    

end

