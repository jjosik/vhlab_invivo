function [ output_args ] = compare_collVF( plot_cells,save_it )
%UNTITLED2 Summary of this function goes here
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

all_traces = cell(size(directories,1),2);
all_x = cell(size(directories,1),2);
all_slopes = cell(size(directories,1),2);
for i = 1:size(directories,1),
    [ temp_pslope,temp_nslope,p_x,gen_ptrace,n_x,gen_ntrace ] = ...
        batch_general_VFfit( plot_cells,save_it );
    all_traces{i,1} = gen_ptrace;
    all_traces{i,2} = gen_ntrace;
    all_x{i,1} = p_x;
    all_x{i,2} = n_x;
    all_slopes{i,1} = temp_pslope;
    all_slopes{i,1} = temp_nslope;
end

figure;
c_array = {'r','b';'r','b'};
style_array = {':',':';'--','--'};
for i = 1:size(all_traces,1),
    for j = 1:size(all_traces,2),
        ln = plot(all_x{i,j},all_traces{i,j},'Color',c_array{i,j});
        ln.LineStyle = style_array{i,j};
        ln.LineWidth = 2.5;
        hold on;
        xlabel('Vm - est. Vth (mV)');
        ylabel('Firing rate (spikes/sec)');
        legend('Young, PREF','Young, NULL','Old, PREF','Old, NULL','Location','NorthWest');
    end
end

figure;
mean_young_p = mean(all_slopes{1,1});
mean_young_n = mean(all_slopes{1,2});
mean_old_p = mean(all_slopes{2,1});
mean_old_n = mean(all_slopes{2,2});
se_young_p = std(all_slopes{1,1})./sqrt(length(all_slopes{1,1}));
se_young_n = std(all_slopes{1,2})./sqrt(length(all_slopes{1,2}));
se_old_p = std(all_slopes{2,1})./sqrt(length(all_slopes{2,1}));
se_old_n = std(all_slopes{2,2})./sqrt(length(all_slopes{2,2}));
b1 = errorbar(1.0,mean_young_p,se_young_p);
b2 = errorbar(2.0,mean_young_n,se_young_n);
b3 = errorbar(3.0,mean_old_p,se_old_p);
b4 = errorbar(4.0,mean_old_n,se_old_n);


        
    
    


end

