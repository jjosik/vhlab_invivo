function [ temp_pslope,temp_nslope,p_x,gen_ptrace,n_x,gen_ntrace ] = ...
    batch_general_VFfit( plot_cells,save_it )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%machine_path = input(['Enter ''1'' for Alectryon path, ''2''',...
%    ' for osik_mac path: ']);

%switch machine_path
%    
%    case 1
%       %Alectryon path:
%        top = ('/Users/vhlab/Documents/intracellular_data/');
%        
%    case 2
%        %osik_mac path:
%        top = ('Users/osik_mac/Documents/MATLAB/data');
%end
%p = addpath(genpath(top));

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
    set(0,'currentfigure',fp);
    gcf;
    [trace,params,gof,fitinfo] = tanhfit(bin_centers,temp_ymeans);
    plot(bin_centers,trace,'Color',[0.3 0 0]);
    ylim([0 180]);
    temp_pslope(i,1) = params.b/params.c;
    hold on;
    
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
    set(0,'currentfigure',fn);
    gcf;
    [ntrace,nparams,ngof,nfitinfo] = tanhfit(nbin_centers,ntemp_ymeans);
    plot(nbin_centers,ntrace,'Color',[0 0 0.3]);
    ylim([0 180]);
    temp_nslope(i,1) = nparams.b/nparams.c;
    hold on;
    if i == size(directories,1),
        set(0,'currentfigure',fn);
        gcf;
        [gen_ntrace,gen_nparams,gen_ngof,gen_nfitinfo] = ...
            tanhfit(n_collated_bin_centers,n_collated_bin_ymean);
        n_x = floor(min(n_collated_bin_centers)):1:floor(max(n_collated_bin_centers));
        fullnt = plot(n_x,gen_ntrace,'b');
        fullnt.LineWidth = 3.0;
        xlabel('Vm - est. Vth (mV)');
        ylabel('Firing rate (spikes/sec)');
        ylim([0 180]);
        n_c50_slope = gen_nparams.b/gen_nparams.c;
        [~,nd] = fileparts(sub_folder);
        title([nd,', NULL direction, ',num2str(gen_nparams.a),...
            ' + ',num2str(gen_nparams.b),'.*tanh((x - ',...
            num2str(gen_nparams.c),')./',num2str(gen_nparams.d),'), Slope: ',num2str(n_c50_slope)]);
        if save_it == 1,
            saveas(gcf,[sub_folder filesep directories(i,1),'NULL_collVF.fig']);
        else
        end
        
        set(0,'currentfigure',fp);
        gcf;
        [gen_ptrace,gen_pparams,gen_pgof,gen_pfitinfo] = ...
            tanhfit(p_collated_bin_centers,p_collated_bin_ymean);
        p_x = floor(min(p_collated_bin_centers)):1:floor(max(p_collated_bin_centers));
        fullpt = plot(p_x,gen_ptrace,'r');
        fullpt.LineWidth = 3.0;
        xlabel('Vm - est. Vth (mV)');
        ylabel('Firing rate (spikes/sec)');
        ylim([0 180]);
        [~,pd] = fileparts(sub_folder);
        p_c50_slope = gen_pparams.b/gen_pparams.c;
        title([pd,', PREF direction, ',num2str(gen_pparams.a),...
            ' + ',num2str(gen_pparams.b),'.*tanh((x - ',...
            num2str(gen_pparams.c),')./',num2str(gen_pparams.d),'), Slope: ',num2str(p_c50_slope)]);
        if save_it == 1,
            saveas(gcf,[sub_folder filesep directories(i,1),'PREF_collVF.fig']);
        else
        end
    else
    end
    
    cd ..
end
collated_bin_ymean = [p_collated_bin_ymean;n_collated_bin_ymean];
collated_bin_centers = [p_collated_bin_centers;n_collated_bin_centers];
  
s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-Inf,0,-100,1e-12],...
    'Upper',[Inf,Inf,100,1e12],...
    'Startpoint',[0,1,0,1]);
ft = fittype('a + b.*tanh((x - c)./d)','options',s);
[params,gof,fitinfo] = fit(collated_bin_centers,collated_bin_ymean,ft);
s_x = [floor(min(collated_bin_centers)):1:ceil(max(collated_bin_centers))];
trace = params.a + params.b.*tanh((s_x - params.c)./params.d);

figure;
scatter(p_collated_bin_centers,p_collated_bin_ymean,'r');
hold on;
scatter(n_collated_bin_centers,n_collated_bin_ymean,'b');
htline = plot(s_x,trace,'k');
htline.LineWidth = 2.5;
xlabel('Vm - est. Vth (mV)');
ylabel('Firing rate (spikes/sec)');
legend('Pref','Null','Location','NorthWest');
ylim([0 180]);
[~,pd] = fileparts(pwd);
c50_slope = params.b/params.c;
title(['General fit, ' pd,': ',num2str(params.a),' + ',num2str(params.b),'.*tanh((x - ',...
    num2str(params.c),')./',num2str(params.d),'), Slope: ',num2str(c50_slope)]);
if save_it == 1,
    saveas(gcf,[sub_folder filesep 'VFgeneral.fig']);
else
end


%alt. exponential fit
%s_e = fitoptions('Method','NonlinearLeastSquares',...
%    'Lower',[-100,0,1],...
%    'Upper',[100,Inf,50],...
%    'Startpoint',[0,1,2]);
%fte = fittype('a + b.*exp(x.*c)','options',s_e);
%[eparams,egof,efitinfo] = fit(collated_bin_centers,collated_bin_ymean,fte);
%e_x = [floor(min(collated_bin_centers)):1:ceil(max(collated_bin_centers))];
%etrace = eparams.a.*exp(eparams.b.*e_x);

%figure;
%scatter(p_collated_bin_centers,p_collated_bin_ymean,'r');
%hold on;
%%cg = [0.6 0.6 0.6];
%%scatter(n_collated_bin_centers,n_collated_bin_ymean,[],cg);
%scatter(n_collated_bin_centers,n_collated_bin_ymean,'b');
%eline = plot(e_x,etrace,'k');
%eline.LineWidth = 2.5;
%xlabel('Vm - est. Vth (mV)');
%ylabel('Firing rate (spikes/sec)');
%legend('Pref','Null','Location','NorthWest');
%[~,pd] = fileparts(pwd);
%title(['General fit, ' pd,': ',num2str(eparams.a),' + ',num2str(eparams.b),...
%    '.*e^(',num2str(eparams.c),'.*x)']);

            
        
    
    

end

