function [ f,h,bin_edges ] = ...
    VF_densitymap_C( coll_Vm,coll_FR,adapt_bins,stimset,scale_z,...
    match_PN,h_margins)


    f = figure;
    if adapt_bins == 1,
        if match_PN == 1,
            loc_min = ...
                min(coll_Vm);
            loc_max = ...
                max(coll_Vm);
            low = round(loc_min/10)*10;
            high = low+50;
            if loc_max>high,
                high = ceil(loc_max/10)*10;
            else
            end
        else
            l_est = min(coll_Vm);
            h_est = max(coll_Vm);
            low = round(l_est/10)*10;
            high = low+50;
            if h_est>high,
                high = ceil(h_est/10)*10;
            else
            end
        end
        v_bin_edges = low:1:high;
        if length(v_bin_edges)==51,
            f_bin_edges = 0:5:250;
        else
            if mod(250,length(v_bin_edges)-1)==0
                f_bin_edges = 0:250/(length(v_bin_edges)-1):250;
            else
                link_div = ceil(250/(length(v_bin_edges)-1));
                f_bin_edges = 0:link_div:(link_div*(length(v_bin_edges)-1));
            end
        end
        bin_edges = {v_bin_edges f_bin_edges};
        if iscell(coll_Vm),
            hc = hist3([coll_Vm{stimset,1}(1:end-1,1),coll_FR{stimset,1}(1:end,1)],'Edges',bin_edges);
        else
            hc = hist3([coll_Vm(1:end-1,1),coll_FR(1:end-1)],'Edges',bin_edges);
        end
    else
        N_bins = [50,50];
        if iscell(coll_Vm),
            hc = hist3([coll_Vm{stimset,1}(1:end-1,1),coll_FR{stimset,1}(1:end,1)],N_bins);
        else
            hc = hist3([coll_Vm(1:end-1,1),coll_FR(1:end-1,1)],N_bins);
        end
    end
    %n = (sqrt(10)-1)/median(reshape(hc,(size(hc,1)*size(hc,2)),1));  %necessary iff line 39 in effect
    if scale_z == 1,
        hc1 = 1+log10(hc');  %remove log10() and offset and multiplier to get back to linear map
        %hc1 = log10(1+n.*hc');
    else
        hc1 = hc';
        %hc1 = 3*(hc' > 0);
    end
    hc1(size(hc,1)+1,size(hc,2)+1) = 0;
    if iscell(coll_Vm),
        xc = linspace(min(coll_Vm{stimset,1}(1:end-1,1)),max(coll_Vm{stimset,1}(1:end-1,1)),size(hc,1)+1);
        yc = linspace(min(coll_FR{stimset,1}(1:end,1)),max(coll_FR{stimset,1}(1:end,1)),size(hc,1)+1);
    else
        coll_FR = reshape(coll_FR,1,length(coll_FR));
        xc = linspace(min(coll_Vm(1:end-1,1)),max(coll_Vm(1:end-1,1)),size(hc,1)+1);
        yc = linspace(min(coll_FR(1,1:end-1)),max(coll_FR(1,1:end-1)),size(hc,1)+1);
    end
    %set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    h = pcolor(xc,yc,hc1);
    %newmap = parula(256);
    %newmap(1,:) = [1 1 1];
    %colormap(newmap);
    %h.ZData = ones(size(hc1))*(-max(max(hc)));
    %caxis(log([h(1) h(length(h))]));
    xlabel('Membrane potential (mV)');
    ylabel('Firing rate');
    axis square;
    
    hold on;


end



