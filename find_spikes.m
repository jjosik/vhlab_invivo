function [ spike_locations,h_margins ] = find_spikes( spike_trace,t_vec,...
    sampling_rate,display_spikes,auto_detect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h_margins = 0.5;
    if auto_detect == 1,
        th_ = -20.0;   %in mV
    else
        f = figure;
        plot(t_vec,spike_trace,'k');
        hold on;
        display(['Use cursor to select coarse detection threshold for the experiment trace.  Press ''Return'' when finished']);
        [~,th_] = ginput;
    end

    spike_locations = [];
    spike_stepLocations = find(spike_trace>th_);
    if isempty(spike_stepLocations),
        msgID = 'MATLAB:th_test';
        msgtext = 'No spikes were detected in this dataset.';
        ME = MException(msgID,msgtext);
        throw(ME)
    else
    end
    onsets = diff(spike_stepLocations);
    sp_edge = find(onsets>1);
    ind_edge = [spike_stepLocations(1,1);spike_stepLocations(sp_edge(:,1)+1,1)];

    for ii = 1:length(ind_edge),
        if ii<length(ind_edge),
            subset = find((spike_stepLocations>ind_edge(ii,1))&...
                (spike_stepLocations<ind_edge(ii+1,1)));
        else
            subset = find(spike_stepLocations>ind_edge(ii,1));
        end
        total_subset = [ind_edge(ii,1);spike_stepLocations(subset,1)];
        [~,sub_loc] = max(spike_trace(total_subset(1:end),1));
        spike_locations(ii,1) = ind_edge(ii,1)+(sub_loc-1);
    end

    if display_spikes == 1,
        if ishandle('f'),
            gcf;
            scatter(spike_locations(:,1)./sampling_rate,...
                zeros(length(spike_locations),1),'r');
            title('Scatterpoints show spike peak raster locations');
            if max(spike_trace) < 0,
                ylim([-100 10]);
            else
                ylim([-100 max(spike_trace)+0.1*abs(max(spike_trace))]);
            end
        else
            plot(t_vec,spike_trace,'k');
            hold on;
            scatter(spike_locations(:,1)./sampling_rate,...
                zeros(length(spike_locations),1),'r');
            title('Scatterpoints show spike peak raster locations');
            if max(spike_trace) < 0,
                ylim([-100 10]);
            else
                ylim([-100 (max(spike_trace)+0.1*abs(max(spike_trace)))]);
            end
        end
    else
    end

end

