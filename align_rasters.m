function [ stim_rasters ] = align_rasters( spike_locations,...
    stimvalues,stimorder,nStim_ON,nStim_OFF,sampling_rate )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

stim_rasters = cell(length(stimvalues),1);
    for b = 1:length(stimvalues),
        current_epochs = find(stimorder == b);
        align_on_ind = nStim_ON(current_epochs).*sampling_rate;
        align_off_ind = nStim_OFF(current_epochs).*sampling_rate;
        for a=1:length(current_epochs),
            rasters{a,1} = spike_locations(find(spike_locations>align_on_ind(a,1)&...
                spike_locations<align_off_ind(a,1)))-...
                (nStim_ON(current_epochs(1,a))*sampling_rate);
        end
        new_rasters = [];
        for a=1:length(current_epochs),
            new_rasters = [new_rasters;rasters{a,1}];
        end
        stim_rasters{b,1} = new_rasters;
    end

end

