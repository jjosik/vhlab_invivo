function [ p_bin_sub,n_bin_sub ] = ON_OFF_pixelsort( kernel )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%negative pixels
working_zero = mean(kernel(:));
n_field = find(kernel(:)<working_zero);
size_n = round(length(n_field)*0.8);
[temp_n_set,temp_nloc] = sort(abs(kernel(n_field)),'descend');
n_threshold = -(temp_n_set(size_n));
sub_n_field = find(kernel(:)<n_threshold);
n_bin_sub = zeros(size(kernel));
n_bin_sub(sub_n_field) = 1;


%positive pixels
p_field = find(kernel(:)>=working_zero);
size_p = round(length(p_field)*0.8);
[temp_p_set,temp_ploc] = sort(kernel(p_field),'descend');
p_threshold = temp_p_set(size_p);
sub_p_field = find(kernel(:)>p_threshold);
p_bin_sub = zeros(size(kernel));
p_bin_sub(sub_p_field) = 1;



end

