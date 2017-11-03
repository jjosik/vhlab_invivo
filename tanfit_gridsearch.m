function [ select_s_x,select_trace,opt_params,rmse_star,si ] = ...
tanfit_gridsearch( Vm,FR_inst,reinit_N,use_random,reps,use_trial_baseline )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a_range = [-200 200];
b_range = [0 200];
c_range = [-120 120];
d_range = [0 100];

init_array = cell(reps,1);
sx_array = cell(reps,1);
trace_array = cell(reps,1);
gof_array = zeros(reps,1);
params_array = cell(reps,1);
for i = 1:reps,
    if use_random,
        a_init = randi(a_range,1,1);
        b_init = randi(b_range,1,1);
        c_init = randi(c_range,1,1);
        d_init = (diff(d_range).*rand(1,1))-abs(min(d_range));
        a_inter = round(abs(max(a_range)-min(a_range))/reinit_N,2);
        b_inter = round(abs(max(b_range)-min(b_range))/reinit_N,2);
        c_inter = round(abs(max(c_range)-min(c_range))/reinit_N,2);
        d_inter = round(abs(max(d_range)-min(d_range))/reinit_N,2);
        a_array = min(a_range):a_inter:max(a_range);
        b_array = min(b_range):b_inter:max(b_range);
        c_array = min(c_range):c_inter:max(c_range);
        d_array = min(d_range):d_inter:max(d_range);
        %note this method is incomplete: set use_random = 1 with high N reps;
    end
    init_array{i,1} = {a_init,b_init,c_init,d_init};
    [s_x,trace,params,gof,fitinfo] = alt_tanhfit(Vm,FR_inst,use_trial_baseline,init_array);
    sx_array{i,1} = s_x;
    trace_array{i,1} = trace;
    gof_array(i,1) = gof.rmse;
    params_array{i,1} = params;
end

[rmse_star,si] = sort(gof_array);
select_s_x = sx_array{si,1};
select_trace = trace_array{si,1};
opt_params = init_array{si,1};

    
    
end

