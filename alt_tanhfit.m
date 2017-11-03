function [ s_x,trace,params,gof,fitinfo ] = alt_tanhfit( Vm,FR_inst,use_trial_baseline,varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin > 3,
    a_init = cell2mat(varargin{1,1}{1,1}(1,1));
    b_init = cell2mat(varargin{1,1}{1,1}(1,2));
    c_init = cell2mat(varargin{1,1}{1,1}(1,3));
    d_init = cell2mat(varargin{1,1}{1,1}(1,4));
else
    a_init = max(FR_inst)/2;
    b_init = 20;
    c_init = x_inflect;
    d_init = 1;
end

if use_trial_baseline == 1,
    x_inflect = 15;
else
    x_inflect = -50;
end

FR_inst = reshape(FR_inst,length(FR_inst),1);
Vm = reshape(Vm,length(Vm),1);
s = fitoptions('Method','NonlinearLeastSquare',...
    'Lower',[-Inf,0,-120,1e-12],...
    'Upper',[Inf,Inf,120,1e12],...
    'Startpoint',[a_init,b_init,c_init,d_init]);
ft = fittype('a + b.*tanh((x - c)./d)','options',s);
[params,gof,fitinfo] = fit(Vm,FR_inst,ft);
s_x = [floor(min(Vm))+0.5:1:floor(max(Vm))+0.5];
trace = params.a + params.b.*tanh((s_x - params.c)./params.d);

end

