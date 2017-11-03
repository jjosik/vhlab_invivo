function [ trace,params,gof,fitinfo ] = tanhfit( x_,y_ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x_ = reshape(x_,length(x_),1);
y_ = reshape(y_,length(y_),1);
s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-Inf,0,-100,1e-12],...
    'Upper',[Inf,Inf,100,1e12],...
    'Startpoint',[0,1,0,1]);
ft = fittype('a + b.*tanh((x - c)./d)','options',s);
[params,gof,fitinfo] = fit(x_,y_,ft);
s_x = [floor(min(x_)):1:floor(max(x_))];
trace = params.a + params.b.*tanh((s_x - params.c)./params.d);

end

