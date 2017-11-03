function [ output ] = mvn2d_objfun( params,x,y,kernel,s_kernel,presort)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
global scale;
try
    ccov = [params(5) 0; params(6) params(7)];  %lower left Cholesky factor of input cov. matrix
    %cinv = inv(ccov);
catch
    ccov = [5 0; 5 5];       %backup lower left Cholesky "guess"
    %cinv = inv(ccov);
end
temp = [x(:) - params(2),y(:) - params(3)]*ccov;
pred = params(1) + params(4)*exp(-sum(temp.*temp,2)/2); 
if presort == 1,
    output = pred - s_kernel(:);
else
    output = pred - kernel(:);
end
    
    

end

