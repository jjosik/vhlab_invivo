function [ I_ex,J_ex,radius ] = exclude_peaks( resampled_kernel,X,Y,Mx,My )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%gtest = gradient(resampled_kernel);
%original_sign = gtest(M_sx(1),M_ty(1))/abs(gtest(M_sx(1),M_ty(1)));
%%negative y
%loc_y_minus = M_ty(1);
%current_sign = original_sign;
%while current_sign == original_sign,
%    loc_y_minus = loc_y_minus-1;
%    current_sign = gtest(M_sx(1),loc_y_minus)/abs(gtest(M_sx(1),loc_y_minus));
%end
%%positive y
%loc_y_plus = M_ty(1);
%current_sign = original_sign;
%while current_sign == original_sign,
%    loc_y_plus = loc_y_plus+1;
%    current_sign = gtest(M_sx(1),loc_y_plus)/abs(gtest(M_sx(1),loc_y_plus));
%end
%%negative x
%loc_x_minus = M_sx(1);
%current_sign = original_sign;
%while current_sign == original_sign,
%    loc_x_minus = loc_x_minus-1;
%    current_sign = gtest(loc_x_minus,M_ty(1))/abs(gtest(loc_x_minus,M_ty(1)));
%end
%%positive x
%loc_x_plus = M_sx(1);
%current_sign = original_sign;
%while current_sign == original_sign,
%    loc_x_plus = loc_x_plus+1;
%    current_sign = gtest(loc_x_plus,M_ty(1))/abs(gtest(loc_x_plus,M_ty(1)));
%end

%ALT. VERSION - CIRCULAR ROI
gtest = gradient(resampled_kernel);
original_sign = gtest(Mx,My)/abs(gtest(Mx,My));
radius = 2;
current_sign = original_sign;
while ~(any(current_sign == -(original_sign))),
    radius = radius + 1;
    circle = (X-Mx).^2+(Y-My).^2 < radius;
    current_sign = gtest(circle)./abs(gtest(circle));
end
ex_coords = find(circle);
[I_ex,J_ex] = ind2sub(size(resampled_kernel),ex_coords);



end

