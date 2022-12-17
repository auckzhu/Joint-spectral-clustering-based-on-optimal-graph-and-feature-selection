function [ M ] = update_M(L_a, Y , G, W,X,b ,para_miu)
%UPDATE_M Summary of this function goes here
%   Detailed explanation goes here

M=L_a*Y+para_miu*G*Y+para_miu*G*X'*W+para_miu*G*b;



end

