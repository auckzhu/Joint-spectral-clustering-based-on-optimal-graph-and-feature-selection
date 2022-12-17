function [ M ] = update_MH(L_a, Y, Q, W , X, paraH_miu)
%UPDATE_MH Summary of this function goes here
%   Detailed explanation goes here

M=L_a*Y+paraH_miu*H'*Q*H*X'*W+H'*Q*H*Y;

end

