function [ b ] = update_b( e,Q,Y,X,W )
%UPDATE_B Summary of this function goes here
%   Detailed explanation goes here

b=inv((e*Q*e'))*(e*Q*(Y-X'*W));

end

