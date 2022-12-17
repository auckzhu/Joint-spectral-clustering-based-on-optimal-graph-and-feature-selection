function [ Y ] = update_Y(M)
%UPDATE_Y Summary of this function goes here
%   Detailed explanation goes here


[U,S,V]=svd(M);

Y=UV';


end

