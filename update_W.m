function [W] = update_W(Y, X,H,para_beta, para_alpha,D,Q)
%UPDATE_W Summary of this function goes here
%   Detailed explanation goes here

W=inv(X*H'*Q*H*X'+para_beta/para_alpha*D)*X*H'*Q*Y

end

