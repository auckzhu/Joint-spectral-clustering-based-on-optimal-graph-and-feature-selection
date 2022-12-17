function [Y,W,b,Term_one,Fterm_21,Sterm_21 ] = joint_clustering(X,para_a,para_b,c,L)
%JO Summary of this function goes here
%   Detailed explanation goes here

% Objective function
%                 min_{Y，W，b}   tr(Y^T*L_a*Y) + para.a*||X^T*W+e^T*b-Y||_2_1
%                                + para.b * ||W||_2_1
%                   
% solution:
%             Y = tr(Y'*L*Y+para.a*Y'*H'*Q*H*Y-2*Y'*H'*Q*H*X'*W)，
%             min {Y,Y'*Y=I} L*Y+para.a*H'*G*H*X'*W+H'*G*H*Y
%            
%             where Q = 0.5*||(H*X^T*W-H*Y)^i||_2, Y'*Y=I,
%                   H = I-inv(e*G*e)*e'*e*G,
%                   G_ii = 0.5*||(X'*W+e'*b-Y)^i||_2;
%
%             W = inv(X*H'*Q*H*X'+para.b/para.a*D*X*H'*D*Y, 
%             where Q_ii = 0.5*||(H*X^T*W-H*Y)^i||_2;
%                   D_ii = 0.5*||W^i||_2;
%
%             b = inv(e*G*e^T)(e*G*(Y-X'*W)), 
%             where Q_ii = 0.5*||(H*X'*W-H*Y)^i||_2;
% Y: n * c
% L: n * n symmetric matrix
% X: d * n
% W: d * c
% H: n * n
% e: 1 * n
% b: 1 * c
% D: c * c
% G: c * c
% Q: n * n
%
% input: clear;clc;joint_example
%
%
%  example
%  clear;clc; [Y,W,b] = joint_clustering(rand(50,500),1,3,20);

%
%% =================== initialize ======================

if ~exist('L','var')
    options = [];
    options.NeighborMode = 'KNN';
    options.k = 7;
    options.WeightMode = 'HeatKernel';
    options.t = 4;
    tempW = constructW(X',options);
    tempD = diag(sum(tempW));
    L = tempD - tempW;
else
   ;
end
% 
% c=length(unique(true_label));

[d,n] = size(X);

b = rand(1,c);
e = ones(1,n);

W      = rand(d,c);
D_temp = sqrt(sum(W.*W,2)+eps);
D_d2   = 0.5./(D_temp);
D      = diag(D_d2);

I = eye(n);
Y = rand(n,c);
H = rand(n,n);

Q = Q_2_1(H,X,W,Y);
H = I-inv(e*Q*e')*e'*e*Q;

% W = inv(X*H'*Q*H*X'+para_b/para_a*D)*X*H'*Q*Y;
    


%% =========iterate============

iter = 1;
obji = 1;
while 1
   
    % update Q
    HXtWHY = H*X'*W-H*Y;
    Q_temp = sqrt(sum(HXtWHY.*HXtWHY,2)+eps);
    Q_d2   = 0.5./(Q_temp);
    Q      = diag(Q_d2);  
    
    % update G
    XtWEtBY= X'*W+e'*b-Y;
    G_temp = sqrt(sum(XtWEtBY.*XtWEtBY,2)+eps);
    G_d2 = 0.5./(G_temp);
    G = diag(G_d2);
    b      = inv((e*G*e'))*(e*G*(Y-X'*W));    
    H      = I-inv(e*G*e')*e'*e*G;
    
    
    % update W
    W      = inv(X*H'*Q*H*X'+para_b/para_a*D)*X*H'*Q*Y;
    D_temp = sqrt(sum(W.*W,2)+eps);
    D_d2   = 0.5./(D_temp);
    D      = diag(D_d2);
  
    % update Y
    % min_{Y'Y=I} Tr(Y'LY-2*para_a*Y'(H'*Q*H*X'*W-H'*Q*H*Y))
    B = para_a*H'*Q*H*X'*W-H'*Q*H*Y;
    Y = GPI(L,B,1);
    
    %obj
    Term_one=trace(Y'*L*Y);
    
    Fterm_i2 = sqrt(sum((X'*W+e'*b-Y).*(X'*W+e'*b-Y))+eps);
    Fterm_21 = sum(Fterm_i2);
     
    Sterm_i2 = sqrt(sum((W.*W))+eps);
    Sterm_21 = sum(Sterm_i2);
    
    
    obj(iter) = Term_one + para_a*Fterm_21 + para_b*Sterm_21;
    cver = abs((obj(iter)-obji)/obji);
    obji = obj(iter); 
    iter = iter+1;
    if (cver < 10^-9 && iter > 2) || iter == 100;
    break, end
end
 plot(obj(2:end))
end


%% ==================== Q_2_1=1/2*||H*X'*W-H*Y||_2) matrix ===================
function Q = Q_2_1(H,X,W,Y)
    Q_temp = sqrt(sum((H*X'*W-H*Y).*(H*X'*W-H*Y),2)+eps);
    Q_d2 = 0.5./(Q_temp);
    Q = diag(Q_d2);
end
    

%% ================== G_2_1=1/2*||X'*W+e'*b-Y||_2 matrix =========================
function G = G_2_1(X,W,e,b,Y)
    G_temp = sqrt(sum((X'*W+e'*b-Y).*(X'*W+e'*b-Y),2)+eps);
    G_d2 = 0.5./(G_temp);
    G = diag(G_d2);
end

%% =================== D_2_1=1/2*||W^i||_2  matrix ==========================

function D = D_2_1(W)
    D_temp = sqrt(sum(W.*W,2)+eps);
    D_d2 = 0.5./(D_temp);
    D = diag(D_d2);
end


