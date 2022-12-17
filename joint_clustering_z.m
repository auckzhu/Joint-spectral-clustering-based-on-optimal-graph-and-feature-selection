function [Y,W,b] = joint_clustering(X,para_a,para_b)
%JO Summary of this function goes here
%   Detailed explanation goes here

% Objective function
%                 min_{Y，W，b}   tr(Y^T*L_a*Y) + para.a*||X^T*W+e^T*b-Y||_2_1
%                                + para.b * ||W||_2_1
%                   
% solution:
%             Y = tr(Y'*L*Y+para.a*Y'*H'*Q*H*Y-2*Y'*H'*Q*H*X'*W)，
%           min {Y,Y'*Y=I} L*Y+para.a*H'*G*H*X'*W+H'*G*H*Y
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
% L: n * n
% X: c * n
% W: c * n
% H: n * n
% e: 1 * n
% b: 1 * c
% D: c * c
% G: c * c
% Q: n * n
%
% clear;clc;load psMCIm
%
%
%


%% =================== initialize ======================
    [c,n] = size(X);

    
    L_a = rand(n);
    
    L = triu(L_a,0)+tril(L_a',-1);
    A = L;
  
    Y=rand(n,c);

    b=rand(1,c);

    e=ones(1,n);
    
    W=rand(c,c);
    
    I=eye(n);

    G_temp = sqrt(sum((X'*W+e'*b-Y).*(X'*W+e'*b-Y),2)+eps);
    G_d2 = 0.5./(G_temp);
    G = diag(G_d2);
    
    H = I-inv(e*G*e')*e'*e*G;
 
    D_temp = sqrt(sum(W.*W,2)+eps);
    D_d2 = 0.5./(D_temp);
    D = diag(D_d2);
    
    Q_temp = sqrt(sum((H*X'*W-H*Y).*(H*X'*W-H*Y),2)+eps);
    Q_d2 = 0.5./(Q_temp);
    Q = diag(Q_d2);
    
    D=eye(c);

    Q=eye(n);

    W = inv(X*H'*Q*H*X'+para_b/para_a*D)*X*H'*Q*Y;





%% ==================== Q_2_1=1/2*||H*X'*W-H*Y||_2) matrix ===================
function Q = Q_2_1(H,X,W,Y)
    Q_temp = sqrt(sum((H*X'*W-H*Y).*(H*X'*W-H*Y),2)+eps);
    Q_d2 = 0.5./(Q_temp);
    Q = diag(Q_d2);
end
    
%     add_column=[];
%     temp_Q=H*X'*W-H*Y; 
% for i=size(temp_Q,2)
%     each_column=1/2*norm(temp_Q(:,i));
%     add_column=[add_column each_column];
% end
%     Q=diag(add_column);

%% ================== G_2_1=1/2*||X'*W+e'*b-Y||_2 matrix =========================
function G = G_2_1(X,W,e,b,Y)
    G_temp = sqrt(sum((X'*W+e'*b-Y).*(X'*W+e'*b-Y),2)+eps);
    G_d2 = 0.5./(G_temp);
    G = diag(G_d2);
end
%     add2_column=[];
%     temp_G=X'*W+e'*b-Y; 
% for i=size(temp_G,2)
%     each_column=1/2*norm(temp_D(:,i));
%     add2_column=[add2_column each_column];
% end
%     G=diag(add2_column);
%% =================== D_2_1=1/2*||W^i||_2  matrix ==========================

function D = D_2_1(W)
    D_temp = sqrt(sum(W.*W,2)+eps);
    D_d2 = 0.5./(D_temp);
    D = diag(D_d2);
end

%     add3_column=[];
%     temp_D=W;
% for i=size(temp_D,2)
%     each_column = 1/2*norm(temp_D(:,i));
%     add3_column = [add3_column each_column];
% end
%     D=diag(add3_column);
%% =========iterate============
iter = 1;
obji = 1;
while 1 

    Q = Q_2_1(H,X,W,Y);

    G = G_2_1(X,W,e,b,Y);

    D = D_2_1(W);
    
    b=inv((e*G*e'))*(e*G*(Y-X'*W));
    
    W=inv(X*H'*Q*H*X'+para_b/para_a*D)*X*H'*Q*Y;
   
    % min_{Y'Y=I} Tr(Y'LY-2Y'(H'*G*H*X'*W+H'*G*H*Y))
    B = H'*G*H*X'*W+H'*G*H*Y;
    Y = GPI(A,B);
    
    Fterm_i2 = sqrt(sum((X'*W+e'*b-Y).*(X'*W+e'*b-Y))+eps);
    Fterm_21 = sum(Fterm_i2);
     
    Sterm_i2 = sqrt(sum((W.*W))+eps);
    Sterm_21 = sum(Sterm_i2);
    
    obj(iter) = trace(Y'*L*Y)+para_a*Fterm_21+para_b*Sterm_21;
    cver = abs((obj(iter)-obji)/obji);
    obji = obj(iter); 
    iter = iter+1;
%     if (cver < 10^-9 && iter > 2) || iter == 20;
%     cver = abs((obj(iter)-obji)/obji);
%     obji = obj(iter); 
%     iter = iter+1;
%     if (cver < 10^-9 && iter > 2) || iter == 70;
end
 plot(obj);
end



