%% sample code for Low-rank Multi-modality (LRMM) feature selection 
% Objective function
%                    |Y1 - W1'X1|_F^1 + |Y2 - W2X2|_F^2 + alpha * |W1|_2p +   beta * |W2|_2p + gamma * |[W1 W1]|_*
% solution:
%             W1 = (X1'X1 + alpha*D1 + gamma*D)\X'Y, where D1_ii =(W1)^i; D_ii = 0.5 * (WW')^(-0.5), W = [W1 W2];
%             W2 = (X2'X2 + alpha*D2 + gamma*D)\X'Y, where D2_ii =(W2)^i; D_ii = 0.5 * (WW')^(-0.5), W = [W1 W2];
% X1,X2: ins * fea
% Y1,Y2: ins * class
% W1,W2: fea * class


% clear;clc;load psMCIm
% X1 = Data(:,94:end);Y1 = Y; X2 = Data(:,1:93);Y2 = Y;
% para.alpha = 1; para.beta = 2;para.gamma = 2;
% [W1,W2,obj] = L21TraceFS(X1,X2,Y1,Y2,para);
%%
function [W1,W2,obj] = L21TraceFS(X1,X2,Y1,Y2,para)
if isfield(para, 'alpha')
    alpha = para.alpha;
else
    alpha = 1;
end
if isfield(para, 'beta')
    beta = para.beta;
else
    beta = 1;
end
if isfield(para, 'gamma')
    gamma = para.gamma;
else
    gamma = 1;
end
if isfield(para, 'p')
    p = para.p;
else
    p = 0.5;
end
if isfield(para, 'flag')
    flag = para.flag;
else
    flag = 0;
end

[n1,d1] = size(X1);
[n2,d2] = size(X2);
c1 = size(Y1,2);
c2 = size(Y2,2);

I1 = eye(n1);
I2 = eye(n2);
e1 = ones(n1,1);
e2 = ones(n2,1);
H1 = I1-e1*e1'/n1;
H2 = I2-e2*e2'/n2;

W1 = rand(d1,c1);
W2 = rand(d2,c2);
W = [W1,W2];

X1 = H1 * X1;
X2 = H2 * X2;
Y1 = H1 * Y1;

X1tX1 = X1'*X1;
X2tX2 = X2'*X2;
X1tY1 = X1'*Y1;
X2tY2 = X2'*Y2;

iter = 1;
obji = 1;
while 1 
    D = 0.5*((W*W'+ eps)^(-0.5)); %appearing complex number
    
    dd1 = (p/2)./(sqrt(sum(W1.*W1,2)+eps).^(2-p));
    DD1 = diag(dd1);        
    W1 = (X1tX1 + alpha * DD1 + beta* D)\X1tY1;
    b1 = (e1'*Y1 - e1'*X1*W1)/n1; 
    
    dd2 = (p/2)./(sqrt(sum(W2.*W2,2)+eps).^(2-p));
    DD2 = diag(dd2);        
    W2 = (X2tX2 + alpha * DD2 + beta* D)\X2tY2;
    b2 = (e2'*Y2 - e2'*X2*W2)/n2;

    W = [W1,W2];

    obj(iter) = sum(sqrt(sum((X1*W1 + e1*b1 - Y1).*(X1*W1 + e1*b1 - Y1),2)+eps)) ...
              + sum(sqrt(sum((X2*W2 + e2*b2 - Y2).*(X2*W2 + e2*b2 - Y2),2)+eps)) ...
              + alpha * sum(sqrt(sum(W1.*W1,2)+eps).^(2-p))...
              + beta * sum(sqrt(sum(W2.*W2,2)+eps).^(2-p))...
              + gamma * sum(svd(W));
    cver = abs((obj(iter)-obji)/obji);
    obji = obj(iter); 
    iter = iter+1;
    if (cver < 10^-9 && iter > 2) || iter == 20;
    break, end
end
if flag == 1
    plot(obj)
else
    
end












