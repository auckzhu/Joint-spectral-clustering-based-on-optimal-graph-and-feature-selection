    clc,clear
    X=rand(80,150);
    [c,n] = size(X);
    
    L_a = rand(n);    
    L = triu(L_a,0)+tril(L_a',-1);
    
    para_a=n/c;
    para_b=n/c;
    
   save joint_example L X para_a para_b
    
    