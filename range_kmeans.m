function [aveAcc,aveMIhat,aveNMI,avePurity,aveARI,W,b, Term_one,Fterm_21,Sterm_21]= range_kmeans(fileName,para_a,para_b)
% example:   clear;clc;[Acc,NMI,MIhat,Purity, ARI]= KMeans_v1(1);

%FileName(31,:)="";


    switch fileName
        case 1
            fileName = 'balance_uni';
            
        case 2
            fileName = 'banknote_authentication';
            
        case 3
            fileName = 'binalpha_uni';
            
        case 4
            fileName = 'Breast_Cancer_Wisconsin_Original';
            
        case 5
            fileName = 'Cardiotocography';
            
        case 6
            fileName = 'chess_uni';
            
        case 7
            fileName = 'Cora_ML_uni';
            
        case 8
            fileName = 'crx_uni';
            
        case 9
            fileName = 'Diabetic_Retinopathy_Debrecen';
            
        case 10
            fileName = 'dig1-10_uni';
            
        case 11
            fileName = 'german_uni';
            
        case 12
            fileName = 'letter_uni';
            
        case 13
            
            fileName = 'Wine_Quality_Red';
        case 14
            fileName = 'MSRA25_uni';
            
        case 15
            fileName = 'Musk_v2';
            
        case 16
            fileName = 'Parkinson_Speech_Dataset';
            
        case 17
            fileName = 'pima_uni';
            
        case 18
            fileName = 'Wireless_Indoor_Localization';
            
        case 19
            fileName = 'segment_uni';
            
        case 20
            fileName = 'solar_uni';
            
            
        case 21
            fileName = 'Spambase';
            
        case 22
            fileName = 'Yale_32x32';
        case 23
            
            fileName = 'yeast_uni';
        case 24
            fileName = 'USPSdata_20_uni';
            
        case 25
            fileName = 'uspst_uni';
            
        case 26
            fileName = 'vehicle_uni';
            
        case 27
            fileName = 'Waveform';
            
        case 28
            fileName = 'waveform-21_uni';
            
        case 29
            fileName = 'WebKB_wisconsin_uni';
            
        case 30
            fileName = 'Website_Phishing';
            
        case 31
            fileName = 'wine';
                 case 32
            fileName = 'ecoli_processed.mat';
        case 33
            fileName = 'glass_processed.mat';
        case 34
            fileName = 'iris_processed.mat';
        case 35
            fileName='yeast_processed.mat';
            
            
    end
    % which(fileName)
    load (fileName);
    c=length(unique(Y));
    X=Data';
    
    [pre_data,W,b,Term_one,Fterm_21,Sterm_21 ] = joint_clustering(X,para_a,para_b,c);
   
    
    
    [Acc,NMI,MIhat,Purity, ARI,aveAcc,aveMIhat,aveNMI,avePurity,aveARI] = Evaluate(pre_data,Y);
    
    %res{:} =res;  
 %   MyFilename=sprintf('%d.%s', Filename, fileName);
 %   save(['C:\Users\tliu\Desktop\PHD_2018\Code_Tong\Results\KMeans\' MyFilename,'kMeans_Results'],'Acc','NMI','MIhat','Purity','ARI','aveAcc','aveNMI', 'aveMIhat','avePurity','aveARI');
 save([fileName,'kMeans_Results'],'Acc','NMI','MIhat','Purity','ARI','aveAcc','aveNMI', 'aveMIhat','avePurity','aveARI','W');
end

function [Acc,NMI,MIhat,Purity,ARI,aveAcc,aveMIhat,aveNMI,avePurity,aveARI, res] = Evaluate(X,Y)
gnd= Y;


zero_index=find(gnd==0,1);

if ~isempty(zero_index)
gnd=gnd+1;
end

k= length(unique(gnd));

for iter = 1:20
b= kmeans(X,k);

res    = bestMap(gnd,b);
Acc(iter) = length(find(gnd == res))/length(gnd);
MIhat(iter)    = MutualInfo(gnd,res);
NMI(iter)   = nmi(gnd,res);
Purity(iter)   = Calculate_purity(gnd,res);
ARI(iter)   = clustereval(gnd,res, 'ari') ;
clear b res
end

aveAcc = mean(Acc(:));
aveMIhat = mean(MIhat(:));
aveNMI = mean(NMI(:));
avePurity = mean(Purity(:));
aveARI = mean(ARI(:));

%aveACC(i) = mean(Acc(:));
%aveMIhat(i) = mean(MIhat(:));
%aveNMI(i) = mean(NMI(:));
%avePurity(i) = mean(Purity(:));
%aveARI(i) = mean(ARI(:));


end