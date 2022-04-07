clear;
close all;
load('Sims2020_12020-06-13-T135211.mat');

DataAccuracy
DataLatency

NCond=1;

for ct=1:NCond
    for age=1:3
        rmseA(ct,age)=sqrt(mean((data(:,age)-result(ct).latenciesCorrect(:,age)).^2,1));
    end
    for exp=1:2
        if exp == 1
            sum_of_squares_rmseE_IOWA=sum(sum((dataIOWA-result(ct).latenciesCorrect(1:5,:)).^2));
            rmseE(ct,exp)=sqrt(sum_of_squares_rmseE_IOWA/15);
            %rmseE(ct,exp)=sqrt(mean((dataIOWA(1:5,age)-result(ct).latenciesCorrect(1:5,age)).^2,1));
        else
            sum_of_squares_rmseE_IOWAC=sum(sum((dataIOWAC-result(ct).latenciesCorrect(6:10,:)).^2));
            rmseE(ct,exp)=sqrt(sum_of_squares_rmseE_IOWAC/15);
            %rmseE(ct,exp)=sqrt(mean((dataIOWAC(1:5,age)-result(ct).latenciesCorrect(6:10,age)).^2,1));
        end
    end
    sum_of_squares=sum(sum((data-result(ct).latenciesCorrect).^2));
    rmseO(ct)=sqrt(sum_of_squares/30);
    %rmseO(ct)=sqrt(mean(mean((data(:,:)-result(ct).latenciesCorrect(:,:)).^2,1)),2);
end

for age=1:3
    %[rmseA_best(age),rmseA_expnum(age)] = min(rmseA(:,age),1);
    %[rmse_best, rmse_expnum] = min(rmseA(:,age),1);
    [rmseA_best_5] = min(rmseA(:,1));
    [rmseA_index_5] = find(any(rmseA(:,1) == rmseA_best_5,2));
    [rmseA_best_7] = min(rmseA(:,2));
    [rmseA_index_7] = find(any(rmseA(:,2) == rmseA_best_7,2));
    [rmseA_best_10] = min(rmseA(:,3));
    [rmseA_index_10] = find(any(rmseA(:,3) == rmseA_best_10,2));
   
end
for exp=1
    %[rmseE_best(exp),rmseE_expnum(exp)] = min(rmseE(:,exp),1);
    [rmseE_best_1] = min(rmseE(:,exp));
    [rmseE_index_Exp1] = find(any(rmseE(:,1) == rmseE_best_1,2));
end
for exp=2
    %[rmseE_best(exp),rmseE_expnum(exp)] = min(rmseE(:,exp),1);
    [rmseE_best_2] = min(rmseE(:,exp));
    [rmseE_index_Exp2] = find(any(rmseE(:,2) == rmseE_best_2,2));
end
rmseO = rmseO.';
[rmseO_best] = min(rmseO)
[rmse_index_O] = find(any(rmseO(:,1) == rmseO_best,2))

for ct=1:NCond
    for age=1:3
        rmseAccA(ct,age)=sqrt(mean((dataAcc(:,age)-result(ct).ratesCorrect(:,age)).^2,1));
    end
    for exp=1:2
        if exp == 1
            sum_of_squares_rmseEAcc=sum(sum((dataAccIOWA-result(ct).ratesCorrect(1:5,:)).^2));
            rmseAccE(ct,exp)=sqrt(sum_of_squares_rmseEAcc/15);
            %rmseAccE(ct,exp)=sqrt(mean((dataAccIOWA(1:5,age)-result(ct).ratesCorrect(1:5,age)).^2,1));
        else
            sum_of_squares_rmseEAcc=sum(sum((dataAccIOWAC-result(ct).ratesCorrect(6:10,:)).^2));
            rmseAccE(ct,exp)=sqrt(sum_of_squares_rmseEAcc/15);
            %rmseAccE(ct,exp)=sqrt(mean((dataAccIOWAC(1:5,age)-result(ct).ratesCorrect(6:10,age)).^2,1));
        end
    end
    sum_of_squares_Acc=sum(sum((dataAcc-result(ct).ratesCorrect).^2));
    rmseAccO(ct)=sqrt(sum_of_squares_Acc/30);
    %rmseAccO(ct)=sqrt(mean((dataAcc(:,age)-result(ct).ratesCorrect(:,age)).^2,1));
    %rmseO(ct)=sqrt(mean(mean((data(:,:)-result(ct).latenciesCorrect(:,:)).^2,1)),2);
end

for age=1:3
    %[rmseA_best(age),rmseA_expnum(age)] = min(rmseA(:,age),1);
    %[rmse_best, rmse_expnum] = min(rmseA(:,age),1);
    [rmseAccA_best_5] = min(rmseAccA(:,1));
    [rmseAccA_index_5] = find(any(rmseAccA(:,1) == rmseAccA_best_5,2));
    [rmseAccA_best_7] = min(rmseAccA(:,2));
    [rmseAccA_index_7] = find(any(rmseAccA(:,2) == rmseAccA_best_7,2));
    [rmseAccA_best_10] = min(rmseAccA(:,3));
    [rmseAccA_index_10] = find(any(rmseAccA(:,3) == rmseAccA_best_10,2));
   
end
for exp=1
    %[rmseE_best(exp),rmseE_expnum(exp)] = min(rmseE(:,exp),1);
    [rmseAccE_best_1] = min(rmseAccE(:,exp));
    [rmseAccE_index_Exp1] = find(any(rmseAccE(:,1) == rmseAccE_best_1,2));
end
for exp=2
    %[rmseE_best(exp),rmseE_expnum(exp)] = min(rmseE(:,exp),1);
    [rmseAccE_best_2] = min(rmseAccE(:,exp));
    [rmseAccE_index_Exp2] = find(any(rmseAccE(:,2) == rmseAccE_best_2,2));
end
rmseAccO = rmseAccO.';
[rmseAccO_best] = min(rmseAccO)
[rmseAcc_index_O] = find(any(rmseAccO(:,1) == rmseAccO_best,2))

bestRTrmseACC = rmseAccO(rmse_index_O)
bestACCrmseRT = rmseO(rmseAcc_index_O)

figure;
plot(rmseAccO);

figure;
plot(rmseO);

%%%% plot best results
Best=1; %rmseAcc_index_O;
ttext={'IOWA Acc: 5mo' 'IOWA Acc: 7mo' 'IOWA Acc: 10mo'}; 
ttextC={'IOWAc Acc: 5mo' 'IOWAc Acc: 7mo' 'IOWAc Acc: 10mo'}; 
ttextRT={'IOWA RT: 5mo' 'IOWA RT: 7mo' 'IOWA RT: 10mo'}; 
ttextRTC={'IOWAc RT: 5mo' 'IOWAc RT: 7mo' 'IOWAc RT: 10mo'}; 
subs=[1 3 5];
subsC=[2 4 6];

rmseAccO(Best)
rmseO(Best)

figure;
for act=1:3
    subplot(3,2,subs(act))
    plot(dataAccIOWA(:,act));
    title(ttext(act));
    hold
    plot(result(Best).ratesCorrect(1:5,act),'--');
end
for act=1:3
    subplot(3,2,subsC(act))
    plot(dataAccIOWAC(:,act));
    title(ttextC(act));
    hold
    plot(result(Best).ratesCorrect(6:10,act),'--');
end

figure;
for act=1:3
    subplot(3,2,subs(act))
    plot(dataIOWA(:,act));
    title(ttextRT(act));
    hold
    plot(result(Best).latenciesCorrect(1:5,act),'--');
end
for act=1:3
    subplot(3,2,subsC(act))
    plot(dataIOWAC(:,act));
    title(ttextRTC(act));
    hold
    plot(result(Best).latenciesCorrect(6:10,act),'--');
end