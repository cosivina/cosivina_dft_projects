%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis for Suanda_Mugwnya_Namy_2014
clear all; close all; % This script generates an output file
%Global Variables
scale_factor=8; nFeatures=2;
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);
%% Raw Data File Name 
simName = 'WPPR_3100_Suanda_Mugwanya_Namy_2014_results'
% loading the file
xsit_result = load ([simName '.mat']);

%Experiment Variables
new_sim_formatting = true; % new format sims with task_variables saved at simulation
if (new_sim_formatting == true)
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrials = xsit_result.train(1).maxTRt;
else  
    nObjects = 8; %%total number of objects
    nTrials=8;
end

%% TEST CONDITION ANALYSIS
numSubjects=size(xsit_result.test,1);
targLookTime=zeros(numSubjects,nTrials);
dstrLookTime=zeros(numSubjects,nTrials);
HCDtargTime=[];MCDtargTime=[];LCDtargTime=[];HCDdistTime=[];MCDdistTime=[];LCDdistTime=[];iH=0;iM=0;iL=0;winnersH=0;winnersM=0;winnersL=0;
for subject=1:numSubjects
    corResp=0;
    for trt=1:nTrials
        Look(1)= sum( xsit_result.test(subject).historyLt(trt,:));%left extreme
        Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,:));%left middle
        Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,:));%right middle
        Look(4)= sum( xsit_result.test(subject).historyRt(trt,:));%right extreme
        
        targetPos = xsit_result.test(subject).targetPosition(trt,:);
        trgL = find(targetPos == 1);%the spatial location of target 
        for i=1:4
            if i==trgL
                targLookTime(subject,trt)=targLookTime(subject,trt)+Look(i);
            else
                dstrLookTime(subject,trt)=dstrLookTime(subject,trt)+Look(i);%time to all 3 distracrors added up
            end
        end
       if  targLookTime(subject,trt) > (dstrLookTime(subject,trt))
           corResp=corResp+1;%number of trials answered/looked correctly
       end
    end
    subj=xsit_result.test(subject).subject; 
    if mod(subj,3)==1
    %% Condition HCD-High
        iH=iH+1;
        HCDtargTime(iH)=mean(targLookTime(subject,:));
        HCDdistTime(iH)=mean(dstrLookTime(subject,:)); 
        %if(HCDtargTime(iH) > HCDdistTime(iH)); winnersH=winnersH+1; end 
        HCDsuc(iH)=corResp/nTrials; if HCDsuc(iH)>0.5; winnersH=winnersH+1; end %proportion of test trials answered correctly
        %H_Cor(iH)=Correct_inTrace(subject);H_Wron(iH)=Wrong_inTrace(subject);
    elseif mod(subj,3)==2
    %% Condition MCD-Moderate
        iM=iM+1;
        MCDtargTime(iM)=mean(targLookTime(subject,:));
        MCDdistTime(iM)=mean(dstrLookTime(subject,:));  
        %if(MCDtargTime(iM) > MCDdistTime(iM)); winnersM=winnersM+1; end 
        MCDsuc(iM)=corResp/nTrials; if MCDsuc(iM)>0.5; winnersM=winnersM+1; end
        %M_Cor(iH)=Correct_inTrace(subject);M_Wron(iH)=Wrong_inTrace(subject);
    elseif mod(subj,3)==0
    %% Condition LCD-Low
        iL=iL+1;
        LCDtargTime(iL)=mean(targLookTime(subject,:));
        LCDdistTime(iL)=mean(dstrLookTime(subject,:)); 
        %if(LCDtargTime(iL) > LCDdistTime(iL)); winnersL=winnersL+1; end 
        LCDsuc(iL)=corResp/nTrials; if LCDsuc(iL)>0.5; winnersL=winnersL+1; end
        %L_Cor(iH)=Correct_inTrace(subject);L_Wron(iH)=Wrong_inTrace(subject);
    end
end

%% Data for Output File Saving 
measurement_i = 'Mean proportion correct looking at test';
empirical_mean_i = [0.48 0.39 0.34];
mean_i = [mean(HCDsuc) mean(MCDsuc) mean(LCDsuc)];
SE_i = [SE(HCDsuc) SE(MCDsuc) SE(LCDsuc)];
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%%
measurement_i = 'Proportion of subjects looking preferentially to target at test';
empirical_mean_i = [0.8 0.66 0.56];
mean_i = [winnersH./iH winnersM./iM winnersL./iL];
SE_i = [0.04 0.04 0.04];
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


figure(311)% Plot Mean proportion correct looking time
blockNames={'High'; 'Moderate'; 'Low'};
sts = [ mean(HCDsuc) 0.48 ;        mean(MCDsuc) 0.39 ;        mean(LCDsuc) 0.34  ];
errY =[ SE(HCDsuc) 0.21/sqrt(28) ; SE(MCDsuc) 0.20/sqrt(28) ; SE(LCDsuc) 0.18 /sqrt(28) ];%28*3 is number of child particpants
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',16)
%legend('Target','Distractor');
ylabel('Mean proportion correct');
legend('WOLVES Model','Suanda, Mugwanya & Namy, 2014');
ylim([0 1]);
grid on

%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])