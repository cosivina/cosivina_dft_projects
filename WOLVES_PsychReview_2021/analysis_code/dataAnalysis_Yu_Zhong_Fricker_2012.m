%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis for Yu Zhong & Fricker (2012)
clear all; close all; % This script generates an output file
%Global Variables
scale_factor=8;nFeatures=2; 

%Plotting Variables
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);
%%  

%% Raw Data File Name %wfChanges_conwmf0_1k15k_hwmc1_fix48_wrdlen2k_new_Yu_Zhong_Fricker_2012_results
simName='WPPR_fam7k_Yu_Zhong_Fricker_2012_results'% 
% loading the file
xsit_result = load( [simName,'.mat']);

numSubjects=size(xsit_result.test,1);

%Experiment (Task) Variables
new_sim_formatting = true; % new format sims with task_variables saved at simulation
if (new_sim_formatting == true)
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrials = xsit_result.train(1).maxTRt; 
    %word_On = xsit_result.test(1).word1_On; word_Off = xsit_result.test(1).word1_Off; 
    %vis_Off = xsit_result.test(1).visuals_Off;
    test_dur = 1000/scale_factor; % mean(word_Off - word_On);
else
    nObjects =  18; %%total number of objects
    nTrials = nObjects; 
    test_dur = 1000/scale_factor;
end
    

%% TEST CONDITION ANALYSIS    
empirical_data_A = [0.9187 0.5812]; %proportion correct word knowledge
lt_1_emp =  [0.41 0.38 0.50 0.59 0.65 0.71]; % proportion target looking strong
lt_0_emp =  [0.35 0.31 0.36 0.41 0.45 0.51]; %average and
lt_11_emp = [0.38 0.30 0.32 0.37 0.35 0.41]; % weak learners 

inA=0;inB=0;inC=0;inD=0; corrMapped=zeros(numSubjects,1);
corrMappedWrd=zeros(numSubjects,2); Look=[];
for subject=1:numSubjects
    Testset = xsit_result.test(subject).Testset; %size
    %Testset, has 3 distractors, 1 target/paired-target and its location
    Words=cell2mat(xsit_result.train(subject).Words);
    for trt=1:nTrials %%
        Look(1)= sum( xsit_result.test(subject).historyLt(trt,1:test_dur));%left extreme
        Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,1:test_dur));%left middle
        Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,1:test_dur));%right middle
        Look(4)= sum( xsit_result.test(subject).historyRt(trt,1:test_dur));%right extreme
        targLoc = Testset(trt,5);%the spatial location of target
        %next, add looking time to all 3 distracrors
        dstrLook=0;for i=1:4;if i~=targLoc; dstrLook= dstrLook + Look(i); end; end
        if (Look(targLoc) > dstrLook) 
            corrMapped(subject)= corrMapped(subject)+1;
            if (Testset(trt,targLoc))== 1 || (Testset(trt,targLoc))==2 || (Testset(trt,targLoc))==3
                corrMappedWrd(subject,1)=corrMappedWrd(subject,1)+1;
            else
                corrMappedWrd(subject,2)=corrMappedWrd(subject,2)+1;
            end
        end   
    end
    if corrMapped(subject) > 13 % based on Yu Zhang Fricker grouping
        learnerType(subject) = 1;
    elseif corrMapped(subject) < 8 % based on Yu Zhang Fricker grouping
        learnerType(subject) = -1;
    else
        learnerType(subject) = 0;
    end
end

figure(171);% proprtion correct response
blockNames={'Pre-trained words';'All other words'};
sts = [empirical_data_A(1) mean(corrMappedWrd(:,1))/3; empirical_data_A(2) mean(corrMappedWrd(:,2))/15 ];
errY =[0.05 SE(corrMappedWrd(:,1))/3;0.05 SE(corrMappedWrd(:,2))/15 ];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'FontSize',14)
ylabel('Proportion Correct');
ylim ([0 1]);
legend('Yu Zhong & Fricker (2012)', 'WOLVES Model');
hold on
a=hline(0.25,'k:');
set(a,'LineWidth',2.0);
grid on


%% Data for Output File Saving 
measurement_i = 'Mean proportion correct words learnt';
empirical_mean_i = empirical_data_A;
mean_i = [mean(corrMappedWrd(:,1))/3  mean(corrMappedWrd(:,2))/15];
SE_i =   [SE(corrMappedWrd(:,1))/3    SE(corrMappedWrd(:,2))/15];
RMSE_i = RMSE(empirical_mean_i, mean_i); MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


%% training analysis

%Experiment (Task) Variables

if (new_sim_formatting == true)
    nTrials = xsit_result.train(1).maxTR; 
    vis_Off = xsit_result.test(1).visuals_Off;
    word_On = xsit_result.train(1).word1_On;
    off_setter = floor(350/scale_factor); word_On = word_On + off_setter;
    word_Off = word_On + floor(1900/scale_factor);%%%1900 word length as specified in paper
    
else   
    nTrials=27; vis_Off=11250/scale_factor;
    word_On  = [2250 4500 6750 9000];    %word_Off = [4500 6750 9000 11250];%
    off_setter = 350;
    word_On = word_On + off_setter;
    word_Off = word_On + 1900;%%%1900 word length as specified in paper
    word_On = floor(word_On/scale_factor); word_Off = floor(word_Off/scale_factor);
end
propLook=NaN(numSubjects,nObjects,6);% 6 occurances of each object
for subject=1:numSubjects
    occurance=zeros(nObjects,1);
    training_pair = xsit_result.train(subject).training_pair; % the objects/words presneted in the trial (actually from the sequence 1-18)
    training_order = xsit_result.train(subject).training_order; %the order in which objects, word were presneted 
    for trt=1:nTrials 
        obj_at_loc = training_order(trt, 1:4);
        word_seq = training_order(trt, 5:8);         
        for w=1:4 % lets take word w  %it was presneted as w temporally.. therefore time on will be word_On(w):word_Off(w)
            look(1) = sum( xsit_result.train(subject).historyL(trt,word_On(w):word_Off(w)));
            look(2) = sum( xsit_result.train(subject).historyLWB(trt,word_On(w):word_Off(w)));
            look(3) = sum( xsit_result.train(subject).historyRWB(trt,word_On(w):word_Off(w)));
            look(4) = sum( xsit_result.train(subject).historyR(trt,word_On(w):word_Off(w)));
            currword = training_pair(trt,w);
            occurance(currword) = occurance(currword)+1;
            for locs =1:4
                if obj_at_loc(locs) == word_seq(w) % if this is the object of this currently played word
                     propLook(subject,currword,occurance(currword)) = look(locs)/((sum(look)));  %looking to target vs mean looking to distractors
                end
            end
        end
    end
end
    
figure(174);%Plot Mean proportion correct looking at every learning instance
errorbar(squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1)),squeeze(SE(nanmean(propLook(learnerType==1,:,:),2)))',plotStyle{1});
hold on
errorbar(squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1)),squeeze(SE(nanmean(propLook(learnerType==0,:,:),2)))',plotStyle{5});
hold on
errorbar(squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1)),squeeze(SE(nanmean(propLook(learnerType==-1,:,:),2)))',plotStyle{3});
hold on
xlim([0.5 6.5]);
ylim([0 1]);
ylabel('Proportion of Time on Target');
legend('Strong','Average','Weak');
set (gca, 'FontSize',14);
grid on

%% Data for Output File Saving 
measurement_i = 'Mean  Proportion of Time on Target';
empirical_mean_i = [lt_1_emp  lt_0_emp  lt_11_emp];
mean_i = [squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1))'  squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1))'  squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1))'];
SE_i =   [squeeze(SE(nanmean(propLook(learnerType==1,:,:),2)))'       squeeze(SE(nanmean(propLook(learnerType==0,:,:),2)))'       squeeze(SE(nanmean(propLook(learnerType==-1,:,:),2)))'];
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])