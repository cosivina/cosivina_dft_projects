%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis for Yurovsky_Yu_Smith, Cognitive Science (2013)
clear all; close all; % This script generates an output file

%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);

%% Raw Data File Name 
simName = 'WPPR_2_Yurovsky_Yu_Smith_2013_results' %WPPR_2_Yurovsky_Yu_Smith_2013_results
% loading the file
xsit_result = load ([simName '.mat']);
numSubjects=size(xsit_result.test,1);
      
% Experiment (Task) Variables
new_sim_formatting = true; % new format sims with task_variables saved at simulation
if (new_sim_formatting == true)
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrials = xsit_result.train(1).maxTRt; 
    %word_On = xsit_result.test(1).word1_On; word_Off = xsit_result.test(1).word1_Off; 
    %vis_Off = xsit_result.test(1).visuals_Off;
    %test_dur = 1000/scale_factor; % mean(word_Off - word_On);
    SINGLE = xsit_result.train(1).SINGLE;
    DOUBLE = xsit_result.train(1).DOUBLE;
    NOISE = xsit_result.train(1).NOISE;
else
    nObjects = 18; %%total number of objects
    nTrials=18;
    SINGLE = [1,5,6,10,14,18]; DOUBLE=[2,4,12,13,15,17]; NOISE=[3,7,8,9,11,16];
end

%% TEST CONDITION ANALYSIS
correctS=zeros(numSubjects,1);eitherD=zeros(numSubjects,1);bothD=zeros(numSubjects,1);
for subject=1:numSubjects

    for trt=1:nTrials
        Look(1)= sum( xsit_result.test(subject).historyLt(trt,:));%left extreme
        Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,:));%left middle
        Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,:));%right middle
        Look(4)= sum( xsit_result.test(subject).historyRt(trt,:));%right extreme
        targetPos = xsit_result.test(subject).test_data(trt,:);
        s1= char(targetPos(5));
        distLook=0;
        if ( strcmp(s1,'S')) %%single word 1 referent
            trgLoc = find(targetPos == 1);%the spatial location of target 
            for i=1:4;if i~=trgLoc; distLook=distLook+Look(i);end; end
            if Look(trgLoc)>distLook
                correctS(subject)=correctS(subject)+1;
            end           
        elseif ( strcmp(s1,'D')) %double word with 2 referents
            trgLoc1 = find(targetPos == 1);%the spatial location of target
            trgLoc2 = find(targetPos == 2);%the spatial location of target
            for i=1:4;if (i~=trgLoc1 && i~=trgLoc2 ) ; distLook=distLook+Look(i);end; end
            if (Look(trgLoc1)>distLook || Look(trgLoc2)>distLook)
                eitherD(subject)=eitherD(subject)+1;
            end
            if (Look(trgLoc1)>distLook && Look(trgLoc2)>distLook)
                bothD(subject)=bothD(subject)+1;
            end
        elseif ( strcmp(s1,'N')) % noise word with no referents
        end        
    end
end

figure(191)% Plot Mean proportion correct looking at words
blockNames={'Single'; 'Either'; 'Both'};
sts = [0.454 mean(correctS)/6 ; 0.698 mean(eitherD)/6 ; 0.301 mean(bothD)/6 ]; %since each 
errY =[0.132/2 SE(correctS)/6   ; 0.105/2 SE(eitherD)/6;    0.073/2 SE(bothD)/6];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',16)
ylabel('Mean proportion correct');
hold on;y = 0.25;line([0.7,1.3],[y,y]);hold on;y = 0.5;line([1.7,2.3],[y,y]);hold on;y = 0.17;line([2.7,3.3],[y,y]);%chance levels 
legend('Yurovsky, Yu & Smith, 2013','WOLVES Model');
grid on

%% Data for Output File Saving 
measurement_i = 'Mean proprotion correct looking at test';
empirical_mean_i = [0.454 0.698 0.301];
mean_i = [mean(correctS)/6  mean(eitherD)/6  mean(bothD)/6];
SE_i =   [SE(correctS)  SE(eitherD)  SE(bothD)];
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])

%% Memory Trace Analysis

Wrong_inTrace= zeros(numSubjects,1); Correct_inTraceS= zeros(numSubjects,1);Correct_inTraceD= zeros(numSubjects,1);
for subject=1:numSubjects    
    inputMapping1=squeeze(xsit_result.train(subject).hwf(1,:,:));
    inputMapping2=squeeze(xsit_result.train(subject).hwf(2,:,:));
    for kk=1:nObjects
        if ismember(kk,SINGLE)
            xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
            xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
            yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
        elseif ismember(kk,DOUBLE)
            xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
            xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
            yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
        elseif ismember(kk,NOISE)
            find(NOISE==kk);
            xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
            xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
            yy(kk)=cell2mat(xsit_result.train(subject).Words(DOUBLE(find(NOISE==kk))));
        end
    end
    C_inTr=0;W_inTr=0;C_inTrS=0;C_inTrD=0;
    for kk=1:nObjects
    %%% calculate the number of associations in the trace for each word 
        as_count1(kk)=0; f_loc=1;
        while f_loc < size(inputMapping1,1)
            if inputMapping1(f_loc,yy(kk))>0.001
                as_count1(kk)=as_count1(kk)+1;
                while (f_loc < size(inputMapping1,1)) && (inputMapping1(f_loc,yy(kk))>=0.001)
                    f_loc=f_loc+1;
                end
            else
                f_loc=f_loc+1;
            end
        end
        as_count2(kk)=0; f_loc=1;
        while f_loc < size(inputMapping2,1)
            if inputMapping2(f_loc,yy(kk))>0.001
                as_count2(kk)=as_count2(kk)+1;
                while (f_loc < size(inputMapping2,1)) && (inputMapping2(f_loc,yy(kk))>=0.001)
                    f_loc=f_loc+1;
                end
            else
                f_loc=f_loc+1;
            end
        end
        %%% calcuate trace strengths
        if ismember(kk, SINGLE)
            a_cvS=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cvS=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
            C_inTrS= C_inTrS+ nanmean([nanmean(a_cvS(a_cvS>0.001))  nanmean(b_cvS(b_cvS>0.001))]);
        else
            a_cvD=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cvD=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
            C_inTrD= C_inTrD+ nanmean([nanmean(a_cvD(a_cvD>0.001))  nanmean(b_cvD(b_cvD>0.001))]);
        end
        inputMapping1(xx1(kk),yy(kk))=0;
        inputMapping2(xx2(kk),yy(kk))=0;
        for jj=1:6
            inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
            inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
        end
        a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
        W_inTr = W_inTr + nanmean([nanmean(a_in(a_in>0.001)) nanmean(b_in(b_in>0.001))]);
    end
    Correct_inTraceS(subject)=C_inTrS/6;
    Correct_inTraceD(subject)=C_inTrD/12;
    Wrong_inTrace(subject)=W_inTr/nObjects;
    InCorr_assocs(subject)=nanmean([as_count1-1 as_count2-1]);
    EntropyTrace(subject)= nanmean( [entropy(inputMapping1) entropy(inputMapping2)] ); 
end

figure(192);%My own Entropy: No of incorrect traces
blockNames={'Single';'Double'};
sts = [nanmean(Correct_inTraceS)*6; nanmean(Correct_inTraceD)*6  ];
errY =[nanstd(Correct_inTraceS)/sqrt(length(Correct_inTraceS)); nanstd(Correct_inTraceD/sqrt(length(Correct_inTraceD)))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'FontSize',16);
ylabel ('Association Trace Strength')
grid on







