%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis for Yu & Smith (2007)
clear all; close all; % This script generates an output file
%Global Variables
scale_factor=8; nFeatures=2; 
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);
%%  

emp_Means = [0.88, 0.76, 0.53]; emp_Errors = [0.05, 0.09, 0.06];
for s=1:3
    %% Raw Data File Names
 if s==1;      simName = 'WPPR_Yu_Smith_Two_2007_results'%'WPPR_10h15k_oneThirdTime_Yu_Smith_Two_2007_results'%
 elseif s==2   simName = 'WPPR_Yu_Smith_Three_2007_results'%'WPPR_10h15k_oneThirdTime_Yu_Smith_Three_2007_results'%
 elseif s==3   simName = 'WPPR_Yu_Smith_Four_2007_results'%'WPPR_10h15k_oneThirdTime_Yu_Smith_Four_2007_results'%
end
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
    nObjects = 18; %%total number of objects
    nTrials = nObjects;
    test_dur = 1000/scale_factor;
end


%% TEST CONDITION ANALYSIS
inA=0;inB=0;inC=0;inD=0; corrMapped=zeros(numSubjects,1);
Wrong_inTrace= zeros(numSubjects,1); Correct_inTrace= zeros(numSubjects,1);
InCorr_assocs= zeros(numSubjects,1); EntropyTrace =  zeros(numSubjects,1);
    for subject=1:numSubjects
        Words =  cell2mat(xsit_result.train(subject).Words);
        %Sequence=  xsit_result.train(subject).Sequence;% 
        Testset = xsit_result.test(subject).Testset; 
        %Testset, has 3 distractors, 1 target/paired-target and its location 
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
            end   
        end

        %% trace analysis
        inputMapping1=squeeze(xsit_result.train(subject).hwf(1,:,:));
        inputMapping2=squeeze(xsit_result.train(subject).hwf(2,:,:));
        for kk=1:nObjects
            xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
            xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
            yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
        end
        C_inTr=0;W_inTr=0;
        for kk=1:nObjects
        %%% calculate the number of associations in the trace for each word 
            as_count1(kk)=0; assoc_c=1;
            while assoc_c < size(inputMapping1,1)
                if inputMapping1(assoc_c,yy(kk))>0.001
                    as_count1(kk)=as_count1(kk)+1;
                    while (assoc_c < size(inputMapping1,1)) && (inputMapping1(assoc_c,yy(kk))>=0.001)
                        assoc_c=assoc_c+1;
                    end
                else
                    assoc_c=assoc_c+1;
                end
            end
            as_count2(kk)=0; assoc_c=1;
            while assoc_c < size(inputMapping2,1)
                if inputMapping2(assoc_c,yy(kk))>0.001
                    as_count2(kk)=as_count2(kk)+1;
                    while (assoc_c < size(inputMapping2,1)) && (inputMapping2(assoc_c,yy(kk))>=0.001)
                        assoc_c=assoc_c+1;
                    end
                else
                    assoc_c=assoc_c+1;
                end
            end
            %%% calcuate trace strengths
            a_cv=inputMapping1(xx1(kk),yy(kk));b_cv=inputMapping2(xx2(kk),yy(kk));
            C_inTr= C_inTr+ nanmean([a_cv b_cv]);
            inputMapping1(xx1(kk),yy(kk))=0;
            inputMapping2(xx2(kk),yy(kk))=0;
            for jj=1:6
                inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
                inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
            end
            a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
            W_inTr = W_inTr + nanmean([nanmean(a_in(a_in>0.001)) nanmean(b_in(b_in>0.001))]);
        end
        Correct_inTrace(subject)=C_inTr/nObjects;
        Wrong_inTrace(subject)=W_inTr/nObjects;
        InCorr_assocs(subject)=nanmean([as_count1-1 as_count2-1]);
        EntropyTrace(subject)= nanmean( [entropy(inputMapping1) entropy(inputMapping2)] ); 

    end
    sim_means(s,1)  =  (mean(corrMapped)/nObjects);
    errY(s,1) =  ((std(corrMapped)/nObjects))/numSubjects;

    meanInCorr(s,1)  =  mean(InCorr_assocs);
    errInCorr(s,1) =  std(InCorr_assocs)/sqrt(length(InCorr_assocs));

    meanCorrStren(s,1)  =  mean(Correct_inTrace); 
    errCorrStren(s,1) =  std(Correct_inTrace)/sqrt(length(Correct_inTrace));
end

%% Data for Output File Saving 
measurement_i = 'Mean proportion correct looking at test';
empirical_mean_i = emp_Means;
mean_i = sim_means(:,1)';
SE_i =   SE(sim_means(:,1)');
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


figure(16)% Plot Mean proportion correct looking at words
blockNames={'2x2'; '3x3'; '4x4'};
sim_means(:,2) =emp_Means'; errY(:,2) = emp_Errors;
b=barwitherr(errY, sim_means);% Plot with errorbars
set(gca,'xticklabel',blockNames)
ylabel('Mean proportion correct');
legend('WOLVES Model','Yu & Smith 2007');
grid on
set(gca,'FontSize',14);
hold on
a=hline(0.25,'k:');
set(a,'LineWidth',2.0);

figure(22);%My own Entropy: No of incorrect traces
subplot(1,2,1);
blockNames={'2x2'; '3x3'; '4x4'};
b=barwitherr(errInCorr, meanInCorr);% Plot with errorbars
set(gca,'xticklabel',blockNames);
ylabel ('Number of incorrect associations');
ylim([0 4.2]);

subplot(1,2,2);%Strength of corrct associations
blockNames={'2x2'; '3x3'; '4x4'};scl = 1.7;
b=barwitherr(errCorrStren, meanCorrStren*scl);% Plot with errorbars
set(gca,'xticklabel',blockNames); set(b, 'FaceColor',[1,0.6,0]);
ylabel ('Association strength');
ylim([0 0.4]);
movegui('north');

%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])
