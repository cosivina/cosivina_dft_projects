%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis
clear all; close all;
%Global Variables
scale_factor=8;MIN_LOOK_DURATION=200/scale_factor; nFeatures=2;
%Plotting Variables
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};
compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);

%% Raw Data File Name             
simName = 'WPPR_Smith_Yu_2013_results'
% loading the file
xsit_result = load ([simName '.mat']); 
numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);%

        
%% TEST CONDITION ANALYSIS
%Experiment Variables

%Experiment (Task) Variables
new_sim_formatting = true; % new format sims with task_variables saved at simulation
if (new_sim_formatting == true)
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrainTrials = xsit_result.train(1).maxTR; 
    nTestTrials = xsit_result.train(1).maxTRt;
    TRAIN_DUR = xsit_result.train(1).t_max - floor(1000/scale_factor);
    TEST_DUR = xsit_result.train(1).t_maxt - floor(1000/scale_factor);
    word_On = xsit_result.test(1).word_On; word_Off = xsit_result.test(1).word_Off; 
    vis_On = xsit_result.test(1).visuals_On;
    vis_Off = xsit_result.test(1).visuals_Off;
    if length(word_On) > 1, words_within_range = length(word_On)-1; else, words_within_range = 1; end
    if length(word_On) > 1, word_effect_Len = floor(1750/scale_factor); else, word_effect_Len=floor(999/scale_factor); end % to map to 1.75 seconds on paper
else
    nObjects=6; nTrainTrials=30; nTestTrials=12;  
    TRAIN_DUR= floor(4000/scale_factor); TEST_DUR= floor(8000/scale_factor);
    
    word_On = floor([8 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
    word_Off = floor((745+[8 1800 3500 5200 6900])/scale_factor);
    word_effect_Len=floor(1750/scale_factor);% to map to 1.75 seconds on paper
    
    % 1-sec test
    %TEST_DUR= floor(1000/scale_factor);
    %word_On = 1;  %% 1 second test 
    %word_Off = floor(1000/scale_factor);
    %
    %
    vis_On = 1;vis_Off = TEST_DUR;
    if length(word_On) > 1, words_within_range = length(word_On)-1; else, words_within_range = 1; end
    if length(word_On) > 1, word_effect_Len = floor(1750/scale_factor); else, word_effect_Len=floor(999/scale_factor); end % to map to 1.75 seconds on paper

end
Look_Smoothening_SY_2AFC; % DATA SMOOTHENING if necessary 

%% 
targLookTime=zeros(numSubjects,nTestTrials);
dstrLookTime=zeros(numSubjects,nTestTrials);
goodLearners=NaN(numSubjects,1);
manIsLearner=NaN(numSubjects,1);
LearntWords= NaN(numSubjects,nObjects);
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
for subject=1:numSubjects
    lcorrect=0;rcorrect=0;targWord=zeros(nObjects,1);dstrWord=zeros(nObjects,1);
    twoSuccess=zeros(nObjects,1);
    for trt=1:nTestTrials
        lLook= sum( xsit_result.test(subject).historyLt(trt,vis_On:vis_Off));%full trial
        rLook= sum( xsit_result.test(subject).historyRt(trt,vis_On:vis_Off));%
        s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
        for kk=1:nObjects    
          if (xsit_result.train(subject).Words{kk} == xsit_result.test(subject).test_pair(trt,2*nFeatures+1))%word index               
               if ( strcmp(s1,'L')) 
                       targWord(kk)=targWord(kk)+lLook;
                       dstrWord(kk)=dstrWord(kk)+rLook;
                       if (lLook > rLook)   %%loking at target on this word trial
                           twoSuccess(kk)=twoSuccess(kk)+1;
                       else
                           twoSuccess(kk)=twoSuccess(kk)+0;
                       end
               elseif ( strcmp(s1,'R'))
                       targWord(kk)=targWord(kk)+rLook;
                       dstrWord(kk)=dstrWord(kk)+lLook;
                       
                       if (lLook < rLook)   %%loking at target on this word trial
                           twoSuccess(kk)=twoSuccess(kk)+1;
                       else
                           twoSuccess(kk)=twoSuccess(kk)+0;
                       end
               else
                       disp('ERROR reading test_pair_char');
               end
          end
        end
        if ( strcmp(s1,'L')) 
            targLookTime(subject,trt)=lLook;
            dstrLookTime(subject,trt)=rLook;
        elseif ( strcmp(s1,'R'))
            targLookTime(subject,trt)=rLook;
            dstrLookTime(subject,trt)=lLook;
        else
              disp('ERROR reading test_pair char');
        end
    end%% trials loop
    
    for kk=1:nObjects
        if (targWord(kk)>dstrWord(kk))
            LearntWords(subject,kk)=1;
        else
            LearntWords(subject,kk)=0;
        end
    end
   if sum(twoSuccess==2) > 3 %% my new criterion (if looked good both times the word was presented)
       manIsLearner(subject)=1;
   else
       manIsLearner(subject)=0;
   end
    if (mean(targLookTime(subject,:)) > mean(dstrLookTime(subject,:)))
        goodLearners(subject)=1;
    else
        goodLearners(subject)=0;
    end    
end
correct_proportion=targLookTime./(targLookTime+dstrLookTime);

% disp('t-test statistics between Target and Distractor Looking');
% [h,p,ci,stats] = ttest(mean(targLookTime,2),mean(dstrLookTime,2),'Tail','right')
% disp('t-test statistics for words learnt');
% [h,p,ci,stats] = ttest(nansum(LearntWords,2),3,'Tail','right')


%% Data for Output File Saving 
%%
measurement_i = 'Preferential looking time ratio (at test)';
empirical_mean_i = [0.54];
mean_i = nanmean(nanmean(correct_proportion,2));
SE_i = SE(nanmean(correct_proportion,2));
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


%%
measurement_i = 'Prop words learnt at test (out of 6)';
%empirical_mean_i = [0.63];
mean_i = mean(sum(LearntWords(),2))/6;
%SE_i = SE(sum(LearntWords(),2));
%RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
%row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
%xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


%%
% The below measure should be averaged over 100 such raw data files /simulations/runs
measurement_i = 'Proportion of Strong Learners (SE) at test';
% empirical_mean_i = [];
 mean_i = sum(goodLearners)/numSubjects; % actually should be mean over this measure
% SE_i = NaN;
% RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
% row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
 xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
%xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);



%% Within Trial Time Measures and Plots
word_Vec=zeros(1,vis_Off);
for wc=1:length(word_On)
    for w = word_On(wc):word_Off(wc)
        word_Vec(w)=1;
    end
end
tLookAtLearnt=zeros(nTestTrials,1);tLookAtNonLrnt=zeros(nTestTrials,1);
for subject=1:numSubjects
    savestate_historyLt = xsit_result.test(subject).historyLt(:,vis_On:vis_Off);
    savestate_historyRt = xsit_result.test(subject).historyRt(:,vis_On:vis_Off);
    % create the off-looking history Vector
    for i=1:nTestTrials
        for j=1:size(savestate_historyLt,2)
           if  (round(savestate_historyLt(i,j)) + round(savestate_historyRt(i,j))) > 0.5
               savestate_historyOt(i,j)=0; 
           else savestate_historyOt(i,j)=1; end
        end
    end
    for tt=1:nTestTrials %% each test trial
        %%% measure sync time between off-looking & word presentation
        wordoff_time=1-savestate_historyOt(tt,:); %sync_time(subject,tt)= sum(wordoff_time.*word_Vec);      
        %%% measure looking ratio at every time instant 
        lookingOn(subject,tt,:)= ((round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:)))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
    end

    % proportion of looking to learnt vs unlearnt words
    for kk=1:nObjects
        show_trial=1;% first time or second time word prsented of 12 trials
        for tt=1:nTestTrials %each test trial 
            if (xsit_result.train(subject).Words{kk}== xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%word prsented in trial
                s1= char(xsit_result.test(subject).test_pair(tt,2*nFeatures+2));
                if (strcmp(s1,'L'))% referent on left
                    % COLLECT TIME DATA for each show
                    if (show_trial==1)  % first presentation of this word of 12 trials  
                        look_at_targ_show1(subject,kk,:)=(round(savestate_historyLt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';%+savestate_historyOt(tt,1:8000)
                        look_at_dist_show1(subject,kk,:)=(round(savestate_historyRt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
                    elseif (show_trial==2) % second presentation of this word
                        look_at_targ_show2(subject,kk,:)=(round(savestate_historyLt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';%add time historyLt
                        look_at_dist_show2(subject,kk,:)=(round(savestate_historyRt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
                   end
                elseif (strcmp(s1,'R')) % referent on right
                    if (show_trial==1)  % first presentation of this word of 12 trials
                        look_at_targ_show1(subject,kk,:)=(round(savestate_historyRt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';%add time historyRt
                        look_at_dist_show1(subject,kk,:)=(round(savestate_historyLt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
                    elseif (show_trial==2) % second of 12 trials with this word
                        look_at_targ_show2(subject,kk,:)=(round(savestate_historyRt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';%add time historyRt
                        look_at_dist_show2(subject,kk,:)=(round(savestate_historyLt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
                    end     
                end
               show_trial=show_trial+1; 
            end
        end
    end
        
    %% AT TEST ...... find total duration of looking to learnt vs unlearnt words
    for kk=1:nObjects
       if (LearntWords(subject,kk)==1)%if the subject has learnt the kk word.. learnt is measured on looking basis.. change to mem_matrix basis if needed
                for tt=1:nTestTrials%each test trial 
                    if(xsit_result.train(subject).Words{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial name1
                        s1=char(xsit_result.test(subject).test_pair(tt,2*nFeatures+2));
                        if (strcmp(s1,'L'))% on left
                            tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyLt(tt,:));%add time historyLt
                        elseif (strcmp(s1,'R'))
                            tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyRt(tt,:));% add time historyRt
                        end
                    end
                end        
        elseif (LearntWords(subject,kk)==0)%if word is NOT learnt
                for tt=1:nTestTrials
                    if(xsit_result.train(subject).Words{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial
                        s1=char(xsit_result.test(subject).test_pair(tt,2*nFeatures+2));
                        if (strcmp(s1,'L'))% on left
                            tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyLt(tt,:));%add time historyL
                        elseif (strcmp(s1,'R'))
                            tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyRt(tt,:));% add time historyR
                        end
                    end
                end
        else
            disp('Error in LearntWords array: supurious data ');
        end
    end
end

all_objLook= NaN (numSubjects,nTestTrials,TEST_DUR);
all_targL= NaN (numSubjects,nObjects,TEST_DUR);
all_distL= NaN (numSubjects,nObjects,TEST_DUR);
larW=NaN (numSubjects,nObjects,TEST_DUR);
ularW=NaN (numSubjects,nObjects,TEST_DUR);
for subject=1:numSubjects
        all_objLook(subject,:,:) = lookingOn(subject,:,:);
        for kk=1:nObjects
            all_targL(subject,kk,:)  = mean([look_at_targ_show1(subject,kk,:) look_at_targ_show2(subject,kk,:)]);
            all_distL(subject,kk,:)  = mean([look_at_dist_show1(subject,kk,:) look_at_dist_show2(subject,kk,:)]);
            if LearntWords (subject,kk) == 1
              larW(subject,kk,:)  = mean([look_at_targ_show1(subject,kk,:) look_at_targ_show2(subject,kk,:)]);
            end
            if LearntWords (subject,kk) == 0
              ularW(subject,kk,:)  = mean([look_at_targ_show1(subject,kk,:) look_at_targ_show2(subject,kk,:)]);
            end
        end
end



%% Plot Learners' proportion looking to T/D following word presentation at test 
twg=[];dwg=[];
targ_avg_word_1=(squeeze(nanmean(nanmean(all_targL(goodLearners()==1,:,:),1),2))); targ_avg_word = reshape(targ_avg_word_1,[],vis_Off);
dist_avg_word_1=(squeeze(nanmean(nanmean(all_distL(goodLearners()==1,:,:),1),2))); dist_avg_word = reshape(dist_avg_word_1,[],vis_Off);
w_gap=word_effect_Len;% to map to 1.75 seconds on paper
for wcc=1:words_within_range
    twg(wcc,:)=targ_avg_word(word_On(wcc):word_On(wcc)+w_gap);
    dwg(wcc,:)=dist_avg_word(word_On(wcc):word_On(wcc)+w_gap);
end
figure (233);%%Plot sum of all bins together
plot (nanmean(twg(:,:),1),'LineWidth',3)
hold on 
plot (nanmean(dwg(:,:),1),'LineWidth',3)
legend('Target',    'Foil');
xlabel('time');
ylabel('proportion of learners looking at target or foil');
ylim([0.2 0.7]);
xticks([0 50 100]);
xticklabels({'','.85','1.7'});
movegui('north');

%% Plot NonLearners' proportion looking to T/D following word presentation at test 
twg=[];dwg=[];
w_gap=word_effect_Len;
targ_avg_word_1=(squeeze(nanmean(nanmean(all_targL(goodLearners()==0,:,:),1),2))); targ_avg_word = reshape(targ_avg_word_1,[],vis_Off);
dist_avg_word_1=(squeeze(nanmean(nanmean(all_distL(goodLearners()==0,:,:),1),2))); dist_avg_word = reshape(dist_avg_word_1,[],vis_Off);
for wcc=1:words_within_range
    twg(wcc,:)=targ_avg_word(word_On(wcc):word_On(wcc)+w_gap);
    dwg(wcc,:)=dist_avg_word(word_On(wcc):word_On(wcc)+w_gap);
end
figure (234);%Plot sum of all bins together
plot (nanmean(twg(:,:),1),'LineWidth',3)
hold on 
plot (nanmean(dwg(:,:),1),'LineWidth',3)
legend('Target',    'Foil');
xlabel('time');
ylabel('proportion of nonlearners looking at target or foil');
ylim([0.2 0.7]);
xticks([0 50 100]);
xticklabels({'','.85','1.7'});
movegui('south');




%% TRAINING CONDITION ANALYSIS trials analysis

word_On = floor([368 1850]/scale_factor);  %% XTRAP    
word_Off = floor((745+[368 1850])/scale_factor);word_Len=floor(745/scale_factor);
C_word_On = word_On+ floor(500/scale_factor);C_word_Off = floor(([1850 4000])/scale_factor);
vis_On=1;vis_Off=(TRAIN_DUR); nFix_limit=10;%max limit 2.5 fixations per sec = 10 per training trial here


%%
corrLookTimeTraining=zeros(numSubjects,nTrainTrials);%number of training trials
incorrLookTimeTraining=zeros(numSubjects,nTrainTrials);
totnlooks=zeros(numSubjects,nTrainTrials);
meanlookdur =zeros(numSubjects,nTrainTrials);
TotalLookTime=zeros(numSubjects,nTrainTrials);
totlonglookdur = zeros(numSubjects,nTrainTrials);
tLookRepeated= zeros(numSubjects,nObjects);tLookVarying=zeros(numSubjects,nObjects);
mLookCorrect= zeros(numSubjects,nObjects);mLookIncorrect=zeros(numSubjects,nObjects);
for subject=1:numSubjects
    savestate_historyL = xsit_result.train(subject).historyL(:,vis_On:vis_Off);
    savestate_historyR = xsit_result.train(subject).historyR(:,vis_On:vis_Off);    
    % create the off-looking history Vector
    for i=1:nTrainTrials
        for j=1:TRAIN_DUR/scale_factor
           if  (round(savestate_historyL(i,j)) + round(savestate_historyR(i,j))) > 0; savestate_historyO(i,j)=0;               
           else savestate_historyO(i,j)=1; end
        end
    end
    
    %%%%% looking during training
    for tr=1:nTrainTrials
        s1=char(xsit_result.train(subject).training_pair(tr,2*nFeatures+3));
        if (strcmp(s1,'P'))% on left
           corrLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyR(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking to target object right side
            
            incorrLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyL(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking wrong way
        elseif (strcmp(s1,'X'))
           corrLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyL(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking to target object right side
            
            incorrLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyR(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking wrong way            
        end
        tLookRepe(subject,tr)=sum(savestate_historyL(tr,:));
        tLookVary(subject,tr)=sum(savestate_historyR(tr,:));
    end
    
    %% no# of looks and fixation-durations calculation
    nlooks=zeros(2,nTrainTrials); %L/R
    longlookdur=zeros(2,nTrainTrials);
    all_look_dur=NaN(2,nTrainTrials,nFix_limit);

    for side=1:2     
        if side == 1
            ldata = savestate_historyL;
        else
            ldata = savestate_historyR;
        end
        for tr=1:size(ldata,1)
            prevlook=0;
            templonglookdur=0;
            for time=1:size(ldata,2)
                if (round(ldata(tr,time)) == 1)
                    if prevlook == 1
                        templonglookdur = templonglookdur+1; 
                    else
                        prevlook = 1;
                        templonglookdur=1;
                    end                    
                else
                    if prevlook == 1
                        if templonglookdur > (MIN_LOOK_DURATION+5)
                            nlooks(side,tr) = nlooks(side,tr)+1;
                            all_look_dur(side,tr,nlooks(side,tr))= templonglookdur;
                            if templonglookdur > longlookdur(side,tr)
                                longlookdur(side,tr) = templonglookdur;
                            end
                            templonglookdur=0;
                        end
                        prevlook = 0;
                    end
                end
            end
            if (round(ldata(tr,time-1)) == 1)
               if templonglookdur > (MIN_LOOK_DURATION+5)
                    nlooks(side,tr) = nlooks(side,tr)+1;
                    all_look_dur(side,tr,nlooks(side,tr))= templonglookdur;
                    if templonglookdur > longlookdur(side,tr)
                        longlookdur(side,tr) = templonglookdur;
                    end
               end
            end
        end   
    end

    for blockz=1:nObjects
        TinA=(nObjects-1)*(blockz-1)+1;
        TinB=(nObjects-1)*(blockz);
        tLookRepeated(subject,blockz)=sum(sum(savestate_historyL(TinA:TinB,:)));
        tLookVarying(subject,blockz)=sum(sum(savestate_historyR(TinA:TinB,:)));

        mLookCorrect(subject,blockz)= mean(corrLookTimeTraining(subject,TinA:TinB));
        mLookIncorrect(subject,blockz)= mean(incorrLookTimeTraining(subject,TinA:TinB));
    end
    
    totnlooks(subject,:)=sum(nlooks,1);
    meanLukhadur(subject,:)=nanmean(nanmean(all_look_dur,3),1);
    totlonglookdur(subject,:)=max(longlookdur,[],1);    
    TotalLookTime(subject,:)=sum(savestate_historyL')+sum(savestate_historyR');    
    meanlookdur(subject,:)= TotalLookTime(subject,:)./totnlooks(subject,:);
     %% calculate entropy in looking on very trial
     total_trialLook_duration=nansum(nansum(all_look_dur,3),1);%1x30 from 2x30x10
     mean_trialLook_duration=nanmean(nanmean(all_look_dur,3),1);%1x30 from 2x30x10
     pdf=NaN(size(all_look_dur));%2x30x10
     variancA=NaN(size(all_look_dur));%2x30x10
     EntropySub(subject,:)= 0;
     VarianceSub(subject,:)=0;
     
     for trial=1:size(nlooks,2) 
        for side=1:size(nlooks,1)
             pdf_side(side,:)=abs(all_look_dur(side,trial,:))./total_trialLook_duration(trial);
             variance_side(side,:)= (all_look_dur(side,trial,:)./total_trialLook_duration(trial)) .*((all_look_dur(side,trial,:)-mean_trialLook_duration(trial)).^2) ;
        end
        pdf=[pdf_side(1,:) pdf_side(2,:)];
        EntropySub(subject,trial)= -1* nansum( pdf(:).*log2(pdf(:)) );  
        variancA=[variance_side(1,:) variance_side(2,:)];
        VarianceSub(subject,trial)= nansum(variancA);
     end     
end


% figure (13);%Plot mean look duration of each fixation
% errorbar(mean(meanLukhadur)*scale_factor/1000,SE(meanLukhadur)*scale_factor/1000, plotStyle{1});% mean look duration % multiped by timing scale factor
% xlabel('training trial');
% ylabel('mean look duration');
% ylim([0.6 1.5]);

%% Plot Habituation over training as proportion of looking Varying and Repeated objects
figure (232);%Plot Mean proportion of looking Varying vs repeated C
rectangle('Position',[0.7,0,4.6,1],'FaceColor',[.1 .5 .5 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[5.7,0,4.6,1],'FaceColor',[.5 .5 .1 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[10.7,0,4.6,1],'FaceColor',[.1 .5 .5 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[15.7,0,4.6,1],'FaceColor',[.5 .5 .1 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[20.7,0,4.6,1],'FaceColor',[.1 .5 .5 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[25.7,0,4.6,1],'FaceColor',[.5 .5 .1 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
errorbar(mean(tLookVary./(tLookVary+tLookRepe)),std(tLookVary./(tLookVary+tLookRepe))./sqrt(length(tLookVary)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepe./(tLookVary+tLookRepe)),std(tLookRepe./(tLookVary+tLookRepe))./sqrt(length(tLookRepe)), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('training trial');
ylim([0.3 0.7]);
xlim([0 31]);
ylabel('Proportion looking time');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('Habituation over training');

%% Plot proportion looking time for Learners to Varying vs Repeated objects
figure (235);%
errorbar(mean(tLookVarying((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:))),std(tLookVarying((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:)))./sqrt(length(tLookVarying((goodLearners()==1),:))),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeated((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:))),std(tLookRepeated((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:)))./sqrt(length(tLookVarying((goodLearners()==1),:))), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('trial block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
ylabel('Proportion looking time Strong Learners');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('Learners');
movegui('west');

%% Plot proportion looking time for Non-learners to Varying vs Repeated objects
figure (236);
errorbar(mean(tLookVarying((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:))),std(tLookVarying((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:)))./sqrt(length(tLookVarying((goodLearners()==0),:))),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeated((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:))),std(tLookRepeated((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:)))./sqrt(length(tLookVarying((goodLearners()==0),:))), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('trial block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
ylabel('Proportion Looking Time Weak Learners');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('NonLearners');
movegui('east');

%% Proportion looking to varying and repeated objects (at training)
empirical_Learners_Varying = [0.61 0.52 0.54 0.47 0.48 0.35]; %empirical values
empirical_nonLearners_Varying = [0.63 0.58 0.62 0.50 0.55 0.44];
empirical_Learners_Repeated = [0.27 0.32 0.26 0.31 0.26 0.39];
empirical_nonLearners_Repeated = [0.24 0.28 0.17 0.22 0.20 0.30];
rough_Errorbars_var_rep = [0.1 0.1 0.1 0.1 0.1 0.1];
%wolves simulation values
prop_var_L=tLookVarying((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:));
prop_var_NL=tLookVarying((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:));
prop_rep_L=tLookRepeated((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:));
prop_rep_NL=tLookRepeated((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:));

measurement_i = 'Proportion looking time varying and repeated objects (at training)';
empirical_mean_i = [empirical_Learners_Varying, empirical_nonLearners_Varying,empirical_Learners_Repeated,empirical_nonLearners_Repeated];
mean_i = [mean(prop_var_L), mean(prop_var_NL), mean(prop_rep_L), mean(prop_rep_NL)];
SE_i = [SE(prop_var_L), SE(prop_var_NL), SE(prop_rep_L), SE(prop_rep_NL)];
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])