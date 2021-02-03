%% Data Analysis 
%xsit_result.sim.saveSettings('test.json')%method to save simulator params
%% LOADING
clear all; close all;
nObjects=6;nFeatures=2;nTrainTrials=30; nTestTrials=12; scale_factor=8; TRAIN_DUR=4000; TEST_DUR=8000;MIN_LOOK_DURATION=160/scale_factor; %check with auto file 160
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};;compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
TASK= 'XSIT';% 'XTRAP';%
% xsit_result.Names={[1] [2] [3] [4] [5] [6]};
%paramValues=[0.5 1.5 2.5 3.5 4.5 5.5];%[4.3 4.5 4.7 5.1 5.5 5.9];% 0 2 4 6  %31 38 40 42 %
%59 %  % 0 1.5 3.1 4.2 5.5 7.1 8.5
%paramValues=[4.85 5.05 5.25 5.45 5.65 5.85];
paramValues=[5 10 15];
%paramValues=[0 1 2 3 4 5 6 7 8 9 10];%[125 500 1000 2000 4000];
%xsit_result.train=[];
%xsit_result.test=[];
for param=1:length(paramValues)
parVal= paramValues(param);
legendInfo{param}= [plotStyle{param} ' '  num2str(parVal)];%2*param -1
%legendInfo{2*param}= [plotStyle{param+compStyle} num2str(parVal)];
simName = ['wfChanges_conwmf15_12h5k_sytst_hmwc' num2str(parVal) '_fix49_Smith_Yu_2013_']; %%simName = 'B_xsitNew_16_sub_';
%  if param==1; simName = 'goodlooking_wordwf42_Smith_Yu_2008_';
%  elseif param==2; simName = 'wfChanges_wordwf42_Smith_Yu_2008_';
% % elseif param==3; simName = 'wfChanges2_38_hwf4_Smith_Yu_2008_';
%  end 
%simName = ['vary_mem_sy_base04-Oct-2018M'  num2str(parVal*1) '_test_only_learners_smithYu_'];
%simName = 'goodlooking_3OrdersTest_Smith_Yu_2008_';
%simName = 'muddled_mem_6_pairs_sy_base05-Oct-2018_test_only_learners_smithYu_'%unified_param_run8_02-Oct-2018_Smith_Yu_2013_'%'bbc_H50F30Asc50_Tzero_Smith_Yu_2013_'; %test4_Smith_Yu_2008_
OutName = [simName,'results.mat']
xsit_result = load (OutName);
%numSubjects=size(xsit_result1.test,1);
%xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 
%xsit_result.train = [xsit_result.train ; xsit_result1.train];  
%xsit_result.test = [xsit_result.test ; xsit_result1.test];
%end
numSubjects=size(xsit_result.test,1);
xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 
%xsit_result.sim.saveSettings('recoverParams.json'); visdiff ('recoverParams1.json', 'recoverParams.json')
%% DATA SMOOTHENING
for subject=1:numSubjects
    for trt=1:nTestTrials
        for side=1:2      
         if side == 1
             ldata = round(xsit_result.test(subject).historyLt(trt,:));
         else
             ldata = round(xsit_result.test(subject).historyRt(trt,:));
         end
            prevlooking=0;
            gaplen=0;
            for time=1:length(ldata)
                if (round(ldata(time)) == 0)
                    if prevlooking==0 %%continue adding off looking gap length
                        gaplen=gaplen+1;
                    else
                        prevlooking=0;
                        gaplen=0;
                    end
                else
                    if gaplen <= MIN_LOOK_DURATION
                        ldata(time-gaplen:time)=1;
                    end
                    gaplen=0;
                end
            end
            if side == 1
                xsit_result.test(subject).historyLt(trt,:)=ldata;
            else
                xsit_result.test(subject).historyRt(trt,:)=ldata;
            end    
        end
    end
    for tr=1:nTrainTrials
        for side=1:2      
         if side == 1
             ldata = round(xsit_result.train(subject).historyL(tr,:));
         else
             ldata = round(xsit_result.train(subject).historyR(tr,:));
         end
            prevlooking=0;
            gaplen=0;
            for time=1:length(ldata)
                if (round(ldata(time)) == 0)
                    if prevlooking==0 %%continue adding off looking gap length
                        gaplen=gaplen+1;
                    else
                        prevlooking=0;
                        gaplen=0;
                    end
                else
                    if gaplen <= MIN_LOOK_DURATION
                        ldata(time-gaplen:time)=1;
                    end
                    gaplen=0;
                end
            end
            if side == 1
                xsit_result.train(subject).historyL(tr,:)= ldata;
            else
                xsit_result.train(subject).historyR(tr,:)= ldata;
            end
        end
    end
end

%% TEST CONDITION ANALYSIS
if strcmp(TASK, 'XSIT')
    word_On = floor([500 2000 4500 6000]/scale_factor);%XSIT    
    word_Off = floor([1500 3000 5500 7000]/scale_factor);word_Len=floor(1000/scale_factor);
elseif strcmp(TASK, 'XTRAP')
    word_On = floor([8 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
    word_Off = floor((745+[8 1800 3500 5200 6900])/scale_factor);word_Len=floor(745/scale_factor);
end
vis_On = 1;vis_Off = floor(TEST_DUR/scale_factor);

correct_proportion=zeros(numSubjects,2);
targLookTime=zeros(numSubjects,nTestTrials);
dstrLookTime=zeros(numSubjects,nTestTrials);
goodLearners=NaN(numSubjects,1);
LearntWords= NaN(numSubjects,nObjects);
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;

for subject=1:numSubjects
    lcorrect=0;rcorrect=0;targWord=zeros(nObjects,1);dstrWord=zeros(nObjects,1);
    for trt=1:nTestTrials
        lLook= sum( xsit_result.test(subject).historyLt(trt,vis_On:vis_Off));%full trial
        rLook= sum( xsit_result.test(subject).historyRt(trt,vis_On:vis_Off));%
%         lLook=0; for wint=1:length(word_On); lLook= lLook+sum(xsit_result.test(subject).word_On(wint):word_Off(wint));end %only word-on time-windows
%         rLook=0; for wint=1:length(word_On); rLook= rLook+sum(xsit_result.test(subject).word_On(wint):word_Off(wint));end %only word-on time-windows

        s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
        for kk=1:nObjects    
          if (xsit_result.train(subject).Words{kk} == xsit_result.test(subject).test_pair(trt,2*nFeatures+1))%word index
                             
               if ( strcmp(s1,'L')) 
                       targWord(kk)=targWord(kk)+lLook;
                       dstrWord(kk)=dstrWord(kk)+rLook;
               elseif ( strcmp(s1,'R'))
                       targWord(kk)=targWord(kk)+rLook;
                       dstrWord(kk)=dstrWord(kk)+lLook;
               else
                       disp('ERROR reading test_pair_char');
               end
          end
        end
        
        if ( strcmp(s1,'L')) 
            targLookTime(subject,trt)=lLook;
            dstrLookTime(subject,trt)=rLook;
            lcorrect=lcorrect+lLook/(lLook+rLook);
        elseif ( strcmp(s1,'R'))
              targLookTime(subject,trt)=rLook;
              dstrLookTime(subject,trt)=lLook;
              rcorrect=rcorrect+rLook/(lLook+rLook);
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
    lcorrect = lcorrect/(nTestTrials);     rcorrect = rcorrect/(nTestTrials);
    correct_proportion(subject,1)=lcorrect;    correct_proportion(subject,2)=rcorrect;
    if (mean(targLookTime(subject,:)) > mean(dstrLookTime(subject,:)))
        goodLearners(subject)=1;
    else
        goodLearners(subject)=0;
    end
end

%goodLearnerlist(param,:)= goodLearners; 
%%%%Fig 1
test_TotalLookingTime_Mean(param)=(mean(mean(targLookTime+dstrLookTime))*(scale_factor))/1000;
test_TotalLookingTime_Error(param)=((std(mean(targLookTime+dstrLookTime,2))*(scale_factor))/1000)/ sqrt(length(targLookTime+dstrLookTime));

%%% fig 2
test_TargetLookingTime_Mean(param)=(mean(mean(targLookTime,2))*(scale_factor))/1000;
test_TargetLookingTime_Error(param)=((std(mean(targLookTime,2))*(scale_factor))/1000)/ sqrt(length(targLookTime));
test_DistractorLookingTime_Mean(param)=(mean(mean(dstrLookTime,2))*(scale_factor))/1000;
test_DistractorLookingTime_Error(param)=((std(mean(dstrLookTime,2))*(scale_factor))/1000)/ sqrt(length(dstrLookTime));

%% fig 21
test_correct_proportion(param)=mean(mean(targLookTime./(targLookTime+dstrLookTime)));

%%% fig 31
test_STargetLookingTimePerTrial_Mean(param,:)=mean(targLookTime((goodLearners()==1),:))*scale_factor/1000;
test_STargetLookingTimePerTrial_Error(param,:)=(std(targLookTime((goodLearners()==1),:))*scale_factor/1000)/ sqrt(length(targLookTime((goodLearners()==1))));
test_SDistractorLookingTimePerTrial_Mean(param,:)=mean(dstrLookTime((goodLearners()==1),:))*scale_factor/1000;
test_SDistractorLookingTimePerTrial_Error(param,:)=(std(dstrLookTime((goodLearners()==1),:))*scale_factor/1000)/ sqrt(length(dstrLookTime((goodLearners()==1))));

%%% fig 32
test_WTargetLookingTimePerTrial_Mean(param,:)=mean(targLookTime((goodLearners()==0),:))*scale_factor/1000;
test_WTargetLookingTimePerTrial_Error(param,:)=(std(targLookTime((goodLearners()==0),:))*scale_factor/1000)/ sqrt(length(targLookTime((goodLearners()==0))));
test_WDistractorLookingTimePerTrial_Mean(param,:)=mean(dstrLookTime((goodLearners()==0),:))*scale_factor/1000;
test_WDistractorLookingTimePerTrial_Error(param,:)=(std(dstrLookTime((goodLearners()==0),:))*scale_factor/1000)/ sqrt(length(dstrLookTime((goodLearners()==0))));


%%%%Fig 4
test_ProportionStrong_Mean(param)= sum(goodLearners)/numSubjects;
test_ProportionStrong_Error(param)=0;
test_ProportionWeak_Mean(param)=(numSubjects-sum(goodLearners))/numSubjects;
test_ProportionWeak_Error(param)=0;
%sqrt((((sum(goodLearners)/numSubjects) - (12/18))^2 + (((numSubjects-sum(goodLearners))/numSubjects) - (6/18))^2)./2)

%%%%Fig 5
test_WordsLearnt_Mean(param)=mean(sum(LearntWords(),2));
test_WordsLearnt_Error(param)=0;


word_Vec=zeros(1,TEST_DUR/scale_factor);
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
        wordoff_time=savestate_historyOt(tt,:); sync_time(subject,tt)= sum(wordoff_time.*word_Vec);      
        %%% measure looking ratio at every time instant 
        lookingOn(subject,tt,:)= ((round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:)))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
    end

    % proportion of subjects looking to learnt vs unlearnt words
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
        
%    %%find total duration of looking to learnt vs unlearnt words
%     for kk=1:nObjects
%        if (LearntWords(subject,kk)==1)%if the subject has learnt the kk word.. learnt is measured on looking basis.. change to mem_matrix basis if needed
%                 for tt=1:nTestTrials%each test trial 
%                     if(xsit_result.train(subject).Words{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial name1
%                         s1=char(xsit_result.test(subject).test_pair(tt,2*nFeatures+2));
%                         if (strcmp(s1,'L'))% on left
%                             tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyLt(tt,:));%add time historyLt
%                         elseif (strcmp(s1,'R'))
%                             tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyRt(tt,:));% add time historyRt
%                         end
%                     end
%                 end        
%         elseif (LearntWords(subject,kk)==0)%if word is NOT learnt
%                 for tt=1:nTestTrials
%                     if(xsit_result.train(subject).Words{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial
%                         s1=char(xsit_result.test(subject).test_pair(tt,2*nFeatures+2));
%                         if (strcmp(s1,'L'))% on left
%                             tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyLt(tt,:));%add time historyL
%                         elseif (strcmp(s1,'R'))
%                             tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyRt(tt,:));% add time historyR
%                         end
%                     end
%                 end
%         else
%             disp('Error in LearntWords array: supurious data ');
%         end
%     end
end

%% fig 51
test_sync_time_Mean(param)=mean(mean(sync_time,2));
test_sync_time_Error(param)=std(mean(sync_time,2))./sqrt(length(sync_time));

% 
% figure (109);%Plot duration of looking to learnt vs unlearnt words
% plot((tLookAtLearnt/numSubjects)*scale_factor/1000); 
% hold on
% plot((tLookAtNonLrnt/numSubjects)*scale_factor/1000); 
% xlabel('per test trial');
% ylabel('duration of looking to words');
% legend('Learned words','Non-learned words');
% ylim([0 4]);


all_objLook= NaN (numSubjects,nTestTrials,TEST_DUR/scale_factor);
all_targL= NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
all_distL= NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
larW=NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
ularW=NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
for subject=1:numSubjects
        all_objLook(subject,:,:) = lookingOn(subject,:,:);
        for kk=1:6
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

%%% fig 60 %% Temporal Analysis Plots
test_PropLookAtObjs_Mean(param,:)=squeeze(nanmean(nanmean(all_objLook,1),2));
test_PropLookAtTarg_Mean(param,:)=squeeze(nanmean(nanmean(all_targL,1),2));
test_PropLookAtDist_Mean(param,:)=squeeze(nanmean(nanmean(all_distL,1),2));

%%% fig 61 %% Temporal Analysis Plots
test_SPropLookAtObjs_Mean(param,:)=squeeze(nanmean(nanmean(all_objLook((goodLearners()==1),:,:),1),2));
test_SPropLookAtTarg_Mean(param,:)=squeeze(nanmean(nanmean(all_targL((goodLearners()==1),:,:),1),2));
test_SPropLookAtDist_Mean(param,:)=squeeze(nanmean(nanmean(all_distL((goodLearners()==1),:,:),1),2));

%%% fig 62 %% Temporal Analysis Plots
test_WPropLookAtObjs_Mean(param,:)=squeeze(nanmean(nanmean(all_objLook((goodLearners()==0),:,:),1),2));
test_WPropLookAtTarg_Mean(param,:)=squeeze(nanmean(nanmean(all_targL((goodLearners()==0),:,:),1),2));
test_WPropLookAtDist_Mean(param,:)=squeeze(nanmean(nanmean(all_distL((goodLearners()==0),:,:),1),2));


xx=['Avg Looking time per test trial is ',num2str(mean(mean(targLookTime+dstrLookTime))*(scale_factor/1000))]; disp(xx);
xx=['Looking time to Target per test trial  is ',num2str(mean(mean(targLookTime))*(scale_factor/1000))]; disp(xx);
xx=['Looking time to Distractor per test trial per is ',num2str(mean(mean(dstrLookTime))*(scale_factor/1000))]; disp(xx);
xx=['Proportion of time looking correctly (Target/Total) is ',num2str(mean(sum (correct_proportion,2)))]; disp(xx);
xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);% 
xx=['Avg # of Words Learnt by Strong Learners is ',num2str(mean(sum(LearntWords(goodLearners()==1,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt by Weak Learners is ',num2str(mean(sum(LearntWords(goodLearners()==0,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt is ',num2str(mean(sum(LearntWords(),2)))]; disp(xx);%
% stanard deviation std(sum(LearntWords(goodLearners()==1,:),2))

%% TRAINING CONDITION ANALYSIS trials analysis NEW
if strcmp(TASK, 'XSIT')
    word_On = floor([500 2000]/scale_factor);%XSIT            
    word_Off = floor([1500 3000 ]/scale_factor);word_Len=floor(1000/scale_factor);
elseif strcmp(TASK, 'XTRAP')
    word_On = floor([368 1850]/scale_factor);  %% XTRAP    
    word_Off = floor((745+[368 1850])/scale_factor);word_Len=floor(745/scale_factor);
end
vis_On=1;vis_Off=(TRAIN_DUR/scale_factor); nFix_limit=10;

for subject=1:numSubjects
    savestate_historyL = xsit_result.train(subject).historyL(:,vis_On:vis_Off);
    savestate_historyR = xsit_result.train(subject).historyR(:,vis_On:vis_Off);    
    nHasLWord=zeros(nTrainTrials,1); tLookAtLearnt=zeros(nTrainTrials,1);
    nHasUWord=zeros(nTrainTrials,1); tLookAtNonLrnt=zeros(nTrainTrials,1);
    for kk=1:nObjects
        if (LearntWords(subject,kk)==1)%if the subject has learnt the kk word      
            for tt=1:nTrainTrials
                if(xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+1))%object exists in trial name1
                    nHasLWord(tt)=1;
                    s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                    if (strcmp(s1,'P'))% on left
                        tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
                    elseif (strcmp(s1,'X'))
                        tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
                    end
                elseif (xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+2))%object exists in trial
                    nHasLWord(tt)=1;
                    s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                    if (strcmp(s1,'X'))% on left
                         tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
                    elseif (strcmp(s1,'P'))
                        tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
                    end
                end
            end
           
        elseif (LearntWords(subject,kk)==0)%if word is NOT learnt     
            for tt=1:nTrainTrials
                if(xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+1))%object exists in trial
                    nHasUWord(tt)=1;
                    s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                    if (strcmp(s1,'P'))% on left
                        tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
                    elseif (strcmp(s1,'X'))
                        tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
                    end
                elseif (xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+2))%object exists in trial
                    nHasUWord(tt)=1;
                    s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                    if (strcmp(s1,'X'))% on left
                         tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
                    elseif (strcmp(s1,'P'))
                        tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
                    end
                end
            end
        else
            disp('Error in LearntWords array: supurious data ');
        end
     end
 
    for tt=1:nTrainTrials
        if nHasLWord(tt)==0; tLookAtLearnt(tt)=NaN;end % take only those trials wherein the word existed
        if nHasUWord(tt)==0; tLookAtNonLrnt(tt)=NaN;end
        tLookAtLearntWords(subject,tt)=tLookAtLearnt(tt);
        tLookAtNonLrntWords(subject,tt)=tLookAtNonLrnt(tt);
    end
end

train_Slook_to_Learnt_words_Mean(param,:)= nanmean(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000;
train_Slook_to_Learnt_words_Error(param,:)=(nanstd(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==1,:)));
train_Wlook_to_Learnt_words_Mean(param,:)= nanmean(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000;
train_Wlook_to_Learnt_words_Mean(param,:)= (nanstd(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==0,:)));


train_Slook_to_NonLearnt_words_Mean(param,:)= nanmean(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000;
train_Slook_to_NonLearnt_words_Error(param,:)=(nanstd(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==1,:)));
train_Wlook_to_NonLearnt_words_Mean(param,:)= nanmean(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000;
train_Wlook_to_NonLearnt_words_Mean(param,:)= (nanstd(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==0,:)));

% figure(8);%Plot looking to learnt vs unlearnt 
% errorbar(,,plotStyle{1});%
% hold on
% errorbar(nanmean(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000,(nanstd(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==1,:))),plotStyle{2+compStyle});%
% xlabel('per training trial');
% ylabel('Looking time (s)');
% legend('Learned Words','Non-Learned Words');
% title('Strong Learners')
% 
% figure(9);%Plot looking to learnt vs unlearnt 
% errorbar(nanmean(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000,(nanstd(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==0,:))),plotStyle{1});%
% hold on
% errorbar(nanmean(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000,(nanstd(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==0,:))),plotStyle{2+compStyle});%
% xlabel('per training trial');
% ylabel('Looking time (s)');
% legend('Learned Words','Non-Learned Words');
% title('Weak Learners')
% %%% above nan lengths are NOT CORRECT... length(x(~isnan(x)))
%% 
targLookTimeTraining=zeros(numSubjects,nTrainTrials);%number of training trials
dstrLookTimeTraining=zeros(numSubjects,nTrainTrials);
totnlooks=zeros(numSubjects,nTrainTrials);
meanlookdur =zeros(numSubjects,nTrainTrials);
TotalLookTime=zeros(numSubjects,nTrainTrials);
totlonglookdur = zeros(numSubjects,nTrainTrials);
tLookRepeated= zeros(numSubjects,nObjects);
tLookVarying=zeros(numSubjects,nObjects);
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
           targLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,word_On(1):word_Off(1))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyR(tr,word_On(2):word_Off(2)));% 2nd audio presentation, Looking to target object right side
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,word_On(1):word_Off(1))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyL(tr,word_On(2):word_Off(2)));% 2nd audio presentation, Looking wrong way
        elseif (strcmp(s1,'X'))
           targLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,word_On(1):word_Off(1))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyL(tr,word_On(2):word_Off(2)));% 2nd audio presentation, Looking to target object right side
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,word_On(1):word_Off(1))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyR(tr,word_On(2):word_Off(2)));% 2nd audio presentation, Looking wrong way            
        end
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

%%%%Fig 11
train_TotalLookingTime_Mean(param)=mean(mean([TotalLookTime(goodLearners()==1,:); TotalLookTime(goodLearners()==0,:)])) *(scale_factor/1000);
train_TotalLookingTime_Error(param)=(std(mean([TotalLookTime(goodLearners()==1,:); TotalLookTime(goodLearners()==0,:)])) *(scale_factor/1000))/ sqrt(length([TotalLookTime(goodLearners()==1,:); TotalLookTime(goodLearners()==0,:)]));

%sqrt(((mean(mean([TotalLookTimeS; TotalLookTimeW])) *(scale_factor/1000) - 3.04)^2 + (mean(mean([TotalLookTimeS; TotalLookTimeW])) *(scale_factor/1000) - 2.99)^2)./2)

%%% Fig 16
train_SLookEntropy_Mean(param,:)=   mean(EntropySub((goodLearners()==1),:));
train_WLookEntropy_Mean(param,:)=   mean(EntropySub((goodLearners()==0),:));
train_SLookEntropy_Error(param,:)=  std(EntropySub((goodLearners()==1),:))/sqrt(length(EntropySub((goodLearners()==1),:)));
train_WLookEntropy_Error(param,:)=  std(EntropySub((goodLearners()==0),:))/sqrt(length(EntropySub((goodLearners()==0),:)));


%%% Fig 107
train_SLookVariance_Mean(param,:)=   mean(VarianceSub((goodLearners()==1),:));
train_WLookVariance_Mean(param,:)=   mean(VarianceSub((goodLearners()==0),:));
train_SLookVariance_Error(param,:)=  std(VarianceSub((goodLearners()==1),:))/sqrt(length(VarianceSub((goodLearners()==1),:)));
train_WLookVariance_Error(param,:)=  std(VarianceSub((goodLearners()==0),:))/sqrt(length(VarianceSub((goodLearners()==0),:)));



%%% fig 7 %%  Plots
train_STotalLookTime_Mean(param,:)=mean(TotalLookTime(goodLearners()==1,:))*scale_factor/1000;
train_WTotalLookTime_Mean(param,:)=mean(TotalLookTime(goodLearners()==0,:))*scale_factor/1000;

%%% fig 8 %% 
train_Stotnlooks_Mean(param,:)=mean(totnlooks(goodLearners()==1,:));
train_Wtotnlooks_Mean(param,:)=mean(totnlooks(goodLearners()==0,:));


summaF=mean(mean(totnlooks(goodLearners()==1,:)));xx=['number of fixations/looks Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(totnlooks(goodLearners()==0,:)));xx=['number of fixations/looks Weak learners ',num2str(summaF)]; disp(xx);

%%% fig 9 %% 
train_Smeanlookdur_Mean(param,:)=mean(meanlookdur(goodLearners()==1,:))*scale_factor/1000;
train_Wmeanlookdur_Mean(param,:)=mean(meanlookdur(goodLearners()==0,:))*scale_factor/1000;

summaF=mean(mean(meanlookdur(goodLearners()==1,:))*scale_factor)/1000;xx=['mean look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanlookdur(goodLearners()==0,:))*scale_factor)/1000;xx=['mean look duration Strong weak learners ',num2str(summaF)]; disp(xx);

%%% fig 10 %% 
train_Stotlonglookdur_Mean(param,:)=mean(totlonglookdur(goodLearners()==1,:))*scale_factor/1000;
train_Wtotlonglookdur_Mean(param,:)=mean(totlonglookdur(goodLearners()==0,:))*scale_factor/1000;

summaF=mean(mean(totlonglookdur(goodLearners()==1,:))*scale_factor)/1000;xx=['Longest look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totlonglookdur(goodLearners()==0,:))*scale_factor)/1000;xx=['Longest look duration Strong weak learners ',num2str(summaF)]; disp(xx);


summaF=mean(mean([TotalLookTime(goodLearners()==1,:); TotalLookTime(goodLearners()==0,:)])) *(scale_factor/1000);xx=['Avg looking time per training trial is ',num2str(summaF)]; disp(xx);
var1= mean(mean(TotalLookTime(goodLearners()==1,:)))*scale_factor/1000;xx=['Avg looking time per training trial per Strong learner is ',num2str(var1)]; disp(xx);
var2= mean(mean(TotalLookTime(goodLearners()==0,:)))*scale_factor/1000;xx=['Avg looking time per training trial per Weak learner is ',num2str(var2)]; disp(xx);


%% trace analysis
Wrong_inTrace= zeros(numSubjects,1); Correct_inTrace= zeros(numSubjects,1);
count_dominant_word_traces=zeros(numSubjects,1);
for subject=1:numSubjects    
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
        a_cv=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cv=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
        C_inTr= C_inTr+ mean([mean(a_cv(a_cv>0.001))  mean(b_cv(b_cv>0.001))]);
        max_corr_strength = mean ([max(a_cv) max(b_cv)]);
        inputMapping1(xx1(kk),yy(kk))=0;
        inputMapping2(xx2(kk),yy(kk))=0;
        for jj=1:6
            inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
            inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
        end
        a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
        W_inTr = W_inTr + mean([mean(a_in(a_in>0.001)) mean(b_in(b_in>0.001))]);
        max_incorr_strength = mean ([max(a_in) max(b_in)]);
        if (max_corr_strength > max_incorr_strength )
            count_dominant_word_traces(subject) = count_dominant_word_traces(subject) + 1;
        end
    end
    Correct_inTrace(subject)=C_inTr/nObjects;
    Wrong_inTrace(subject)=W_inTr/nObjects;
    InCorr_assocs(subject)=mean([as_count1-1 as_count2-1])/nObjects;
    EntropyTrace(subject)= mean( [entropy(inputMapping1) entropy(inputMapping2)] ); 
end

%%% Fig 141
train_InCorr_assocs_Mean(param)=   mean(InCorr_assocs);
train_InCorr_assocs_Error(param)=  std(InCorr_assocs)/sqrt(length(InCorr_assocs));
train_SInCorr_assocs_Mean(param)=   mean(InCorr_assocs((goodLearners()==1)));
train_WInCorr_assocs_Mean(param)=   mean(InCorr_assocs((goodLearners()==0)));
train_SInCorr_assocs_Error(param)=  std(InCorr_assocs((goodLearners()==1)))/sqrt(length(InCorr_assocs((goodLearners()==1))));
train_WInCorr_assocs_Error(param)=  std(InCorr_assocs((goodLearners()==0)))/sqrt(length(InCorr_assocs((goodLearners()==0))));


%%% Fig 142
train_count_dominant_word_traces_Mean(param)=   mean(count_dominant_word_traces);
train_count_dominant_word_traces_Error(param)=   std(count_dominant_word_traces)/sqrt(length(count_dominant_word_traces));

%%% Fig 15
train_TraceEntropy_Mean(param)=   mean(EntropyTrace);
train_TraceEntropy_Error(param)=  std(EntropyTrace)/sqrt(length(EntropyTrace));
train_STraceEntropy_Mean(param)=   mean(EntropyTrace((goodLearners()==1)));
train_WTraceEntropy_Mean(param)=   mean(EntropyTrace((goodLearners()==0)));
train_STraceEntropy_Error(param)=  std(EntropyTrace((goodLearners()==1)))/sqrt(length(EntropyTrace((goodLearners()==1))));
train_WTraceEntropy_Error(param)=  std(EntropyTrace((goodLearners()==0)))/sqrt(length(EntropyTrace((goodLearners()==0))));

%%fig 17
train_Correct_inTrace_Strength_Mean(param)= nanmean(Correct_inTrace);
train_Correct_inTrace_Strength_Error(param)=  nanstd(Correct_inTrace)/sqrt(length(Correct_inTrace));

train_Wrong_inTrace_Strength_Mean(param)= nanmean(Wrong_inTrace);
train_Wrong_inTrace_Strength_Error(param)=  nanstd(Wrong_inTrace)/sqrt(length(Wrong_inTrace));


%%%%%%% Association hwf trace analysis
corrAsocn=zeros(numSubjects,nObjects);
cS=1;cW=1;
for subject=1:numSubjects  
    
    AsocMat=squeeze(xsit_result.train(subject).hwf(1,:,:));% 1 for first feature (shape) only
    inputMapping=zeros(size(AsocMat));
     for kk=1:nObjects
         inputMapping(cell2mat(xsit_result.train(subject).Feature1(kk)),cell2mat(xsit_result.train(subject).Words(kk)))=1;
     end    
    for kk=1:nObjects        
        temp=[];
        temp=AsocMat(:,cell2mat(xsit_result.train(subject).Words(kk)));  
        maxAsocnVal(subject,kk) = max(temp);
        [temp2 in2] = max(temp); 
        NxtmaxAsocnVal(subject,kk)= max(max (temp(1:max(1,in2-5))), max (temp(min(size(temp),in2+5):size(temp))));
        %NxtmaxAsocnVal(subject,kk) = max(temp(temp<max(temp)));
        ratioMax(subject,kk)= maxAsocnVal(subject,kk)./NxtmaxAsocnVal(subject,kk);
        prodtMR(subject,kk)=ratioMax(subject,kk).*maxAsocnVal(subject,kk);
        
        
        [maxIn(kk) indIn(kk)] = max(inputMapping(:,cell2mat(xsit_result.train(subject).Words(kk))));
        [maxAs(kk) indAs(kk)] = max(AsocMat(:,cell2mat(xsit_result.train(subject).Words(kk))));       
        if (abs(indIn(kk)-indAs(kk)) <= 2)%if association is correct i..e same as input?
           corrAsocn(subject, kk)=1; % wrongAssocn = 6-corrAsocn
        end
    end 
end

% 
SLer=[];WLer=[];SNon=[];WNon=[];SLer2=[];WLer2=[];SNon2=[];WNon2=[]; %totM_L=[];totM_N=[];
for subject=1:numSubjects 
    %totM_S=[totM_S maxAsocnVal(subject,LearntWords(subject,:)==1)];
    %totM_W=[totM_W maxAsocnVal(subject,LearntWords(subject,:)==0)];
     
    if(goodLearners(subject)==1)
        SLer=[SLer maxAsocnVal(subject,LearntWords(subject,:)==1)];
        SNon=[SNon maxAsocnVal(subject,LearntWords(subject,:)==0)];
        
        SLer2=[SLer2 ratioMax(subject,LearntWords(subject,:)==1)];
        SNon2=[SNon2 ratioMax(subject,LearntWords(subject,:)==0)];
     elseif (goodLearners(subject)==0)
         WLer=[WLer maxAsocnVal(subject,LearntWords(subject,:)==1)];
         WNon=[WNon maxAsocnVal(subject,LearntWords(subject,:)==0)];
        
         WLer2=[WLer2 ratioMax(subject,LearntWords(subject,:)==1)];
         WNon2=[WNon2 ratioMax(subject,LearntWords(subject,:)==0)];
    end
end
if ((size(SLer,2)+size(WLer,2)+size(SNon,2)+size(WNon,2))./numSubjects ~= nObjects), disp('ERROR ERROR ERROR ERROR'), end

%%% fig 13 %% 
totML_Mean(param)=mean([SLer WLer]);
totMN_Mean(param)=mean([SNon WNon]);
totML_Error(param)=(std(mean([SLer WLer])))/ sqrt(length([SLer WLer]));
totMN_Error(param)=(std(mean([SNon WNon])))/ sqrt(length([SNon WNon]));


totMS_Mean(param)=mean([SLer SNon]);
totMW_Mean(param)=mean([WLer WNon]);
totMS_Error(param)=(std(mean([SLer SNon])))/ sqrt(length([SLer SNon]));
totMW_Error(param)=(std(mean([WLer WNon])))/ sqrt(length([WLer WNon]));


 varb=mean(SLer); xx=['Avg association strength for Learnt words in Strong learners ',num2str(varb)]; disp(xx);
 varb=mean(WLer); xx=['Avg association strength for Learnt words in Weak learners ',num2str(varb)]; disp(xx);
 varb=mean(SNon); xx=['Avg association strength for NonLearnt words in Strong learners ',num2str(varb)]; disp(xx);
 varb=mean(WNon); xx=['Avg association strength for NonLearnt words in Weak learners ',num2str(varb)]; disp(xx);
  varb=mean(SLer2); xx=['Avg Ratio of 2Maximums for Learnt words in Strong learners ',num2str(varb)]; disp(xx);
 varb=mean(WLer2); xx=['Avg Ratio of 2Maximums for Learnt words in Weak learners ',num2str(varb)]; disp(xx);
 varb=mean(SNon2); xx=['Avg Ratio of 2Maximums for NonLearnt words in Strong learners ',num2str(varb)]; disp(xx);
 varb=mean(WNon2); xx=['Avg Ratio of 2Maximums for NonLearnt words in Weak learners ',num2str(varb)]; disp(xx);
% %wordsAssoc=sum(corrAsocn,2);
 varb=mean(sum(corrAsocn,2))/nObjects; xx=['Avg proportion of correctly associated words in Memory ',num2str(varb)]; disp(xx);
varb=mean(sum(corrAsocn((goodLearners()==1),:),2)); xx=['Avg # of correctly associated words for Strong in Memory ',num2str(varb)]; disp(xx);
varb=mean(sum(corrAsocn((goodLearners()==0),:),2)); xx=['Avg # of correctly associated words for Weak in Memory ',num2str(varb)]; disp(xx);
% %Learnt non learnt by strong weak 
% varb=mean(corrWordsStrong); xx=['Avg # of correctly associated words in Memory for Strong counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% std(corrWordsStrong);
% varb=mean(corrWordsWeak); xx=['Avg # of correctly associated words in Memory for Weak counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% std(corrWordsWeak);
% varb=sum((corrAsocn(:,1:6).*LearntWords(:,1:6)));
% xx=['# of subjects with correctly associated word in Memory counted as Learnt thru looking for ',num2str(numSubjects),' subjects is ' num2str(varb)]; disp(xx);

% % % for subject=1:numSubjects
% % % 
% % % inputMapping1=zeros(306,20);
% % % inputMapping2=zeros(306,20);
% % % muddled_pairs=0;
% % % muddled_mem_val=1;
% % %     for kk=1:nObjects
% % %         xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
% % %         xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
% % %         yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
% % %     end
% % %     for kk=1:nObjects       
% % %         inputMapping1(xx1(kk),yy(kk))=1;
% % %         inputMapping2(xx2(kk),yy(kk))=1;
% % %         for jj=1:8
% % %             inputMapping1(xx1(kk)+jj-4,yy(kk))=1;
% % %             inputMapping2(xx2(kk)+jj-4,yy(kk))=1;
% % %         end   
% % %     end
% % %     
% % %     figure(21)
% % %     
% % %     lsp=subplot(1,2,1);
% % %     %surface(inputMapping1)
% % %     %hold on
% % %     [mA, iA] = max(squeeze(xsit_result.train(subject).hwf(1,:,:)));
% % %     surface(squeeze(xsit_result.train(subject).hwf(1,:,:)));
% % %     title([num2str(mA)]);
% % %     shading flat; 
% % %     
% % %     rsp= subplot(1,2,2);
% % %     hold on
% % %     
% % %     title('Input');
% % %     
% % %     surface(squeeze(xsit_result.test(subject).hwft(1,:,:)));
% % %     title('Test');
% % %     shading flat;
% % %     
% % %     if (goodLearners(subject)==1), suptitle(['Strong # ' num2str(subject)]),
% % %     elseif (goodLearners(subject)==0), suptitle(['Weak # ' num2str(subject)]), end
% % %     pause(7);
% % %     clf;
% % % end



% % % lookR = savestate_historyR';
% % % lookL = savestate_historyL';
% % % vecL = lookL(:);
% % % vecR = lookR(:);
% % % c = [vecL, vecR];
% % % size(c)
% % % area(c)


% % xA=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
% % %yA=[22/72 9/72 15/72 8/72 28/72 52/72 66/72 71/72 72/72 72/72 72/72];
% % %yA =[3.68 3.3 3.3 3.25 3.9 3.9 4.04 4.74 5.2 5.4 5.5]
% % %yA =[2.02 1.7 1.7 2.1 2.3 2.85 3.16 3 0 0 0]
% % yA= [2.57 1.9 2.08 2.18 2.98 3.6 3.97 4.7 5.2 5.4 5.5]
% % plot (xA,yA)
% % xlabel('memory trace strength');
% % ylabel('Avg words learnt');
% % %ylabel('proportion of Strong Learners');
% % grid on
end














