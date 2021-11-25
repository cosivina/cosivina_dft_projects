%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis 
clear all; close all; % This script generates an output file
%Global Variables
scale_factor=8;nFeatures=2;MIN_LOOK_DURATION=200/scale_factor; 
%Plotting Variables
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);

%% Raw Data File Name             
simName = 'Test2_Smith_Yu_2008_2011_results'
xsit_result = load ([simName '.mat']); % loading the file
numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 
%xsit_result.sim.saveSettings('tesw7.json'); %visdiff ('tesw2.json', 'tesw7.json')

%Experiment (Task) Variables
TASK= 'XSIT';% 'XTRAP';% 
new_sim_formatting = true; % new format sims with task_variables saved at simulation

if (new_sim_formatting == true)
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrainTrials = xsit_result.train(1).maxTR; 
    nTestTrials = xsit_result.train(1).maxTRt;
    TRAIN_DUR = xsit_result.train(1).t_max - floor(1000/scale_factor);
    TEST_DUR = xsit_result.train(1).t_maxt - floor(1000/scale_factor);
    word_On = xsit_result.test(1).word_On; word_Off = xsit_result.test(1).word_Off; 
    word_Len=floor(1000/scale_factor);
    vis_On = xsit_result.test(1).visuals_On;
    vis_Off = xsit_result.test(1).visuals_Off;
    if length(word_On) > 1, words_within_range = length(word_On)-1; else, words_within_range = 1; end
    if length(word_On) > 1, word_effect_Len = floor(1750/scale_factor); else, word_effect_Len=floor(999/scale_factor); end % to map to 1.75 seconds on paper
else
    nObjects=6; nTrainTrials=30; nTestTrials=12;  
    TRAIN_DUR= floor(4000/scale_factor); TEST_DUR= floor(8000/scale_factor);
    
    %%% XSIT test
    word_On = floor([500 2000 4500 6000]/scale_factor);    
    word_Off = floor([1500 3000 5500 7000]/scale_factor);
    word_Len=floor(1000/scale_factor);
    
%%% 1-sec test
     %TEST_DUR= floor(1000/scale_factor);
     %word_On = 1;  %% 1 second test 
     %word_Off = floor(1000/scale_factor);

%%% XTRAP test
    %word_On = floor([8 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
    %word_Off = floor((745+[8 1800 3500 5200 6900])/scale_factor);
    %word_Len=floor(745/scale_factor);
    
    vis_On = 1;vis_Off = TEST_DUR;
    if length(word_On) > 1, words_within_range = length(word_On)-1; else, words_within_range = 1; end
    if length(word_On) > 1, word_effect_Len = floor(1750/scale_factor); else, word_effect_Len=floor(999/scale_factor); end % to map to 1.75 seconds on paper

end

Look_Smoothening_SY_2AFC; % DATA SMOOTHENING if necessary 
        
%% TEST CONDITION ANALYSIS
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
measurement_i = 'Mean looking time per 8s test trial';
empirical_mean_i = [6.10 5.92];
mean_i = mean(mean(targLookTime+dstrLookTime))*(scale_factor/1000);
SE_i = SE(mean(targLookTime+dstrLookTime))*(scale_factor/1000);
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%%
measurement_i = 'Preferential looking time ratio (at test)';
empirical_mean_i = [0.60 0.54];
mean_i = nanmean(nanmean(correct_proportion,2));
SE_i = SE(nanmean(correct_proportion,2));
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%%
measurement_i = 'Mean words learnt at test (out of 6)';
empirical_mean_i = [4.0 3.5];
mean_i = mean(sum(LearntWords(),2));
SE_i = SE(sum(LearntWords(),2));
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%%
% The below measure should be averaged over 100 such raw data files /simulations/runs
% measurement_i = 'Proportion of Strong Learners (SE) at test';
% empirical_mean_i = [0.67];
 mean_i = sum(goodLearners)/numSubjects; % actually should be mean over this measure
% SE_i = NaN;
% RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
% row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
 xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
%xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%%
measurement_i = 'Mean looking time to target per test trial';
empirical_mean_i = [3.6 3.25];
mean_i = mean(mean(targLookTime))*(scale_factor/1000);
SE_i = SE(mean(targLookTime))*(scale_factor/1000);
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%%
measurement_i = 'Mean looking time to distractor per test trial';
empirical_mean_i = [2.5 2.67];
mean_i = mean(mean(dstrLookTime))*(scale_factor/1000);
SE_i = SE(mean(dstrLookTime))*(scale_factor/1000);
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


%% other useful measures
xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);% 
xx=['Looking time to Target per by strong test trial  is ',num2str(mean(mean(targLookTime(goodLearners()==1,:)))*(scale_factor/1000))]; disp(xx);
xx=['Looking time to Distractor by strong per test trial per is ',num2str(mean(mean(dstrLookTime(goodLearners()==1,:)))*(scale_factor/1000))]; disp(xx);
xx=['Looking time to Target per by Weak test trial  is ',num2str(mean(mean(targLookTime(goodLearners()==0,:)))*(scale_factor/1000))]; disp(xx);
xx=['Looking time to Distractor by Weak per test trial per is ',num2str(mean(mean(dstrLookTime(goodLearners()==0,:)))*(scale_factor/1000))]; disp(xx);
xx=['Avg # of Words Learnt by Strong Learners is ',num2str(mean(sum(LearntWords(goodLearners()==1,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt by Weak Learners is ',num2str(mean(sum(LearntWords(goodLearners()==0,:),2)))]; disp(xx);%  
xx=['Avg Looking time per test trial by Strong is ',num2str(mean(mean(targLookTime(goodLearners()==1,:)+dstrLookTime(goodLearners()==1,:)))*(scale_factor/1000))]; disp(xx);
xx=['Avg Looking time per test trial by Weak is ',num2str(mean(mean(targLookTime(goodLearners()==0,:)+dstrLookTime(goodLearners()==0,:)))*(scale_factor/1000))]; disp(xx);


%% Plot total looking time during test trial
barColorMap(1,:) = [.2 .71 .3];	% Green Color for segment 1.
barColorMap(2,:) = [.25 .55 .79];	% Blue Color for segment 2.
barColorMap(3,:) = [.9 .1 .14];	% Red Color for segment 3.
barColorMap(4,:) = [.9 .9 .14];	% Yellow Color for segment 4.
   
figure(1)
sts2 = [6.1;     5.92; (mean(mean(targLookTime+dstrLookTime))*(scale_factor))/1000];
errY2=[0.05;     0.05; SE(mean(targLookTime+dstrLookTime,2))*(scale_factor/1000)];
ylabel('Looking time (seconds)');
ylim([0 8]);
hold on;
barFontSize = 11;
for b = 1 : length(sts2)
	% Plot one single bar as a separate bar series.
	handleToThisBarSeries(b) = bar(b, sts2(b), 'BarWidth', 0.6);
	% Apply the color to this bar series.
	set(handleToThisBarSeries(b), 'FaceColor', barColorMap(b,:));
	% Place text atop the bar
	barTopper = sprintf('%.3f', sts2(b));
	text(b-0.1, sts2(b)+0.3, barTopper, 'FontSize', barFontSize);
	hold on;
end
legend('Smith & Yu 2008 (14 m)', 'Yu & Smith 2011','WOLVES Model');



%% Plot Target vs Distractor looking time during test trial
mT=(mean(mean(targLookTime,2))*(scale_factor))/1000;mD=(mean(mean(dstrLookTime,2))*(scale_factor))/1000;
eT= SE(mean(targLookTime,2))*(scale_factor/1000);eD= SE(mean(dstrLookTime,2))*(scale_factor/1000);
figure(2)
blockNames={'Target'; 'Distractor'};
sts = [ 3.6/(3.6+2.5)    3.25/(3.25+2.67) mT/(mT+mD) ;  2.5/(3.6+2.5)  2.67/(3.25+2.67) mD/(mT+mD)];
errY =[ 0.2/(3.6+2.5)    0.49/(3.25+2.67) eT/(mT+mD) ; 0.25/(3.6+2.5)  0.38/(3.25+2.67) eD/(mT+mD)];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',18);
legend('Smith & Yu 2008 14 month olds', 'Yu & Smith 2011','WOLVES Model');
ylabel('Looking time(seconds)');
ylim([0 1]);


%% Within Trial Time Measures and Plots
word_Vec=zeros(1,TEST_DUR);
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

figure (3);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
for wc=1:length(word_On)
 rectangle('Position',[word_On(wc),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
end
plot((squeeze(nanmean(nanmean(all_objLook(:,:,:),1),2))),'LineWidth',3)
hold on
plot((squeeze(nanmean(nanmean(all_targL(:,:,:),1),2))),'LineWidth',3)
hold on
plot((squeeze(nanmean(nanmean(all_distL(:,:,:),1),2))),'LineWidth',3)
hold on
hline(0.5);
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Both Objects','Target Object','Distractor');
xlabel('time');
ylabel('proportion of looking');
set (gca, 'FontSize',16);
ylim([0 1]);


%% TRAINING CONDITION ANALYSIS trials analysis

if (new_sim_formatting == true)
    word_On = xsit_result.train(1).word1_On; 
    word_Off = xsit_result.train(1).word1_Off; 
    C_word_On = word_On+ floor(500/scale_factor); C_word_Off = floor([2000 4000]/scale_factor);
else

    if strcmp(TASK, 'XSIT')
        word_On = floor([500 2000]/scale_factor);%XSIT            
        word_Off = floor([1500 3000 ]/scale_factor);word_Len=floor(1000/scale_factor);
        C_word_On = word_On+ floor(500/scale_factor); C_word_Off = floor([2000 4000]/scale_factor);

    elseif strcmp(TASK, 'XTRAP')
        word_On = floor([368 1850]/scale_factor);  %% XTRAP    
        word_Off = floor((745+[368 1850])/scale_factor);word_Len=floor(745/scale_factor);
        C_word_On = word_On+ floor(500/scale_factor);C_word_Off = floor(([1850 4000])/scale_factor);
    end
end

vis_On=1;vis_Off=TRAIN_DUR; nFix_limit=10;


%%
empirical_Learners_Varying = [0.61 0.52 0.54 0.47 0.48 0.35];
empirical_nonLearners_Varying = [0.63 0.58 0.62 0.50 0.55 0.44];
empirical_Learners_Repeated = [0.27 0.32 0.26 0.31 0.26 0.39];
empirical_nonLearners_Repeated = [0.24 0.28 0.17 0.22 0.20 0.30];
rough_Errorbars_var_rep = [0.1 0.1 0.1 0.1 0.1 0.1];
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
        for j=1:TRAIN_DUR
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

%%
measurement_i = 'Mean looking time per 4s training trial';
empirical_mean_i = [3.04*0.67 + 2.99*0.33];
mean_i = mean(mean(TotalLookTime))*scale_factor/1000;
SE_i =  SE(mean(TotalLookTime))*scale_factor/1000;
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%%
measurement_i = 'number of fixations(looks) per training trial';
empirical_mean_i = [2.75*0.67 + 3.82*0.33];
mean_i = mean(mean(totnlooks));
SE_i =  SE(mean(totnlooks));
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


%%
measurement_i = 'Mean fixation (longest look) duration';
empirical_mean_i = [1.69*0.67 + 1.21*0.33];
mean_i = mean(mean(totlonglookdur)*scale_factor)/1000;
SE_i =  SE(mean(totlonglookdur)*scale_factor)/1000;
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%% Other training data measures
summaF=mean(mean(VarianceSub((goodLearners()==1),:)));
xx=['Variance Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(VarianceSub((goodLearners()==0),:)));
xx=['Variance Weak learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(EntropySub((goodLearners()==1),:)));
xx=['Entropy Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(EntropySub((goodLearners()==0),:)));
xx=['Entropy Weak learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totnlooks(goodLearners()==1,:)));
xx=['number of fixations/looks Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(totnlooks(goodLearners()==0,:)));
xx=['number of fixations/looks Weak learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanLukhadur(goodLearners()==1,:))*scale_factor)/1000;
xx=['mean look duration Strong learners (indiv calcs) ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanLukhadur(goodLearners()==0,:))*scale_factor)/1000;
xx=['mean look duration Weak learners (indiv calcs) ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totlonglookdur(goodLearners()==1,:))*scale_factor)/1000;
xx=['Longest look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totlonglookdur(goodLearners()==0,:))*scale_factor)/1000;
xx=['Longest look duration Weak learners ',num2str(summaF)]; disp(xx);
    

%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])