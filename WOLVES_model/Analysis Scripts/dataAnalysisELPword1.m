%% Data Analysis File CHANGE MADE...  : = 1:word_Off(1)
clear all; close all; %% 
nObjects=6;nFeatures=2;nTrainTrials=30; nTestTrials=12; scale_factor=8; TRAIN_DUR=4000; TEST_DUR=8000;MIN_LOOK_DURATION=160/scale_factor; %check with auto file 160
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];

%simName = 'muddled_mem_6_pairs_sy_base05-Oct-2018_test_only_learners_smithYu_'%unified_param_run8_02-Oct-2018_Smith_Yu_2013_'%'bbc_H50F30Asc50_Tzero_Smith_Yu_2013_'; %test4_Smith_Yu_2008_
simName = 'wfChanges_conwmf15_hwf20_Smith_Yu_2008_';OutName = [simName,'results.mat']
xsit_result = load (OutName);

% add another simulation
% simName = 'goodlooking1-07-Nov-18_Smith_Yu_2013_';
% OutName1 = [simName,'results.mat']
% xsit_result1 = load (OutName1);
% xsit_result.train = [xsit_result.train ; xsit_result1.train];  
% xsit_result.test = [xsit_result.test ; xsit_result1.test];  
% %% 

%xsit_result.sim.saveSettings('tesw26.json'); %visdiff ('tesw2.json', 'tesw26.json')
numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 

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
word_On = floor([500 2000 4500 6000]/scale_factor);%XSIT    
word_Off = floor([1500 3000 5500 7000]/scale_factor);word_Len=floor(1000/scale_factor);
% word_On = floor([8 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
% word_Off = floor((745+[8 1800 3500 5200 6900])/scale_factor);word_Len=floor(745/scale_factor);
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
        lLook= sum( xsit_result.test(subject).historyLt(trt,vis_On:word_Off(1)));%full trial
        rLook= sum( xsit_result.test(subject).historyRt(trt,vis_On:word_Off(1)));%
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

barColorMap(1,:) = [.2 .71 .3];	% Green Color for segment 1.
barColorMap(2,:) = [.25 .55 .79];	% Blue Color for segment 2.
barColorMap(3,:) = [.9 .1 .14];	% Red Color for segment 3.
barColorMap(4,:) = [.9 .9 .14];	% Yellow Color for segment 4.
   
figure(1)% Plot total looking time during test trial
sts = [6.1;     5.92; (mean(mean(targLookTime+dstrLookTime))*(scale_factor))/1000];
errY=[0.05;     0.05; ((std(mean(targLookTime+dstrLookTime,2))*(scale_factor))/1000)./sqrt(length(targLookTime))];
%hb=barwitherr(errY,sts,0.6);% Plot with errorbars
%title ('total looking time per test trial');
ylabel('Looking time (seconds)');
ylim([0 8]);
%set(gca,'xticklabel',,'fontsize',11);
%sqrt((((mean(mean(targLookTime+dstrLookTime))*(scale_factor*slice_wordOnset4))/1000 - 6.1)^2 + ((mean(mean(targLookTime+dstrLookTime))*(scale_factor*slice_wordOnset4))/1000 - 5.92)^2)./2)
hold on;
barFontSize = 11;
for b = 1 : length(sts)
	% Plot one single bar as a separate bar series.
	handleToThisBarSeries(b) = bar(b, sts(b), 'BarWidth', 0.6);
	% Apply the color to this bar series.
	set(handleToThisBarSeries(b), 'FaceColor', barColorMap(b,:));
	% Place text atop the bar
	barTopper = sprintf('%.3f', sts(b));
	text(b-0.1, sts(b)+0.3, barTopper, 'FontSize', barFontSize);
	hold on;
end
legend('Smith & Yu 2008 (14 m)', 'Yu & Smith 2011','WOLVES Model');
%errorbar(sts,errY,'o');


figure(101)% Plot total looking time during test trial
blockNames={'WOLVES Model XTRAP';'WOLVES in XSIT'};
sts = [(mean(mean(targLookTime+dstrLookTime))*(scale_factor))/1000;    5.89];
errY=[((std(mean(targLookTime+dstrLookTime,2))*(scale_factor))/1000)./sqrt(length(targLookTime));   0.1];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Avg Looking Time during a test trial');
ylabel('Looking time (ms)');
ylim([0 8]);
%legend('WOLVES Model','Smith & Yu 2008 14 month olds', 'Yu & Smith 2011');
grid on

figure(102)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ (mean(mean(targLookTime))*(scale_factor))/1000  3.637;  (mean(mean(dstrLookTime))*(scale_factor))/1000 2.559];
errY =[((std(mean(targLookTime,2))*(scale_factor))/1000)./sqrt(length(targLookTime))  0.33; ((std(mean(dstrLookTime,2))*(scale_factor))/1000)./sqrt(length(dstrLookTime))  0.32];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
legend('WOLVES XTRAP', 'WOLVES XSIT');
title ('Target vs Distractor Looking Time');
ylabel('Looking time per test trial');
ylim([0 4]);

figure(2)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ 3.6    3.25  (mean(mean(targLookTime,2))*(scale_factor))/1000;  2.5  2.67 (mean(mean(dstrLookTime,2))*(scale_factor))/1000];
errY =[ 0.2    0.49  ((std(mean(targLookTime,2))*(scale_factor))/1000)./sqrt(length(targLookTime)); 0.25  0.38 ((std(mean(dstrLookTime,2))*(scale_factor))/1000)./sqrt(length(dstrLookTime))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',12);
legend('Smith & Yu 2008 14 month olds', 'Yu & Smith 2011','WOLVES Model');
%title ('Target vs Distractor Looking Time');
ylabel('Looking time(seconds)');
ylim([0 4]);
%sqrt((((mean(mean(targLookTime))*(scale_factor*slice_wordOnset4))/1000 - 3.6)^2 + ((mean(mean(targLookTime))*(scale_factor*slice_wordOnset4))/1000 - 3.25)^2 + ((mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4))/1000 - 2.5)^2 + ((mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4))/1000 - 2.67)^2 )./4)

figure(201)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ (mean(mean(targLookTime))*(scale_factor))/1000    3.25;  (mean(mean(dstrLookTime))*(scale_factor))/1000  2.67];
errY =[(std(mean(targLookTime,2))*(scale_factor))/1000./sqrt(length(targLookTime))    0.49; (std(mean(dstrLookTime,2))*(scale_factor))/1000./sqrt(length(dstrLookTime))    0.38];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
legend('WOLVES Model', 'Yu & Smith 2011');
title ('Target vs Distractor Looking Time');
ylabel('Looking time per test trial');
ylim([0 4]);

figure(202);%Plot Target vs Distractor looking time during testing
errorbar(mean(targLookTime((goodLearners()==1),:))*scale_factor/1000,(std(targLookTime((goodLearners()==1),:))*scale_factor/1000)./sqrt(length(targLookTime((goodLearners()==1),:))) ,plotStyle{1});%
hold on
errorbar(mean(dstrLookTime((goodLearners()==1),:))*scale_factor/1000,(std(dstrLookTime((goodLearners()==1),:))*scale_factor/1000)./sqrt(length(dstrLookTime((goodLearners()==1),:))),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Testing Trial');
ylabel('Strong Learners: Looking Time at test');
ylim([0 8]);

figure(203);%Plot Target vs Distractor looking time during testing
errorbar(mean(targLookTime((goodLearners()==0),:))*scale_factor/1000,(std(targLookTime((goodLearners()==0),:))*scale_factor/1000)./sqrt(length(targLookTime((goodLearners()==0),:))) ,plotStyle{1});%
hold on
errorbar(mean(dstrLookTime((goodLearners()==0),:))*scale_factor/1000,(std(dstrLookTime((goodLearners()==0),:))*scale_factor/1000)./sqrt(length(dstrLookTime((goodLearners()==0),:))),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Testing Trial');
ylabel('Weak Learners: Looking Time at test');
ylim([0 8]);


figure(3)%Plot Targ vs Dist looking time for Strong learners
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime(goodLearners()==1,:)))*scale_factor/1000    mean(mean(dstrLookTime(goodLearners()==1,:)))*scale_factor/1000];
errY=[ (std(mean(targLookTime(goodLearners()==1,:),2))*scale_factor/1000)./sqrt(length(targLookTime((goodLearners()==1),:)))     (std(mean(dstrLookTime(goodLearners()==1,:),2))*scale_factor/1000)./sqrt(length(dstrLookTime((goodLearners()==1),:)))];
%b=bar(sts)
b = barwitherr(errY, sts, 0.6);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',11);
%legend('Target ','Distractor ');
title ('Strong Learners');
ylabel('Looking time (ms)');
ylim([0 4]);

figure(4)%Plot Targ vs Dist looking time for Strong learners
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime(goodLearners()==0,:)))*scale_factor/1000    mean(mean(dstrLookTime(goodLearners()==0,:)))*scale_factor/1000];
errY=[ (std(mean(targLookTime(goodLearners()==0,:),2))*scale_factor/1000)./sqrt(length(targLookTime((goodLearners()==0),:)))     (std(mean(dstrLookTime(goodLearners()==0,:),2))*scale_factor/1000)./sqrt(length(dstrLookTime((goodLearners()==0),:)))];
%b=bar(sts)
b = barwitherr(errY, sts, 0.6);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',11);
%legend('Target ','Distractor ');
title ('Weak Learners');
ylabel('Looking time (ms)');
ylim([0 4]);


figure(5);%Plot Strong vs Weak looking time during over TEST trials
errorbar(mean(targLookTime,1)*scale_factor/1000,(std(targLookTime,1)*scale_factor/1000)./sqrt(length(targLookTime)),plotStyle{1});%
hold on
errorbar(mean(dstrLookTime,1)*scale_factor/1000,(std(dstrLookTime,1)*scale_factor/1000)./sqrt(length(dstrLookTime)),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('test trial');
ylabel('total looking time Target vs Distractor');
ylim([0 8]);


figure(6)%Plot Proportion of Strong/ Weak learners
blockNames={'Strong'; 'Weak'};
sts = [12/18 sum(goodLearners)/numSubjects; 6/18  (numSubjects-sum(goodLearners))/numSubjects];
errY =[ 0 0; 0 0];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',12);
legend('Yu & Smith 2011','WOLVES Model');
%title ('Target vs Distractor Looking Time');
ylabel('Proportion of Learners');
grid on
ylim([0 1]);
%sqrt((((sum(goodLearners)/numSubjects) - (12/18))^2 + (((numSubjects-sum(goodLearners))/numSubjects) - (6/18))^2)./2)

figure(7)% Avg # of Words Learnt
sts = [4; 3.5; mean(sum(LearntWords(),2))];
errY= [0; 0; 0];
hb=barwitherr(errY,sts,0.6);% Plot with errorbars
set(gca,'xticklabel',{'Smith & Yu 2008 (14 m)'; 'Yu & Smith 2011';'WOLVES Model'},'fontsize',11);
ylabel('Average number of words learnt');
ylim([0 6]);
%sqrt(((mean(sum(LearntWords(),2))- 4)^2 + (mean(sum(LearntWords(),2)) - 3.5)^2)./2)
hold on
barFontSize = 11;
for b = 1 : length(sts)
	% Plot one single bar as a separate bar series.
	handleToThisBarSeries(b) = bar(b, sts(b), 'BarWidth', 0.6);
	% Apply the color to this bar series.
	set(handleToThisBarSeries(b), 'FaceColor', barColorMap(b,:));
	% Place text atop the bar
	barTopper = sprintf('%.3f', sts(b));
	text(b-0.1, sts(b)+0.15, barTopper, 'FontSize', barFontSize);
	hold on;
end


word_Vec=zeros(1,TEST_DUR/scale_factor);
for wc=1:length(word_On)
    for w = word_On(wc):word_Off(wc)
        word_Vec(w)=1;
    end
end
tLookAtLearnt=zeros(nTestTrials,1);tLookAtNonLrnt=zeros(nTestTrials,1);
for subject=1:numSubjects
    savestate_historyLt = xsit_result.test(subject).historyLt(:,vis_On:word_Off(1));
    savestate_historyRt = xsit_result.test(subject).historyRt(:,vis_On:word_Off(1));
    % create the off-looking history Vector
    for i=1:nTestTrials
        for j=1:word_Off(1)
           if  (round(savestate_historyLt(i,j)) + round(savestate_historyRt(i,j))) > 0.5
               savestate_historyOt(i,j)=0; 
           else savestate_historyOt(i,j)=1; end
        end
    end
    for tt=1:nTestTrials %% each test trial
        %%% measure sync time between off-looking & word presentation
        wordoff_time=savestate_historyOt(tt,1:word_Off(1)); sync_time(subject,tt)= sum(wordoff_time.*word_Vec(1:word_Off(1)));      
        %%% measure looking ratio at every time instant 
        lookingOn(subject,tt,1:word_Off(1))= ((round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1))))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1)))+savestate_historyOt(tt,1:word_Off(1))))';
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
                        look_at_targ_show1(subject,kk,:)=(round(savestate_historyLt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1)))+savestate_historyOt(tt,1:word_Off(1))))';%+savestate_historyOt(tt,1:8000)
                        look_at_dist_show1(subject,kk,:)=(round(savestate_historyRt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1)))+savestate_historyOt(tt,1:word_Off(1))))';
                    elseif (show_trial==2) % second presentation of this word
                        look_at_targ_show2(subject,kk,:)=(round(savestate_historyLt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1)))+savestate_historyOt(tt,1:word_Off(1))))';%add time historyLt
                        look_at_dist_show2(subject,kk,:)=(round(savestate_historyRt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1)))+savestate_historyOt(tt,1:word_Off(1))))';
                   end
                elseif (strcmp(s1,'R')) % referent on right
                    if (show_trial==1)  % first presentation of this word of 12 trials
                        look_at_targ_show1(subject,kk,:)=(round(savestate_historyRt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';%add time historyRt
                        look_at_dist_show1(subject,kk,:)=(round(savestate_historyLt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
                    elseif (show_trial==2) % second of 12 trials with this word
                        look_at_targ_show2(subject,kk,:)=(round(savestate_historyRt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1)))+savestate_historyOt(tt,1:word_Off(1))))';%add time historyRt
                        look_at_dist_show2(subject,kk,:)=(round(savestate_historyLt(tt,1:word_Off(1)))./(round(savestate_historyLt(tt,1:word_Off(1)))+round(savestate_historyRt(tt,1:word_Off(1)))+savestate_historyOt(tt,1:word_Off(1))))';
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
                            tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyLt(tt,1:word_Off(1)));%add time historyLt
                        elseif (strcmp(s1,'R'))
                            tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyRt(tt,1:word_Off(1)));% add time historyRt
                        end
                    end
                end        
        elseif (LearntWords(subject,kk)==0)%if word is NOT learnt
                for tt=1:nTestTrials
                    if(xsit_result.train(subject).Words{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial
                        s1=char(xsit_result.test(subject).test_pair(tt,2*nFeatures+2));
                        if (strcmp(s1,'L'))% on left
                            tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyLt(tt,1:word_Off(1)));%add time historyL
                        elseif (strcmp(s1,'R'))
                            tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyRt(tt,1:word_Off(1)));% add time historyR
                        end
                    end
                end
        else
            disp('Error in LearntWords array: supurious data ');
        end
    end
end


figure (1008); % plot sync time vs learning
legend('sync mean');
T1=targLookTime(:); S1=sync_time(:);
scatter(S1,T1);
ylabel('looking to target');
xlabel('sync time');

figure (109);%Plot duration of looking to learnt vs unlearnt words
plot(mean(tLookAtLearnt)*scale_factor/1000); 
hold on
plot((tLookAtNonLrnt)*scale_factor/1000); 
xlabel('per TEST trial');
ylabel('duration of looking to words');
legend('Learned words','Non-learned words');
ylim([0 4]);

all_objLook= NaN (numSubjects,nTestTrials,TEST_DUR/scale_factor);
all_targL= NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
all_distL= NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
larW=NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
ularW=NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
for subject=1:numSubjects
        all_objLook(subject,:,1:word_Off(1)) = lookingOn(subject,:,1:word_Off(1));
        for kk=1:6
            all_targL(subject,kk,1:word_Off(1))  = mean([look_at_targ_show1(subject,kk,:) look_at_targ_show2(subject,kk,:)]);
            all_distL(subject,kk,1:word_Off(1))  = mean([look_at_dist_show1(subject,kk,:) look_at_dist_show2(subject,kk,:)]);
            if LearntWords (subject,kk) == 1
              larW(subject,kk,1:word_Off(1))  = mean([look_at_targ_show1(subject,kk,:) look_at_targ_show2(subject,kk,:)]);
            end
            if LearntWords (subject,kk) == 0
              ularW(subject,kk,1:word_Off(1))  = mean([look_at_targ_show1(subject,kk,:) look_at_targ_show2(subject,kk,:)]);
            end
        end
end

figure (105);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot((squeeze(nanmean(nanmean(all_objLook(goodLearners()==1,:,:),1),2))))
hold on
plot((squeeze(nanmean(nanmean(all_targL(goodLearners()==1,:,:),1),2))))
hold on
plot((squeeze(nanmean(nanmean(all_distL(goodLearners()==1,:,:),1),2))))
hold on
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Both Objects','Target Object','Distractor');
xlabel('time');
ylabel('Strong Learners: proportion of subjects looking');
ylim([0 1]);

figure (111);%Plot timecourse looking learnt vs unlearnt 
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot(squeeze(nanmean(nanmean(larW(goodLearners()==1,:,:),1),2)))
hold on 
plot(squeeze(nanmean(nanmean(ularW(goodLearners()==1,:,:),1),2)))
hold on 
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Learnt words',    'Unlearnt words');
xlabel('time');
ylabel('Strong Learners: proportion of subjects looking');
ylim([0 1]);


figure (106);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot((squeeze(nanmean(nanmean(all_objLook(goodLearners()==0,:,:),1),2))))
hold on
plot((squeeze(nanmean(nanmean(all_targL(goodLearners()==0,:,:),1),2))))
hold on
plot((squeeze(nanmean(nanmean(all_distL(goodLearners()==0,:,:),1),2))))
hold on
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Both Objects','Target Object','Distractor');
xlabel('time');
ylabel('Weak Learners: proportion of subjects looking');
ylim([0 1]);

figure (1110);%Plot timecourse looking learnt vs unlearnt 
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot(squeeze(nanmean(nanmean(larW(goodLearners()==0,:,:),1),2)))
hold on 
plot(squeeze(nanmean(nanmean(ularW(goodLearners()==0,:,:),1),2)))
hold on 
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Learnt words',    'Unlearnt words');
xlabel('time');
ylabel('Weak Learners: proportion of subjects looking');
ylim([0 1]);
 
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



%% TRAINING CONDITION ANALYSIS trials analysis
word_On = floor([500 2000]/scale_factor);%XSIT            
word_Off = floor([1500 3000 ]/scale_factor);word_Len=floor(1000/scale_factor);
% word_On = floor([368 1850]/scale_factor);  %% XTRAP    
% word_Off = floor((745+[368 1850])/scale_factor);word_Len=floor(745/scale_factor);
vis_On=1;vis_Off=(TRAIN_DUR/scale_factor); nFix_limit=10;

tLookAtLearntWords=NaN(numSubjects,nTrainTrials);
tLookAtNonLrntWords=NaN(numSubjects,nTrainTrials);
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



 strongLtime=nanmean(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000;
 SLerr=(nanstd(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==1,:)));
 %
 strongNLtime=nanmean(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000;
 SNLerr=(nanstd(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==1,:)));
%
 weakLtime=nanmean(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000;
 WLerr=(nanstd(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==0,:)));
% 
 weakNLtime=nanmean(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000;
 WNLerr=(nanstd(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==0,:)));
%%% above nan lengths are NOT CORRECT... length(x(~isnan(x)))

figure(8);%Plot looking to learnt vs unlearnt 
errorbar(strongLtime,SLerr,plotStyle{1});%
hold on
errorbar(strongNLtime,SNLerr,plotStyle{2+compStyle});%
xlabel('per training trial');
ylabel('Looking time (s)');
legend('Learned Words','Non-Learned Words');
title('Strong Learners')

figure(801)% Plot looking time at leant vs non-learnt words for Strong learners
blockNames={'1to10','11to20','21to30'};%
sts = [ mean(strongLtime(1:10)) mean(strongNLtime(1:10)); mean(strongLtime(11:20)) mean(strongNLtime(11:20)); mean(strongLtime(21:30)) mean(strongNLtime(21:30))];
errY =[ mean(SLerr(1:10)) mean(SNLerr(1:10)); mean(SLerr(11:20)) mean(SNLerr(11:20)); mean(SLerr(21:30)) mean(SNLerr(21:30))];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Strong Learners');
ylabel('Looking time (ms)');
legend('Avg Learned Words','Avg NonLearned Words');

figure(9);%Plot looking to learnt vs unlearnt 
errorbar(weakLtime,WLerr,plotStyle{1});%
hold on
errorbar(weakNLtime,WNLerr,plotStyle{2+compStyle});%
xlabel('per training trial');
ylabel('Looking time (s)');
legend('Learned Words','Non-Learned Words');
title('Weak Learners')

figure(901)% Plot looking time at leant vs non-learnt words for Weak learners
blockNames={'1to10','11to20','21to30'};%
sts = [ mean(weakLtime(1:10)) mean(weakNLtime(1:10)); mean(weakLtime(11:20)) mean(weakNLtime(11:20)); mean(weakLtime(21:30)) mean(weakNLtime(21:30))];
errY =[ mean(WLerr(1:10)) mean(WNLerr(1:10)); mean(WLerr(11:20)) mean(WNLerr(11:20)); mean(WLerr(21:30)) mean(WNLerr(21:30))];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Weak Learners');
ylabel('Looking time (ms)');
legend('Avg Learnt Words','Avg Non-Learnt Words');


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

figure (10);% Plot entropy in looking fixation durations
errorbar(mean(VarianceSub((goodLearners()==1),:))*scale_factor/1000,(std(VarianceSub((goodLearners()==1),:))*scale_factor/1000)./sqrt( length( VarianceSub((goodLearners()==1),:) )),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(mean(VarianceSub((goodLearners()==0),:))*scale_factor/1000,(std(VarianceSub((goodLearners()==0),:))*scale_factor/1000)./sqrt( length( VarianceSub((goodLearners()==0),:) )),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('Variance Strong vs Weak learners');
%ylabel('number of fixations/looks Strong learners');
legend('Strong','Weak');
%ylim([0 2.5]);
%ylim([1.5 3.5])
%hold off
summaF=mean(mean(VarianceSub((goodLearners()==1),:)));
xx=['Variance Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(VarianceSub((goodLearners()==0),:)));
xx=['Variance Weak learners ',num2str(summaF)]; disp(xx);

figure (1001);% Plot entropy in looking fixation durations
errorbar(mean(EntropySub((goodLearners()==1),:)),std(EntropySub((goodLearners()==1),:))./sqrt( length( EntropySub((goodLearners()==1),:) )),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(mean(EntropySub((goodLearners()==0),:)),std(EntropySub((goodLearners()==0),:))./sqrt( length( EntropySub((goodLearners()==0),:) )),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('Entropy Strong vs Weak learners');
%ylabel('number of fixations/looks Strong learners');
legend('Strong','Weak');
summaF=mean(mean(EntropySub((goodLearners()==1),:)));
xx=['Entropy Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(EntropySub((goodLearners()==0),:)));
xx=['Entropy Weak learners ',num2str(summaF)]; disp(xx);

figure(11);%Plot Strong vs Weak looking time during a training trial
errorbar(mean(TotalLookTime(goodLearners()==1,:))*scale_factor/1000, (std(TotalLookTime(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(TotalLookTime(goodLearners()==1,:))),plotStyle{1});%
hold on
errorbar(mean(TotalLookTime(goodLearners()==0,:))*scale_factor/1000,(std(TotalLookTime(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(TotalLookTime(goodLearners()==0,:))),plotStyle{2+compStyle});%
legend('Strong','Weak');
xlabel('per training trial');
ylabel('total looking time Strong vs Weak learners');
ylim([2 4]);
%ylabel('total looking time Strong learners');
%hold off 

figure (12);% Plot number of fixations/looks over training trials
errorbar(mean(totnlooks(goodLearners()==1,:)),std(totnlooks(goodLearners()==1,:))./sqrt(length(totnlooks(goodLearners()==1,:))),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(mean(totnlooks(goodLearners()==0,:)),std(totnlooks(goodLearners()==0,:))./sqrt(length(totnlooks(goodLearners()==0,:))),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('number of fixations/looks Strong vs Weak learners');
%ylabel('number of fixations/looks Strong learners');
legend('Strong','Weak');
ylim([1.5 4.5])
%hold off
summaF=mean(mean(totnlooks(goodLearners()==1,:)));
xx=['number of fixations/looks Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(totnlooks(goodLearners()==0,:)));
xx=['number of fixations/looks Weak learners ',num2str(summaF)]; disp(xx);

figure (13);%Plot mean look duration of each fixation
errorbar(mean(meanlookdur(goodLearners()==1,:))*scale_factor/1000,(std(meanlookdur(goodLearners()==1,:))*scale_factor/1000)./sqrt( length( meanlookdur((goodLearners()==1),:) )), plotStyle{1});% mean look duration % multiped by timing scale factor
hold on
errorbar(mean(meanlookdur(goodLearners()==0,:))*scale_factor/1000,(std(meanlookdur(goodLearners()==0,:))*scale_factor/1000)./sqrt( length( meanlookdur((goodLearners()==0),:))), plotStyle{2+compStyle});% mean look duration % multiped by timing scale factor
xlabel('per training trial');
ylabel('mean look duration Strong vs Weak learners');
%ylabel('mean look duration Strong learners');
legend('Strong','Weak');
ylim([0 2.5]);
%hold off
summaF=mean(mean(meanlookdur(goodLearners()==1,:))*scale_factor)/1000;
xx=['mean look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanlookdur(goodLearners()==0,:))*scale_factor)/1000;
xx=['mean look duration Strong weak learners ',num2str(summaF)]; disp(xx);


figure (1301);%Plot mean look duration of each fixation
errorbar(mean(meanLukhadur(goodLearners()==1,:))*scale_factor/1000,(std(meanLukhadur(goodLearners()==1,:))*scale_factor/1000)./sqrt( length( meanLukhadur((goodLearners()==1),:) )), plotStyle{1});% mean look duration % multiped by timing scale factor
hold on
errorbar(mean(meanLukhadur(goodLearners()==0,:))*scale_factor/1000,(std(meanLukhadur(goodLearners()==0,:))*scale_factor/1000)./sqrt( length( meanLukhadur((goodLearners()==0),:) )), plotStyle{2+compStyle});% mean look duration % multiped by timing scale factor
xlabel('per training trial (from durations)');
ylabel('mean look duration Strong vs Weak learners');
%ylabel('mean look duration Strong learners');
legend('Strong','Weak');
ylim([0 2.5]);
%hold off
summaF=mean(mean(meanLukhadur(goodLearners()==1,:))*scale_factor)/1000;
xx=['mean look duration Strong learners (indiv calcs) ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanLukhadur(goodLearners()==0,:))*scale_factor)/1000;
xx=['mean look duration Strong weak learners (indiv calcs) ',num2str(summaF)]; disp(xx);

figure (14);%Plot duration of longest look per trial
errorbar(mean(totlonglookdur(goodLearners()==1,:))*scale_factor/1000,(std(totlonglookdur(goodLearners()==1,:))*scale_factor/1000)./sqrt( length( totlonglookdur((goodLearners()==1),:) )),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(totlonglookdur(goodLearners()==0,:))*scale_factor/1000,(std(totlonglookdur(goodLearners()==0,:))*scale_factor/1000)./sqrt( length( totlonglookdur((goodLearners()==0),:) )), plotStyle{2+compStyle}); %duration of longest look per trial
xlabel('per training trial');
ylabel('duration of longest look Strong vs Weak learners');
ylim([0 4]);
%ylabel('duration of longest look Strong learners');
legend('Strong','Weak');
%hold off
summaF=mean(mean(totlonglookdur(goodLearners()==1,:))*scale_factor)/1000;
xx=['Longest look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totlonglookdur(goodLearners()==0,:))*scale_factor)/1000;
xx=['Longest look duration Strong weak learners ',num2str(summaF)]; disp(xx);

figure(15);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(mean(targLookTimeTraining((goodLearners()==1),:))*scale_factor/1000,(std(targLookTimeTraining((goodLearners()==1),:)) *scale_factor/1000)./sqrt(length(targLookTimeTraining((goodLearners()==1),:))),plotStyle{1});%
hold on
errorbar(mean(dstrLookTimeTraining((goodLearners()==1),:))*scale_factor/1000,(std(dstrLookTimeTraining((goodLearners()==1),:)) *scale_factor/1000)./sqrt(length(dstrLookTimeTraining((goodLearners()==1),:))),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Training Trial');
ylabel('Strong Learners: Looking Time when words are ON');
ylim([0.2 1.75]);

figure(16);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(mean(targLookTimeTraining((goodLearners()==0),:))*scale_factor/1000,(std(targLookTimeTraining((goodLearners()==0),:)) *scale_factor/1000)./sqrt(length(targLookTimeTraining((goodLearners()==0),:))),plotStyle{1});%
hold on
errorbar(mean(dstrLookTimeTraining((goodLearners()==0),:))*scale_factor/1000,(std(dstrLookTimeTraining((goodLearners()==0),:)) *scale_factor/1000)./sqrt(length(dstrLookTimeTraining((goodLearners()==0),:))),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Training Trial');
ylabel('Weak Learners: Looking Time when words are ON');
ylim([0.2 1.75]);

figure (17);%Plot Mean proportion of looking Varying vs repeated C
errorbar(mean(tLookVarying./(tLookVarying+tLookRepeated)),std(tLookVarying./(tLookVarying+tLookRepeated))./sqrt(length(tLookVarying)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeated./(tLookVarying+tLookRepeated)),std(tLookRepeated./(tLookVarying+tLookRepeated))./sqrt(length(tLookRepeated)), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('trial block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
ylabel('Proportion looking time');
%ylabel('duration of longest look Strong learners');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('Habituation over training');

figure (1711);%Plot Mean proportion of looking Varying vs repeated C
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
%ylabel('duration of longest look Strong learners');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('Habituation over training');

figure (18);%Plot Mean propertion of looking Varying vs repeated C
errorbar(mean(tLookVarying((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:))),std(tLookVarying((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:)))./sqrt(length(tLookVarying((goodLearners()==1),:))),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeated((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:))),std(tLookRepeated((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:)))./sqrt(length(tLookVarying((goodLearners()==1),:))), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('trial block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
ylabel('Proportion looking time Strong Learners');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('Habituation over training');

figure (19);%Plot Mean propertion of looking Varying vs repeated C
errorbar(mean(tLookVarying((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:))),std(tLookVarying((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:)))./sqrt(length(tLookVarying((goodLearners()==0),:))),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeated((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:))),std(tLookRepeated((goodLearners()==0),:)./(tLookVarying((goodLearners()==0),:)+tLookRepeated((goodLearners()==0),:)))./sqrt(length(tLookVarying((goodLearners()==0),:))), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('trial block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
ylabel('Proportion Looking Time Weak Learners');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('Habituation over training');

summaF=mean(mean([TotalLookTime; TotalLookTime])) *(scale_factor/1000);
xx=['Avg looking time per training trial is ',num2str(summaF)]; disp(xx);
var1= mean(mean(TotalLookTime))*scale_factor/1000;
xx=['Avg looking time per training trial per Strong learner is ',num2str(var1)]; disp(xx);
var2= mean(mean(TotalLookTime))*scale_factor/1000;
xx=['Avg looking time per training trial per Weak learner is ',num2str(var2)]; disp(xx);

figure(20);%Plot mean looking time per TEST trial
sts = [3.04;     2.99; mean(mean([TotalLookTime((goodLearners()==1),:); TotalLookTime((goodLearners()==0),:)])) *(scale_factor/1000)];
errY=[0.0;        0.0;  std(mean([TotalLookTime((goodLearners()==1),:); TotalLookTime((goodLearners()==0),:)])) *(scale_factor/1000)];
hb=barwitherr(errY,sts,0.6);% Plot with errorbars
set(gca,'xticklabel',{'Smith & Yu 2008 (14 m)'; 'Yu & Smith 2011';'WOLVES Model'},'fontsize',11)
%title ('total looking time per test trial');
ylabel('Looking time (seconds) at training');
ylim([0 4]);
%sqrt(((mean(mean([TotalLookTimeS; TotalLookTimeW])) *(scale_factor/1000) - 3.04)^2 + (mean(mean([TotalLookTimeS; TotalLookTimeW])) *(scale_factor/1000) - 2.99)^2)./2)
hold on
barFontSize = 11;
for b = 1 : length(sts)
	% Plot one single bar as a separate bar series.
	handleToThisBarSeries(b) = bar(b, sts(b), 'BarWidth', 0.6);
	% Apply the color to this bar series.
	set(handleToThisBarSeries(b), 'FaceColor', barColorMap(b,:));
	% Place text atop the bar
	barTopper = sprintf('%.3f', sts(b));
	text(b-0.1, sts(b)+0.15, barTopper, 'FontSize', barFontSize);
	hold on;
end




%% trace analysis
Wrong_inTrace= zeros(numSubjects,1); Correct_inTrace= zeros(numSubjects,1);
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
        inputMapping1(xx1(kk),yy(kk))=0;
        inputMapping2(xx2(kk),yy(kk))=0;
        for jj=1:6
            inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
            inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
        end
        a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
        W_inTr = W_inTr + mean([mean(a_in(a_in>0.001)) mean(b_in(b_in>0.001))]);
    end
    Correct_inTrace(subject)=C_inTr/nObjects;
    Wrong_inTrace(subject)=W_inTr/nObjects;
    InCorr_assocs(subject)=mean([as_count1-1 as_count2-1])/nObjects;
    EntropyTrace(subject)= mean( [entropy(inputMapping1) entropy(inputMapping2)] ); 
end

figure(21);%Entropy
blockNames={'Strong';'Weak'};
sts = [mean(EntropyTrace((goodLearners()==1))); mean(EntropyTrace((goodLearners()==0)))  ];
errY =[std(EntropyTrace((goodLearners()==1)))/sqrt(length(EntropyTrace((goodLearners()==1)))); std(EntropyTrace((goodLearners()==0)))/sqrt(length(EntropyTrace((goodLearners()==0))))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Entropy in the traces');
%ylabel('Looking time per test trial');
%ylim([0 4]);

figure(22);%My own Entropy: No of incorrect traces
blockNames={'Strong';'Weak'};
sts = [mean(InCorr_assocs((goodLearners()==1))); mean(InCorr_assocs((goodLearners()==0)))  ];
errY =[std(InCorr_assocs((goodLearners()==1)))/sqrt(length(InCorr_assocs((goodLearners()==1)))); std(InCorr_assocs((goodLearners()==0)))/sqrt(length(InCorr_assocs((goodLearners()==0))))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Ptoportion of incorrect assocs in the traces');
ylabel('# of incorrect assocs/# of words (6)');

figure(221);%My own Entropy: No of incorrect traces
blockNames={'Strong';'Weak'};
sts = [nanmean(Correct_inTrace((goodLearners()==1))); nanmean(Correct_inTrace((goodLearners()==0)))  ];
errY =[nanstd(Correct_inTrace((goodLearners()==1)))/sqrt(length(Correct_inTrace((goodLearners()==1)))); nanstd(Correct_inTrace((goodLearners()==0)))/sqrt(length(Correct_inTrace((goodLearners()==0))))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Strength of correct assocs in the traces');

figure(222);%My own Entropy: No of incorrect traces
blockNames={'Strong';'Weak'};
sts = [nanmean(Wrong_inTrace((goodLearners()==1))); nanmean(Wrong_inTrace((goodLearners()==0)))  ];
errY =[nanstd(Wrong_inTrace((goodLearners()==1)))/sqrt(length(Wrong_inTrace((goodLearners()==1)))); nanstd(Wrong_inTrace((goodLearners()==0)))/sqrt(length(Wrong_inTrace((goodLearners()==0))))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Strength of INCORRECT assocs in the traces');

%%%%%%% Association hwf trace analysis
corrAsocn=zeros(numSubjects,nObjects);
cS=1;cW=1;
for subject=1:numSubjects  
    
    AsocMat=squeeze(xsit_result.train(subject).hwf(1,:,:));
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
        
        
        [maxIn(kk), indIn(kk)] = max(inputMapping(:,cell2mat(xsit_result.train(subject).Words(kk))));
        [maxAs(kk) indAs(kk)] = max(AsocMat(:,cell2mat(xsit_result.train(subject).Words(kk))));       
        if (abs(indIn(kk)-indAs(kk)) <= 2)%if association is correct i..e same as input?
           corrAsocn(subject, kk)=1; % wrongAssocn = 6-corrAsocn
        end
    end 
end

% 
SLer=[];WLer=[];SNon=[];WNon=[];SLer2=[];WLer2=[];SNon2=[];WNon2=[];
for subject=1:numSubjects 
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

% for subject=1:numSubjects
% 
% inputMapping1=zeros(306,20);
% inputMapping2=zeros(306,20);
% muddled_pairs=0;
% muddled_mem_val=1;
%     for kk=1:nObjects
%         xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
%         xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
%         yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
%     end
%     for kk=1:nObjects       
%         inputMapping1(xx1(kk),yy(kk))=1;
%         inputMapping2(xx2(kk),yy(kk))=1;
%         for jj=1:8
%             inputMapping1(xx1(kk)+jj-4,yy(kk))=1;
%             inputMapping2(xx2(kk)+jj-4,yy(kk))=1;
%         end   
%     end
%     
%     figure(23)
%     
%     lsp=subplot(1,2,1);
%     %surface(inputMapping1)
%     %hold on
%     [mA, iA] = max(squeeze(xsit_result.train(subject).hwf(1,:,:)));
%     surface(squeeze(xsit_result.train(subject).hwf(1,:,:)));
%     title([num2str(mA)]);
%     shading flat;     
%     
%     rsp= subplot(1,2,2);
%     hold on
%     
%     title('Input');
%     
%     surface(squeeze(xsit_result.test(subject).hwft(1,:,:)));
%     title('Test');
%     shading flat;
%     
%     if (goodLearners(subject)==1), suptitle(['Strong # ' num2str(subject)]),
%     elseif (goodLearners(subject)==0), suptitle(['Weak # ' num2str(subject)]), end
%     pause(7);
%     clf;
% end



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

%%% fixations agianst other measures
fix_count_All=mean(totnlooks,2);
TargDis_All=mean((targLookTime./(targLookTime+dstrLookTime)),2);
Words_All=mean(LearntWords,2);
T1=fix_count_All(:);
S1=InCorr_assocs(:); S2=TargDis_All(:); S3=Words_All(:);
figure (31)
scatter(T1,S1);
yb = scatstat1(T1,S1,0,@mean);
plot(T1,yb,'bo')
xlabel('fixation count')
ylabel('# incorrect associations');

figure (32)
scatter(T1,S2);
%hold on
yb = scatstat1(T1,S2,0,@mean);
plot(T1,yb,'bo')
xlabel('fixation count')
ylabel('Prop time looking to target');

figure (33)
scatter(T1,S3);
%hold on
yb = scatstat1(T1,S3,0,@mean);
plot(T1,yb,'bo')
xlabel('fixation count')
ylabel('words learnt');




%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%old fixation code 2
% %     nlooks=zeros(2,size(savestate_historyL,1)); %L/R %%SAVE US!! 
% %     longlookdur=zeros(2,size(savestate_historyL,1));
% %     for side=1:2      
% %         if side == 1
% %             ldata = savestate_historyL;
% %             cdata = savestate_historyR;
% %         else
% %             ldata = savestate_historyR;
% %             cdata = savestate_historyL;
% %         end                      
% %         for tr=1:size(ldata,1) %% or size(rdata,1) 
% %             looking=0;
% %             templonglookdur=0;
% %             lastLookEndTime=1;
% %             for time=1:size(ldata,2) %% or size(rdata,2) 
% %                 if (round(ldata(tr,time)) == 1)%% if model looking to left
% %                     if looking == 0 %% if the model is not already looking 
% %                         if sum(round(cdata(tr,lastLookEndTime:time))) > 10 ||  sum(nlooks(side,tr))== 0 %%%% if it was looing to some other object before or no looks so far
% %                             nlooks(side,tr) = nlooks(side,tr)+1;
% %                             if templonglookdur > longlookdur(side,tr)
% %                                 longlookdur(side,tr) = templonglookdur;
% %                             end
% %                             templonglookdur=0;
% %                         else
% %                             templonglookdur= templonglookdur + (time-lastLookEndTime); %%% add the gap between consecutive looks to same side
% %                         end
% %                     end
% %                     looking = 1;
% %                     templonglookdur = templonglookdur+1;
% %                 else
% %                     if looking == 1; lastLookEndTime=time;end
% %                     looking = 0;
% %                 end
% %             end
% % %             if (round(ldata(tr,time-1)) == 1)
% % %                 nlooks(side,tr) = nlooks(side,tr)+1;
% % %             end
% %             if templonglookdur > longlookdur(side,tr)
% %                 longlookdur(side,tr) = templonglookdur;
% %             end
% %         end   
% %     end


%%%%% olf fixation code 1
% %     for side=1:2      
% %         if side == 1
% %             ldata = savestate_historyL;
% %             cdata = savestate_historyR;
% %             odata = savestate_historyO;
% %         else
% %             ldata = savestate_historyR;
% %             cdata = savestate_historyL;
% %             odata = savestate_historyO;
% %         end                      
% %         for tr=1:size(ldata,1) %% for each training trial
% %             prevlooking=0;
% %             templonglookdur=0;
% %             dummylookdur=0;
% %             lastLookEndTime=1;
% %             for time=1:size(ldata,2) %% for each time-event/millisec  %%%sum(round(odata(tr,lastLookEndTime:time))) > prev_look_threshold ||
% %                 if (round(ldata(tr,time)) == 1)%% if child is looking now at chosen side 
% %                     if prevlooking == 0 %% if child was not looking to chosen side in last millisec
% %                         if (sum(round(odata(tr,lastLookEndTime:time))) > (MIN_LOOK_DURATION+2)) || (sum(round(cdata(tr,lastLookEndTime:time))) > (MIN_LOOK_DURATION+2)) || ( nlooks(side,tr) == 0) %%%% if it was looking off or to some other object before this look or away or no looks so far                          
% %                             nlooks(side,tr) = nlooks(side,tr)+1;
% %                             if templonglookdur > longlookdur(side,tr) % save longest looks
% %                                 longlookdur(side,tr) = templonglookdur;
% %                             end
% %                             templonglookdur=0;% reset look-length counter
% %                         else
% %                            templonglookdur= templonglookdur + (time-lastLookEndTime); %add off-looking flicker bit to current look length
% %                         end
% %                         dummylookdur=0; % result this possible look's length
% %                     end
% %                     prevlooking = 1; 
% %                     dummylookdur=dummylookdur+1; % add current look length
% %                 else
% %                     if prevlooking == 1 % checks if model was looking in the last time-event just before this
% %                         templonglookdur= templonglookdur+dummylookdur; %add this looks possible length to the look-duration
% %                         lastLookEndTime=time-1; % reset nd of last look to previos time-event/millisec
% %                         if (templonglookdur < MIN_LOOK_DURATION) % checks if this look was too short, the remove it
% %                             nlooks(side,tr) = nlooks(side,tr)-1;
% %                             templonglookdur
% %                             dummylookdur
% %                         end
% %                     end
% %                     prevlooking = 0;
% %                 end
% %             end
% %             if templonglookdur > longlookdur(side,tr)
% %                 longlookdur(side,tr) = templonglookdur;
% %             end
% %         end   
% %     end

%%%%
%%% old fixation code 0 primary
% % % % % %     for side=1:2
% % % % % %         
% % % % % %         if side == 1
% % % % % %             ldata = savestate_historyL;
% % % % % %         else
% % % % % %             ldata = savestate_historyR;
% % % % % %         end
% % % % % %         for tr=1:size(ldata,1)
% % % % % %             look=0;
% % % % % %             templonglookdur=0;
% % % % % %             for time=1:size(ldata,2)
% % % % % %                 if (round(ldata(tr,time)) == 1)
% % % % % %                     if look == 0
% % % % % %                         nlooks(side,tr) = nlooks(side,tr)+1;
% % % % % %                         if templonglookdur > longlookdur(side,tr)
% % % % % %                             longlookdur(side,tr) = templonglookdur;
% % % % % %                         end
% % % % % %                         templonglookdur=0;
% % % % % %                     end
% % % % % %                     look = 1;
% % % % % %                     templonglookdur = templonglookdur+1;
% % % % % %                 else
% % % % % %                     look = 0;
% % % % % %                 end
% % % % % %             end
% % % % % %             if (round(ldata(tr,time-1)) == 1)
% % % % % %                 nlooks(side,tr) = nlooks(side,tr)+1;
% % % % % %             end
% % % % % %             if templonglookdur > longlookdur(side,tr)
% % % % % %                 longlookdur(side,tr) = templonglookdur;
% % % % % %             end
% % % % % %         end   
% % % % % %     end