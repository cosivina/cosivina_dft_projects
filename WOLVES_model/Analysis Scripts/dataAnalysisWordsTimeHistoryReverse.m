%% Data Analysis 
%%%sim.saveSettings('bbc.json')%method to save simulator
%% LOADING
clear all; close all;
nObjects=6;nFeatures=2; scale_factor=8; TRAIN_DUR=4000; TEST_DUR=8000; %check with auto file
plotStyle = {'k-','b-','g-','c-','r-','m-','y-','k-.','b-.','g-.','c-.','r-.','m-.','y-.','b:','w.'};
compStyle=7;
blockNames{10}=[];sts{10}=[];errY{10}=[];
% xsit_result.Names={[1] [2] [3] [4] [5] [6]};
%xsit_result.train =load('5atn_13wf_xsitNF_surprise3_1_train.mat');
%xsit_result.test =load('5atn_13wf_xsitNF_surprise3_1_test.mat');

simName = 'less_offlooking5_noWord_hwmf_05-28-Oct-18_Smith_Yu_2013_';
%simName = 'muddled_mem_6_pairs_sy_base05-Oct-2018_test_only_learners_smithYu_'%unified_param_run8_02-Oct-2018_Smith_Yu_2013_'%'bbc_H50F30Asc50_Tzero_Smith_Yu_2013_'; %test4_Smith_Yu_2008_
OutName = [simName,'results.mat']
xsit_result = load (OutName);

%% add another simulation
% simName = 'less_offlooking5b-27-Oct-18_Smith_Yu_2013_';
% OutName1 = [simName,'results.mat']
% xsit_result1 = load (OutName1);
% xsit_result.train = [xsit_result.train ; xsit_result1.train];  
% xsit_result.test = [xsit_result.test ; xsit_result1.test];  
%% 
numSubjects=size(xsit_result.test,1);
xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 
%xsit_result.sim.saveSettings('recoverParams.json'); visdiff ('recoverParams1.json', 'recoverParams.json')

%% TEST CONDITION ANALYSIS
correct_proportion=zeros(numSubjects,2);
targLookTime=zeros(numSubjects,2*nObjects);
dstrLookTime=zeros(numSubjects,2*nObjects);
goodLearners=999*ones(numSubjects,1);
LearntWords= 999*ones(numSubjects,nObjects);
%word_On = floor([500 2000 4500 6000]/scale_factor);%[0 450 875 1300 ];    
%word_Off = floor([1500 3000 5500 7000]/scale_factor);%[375 650 1075 1500 ];
slice_wordOnset4=1;%6/8;
vis_On = 1;
vis_Off = floor(TEST_DUR/scale_factor);
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
for subject=1:numSubjects
    lcorrect=0;rcorrect=0;targWord=zeros(nObjects,1);dstrWord=zeros(nObjects,1);
    for trt=1:2*nObjects
          lLook= sum( xsit_result.test(subject).historyLt(trt,vis_On:vis_Off));%250:2000 
          rLook= sum( xsit_result.test(subject).historyRt(trt,vis_On:vis_Off));%250:2000
%        lLook= sum( xsit_result.test(subject).historyLt(trt,floor(500/scale_factor):floor(1750/scale_factor)));%250:2000
%        lLook=lLook+sum( xsit_result.test(subject).historyLt(trt,floor(2000/scale_factor):floor(3250/scale_factor)));
%        lLook=lLook+sum( xsit_result.test(subject).historyLt(trt,floor(4500/scale_factor):floor(5750/scale_factor)));
%        lLook=lLook+sum( xsit_result.test(subject).historyLt(trt,floor(6000/scale_factor):floor(7250/scale_factor)));   
%         
%        rLook= sum( xsit_result.test(subject).historyRt(trt,floor(500/scale_factor):floor(1750/scale_factor)));%250:2000
%        rLook=rLook+sum( xsit_result.test(subject).historyRt(trt,floor(2000/scale_factor):floor(3250/scale_factor)));
%        rLook=rLook+sum( xsit_result.test(subject).historyRt(trt,floor(4500/scale_factor):floor(5750/scale_factor)));
%        rLook=rLook+sum( xsit_result.test(subject).historyRt(trt,floor(6000/scale_factor):floor(7250/scale_factor)));   
       %%%%%%NEW
       for kk=1:nObjects
           %xsit_result.Names{kk}      
          if (xsit_result.train(subject).Words{kk}== xsit_result.test(subject).test_pair(trt,2*nFeatures+1))%word index
               s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
               
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
%%%%%
       s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
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
    end
    for kk=1:nObjects
        if (targWord(kk)>dstrWord(kk))
            LearntWords(subject,kk)=1;
        else
            LearntWords(subject,kk)=0;
        end
    end
    
    lcorrect = lcorrect/(2*nObjects);     rcorrect = rcorrect/(2*nObjects);
    correct_proportion(subject,1)=lcorrect;
    correct_proportion(subject,2)=rcorrect;
    %targLookTime(subject) =targLookTime(subject)./(2*nObjects);
    %dstrLookTime(subject) =dstrLookTime(subject)./(2*nObjects);
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
sts = [6.1;     5.92; (mean(mean(targLookTime+dstrLookTime))*(scale_factor*slice_wordOnset4))/1000];
errY=[0.05;     0.05; (std(mean(targLookTime+dstrLookTime,2))*(scale_factor*slice_wordOnset4))/1000];
hb=barwitherr(errY,sts,0.6);% Plot with errorbars
set(gca,'xticklabel',{'Smith & Yu 2008 (14 m)'; 'Yu & Smith 2011';'WOLVES Model'},'fontsize',11);
%title ('total looking time per test trial');
ylabel('Looking time (seconds)');
ylim([0 8]);
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
errorbar(sts,errY,'o');


figure(1000)% Plot total looking time during test trial
blockNames={'WOLVES Model XTRAP';'WOLVES in XSIT'};
sts = [(mean(mean(targLookTime+dstrLookTime))*(scale_factor*slice_wordOnset4))/1000;    5.89];
errY=[(std(mean(targLookTime+dstrLookTime,2))*(scale_factor*slice_wordOnset4))/1000;   0.1];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Avg Looking Time during a test trial');
ylabel('Looking time (ms)');

%legend('WOLVES Model','Smith & Yu 2008 14 month olds', 'Yu & Smith 2011');
grid on

figure(200)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ (mean(mean(targLookTime))*(scale_factor*slice_wordOnset4))/1000  3.43;  (mean(mean(dstrLookTime))*(scale_factor*slice_wordOnset4))/1000 2.67];
errY =[(std(mean(targLookTime,2))*(scale_factor*slice_wordOnset4))/1000  0.49; (std(mean(dstrLookTime,2))*(scale_factor*slice_wordOnset4))/1000  0.38];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
legend('WOLVES XTRAP', 'WOLVES XSIT');
title ('Target vs Distractor Looking Time');
ylabel('Looking time per test trial');


figure(2)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ 3.6    3.25  (mean(mean(targLookTime,2))*(scale_factor*slice_wordOnset4))/1000;  2.5  2.67 (mean(mean(dstrLookTime,2))*(scale_factor*slice_wordOnset4))/1000];
errY =[ 0.2    0.49  (std(mean(targLookTime,2))*(scale_factor*slice_wordOnset4))/1000; 0.25  0.38 (std(mean(dstrLookTime,2))*(scale_factor*slice_wordOnset4))/1000];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',12);
legend('Smith & Yu 2008 14 month olds', 'Yu & Smith 2011','WOLVES Model');
%title ('Target vs Distractor Looking Time');
ylabel('Looking time(seconds)');
%sqrt((((mean(mean(targLookTime))*(scale_factor*slice_wordOnset4))/1000 - 3.6)^2 + ((mean(mean(targLookTime))*(scale_factor*slice_wordOnset4))/1000 - 3.25)^2 + ((mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4))/1000 - 2.5)^2 + ((mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4))/1000 - 2.67)^2 )./4)

figure(2121)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ (mean(mean(targLookTime))*(scale_factor*slice_wordOnset4))/1000    3.25;  (mean(mean(dstrLookTime))*(scale_factor*slice_wordOnset4))/1000  2.67];
errY =[(std(mean(targLookTime,2))*(scale_factor*slice_wordOnset4))/1000    0.49; (std(mean(dstrLookTime,2))*(scale_factor*slice_wordOnset4))/1000    0.38];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
legend('WOLVES Model', 'Yu & Smith 2011');
title ('Target vs Distractor Looking Time');
ylabel('Looking time per test trial');

figure(1600);%Plot Target vs Distractor looking time during testing
errorbar(mean(targLookTime((goodLearners()==1),:))*scale_factor,std(targLookTime((goodLearners()==1),:))*scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTime((goodLearners()==1),:))*scale_factor,std(dstrLookTime((goodLearners()==1),:))*scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Testing Trial');
ylabel('Strong Learners: Looking Time at test');

figure(1601);%Plot Target vs Distractor looking time during testing
errorbar(mean(targLookTime((goodLearners()==0),:))*scale_factor,std(targLookTime((goodLearners()==0),:))*scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTime((goodLearners()==0),:))*scale_factor,std(dstrLookTime((goodLearners()==0),:))*scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Testing Trial');
ylabel('Weak Learners: Looking Time at test');


figure(3)%Plot Targ vs Dist looking time for Strong learners
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime(goodLearners()==1,:)))*scale_factor    mean(mean(dstrLookTime(goodLearners()==1,:)))*scale_factor];
errY=[ std(mean(targLookTime(goodLearners()==1,:),2))*scale_factor     std(mean(dstrLookTime(goodLearners()==1,:),2))*scale_factor];
%b=bar(sts)
b = barwitherr(errY, sts, 0.6);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',11);
%legend('Target ','Distractor ');
title ('Strong Learners');
ylabel('Looking time (ms)');

figure(4)%Plot Targ vs Dist looking time for Weak learners
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime(goodLearners()==0,:)))*scale_factor   mean(mean(dstrLookTime(goodLearners()==0,:)))*scale_factor];
errY=[std(mean(targLookTime(goodLearners()==0,:),2))*scale_factor     std(mean(dstrLookTime(goodLearners()==0,:),2))*scale_factor];
b = barwitherr(errY, sts,0.6);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',11);
title ('Weak Learners');
%legend('Target','Distractor');
%xlabel('Fixed order # ');
ylabel('Looking time (ms)');

figure(5);%Plot Strong vs Weak looking time during over TEST trials
errorbar(mean(targLookTime,1)*scale_factor,std(targLookTime,1)*scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTime,1)*scale_factor,std(dstrLookTime,1)*scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('test trial');
ylabel('total looking time Target vs Distractor');

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

tLookAtLearnt=zeros(12,1);
tLookAtNonLrnt=zeros(12,1);
for subject=1:numSubjects
    savestate_historyLt = fliplr(xsit_result.test(subject).historyLt(:,vis_On:vis_Off));
    savestate_historyRt = fliplr(xsit_result.test(subject).historyRt(:,vis_On:vis_Off));
    % create the off-looking history Vector
    for i=1:size(savestate_historyLt,1)
        for j=1:size(savestate_historyLt,2)
           if  (round(savestate_historyLt(i,j)) + round(savestate_historyRt(i,j))) > 0
               savestate_historyOt(i,j)=0; 
           else savestate_historyOt(i,j)=1; end
        end
    end
    for tt=1:size(savestate_historyLt,1) %each test trial 
        lookingOn(subject,tt,:)= ((round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:)))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
    %scale_factor(subject)=1;%TEST_DURATION/size(xsit_result.test(subject).historyOt,2);
    end
    % proportion of subjects looking to learnt vs unlearnt words
    for kk=1:nObjects
        show_trial=1;% first time or second time word prsented of 12 trials
        for tt=1:size(savestate_historyLt,1) %each test trial 
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
        
    %%find total duration of looking to learnt vs unlearnt words
    for kk=1:nObjects
       if (LearntWords(subject,kk)==1)%if the subject has learnt the kk word.. learnt is measured on looking basis.. change to mem_matrix basis if needed
                for tt=1:12%each test trial 
                    if(xsit_result.train(subject).Words{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial name1
                        %nHasLWord(block)=nHasLWord(block)+1;
                        s1=char(xsit_result.test(subject).test_pair(tt,2*nFeatures+2));
                        if (strcmp(s1,'L'))% on left
                            tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyLt(tt,:));%add time historyLt
                        elseif (strcmp(s1,'R'))
                            tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyRt(tt,:));% add time historyRt
                        end
                    end
                end
          
        elseif (LearntWords(subject,kk)==0)%if word is NOT learnt

                for tt=1:12%each trial may be 10 only
                    if(xsit_result.train(subject).Words{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial
                        %nHasUWord(block)=nHasUWord(block)+1;
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

figure (109);%Plot duration of looking to learnt vs unlearnt words
plot((tLookAtLearnt/numSubjects)*scale_factor); %duration of longest look per trial
hold on
plot((tLookAtNonLrnt/numSubjects)*scale_factor); %duration of longest look per trial
xlabel('per test trial');
ylabel('duration of looking to words');
legend('Learned words','Non-learned words');
%hold off
% figure (111);%Plot 
% plot(squeeze(nanmean(nanmean(look_at_targ_show1_bin1)))');
% figure (112)
% vv= (squeeze((nanmean(look_at_targ_show2)))') + (squeeze((nanmean(look_at_targ_show2(LearntWords==1))))');
% %plot(squeeze(nanmean(nanmean(look_at_targ_show1_bin1+look_at_targ_show2_bin1+look_at_targ_show1_bin2+look_at_targ_show2_bin2+ ... 
%   %  look_at_targ_show1_bin3+look_at_targ_show2_bin3+look_at_targ_show1_bin4+look_at_targ_show2_bin4)/8))'); %
% for i=1:6
%     plot ( (vv(:,i) )/2');
%     hold on
% end
% legend('manu',    'bosa',  'gake',   'regli',    'kaki', 'colat');
% %plot(nanmean(look_at_targ_showB)); %duration of longest look per trial

%hold off
all_objLook= NaN (numSubjects,size(savestate_historyLt,1),size(savestate_historyLt,2));
all_targL= NaN (numSubjects,6,size(savestate_historyLt,2));
all_distL= NaN (numSubjects,6,size(savestate_historyLt,2));
larW=NaN (numSubjects,6,size(savestate_historyLt,2));
ularW=NaN (numSubjects,6,size(savestate_historyLt,2));
for subject=1:numSubjects
    if goodLearners(subject)==1
        all_objLook(subject,:,:) = lookingOn(subject,:,:);
        for kk=1:6
            all_targL(subject,kk,:)  = (look_at_targ_show1(subject,kk,:)+look_at_targ_show2(subject,kk,:))./2;
            all_distL(subject,kk,:)  = (look_at_dist_show1(subject,kk,:)+look_at_dist_show2(subject,kk,:))./2;
            if LearntWords (subject,kk) == 1
              larW(subject,kk,:)  = (look_at_targ_show1(subject,kk,:)+look_at_targ_show2(subject,kk,:))./2;
            end
            if LearntWords (subject,kk) == 0
              ularW(subject,kk,:)  = (look_at_targ_show1(subject,kk,:)+look_at_targ_show2(subject,kk,:))./2;
            end
        end
    end
end
figure (105);%Plot XSIT
plot(squeeze(nanmean(nanmean(all_objLook,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_targL,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_distL,1),2)))
hold on
word_On = floor([500 2000 4500 6000]/scale_factor);%XSIT    
word_Off = floor([1500 3000 5500 7000]/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Both Objects','Target Object','Distractor');
title('XSIT');
xlabel('time');
ylabel('Strong Learners: proportion of subjects looking');

figure (106);%Plot XTRAP
plot(squeeze(nanmean(nanmean(all_objLook,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_targL,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_distL,1),2)))
hold on
word_On = floor([0 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
word_Off = floor((745+[0 1800 3500 5200 6900])/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');vline(word_On(5),'r-.','On');vline(word_Off(5),'g:','Off');
legend('Both Objects','Target Object','Distractor');
title('XTRAP');
xlabel('time');
ylabel('Strong Learners: proportion of subjects looking at');


figure (111);%Plot XSIT
plot(squeeze(nanmean(nanmean(larW,1),2)))
hold on 
plot(squeeze(nanmean(nanmean(ularW,1),2)))
hold on 
word_On = floor([500 2000 4500 6000]/scale_factor);%XSIT    
word_Off = floor([1500 3000 5500 7000]/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Learnt words',    'Unlearnt words');
title('XSIT');
xlabel('time');
ylabel('Strong Learners: proportion of subjects looking');

figure (112);%Plot XTRAP
plot(squeeze(nanmean(nanmean(larW,1),2)))
hold on 
plot(squeeze(nanmean(nanmean(ularW,1),2)))
hold on 
word_On = floor([0 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
word_Off = floor((745+[0 1800 3500 5200 6900])/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');vline(word_On(5),'r-.','On');vline(word_Off(5),'g:','Off');
legend('Learnt words',    'Unlearnt words');
title('XTRAP');
xlabel('time');
ylabel('Strong Learners: proportion of subjects looking');
 
all_objLook= NaN (numSubjects,size(savestate_historyLt,1),size(savestate_historyLt,2));
all_targL= NaN (numSubjects,6,size(savestate_historyLt,2));
all_distL= NaN (numSubjects,6,size(savestate_historyLt,2));
larW=NaN (numSubjects,6,size(savestate_historyLt,2));
ularW=NaN (numSubjects,6,size(savestate_historyLt,2));
for subject=1:numSubjects
    if goodLearners(subject)==0
        all_objLook(subject,:,:) = lookingOn(subject,:,:);
        for kk=1:6
            all_targL(subject,kk,:)  = (look_at_targ_show1(subject,kk,:)+look_at_targ_show2(subject,kk,:))./2;
            all_distL(subject,kk,:)  = (look_at_dist_show1(subject,kk,:)+look_at_dist_show2(subject,kk,:))./2;
            if LearntWords (subject,kk) == 1
              larW(subject,kk,:)  = (look_at_targ_show1(subject,kk,:)+look_at_targ_show2(subject,kk,:))./2;
            end
            if LearntWords (subject,kk) == 0
              ularW(subject,kk,:)  = (look_at_targ_show1(subject,kk,:)+look_at_targ_show2(subject,kk,:))./2;
            end
        end
    end
end
figure (1212);%Plot XSIT
plot(squeeze(nanmean(nanmean(all_objLook,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_targL,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_distL,1),2)))
hold on
word_On = floor([500 2000 4500 6000]/scale_factor);%XSIT    
word_Off = floor([1500 3000 5500 7000]/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Both Objects','Target Object','Distractor');
title('XSIT');
xlabel('time');
ylabel('Weak Learners: proportion of subjects looking');

figure (1113);%Plot XTRAP
plot(squeeze(nanmean(nanmean(all_objLook,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_targL,1),2)))
hold on
plot(squeeze(nanmean(nanmean(all_distL,1),2)))
hold on
word_On = floor([0 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
word_Off = floor((745+[0 1800 3500 5200 6900])/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');vline(word_On(5),'r-.','On');vline(word_Off(5),'g:','Off');
legend('Both Objects','Target Object','Distractor');
title('XTRAP');
xlabel('time');
ylabel('Weak Learners: proportion of subjects looking at');

figure (113);%Plot XSIT
plot(squeeze(nanmean(nanmean(larW,1),2)))
hold on 
plot(squeeze(nanmean(nanmean(ularW,1),2)))
hold on 
word_On = floor([500 2000 4500 6000]/scale_factor);%XSIT    
word_Off = floor([1500 3000 5500 7000]/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Learnt words',    'Unlearnt words');
xlabel('time');
title('XSIT');
ylabel('Weak Learners: proportion of subjects looking');

figure (114);%Plot 
plot(squeeze(nanmean(nanmean(larW,1),2)))
hold on 
plot(squeeze(nanmean(nanmean(ularW,1),2)))
hold on 
word_On = floor([0 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
word_Off = floor((745+[0 1800 3500 5200 6900])/scale_factor);
vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');vline(word_On(5),'r-.','On');vline(word_Off(5),'g:','Off');
legend('Learnt words',    'Unlearnt words');
xlabel('time');
title('XTRAP');
ylabel('Weak Learners: proportion of subjects looking');

%%% CHANGE POSSIBLY OFF_LOOKING REMOVE

        %%% structure in word_On = floor([500 2000 4500 6000]% starting from 250ms before it    
        %word_Off = floor([1500 3000 5500 7000]% ending  250ms after it
%         bin1_ON  = floor(250/scale_factor(subject));  bin2_ON  = floor(1750/scale_factor(subject)); bin3_ON  = floor(4250/scale_factor(subject)); bin4_ON  = floor(5750/scale_factor(subject)); %time bins
%         bin1_OFF = floor(1750/scale_factor(subject)); bin2_OFF = floor(3250/scale_factor(subject)); bin3_OFF = floor(5750/scale_factor(subject)); bin4_OFF = floor(7250/scale_factor(subject));

xx=['Avg Looking time per test trial is ',num2str(mean(mean(targLookTime+dstrLookTime))*(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Looking time to Target per test trial  is ',num2str(mean(mean(targLookTime))*(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Looking time to Distractor per test trial per is ',num2str(mean(mean(dstrLookTime))*(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Proportion of time looking correctly (Target/Total) is ',num2str(mean(sum (correct_proportion,2)))]; disp(xx);
xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);% 
xx=['Avg # of Words Learnt by Strong Learners is ',num2str(mean(sum(LearntWords(goodLearners()==1,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt by Weak Learners is ',num2str(mean(sum(LearntWords(goodLearners()==0,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt is ',num2str(mean(sum(LearntWords(),2)))]; disp(xx);%
% stanard deviation std(sum(LearntWords(goodLearners()==1,:),2))



%% TRAINING CONDITION ANALYSIS trials analysis NEW
vis_On=1;
vis_Off=(TRAIN_DUR/scale_factor);
tLookAtLearntWords=zeros(numSubjects,3);%3 blocks
tLookAtNonLrntWords=zeros(numSubjects,3);%3 blocks
for subject=1:numSubjects
    savestate_historyL = fliplr(xsit_result.train(subject).historyL(:,vis_On:vis_Off));
    savestate_historyR = fliplr(xsit_result.train(subject).historyR(:,vis_On:vis_Off));    
    nHasLWord=zeros(3,1); tLookAtLearnt=zeros(3,1);
    nHasUWord=zeros(3,1); tLookAtNonLrnt=zeros(3,1);
    for kk=1:nObjects
        if (LearntWords(subject,kk)==1)%if the subject has learnt the kk word
           for block=1:3
                for tr=1:10%each trial may be 10 only
                    tt=10*(block-1)+tr;
                    if(xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+1))%object exists in trial name1
                        nHasLWord(block)=nHasLWord(block)+1;
                        s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                        if (strcmp(s1,'P'))% on left
                            tLookAtLearnt(block)=tLookAtLearnt(block)+sum(savestate_historyL(tt,:));%add time historyL
                        elseif (strcmp(s1,'X'))
                            tLookAtLearnt(block)=tLookAtLearnt(block)+sum(savestate_historyR(tt,:));% add time historyR
                        end
                    elseif (xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+2))%object exists in trial
                        nHasLWord(block)=nHasLWord(block)+1;
                        s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                        if (strcmp(s1,'X'))% on left
                             tLookAtLearnt(block)=tLookAtLearnt(block)+sum(savestate_historyL(tt,:));%add time historyL
                        elseif (strcmp(s1,'P'))
                            tLookAtLearnt(block)=tLookAtLearnt(block)+sum(savestate_historyR(tt,:));% add time historyR
                        end
                    end
                end
           end
        elseif (LearntWords(subject,kk)==0)%if word is NOT learnt
           for block=1:3
                for tr=1:10%each trial may be 10 only
                    tt=10*(block-1)+tr;
                    if(xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+1))%object exists in trial
                        nHasUWord(block)=nHasUWord(block)+1;
                        s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                        if (strcmp(s1,'P'))% on left
                            tLookAtNonLrnt(block)=tLookAtNonLrnt(block)+sum(savestate_historyL(tt,:));%add time historyL
                        elseif (strcmp(s1,'X'))
                            tLookAtNonLrnt(block)=tLookAtNonLrnt(block)+sum(savestate_historyR(tt,:));% add time historyR
                        end
                    elseif (xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+2))%object exists in trial
                        nHasUWord(block)=nHasUWord(block)+1;
                        s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
                        if (strcmp(s1,'X'))% on left
                             tLookAtNonLrnt(block)=tLookAtNonLrnt(block)+sum(savestate_historyL(tt,:));%add time historyL
                        elseif (strcmp(s1,'P'))
                            tLookAtNonLrnt(block)=tLookAtNonLrnt(block)+sum(savestate_historyR(tt,:));% add time historyR
                        end
                    end
                end
           end
        else
            disp('Error in LearntWords array: supurious data ');
        end
     end
 
    for block=1:3
        tLookAtLearntWords(subject,block)=tLookAtLearnt(block)./nHasLWord(block);%dividi time by occurances
        tLookAtNonLrntWords(subject,block)=tLookAtNonLrnt(block)./nHasUWord(block);%dividi time by occurances
        if isnan(tLookAtLearntWords(subject,block)), tLookAtLearntWords(subject,block)=0;end
        if isnan(tLookAtNonLrntWords(subject,block)), tLookAtNonLrntWords(subject,block)=0;end
    end
end

strongLtime=mean(tLookAtLearntWords((goodLearners()==1),:))*(scale_factor);
SLerr=std(tLookAtLearntWords((goodLearners()==1),:))*(scale_factor);

strongNLtime=mean(tLookAtNonLrntWords((goodLearners()==1),:))*(scale_factor);
SNLerr=std(tLookAtNonLrntWords((goodLearners()==1),:))*(scale_factor);

weakLtime=mean(tLookAtLearntWords((goodLearners()==0),:))*(scale_factor);
WLerr=std(tLookAtLearntWords((goodLearners()==0),:))*(scale_factor);

weakNLtime=mean(tLookAtNonLrntWords((goodLearners()==0),:))*(scale_factor);
WNLerr=std(tLookAtNonLrntWords((goodLearners()==0),:))*(scale_factor);


% figure(8)% Plot looking time at leant vs non-learnt words for Strong learners
% blockNames={'1to10','11to20','21to30'};%
% sts = [ strongLtime(1) strongNLtime(1); strongLtime(2) strongNLtime(2);strongLtime(3) strongNLtime(3)];
% errY =[ SLerr(1) SNLerr(1); SLerr(2) SNLerr(2); SLerr(3) SNLerr(3)];
% barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% title ('Strong Learners');
% ylabel('Looking time (ms)');
% legend('Avg Learned Words','Avg NonLearned Words');
% 
% figure(9)% Plot looking time at leant vs non-learnt words for Weak learners
% blockNames={'1to10','11to20','21to30'};%
% sts = [ weakLtime(1) weakNLtime(1); weakLtime(2) weakNLtime(2);weakLtime(3) weakNLtime(3)];
% errY =[ WLerr(1) WNLerr(1); WLerr(2) WNLerr(2); WLerr(3) WNLerr(3)];
% barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% title ('Weak Learners');
% ylabel('Looking time (ms)');
% legend('Avg Learnt Words','Avg Non-Learnt Words');
%
%% 
targLookTimeTraining=zeros(numSubjects,size(savestate_historyL,1));%number of training trials
dstrLookTimeTraining=zeros(numSubjects,size(savestate_historyL,1));
totnlooksS=zeros(sum(goodLearners),size(savestate_historyL,1));
meanlookdurS =zeros(sum(goodLearners),size(savestate_historyL,1));
totlonglookdurS = zeros(sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeS=zeros(sum(goodLearners),size(savestate_historyL,1));

totnlooksW=zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
meanlookdurW =zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
totlonglookdurW = zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeW=zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
tLookRepeatedS= zeros(sum(goodLearners),6);
tLookVaryingS=zeros(sum(goodLearners),6);
tLookRepeatedW=zeros(numSubjects-sum(goodLearners),6);
tLookVaryingW=zeros(numSubjects-sum(goodLearners),6);
indS=1;indW=1; summa=0;
for subject=1:numSubjects
    savestate_historyL = fliplr(xsit_result.train(subject).historyL(:,vis_On:vis_Off));
    savestate_historyR = fliplr(xsit_result.train(subject).historyR(:,vis_On:vis_Off)); 
    %%%%% looking during training
    for tr=1:size(savestate_historyL,1)
        s1=char(xsit_result.train(subject).training_pair(tr,2*nFeatures+3));
        word1_ON=500;
        word1_OFF=1500;%1500
        word2_ON=2000;
        word2_OFF=3000;%3000
        if (strcmp(s1,'P'))% on left
           targLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,floor(word1_ON/scale_factor):floor(word1_OFF/scale_factor))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyR(tr,floor(word2_ON/scale_factor):floor(word2_OFF/scale_factor)));% 2nd audio presentation, Looking to target object right side
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,floor(word1_ON/scale_factor):floor(word1_OFF/scale_factor))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyL(tr,floor(word2_ON/scale_factor):floor(word2_OFF/scale_factor)));% 2nd audio presentation, Looking wrong way
        elseif (strcmp(s1,'X'))
           targLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,floor(word1_ON/scale_factor):floor(word1_OFF/scale_factor))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyL(tr,floor(word2_ON/scale_factor):floor(word2_OFF/scale_factor)));% 2nd audio presentation, Looking to target object right side
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,floor(word1_ON/scale_factor):floor(word1_OFF/scale_factor))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyR(tr,floor(word2_ON/scale_factor):floor(word2_OFF/scale_factor)));% 2nd audio presentation, Looking wrong way
            
        end
    end
    
    %%% no of looks and fixation duration calculation
    nlooks=zeros(2,size(savestate_historyL,1)); %L/R %%SAVE US!! 
    longlookdur=zeros(2,size(savestate_historyL,1));
    for side=1:2      
        if side == 1
            ldata = savestate_historyL;
            cdata = savestate_historyR;
        else
            ldata = savestate_historyR;
            cdata = savestate_historyL;
        end                      
        for tr=1:size(ldata,1) %% or size(rdata,1) 
            looking=0;
            templonglookdur=0;
            lastLookEndTime=1;
            for time=1:size(ldata,2) %% or size(rdata,2) 
                if (round(ldata(tr,time)) == 1)%% if model looking to left
                    if looking == 0 %% if the model is not already looking 
                        if sum(round(cdata(tr,lastLookEndTime:time))) > 10 ||  sum(nlooks(side,tr))== 0 %%%% if it was looing to some other object before or no looks so far
                            nlooks(side,tr) = nlooks(side,tr)+1;
                            if templonglookdur > longlookdur(side,tr)
                                longlookdur(side,tr) = templonglookdur;
                            end
                            templonglookdur=0;
                        else
                            templonglookdur= templonglookdur + (time-lastLookEndTime); %%% add the gap between consecutive looks to same side
                        end
                    end
                    looking = 1;
                    templonglookdur = templonglookdur+1;
                else
                    if looking == 1; lastLookEndTime=time;end
                    looking = 0;
                end
            end
%             if (round(ldata(tr,time-1)) == 1)
%                 nlooks(side,tr) = nlooks(side,tr)+1;
%             end
            if templonglookdur > longlookdur(side,tr)
                longlookdur(side,tr) = templonglookdur;
            end
        end   
    end
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

     if(goodLearners(subject)==1)
         for blockz=1:6
              TinA=5*(blockz-1)+1;
              TinB=5*(blockz);
            tLookRepeatedS(indS,blockz)=sum(sum(savestate_historyL(TinA:TinB,:)));
            tLookVaryingS(indS,blockz)=sum(sum(savestate_historyR(TinA:TinB,:)));
          end
        totnlooksS(indS,:)=sum(nlooks,1);
        meanlookdurS(indS,:)=(sum(savestate_historyL')+sum(savestate_historyR'))./totnlooksS(indS,:);
        totlonglookdurS(indS,:)=max(longlookdur,[],1);
        TotalLookTimeS(indS,:)=sum(savestate_historyL')+sum(savestate_historyR');
        indS=indS+1;
     elseif (goodLearners(subject)==0)
         for blockz=1:6
              TinA=5*(blockz-1)+1;
              TinB=5*(blockz);
            tLookRepeatedW(indW,blockz)=sum(sum(savestate_historyL(TinA:TinB,:)));
            tLookVaryingW(indW,blockz)=sum(sum(savestate_historyR(TinA:TinB,:)));
          end
        totnlooksW(indW,:)=sum(nlooks,1);
        meanlookdurW(indW,:)=(sum(savestate_historyL')+sum(savestate_historyR'))./totnlooksW(indW,:);
        totlonglookdurW(indW,:)=max(longlookdur,[],1);
        TotalLookTimeW(indW,:)=sum(savestate_historyL')+sum(savestate_historyR');
        indW=indW+1;
     else
         disp('ERROR goodLearners bad data');
     end
               
end
%%multi file code
figure(11);%Plot Strong vs Weak looking time during a training trial
errorbar(mean(TotalLookTimeS)*scale_factor,std(TotalLookTimeS)*scale_factor,plotStyle{1});%
hold on
errorbar(mean(TotalLookTimeW)*scale_factor,std(TotalLookTimeW)*scale_factor,plotStyle{2+compStyle});%
legend('Strong','Weak');
xlabel('per training trial');
ylabel('total looking time Strong vs Weak learners');
%ylabel('total looking time Strong learners');
%hold off 

figure (12);% Plot number of fixations/looks over training trials
errorbar(mean(totnlooksS),std(totnlooksS),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(mean(totnlooksW),std(totnlooksW),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('number of fixations/looks Strong vs Weak learners');
%ylabel('number of fixations/looks Strong learners');
legend('Strong','Weak');
%ylim([1.5 3.5])
%hold off
summaF=mean(mean(totnlooksS));
xx=['number of fixations/looks Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(totnlooksW));
xx=['number of fixations/looks Weak learners ',num2str(summaF)]; disp(xx);

figure (13);%Plot mean look duration of each fixation
errorbar(mean(meanlookdurS)*scale_factor,std(meanlookdurS)*scale_factor, plotStyle{1});% mean look duration % multiped by timing scale factor
hold on
errorbar(mean(meanlookdurW)*scale_factor,std(meanlookdurW)*scale_factor, plotStyle{2+compStyle});% mean look duration % multiped by timing scale factor
xlabel('per training trial');
ylabel('mean look duration Strong vs Weak learners');
%ylabel('mean look duration Strong learners');
legend('Strong','Weak');
%ylim([1000 2000])
%hold off
summaF=mean(mean(meanlookdurS)*scale_factor)/1000;
xx=['mean look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanlookdurW)*scale_factor)/1000;
xx=['mean look duration Strong weak learners ',num2str(summaF)]; disp(xx);

figure (14);%Plot duration of longest look per trial
errorbar(mean(totlonglookdurS)*scale_factor,std(totlonglookdurS)*scale_factor,plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(totlonglookdurW)*scale_factor,std(totlonglookdurW)*scale_factor, plotStyle{2+compStyle}); %duration of longest look per trial
xlabel('per training trial');
ylabel('duration of longest look Strong vs Weak learners');
%ylabel('duration of longest look Strong learners');
legend('Strong','Weak');
%hold off
summaF=mean(mean(totlonglookdurS)*scale_factor)/1000;
xx=['Longest look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totlonglookdurW)*scale_factor)/1000;
xx=['Longest look duration Strong weak learners ',num2str(summaF)]; disp(xx);

figure(15);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(mean(targLookTimeTraining((goodLearners()==1),:))*scale_factor,std(targLookTimeTraining((goodLearners()==1),:))*scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTimeTraining((goodLearners()==1),:))*scale_factor,std(dstrLookTimeTraining((goodLearners()==1),:))*scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Training Trial');
ylabel('Strong Learners: Looking Time when words are ON');

figure(16);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(mean(targLookTimeTraining((goodLearners()==0),:))*scale_factor,std(targLookTimeTraining((goodLearners()==0),:))*scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTimeTraining((goodLearners()==0),:))*scale_factor,std(dstrLookTimeTraining((goodLearners()==0),:))*scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Training Trial');
ylabel('Weak Learners: Looking Time when words are ON');

figure (17);%Plot Mean propertion of looking Varying vs repeated C
errorbar(mean(tLookVaryingS./(tLookVaryingS+tLookRepeatedS)),std(tLookVaryingS./(tLookVaryingS+tLookRepeatedS)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeatedS./(tLookVaryingS+tLookRepeatedS)),std(tLookRepeatedS./(tLookVaryingS+tLookRepeatedS)), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('per block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
ylabel('Proportion Looking Time Strong Learners');
%ylabel('duration of longest look Strong learners');
legend('Varying','Repeated');


figure (18);%Plot Mean propertion of looking Varying vs repeated WEAK
errorbar(mean(tLookVaryingW./(tLookVaryingW+tLookRepeatedW)),std(tLookVaryingW./(tLookVaryingW+tLookRepeatedW)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeatedW./(tLookVaryingW+tLookRepeatedW)),std(tLookRepeatedW./(tLookVaryingW+tLookRepeatedW)), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('per block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
ylabel('Proportion Looking Time Weak Learners');
%ylabel('duration of longest look Strong learners');
legend('Varying','Repeated');


summaF=mean(mean([TotalLookTimeS; TotalLookTimeW])) *(scale_factor/1000);
xx=['Avg looking time per training trial is ',num2str(summaF)]; disp(xx);
var1= mean(mean(TotalLookTimeS))*scale_factor;
xx=['Avg looking time per training trial per Strong learner is ',num2str(var1)]; disp(xx);
var2= mean(mean(TotalLookTimeW))*scale_factor;
xx=['Avg looking time per training trial per Weak learner is ',num2str(var2)]; disp(xx);

figure(51);%Plot mean looking time per TEST trial
sts = [3.04;     2.99; mean(mean([TotalLookTimeS; TotalLookTimeW])) *(scale_factor/1000)];
errY=[0.0;        0.0;  std(mean([TotalLookTimeS; TotalLookTimeW])) *(scale_factor/1000)];
hb=barwitherr(errY,sts,0.6);% Plot with errorbars
set(gca,'xticklabel',{'Smith & Yu 2008 (14 m)'; 'Yu & Smith 2011';'WOLVES Model'},'fontsize',11)
%title ('total looking time per test trial');
ylabel('Looking time (seconds)');
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
        
        
        [maxIn(kk) indIn(kk)] = max(inputMapping(:,cell2mat(xsit_result.train(subject).Words(kk))));
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















