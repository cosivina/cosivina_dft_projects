%% Data Analysis 
clear all; close all;
nObjects=12;nFeatures=2;nTrainTrials=36; nTestTrials=12; scale_factor=8; TRAIN_DUR=4000; TEST_DUR=1000; MIN_LOOK_DURATION=200/scale_factor; %check with auto file 160
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};
compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];

xsit_result.train=[];xsit_result.test=[];
%base_tauM_stim525_wmc9_hwmf15_Vlach_Johnson_2013_
simName = 'wfChanges_conwmf0_hwmc1_fix48_Test1k_12h3000_Vlach_Johnson_2013_'     %goodlooking3-12-Nov-18_Vlach_Johnson_2013_'          % base2011_Vlach_Johnson_2013_   %wmc2_sim_Vlach_Johnson_2013_
OutName = [simName,'results.mat']; xsit_result1 = load (OutName);
xsit_result.train = [xsit_result.train ; xsit_result1.train];  
xsit_result.test = [xsit_result.test ; xsit_result1.test];
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
% word_On = floor([0 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
% word_Off = floor((745+[0 1800 3500 5200 6900])/scale_factor);word_Len=floor(1000/scale_factor);
vis_On = 1;vis_Off = floor(TEST_DUR/scale_factor);

correct_proportion=zeros(numSubjects,2);
targLookTime=zeros(numSubjects,nTestTrials);
dstrLookTime=zeros(numSubjects,nTestTrials);
goodLearners=NaN(numSubjects,1);
LearntWords= NaN(numSubjects,nObjects);
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
MassedTargLookTime=zeros(numSubjects,1);
MassedDistLookTime=zeros(numSubjects,1);
IleavedTargLookTime=zeros(numSubjects,1);
IleavedDistLookTime=zeros(numSubjects,1);
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
        
        for ss= 1:nObjects/2
           if (xsit_result.train(subject).Words{ss}== xsit_result.test(subject).test_pair(trt,2*nFeatures+1))
               s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
               if ( strcmp(s1,'L')) 
                   MassedTargLookTime(subject)=MassedTargLookTime(subject)+lLook;
                   MassedDistLookTime(subject)= MassedDistLookTime(subject)+rLook;
               elseif ( strcmp(s1,'R'))
                   MassedTargLookTime(subject)= MassedTargLookTime(subject) + rLook;
                   MassedDistLookTime(subject)= MassedDistLookTime(subject)+ lLook;
               else
                   disp('ERROR reading test_pair char');
               end
           end
        end

        for ss= (nObjects/2)+1:nObjects
           if (xsit_result.train(subject).Words{ss}== xsit_result.test(subject).test_pair(trt,2*nFeatures+1))
               s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
               if ( strcmp(s1,'L')) 
                   IleavedTargLookTime(subject)= IleavedTargLookTime(subject) + lLook;
                   IleavedDistLookTime(subject)= IleavedDistLookTime(subject) + rLook;
                   lcorrect=lcorrect+lLook/(lLook+rLook);
               elseif ( strcmp(s1,'R'))
                   IleavedTargLookTime(subject)= IleavedTargLookTime(subject) + rLook;
                   IleavedDistLookTime(subject)= IleavedDistLookTime(subject) + lLook;
                   rcorrect=rcorrect+rLook/(lLook+rLook);
               else
                   disp('ERROR reading test_pair char');
               end
           end
        end%%massed/interleaved
    end%% trials loop
    
    for kk=1:nObjects
        if (targWord(kk)>dstrWord(kk))
            LearntWords(subject,kk)=1;
        else
            LearntWords(subject,kk)=0;
        end
    end  
    if (mean(targLookTime(subject,:)) > mean(dstrLookTime(subject,:)))
        goodLearners(subject)=1;
    else
        goodLearners(subject)=0;
    end    
end
correct_proportion=targLookTime./(targLookTime+dstrLookTime);

disp('t-test statistics between Target and Distractor Looking');
[h,p,ci,stats] = ttest(mean(targLookTime,2),mean(dstrLookTime,2),'Tail','right')

disp('t-test statistics for massed target looking');
[h,p,ci,stats] = ttest((MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)),0.5,'Tail','right')

disp('t-test statistics for interleaved target looking');
[h,p,ci,stats] = ttest((IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime)),0.5,'Tail','right')

disp('t-test statistics for words learnt');
[h,p,ci,stats] = ttest(sum(LearntWords,2),6,'Tail','right')

disp('t-test statistics for massed words learnt');
[h,p,ci,stats] = ttest(sum(LearntWords(:,1:6),2),3,'Tail','right')

disp('t-test statistics for interleaved words learnt');
[h,p,ci,stats] = ttest(sum(LearntWords(:,7:12),2),3,'Tail','right')




perform_CSL = sum(LearntWords,2);
freq_CSL=zeros(1,nObjects);
for subject=1:numSubjects
    for trt=1:nObjects
        if trt==perform_CSL(subject)
            freq_CSL(trt)= freq_CSL(trt)+1;
        end
    end
end
mean(perform_CSL)
std(perform_CSL)./sqrt(length(perform_CSL))

xx=['Avg Looking time per test trial is ',num2str(mean(mean(targLookTime+dstrLookTime))*(scale_factor/1000))]; disp(xx);
xx=['Looking time to Target per test trial  is ',num2str(mean(mean(targLookTime))*(scale_factor/1000))]; disp(xx);
xx=['Looking time to Distractor per test trial per is ',num2str(mean(mean(dstrLookTime))*(scale_factor/1000))]; disp(xx);
xx=['Proportion of time looking correctly (Target/Total) is ',num2str(nanmean(nanmean(correct_proportion,2)))]; disp(xx);
xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);% 
xx=['Avg # of Words Learnt by Strong Learners is ',num2str(mean(sum(LearntWords(goodLearners()==1,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt by Weak Learners is ',num2str(mean(sum(LearntWords(goodLearners()==0,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt is ',num2str(mean(sum(LearntWords(),2)))]; disp(xx);%
xx=['Avg # of Massed Words Learnt is ',num2str(mean(sum(LearntWords(:,1:6),2)))]; disp(xx);%
xx=['Avg # of Interleaved Words Learnt is ',num2str(mean(sum(LearntWords(:,7:12),2)))]; disp(xx);%
% stanard deviation std(sum(LearntWords(goodLearners()==1,:),2))

figure(1)% Plot Massed vs Interleaved looking time during test trial
blockNames={'Massed'; 'Interleaved'};
sts = [0.55 mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) ;  0.48 mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))-0.01];
errY =[0.04 std(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime))./sqrt(length(MassedTargLookTime));    0.04 std(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))./sqrt(length(IleavedTargLookTime))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames, 'fontsize',18)
ylabel('Prop. Looking Time to Target');
legend('Vlach & Johnson 2013 (16 m)', 'WOLVES Model');
hold on
hline(0.05,'k:');
grid on
ylim([0 0.6]);

rmse=sqrt(((mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) - 0.55)^2 + (mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime)) - 0.48)^2)./2);
disp(['RMSE error: Massed & Interleaved 16 month is ',num2str(rmse)]);
pe = mean ([abs(mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) - 0.55)/0.55*100 abs(mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime)) - 0.48)/0.48*100 ] )


figure(101)% Plot Massed vs Interleaved looking time during test trial
blockNames={'Massed'; 'Interleaved'};
sts = [0.54 mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime))-0.01 ;  0.55 mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
errY =[0.04 std(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime))./sqrt(length(MassedTargLookTime));    0.04 std(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))./sqrt(length(IleavedTargLookTime))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames, 'fontsize',18)
ylabel('Mean prop of looking time to target');
legend('Vlach & Johnson 2013 (20 m)', 'WOLVES Model');
grid on
%ylim([0 0.6]);

rmse=sqrt(((mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) - 0.54)^2 + (mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime)) - 0.55)^2)./2);
disp(['RMSE error: Massed & Interleaved 20 month is ',num2str(rmse)]);
pe = mean ([abs(mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) - 0.54)/0.54*100 abs(mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime)) - 0.55)/0.55*100 ] )


figure(102)% Plot Massed vs Interleaved looking time during test trial
blockNames={'Massed'; 'Interleaved'};
%sts = [3.55 3.8 ;  4.11 4.3];
%errY =[1.0  0.7 ;  0.9 0.7 ];
sts = [3.556 mean(sum(LearntWords(:,1:6),2)) ;  4.111 mean(sum(LearntWords(:,7:12),2)) ];
errY =[1.042./sqrt(18)  std(sum(LearntWords(:,1:6),2))./sqrt(length(LearntWords(:,1:6))) ;  0.900./sqrt(18)  std(sum(LearntWords(:,7:12),2))./sqrt(length(LearntWords(:,7:12))) ];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames, 'fontsize',16)
ylabel('Mean prop of looking time to target');
legend('Vlach & Johnson 2019 (47-58 m)', 'WOLVES Model');
grid on
hold on
hline(3)
%ylim([0 0.6]);
rmse=sqrt(((3.556-mean(sum(LearntWords(:,1:6),2)))^2 + (4.111-mean(sum(LearntWords(:,7:12),2)))^2)./2);
disp(['RMSE error: Massed & Interleaved (47-58) month is ',num2str(rmse)]);
pe = mean ([abs(3.556-mean(sum(LearntWords(:,1:6),2)))/3.556*100 abs(4.111-mean(sum(LearntWords(:,7:12),2)))/4.111*100 ] )



barColorMap(1,:) = [.2 .71 .3];	% Green Color for segment 1.
barColorMap(2,:) = [.25 .55 .79];	% Blue Color for segment 2.
barColorMap(3,:) = [.9 .1 .14];	% Red Color for segment 3.
barColorMap(4,:) = [.9 .9 .14];	% Yellow Color for segment 4.
   
figure(100)% Plot total looking time during test trial
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
sts = [4/6; 3.5/6; mean(sum(LearntWords(),2))/nObjects];
errY= [0; 0; 0];
hb=barwitherr(errY,sts,0.6);% Plot with errorbars
set(gca,'xticklabel',{'Smith & Yu 2008 (14 m)'; 'Yu & Smith 2011';'WOLVES Model'},'fontsize',11);
ylabel('Proportion of words learnt');
ylim([0 1]);
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





% strongLtime=mean(tLookAtLearntWords((goodLearners()==1),:))/(scale_factor);
% SLerr=std(tLookAtLearntWords((goodLearners()==1),:))/(scale_factor);
% 
% strongNLtime=mean(tLookAtNonLrntWords((goodLearners()==1),:))/(scale_factor);
% SNLerr=std(tLookAtNonLrntWords((goodLearners()==1),:))/(scale_factor);
% 
% weakLtime=mean(tLookAtLearntWords((goodLearners()==0),:))/(scale_factor);
% WLerr=std(tLookAtLearntWords((goodLearners()==0),:))/(scale_factor);
% 
% weakNLtime=mean(tLookAtNonLrntWords((goodLearners()==0),:))/(scale_factor);
% WNLerr=std(tLookAtNonLrntWords((goodLearners()==0),:))/(scale_factor);
% 
% 
% figure(8)% Plot looking time at leant vs non-learnt words for Strong learners
% blockNames={'1to6','7to12','13to18','19to24','25to30','31to36'};%
% sts = [ strongLtime(1) strongNLtime(1); strongLtime(2) strongNLtime(2);strongLtime(3) strongNLtime(3);strongLtime(4) strongNLtime(4); strongLtime(5) strongNLtime(5);strongLtime(6) strongNLtime(6)];
% errY =[ SLerr(1) SNLerr(1); SLerr(2) SNLerr(2); SLerr(3) SNLerr(3); SLerr(4) SNLerr(4); SLerr(5) SNLerr(5); SLerr(6) SNLerr(6)];
% barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% title ('Strong Learners');
% ylabel('Looking time (ms)');
% legend('Avg Learned Words','Avg NonLearned Words');
% 
% figure(9)% Plot looking time at leant vs non-learnt words for Weak learners
% blockNames={'1to6','7to12','13to18','19to24','25to30','31to36'};%
% sts = [ weakLtime(1) weakNLtime(1); weakLtime(2) weakNLtime(2);weakLtime(3) weakNLtime(3);weakLtime(4) weakNLtime(4); weakLtime(5) weakNLtime(5);weakLtime(6) weakNLtime(6)];
% errY =[ WLerr(1) WNLerr(1); WLerr(2) WNLerr(2); WLerr(3) WNLerr(3); WLerr(4) WNLerr(4); WLerr(5) WNLerr(5); WLerr(6) WNLerr(6)];
% barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% title ('Weak Learners');
% ylabel('Looking time (ms)');
% legend('Avg Learnt Words','Avg Non-Learnt Words');
%

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
                    % COLLECT TIME DATA for each group  
                        look_at_targ_obj(subject,kk,:)=(round(savestate_historyLt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';%+savestate_historyOt(tt,1:8000)
                        look_at_dist_obj(subject,kk,:)=(round(savestate_historyRt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';
                elseif (strcmp(s1,'R')) % referent on right
                        look_at_targ_obj(subject,kk,:)=(round(savestate_historyRt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';%add time historyRt
                        look_at_dist_obj(subject,kk,:)=(round(savestate_historyLt(tt,:))./(round(savestate_historyLt(tt,:))+round(savestate_historyRt(tt,:))+savestate_historyOt(tt,:)))';     
                end
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
all_objLook= NaN (numSubjects,nTestTrials,TEST_DUR/scale_factor);
all_targL= NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
all_distL= NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
larW=NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
ularW=NaN (numSubjects,nObjects,TEST_DUR/scale_factor);
for subject=1:numSubjects
        all_objLook(subject,:,:) = lookingOn(subject,:,:);
        for kk=1:nObjects
            all_targL(subject,kk,:)  = look_at_targ_obj(subject,kk,:);
            all_distL(subject,kk,:)  = look_at_dist_obj(subject,kk,:);
            if LearntWords (subject,kk) == 1
              larW(subject,kk,:)  = look_at_targ_obj(subject,kk,:);
            end
            if LearntWords (subject,kk) == 0
              ularW(subject,kk,:)  = look_at_targ_obj(subject,kk,:);
            end
        end
        massedWords_LTarg(subject,:) = squeeze( mean (look_at_targ_obj(subject,1:(nObjects/2),:)));
        inleavedWords_LTarg(subject,:) = mean (look_at_targ_obj(subject,(nObjects/2+1):nObjects,:));
        massedWords_LDist(subject,:) = mean (look_at_dist_obj(subject,1:(nObjects/2),:));
        inleavedWords_LDist(subject,:) = mean (look_at_dist_obj(subject,(nObjects/2+1):nObjects,:));
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

figure (1003);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
%plot((squeeze(nanmean(nanmean(all_objLook(goodLearners()==1,:)),2))));
plot((squeeze(nanmean(massedWords_LTarg+massedWords_LDist))));
hold on
%plot((squeeze(nanmean(massedWords_LTarg(goodLearners()==1,:)))));
plot((squeeze(nanmean(massedWords_LTarg))));
hold on
%plot((squeeze(nanmean(massedWords_LDist(goodLearners()==1,:)))));
plot((squeeze(nanmean(massedWords_LDist))));
hold on
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Both Objects','Target Object','Distractor');
xlabel('time');
ylabel('Proportion of subjects looking to Massed ');
ylim([0 1]);

figure (1013);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
%plot((squeeze(nanmean(nanmean(inleavedWords_LTarg(goodLearners()==1,:)),2))));
plot((squeeze(nanmean(inleavedWords_LTarg+inleavedWords_LDist))));
hold on
%plot((squeeze(nanmean(inleavedWords_LTarg(goodLearners()==1,:)))));
plot((squeeze(nanmean(inleavedWords_LTarg))));
hold on
%plot((squeeze(nanmean(inleavedWords_LDist(goodLearners()==1,:)))));
plot((squeeze(nanmean(inleavedWords_LDist))));
hold on
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Both Objects','Target Object','Distractor');
xlabel('time');
ylabel('Proportion of subjects looking to Interleaved ');
ylim([0 1]);

figure (105);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot((squeeze(nanmean(nanmean(all_objLook(goodLearners()==1,:,:),1),2))));
hold on
plot((squeeze(nanmean(nanmean(all_targL(goodLearners()==1,:,:),1),2))));
hold on
plot((squeeze(nanmean(nanmean(all_distL(goodLearners()==1,:,:),1),2))));
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
ylabel('Strong Learners: proportion of subjects looking');
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
ylabel('Strong Learners: proportion of subjects looking');
ylim([0 1]);
 




%% TRAINING CONDITION ANALYSIS trials analysis
word_On = floor([500 2000]/scale_factor);%XSIT            
word_Off = floor([1500 3000 ]/scale_factor);word_Len=floor(1000/scale_factor);
% word_On = floor([0 1800 3500]/scale_factor);  %% XTRAP    
% word_Off = floor((745+[0 1800 3500])/scale_factor);word_Len=floor(1000/scale_factor);
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
        
%         %%massed interlraved
%         for tt=1:nTrainTrials
%                 if(xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+1))%object exists in trial name1
%                 
%                 end
%         end
                 
     end
 
    for tt=1:nTrainTrials
        if nHasLWord(tt)==0; tLookAtLearnt(tt)=NaN;end % take only those trials wherein the word existed
        if nHasUWord(tt)==0; tLookAtNonLrnt(tt)=NaN;end
        tLookAtLearntWords(subject,tt)=tLookAtLearnt(tt);
        tLookAtNonLrntWords(subject,tt)=tLookAtNonLrnt(tt);
    end
end

figure(8);%Plot looking to learnt vs unlearnt 
errorbar(nanmean(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000,(nanstd(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==1,:))),plotStyle{1});%
hold on
errorbar(nanmean(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000,(nanstd(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==1,:))),plotStyle{2+compStyle});%
xlabel('per training trial');
ylabel('Looking time (s)');
legend('Learned Words','Non-Learned Words');
title('Strong Learners')

figure(9);%Plot looking to learnt vs unlearnt 
errorbar(nanmean(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000,(nanstd(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==0,:))),plotStyle{1});%
hold on
errorbar(nanmean(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000,(nanstd(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==0,:))),plotStyle{2+compStyle});%
xlabel('per training trial');
ylabel('Looking time (s)');
legend('Learned Words','Non-Learned Words');
title('Weak Learners')
%%% above nan lengths are NOT CORRECT... length(x(~isnan(x)))
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
    masList=cell2mat(xsit_result.train(subject).Feature1(xsit_result.train(subject).Sequence(1:6)));
    ilvList=cell2mat(xsit_result.train(subject).Feature1(xsit_result.train(subject).Sequence(7:12)));
    %%%%% looking during training
    for tr=1:nTrainTrials
        trial_obj1=xsit_result.train(subject).training_pair(tr,1);       
        %% finding if object on left is massed or interleaved
        if ismember(trial_obj1,masList)
            %leftobject is massed, right is interleaved
            massLookTimeTraining(subject,tr)= sum(savestate_historyL(tr,vis_On:vis_Off));
            itlvLookTimeTraining(subject,tr)= sum(savestate_historyR(tr,vis_On:vis_Off));
        elseif ismember(trial_obj1,ilvList) % only else may be fine too
            %leftobject is interlraved, right is massed
            itlvLookTimeTraining(subject,tr)= sum(savestate_historyL(tr,vis_On:vis_Off));
            massLookTimeTraining(subject,tr)= sum(savestate_historyR(tr,vis_On:vis_Off));              
        end
        
        % looking when words are ON
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

    for blockz=1:nObjects/2
        TinA=((nObjects/2)*(blockz-1)) + 1;
        TinB=(nObjects/2)*(blockz);
        tLookRepeated(subject,blockz)=sum(massLookTimeTraining(subject,TinA:TinB));
        tLookVarying(subject,blockz)=sum(itlvLookTimeTraining(subject,TinA:TinB));
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
ylim([0 2.5]);

figure(16);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(mean(targLookTimeTraining((goodLearners()==0),:))*scale_factor/1000,(std(targLookTimeTraining((goodLearners()==0),:)) *scale_factor/1000)./sqrt(length(targLookTimeTraining((goodLearners()==0),:))),plotStyle{1});%
hold on
errorbar(mean(dstrLookTimeTraining((goodLearners()==0),:))*scale_factor/1000,(std(dstrLookTimeTraining((goodLearners()==0),:)) *scale_factor/1000)./sqrt(length(dstrLookTimeTraining((goodLearners()==0),:))),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Training Trial');
ylabel('Weak Learners: Looking Time when words are ON');
ylim([0 2.5]);

figure (1711);%Plot Mean proportion of looking Varying vs repeated C
rectangle('Position',[0.7,0,5.6,1],'FaceColor',[.1 .5 .5 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[6.7,0,5.6,1],'FaceColor',[.5 .5 .1 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[12.7,0,5.6,1],'FaceColor',[.1 .5 .5 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[18.7,0,5.6,1],'FaceColor',[.5 .5 .1 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[24.7,0,5.6,1],'FaceColor',[.1 .5 .5 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[30.7,0,5.6,1],'FaceColor',[.5 .5 .1 0.3],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
errorbar(mean(itlvLookTimeTraining./(itlvLookTimeTraining+massLookTimeTraining)),std(itlvLookTimeTraining./(itlvLookTimeTraining+massLookTimeTraining))./sqrt(length(itlvLookTimeTraining)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(massLookTimeTraining./(massLookTimeTraining+itlvLookTimeTraining)),std(massLookTimeTraining./(massLookTimeTraining+itlvLookTimeTraining))./sqrt(length(massLookTimeTraining)), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('training trial');
ylim([0.2 0.8]);
xlim([0 36]);
ylabel('Proportion looking time');
%ylabel('duration of longest look Strong learners');
legend('Interleaved Objects','Massed Objects');
set (gca, 'FontSize',12);
title('Habituation over training');


figure (17);%Plot Mean proportion of looking Varying vs repeated C
errorbar(mean(tLookVarying./(tLookVarying+tLookRepeated)),std(tLookVarying./(tLookVarying+tLookRepeated))./sqrt(length(tLookVarying)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeated./(tLookVarying+tLookRepeated)),std(tLookRepeated./(tLookVarying+tLookRepeated))./sqrt(length(tLookVarying)), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('trial block');
xlim([0.5 6.5]);
ylim([0.3 0.7]);
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
Wrong_inTraceM= zeros(numSubjects,6); Correct_inTraceM= zeros(numSubjects,6);
Wrong_inTraceI= zeros(numSubjects,6); Correct_inTraceI= zeros(numSubjects,6);
for subject=1:numSubjects    
    inputMapping1=squeeze(xsit_result.train(subject).hwf(1,:,:));
    inputMapping2=squeeze(xsit_result.train(subject).hwf(2,:,:));
    for kk=1:nObjects
        xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
        xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
        yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
    end
    C_inTrM=0;W_inTrM=0;
    for kk=1:nObjects/2
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
%         a_cv=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cv=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
%         C_inTr= C_inTr+ nanmean([nanmean(a_cv(a_cv>0.001))  nanmean(b_cv(b_cv>0.001))]);
        a_cv=inputMapping1(xx1(kk),yy(kk));b_cv=inputMapping2(xx2(kk),yy(kk));
        C_inTrM(kk)=nanmean([a_cv b_cv]);
        inputMapping1(xx1(kk),yy(kk))=0;
        inputMapping2(xx2(kk),yy(kk))=0;
        for jj=1:6
            inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
            inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
        end
        a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
        W_inTrM(kk) = nanmean([nanmean(a_in(a_in>0.001)) nanmean(b_in(b_in>0.001))]);
    end
    Correct_inTraceM(subject,:)=C_inTrM;
    Wrong_inTraceM(subject,:)=W_inTrM;
    InCorr_assocsM(subject,:)=nanmean([as_count1(1:6)-1; as_count2(1:6)-1]);
    EntropyTraceM(subject)= nanmean( [entropy(inputMapping1); entropy(inputMapping2)] ); 
        
    C_inTrI=0;W_inTrI=0;
    for kk=7:nObjects
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
%         a_cv=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cv=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
%         C_inTr= C_inTr+ nanmean([nanmean(a_cv(a_cv>0.001))  nanmean(b_cv(b_cv>0.001))]);
        a_cv=inputMapping1(xx1(kk),yy(kk));b_cv=inputMapping2(xx2(kk),yy(kk));
        C_inTrI(kk-6)= nanmean([a_cv b_cv]);
        inputMapping1(xx1(kk),yy(kk))=0;
        inputMapping2(xx2(kk),yy(kk))=0;
        for jj=1:6
            inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
            inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
        end
        a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
        W_inTrI(kk-6) = nanmean([nanmean(a_in(a_in>0.001)) nanmean(b_in(b_in>0.001))]);
    end
    Correct_inTraceI(subject,:)=C_inTrI;
    Wrong_inTraceI(subject,:)=W_inTrI;
    InCorr_assocsI(subject,:)=nanmean([as_count1(7:12)-1; as_count2(7:12)-1]);
    EntropyTraceI(subject)= nanmean( [entropy(inputMapping1) entropy(inputMapping2)] ); 
end

figure(21);%Entropy
blockNames={'Massed';'Interleaved'};
sts = [mean(EntropyTraceM); mean(EntropyTraceI)  ];
errY =[std(EntropyTraceM)/sqrt(length(EntropyTraceM)); std(EntropyTraceI)/sqrt(length(EntropyTraceI))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Entropy in the traces');
%ylabel('Looking time per test trial');
%ylim([0 4]);

figure(212);%My own Entropy: No of incorrect traces
blockNames={'Massed';'Interleaved'};
sts = [mean(mean(InCorr_assocsM)); mean(mean(InCorr_assocsI))  ];
errY =[std(mean(InCorr_assocsM,2))/sqrt(length(InCorr_assocsM)); std(mean(InCorr_assocsI,2))/sqrt(length(InCorr_assocsI))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Proportion of incorrect assocs in the traces');
ylabel('# of incorrect assocs/# of words (6)');


figure(22);%My own Entropy: No of incorrect traces
%blockNames={'Massed';'Interleaved'};
errorbar(nanmean(InCorr_assocsM),std(InCorr_assocsM)/sqrt(length(InCorr_assocsM)))
hold on
errorbar(nanmean(InCorr_assocsI),std(InCorr_assocsI)/sqrt(length(InCorr_assocsI)))
%set(gca,'xticklabel',blockNames);
legend('Massed','Interleaved');
title ('Proportion of incorrect assocs in the traces');
%ylabel('# of incorrect assocs/# of words (6)');

figure(221);%trace strength
blockNames={'Massed';'Interleaved'};
sts = [nanmean(nanmean(Correct_inTraceM,2)); nanmean(nanmean(Correct_inTraceI,2))  ];
errY =[nanstd(nanmean(Correct_inTraceM,2))/sqrt(length(Correct_inTraceM)); nanstd(nanmean(Correct_inTraceI,2))/sqrt(length(Correct_inTraceI))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Strength of correct assocs in the traces');

figure(232)
errorbar(nanmean(Correct_inTraceM),std(Correct_inTraceM)/sqrt(length(Correct_inTraceM)),plotStyle{5})
hold on
errorbar(nanmean(Correct_inTraceI),std(Correct_inTraceI)/sqrt(length(Correct_inTraceI)),plotStyle{2})
%set(gca,'xticklabel',blockNames);
legend('Massed','Interleaved');
%title ('Strength of correct assocs in the traces');
ylim([0 0.25]);
ylabel('Association Trace Strength');
xlabel('Word-Object Pair (Presented in Temporal Order)');
hold on
hline(0.06,'k:');
grid on


% figure(231);%My own Entropy: No of incorrect traces
% legend('Massed','Interleaved');
% sts = [nanmean(Wrong_inTraceM); nanmean(Wrong_inTraceI)  ];
% errY =[nanstd(Wrong_inTraceM)/sqrt(length(Wrong_inTraceM)); nanstd(Wrong_inTraceI)/sqrt(length(Wrong_inTraceI))];
% b=barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% title ('Strength of INCORRECT assocs in the traces');



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
