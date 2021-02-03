%% Data Analysis 
clear all; close all;
nObjects=6;nFeatures=2; scale_factor=1/8;  %check with auto file
plotStyle = {'k-','b-','g-','c-','r-','m-','y-','k-.','b-.','g-.','c-.','r-.','m-.','y-.','b:','w.'};
compStyle=7;
blockNames{50}=[];sts{50}=[];errY{50}=[];
% xsit_result.Names={[1] [2] [3] [4] [5] [6]};
%xsit_result.train =load('5atn_13wf_xsitNF_surprise3_1_train.mat');
%xsit_result.test =load('5atn_13wf_xsitNF_surprise3_1_test.mat');
%paramValues=[0.0 0.5 1.0 2.0 4.0 6.0];%hwf->wf
%paramValues=[1.0 2.0 3.0 5.0 8.0 ];%wf->wf
paramValues=[5.0 5.3 5.6 5.9];%wf->antn_f 0.0 1.0 2.0 
for param=1:1 
parVal= paramValues(param);
legendInfo{2*param -1}= [plotStyle{param}  num2str(parVal)];
legendInfo{2*param}= [plotStyle{param+compStyle} num2str(parVal)];
 simName = ['Hab60Atn53atnc2atnsa'  num2str(parVal*10) '_Smith_Yu_2013_']; %%simName = 'B_xsitNew_16_sub_';
%simName = 'h45a80aa50stim65_Smith_Yu_2008_';
%simName = ['Ahwfwf2_wfatnf1.3_xsitNF_sub_']; %%simName = 'B_xsitNew_16_sub_';
%if param==1, simName = 'weak10hab05_xtrapNF_', elseif param==2, simName = 'Spaceweak10hab_xtrapNF_',end 
OutName = [simName,'results.mat']
xsit_result = load (OutName);
numSubjects=size(xsit_result.test,1);
%% Test trials
correct_proportion=zeros(numSubjects,2);
targLookTime=zeros(numSubjects,2*nObjects);
dstrLookTime=zeros(numSubjects,2*nObjects);
goodLearners=999*ones(numSubjects,1);
LearntWords= 999*ones(numSubjects,nObjects);
%word_On = floor([500 2000 4500 6000]*scale_factor);%[0 450 875 1300 ];    
%word_Off = floor([1500 3000 5500 7000]*scale_factor);%[375 650 1075 1500 ];
slice_wordOnset4=1;%6/8;
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
for subject=1:numSubjects
    lcorrect=0;rcorrect=0;targWord=zeros(nObjects,1);dstrWord=zeros(nObjects,1);
    for trt=1:2*nObjects
          lLook= sum( xsit_result.test(subject).historyLt(trt,:));%250:2000 
          rLook= sum( xsit_result.test(subject).historyRt(trt,:));%250:2000
%        lLook= sum( xsit_result.test(subject).historyLt(trt,floor(500*scale_factor):floor(2000*scale_factor)));%250:2000
%        lLook=lLook+sum( xsit_result.test(subject).historyLt(trt,floor(2000*scale_factor):floor(3500*scale_factor)));
%        lLook=lLook+sum( xsit_result.test(subject).historyLt(trt,floor(4500*scale_factor):floor(6000*scale_factor)));
%        lLook=lLook+sum( xsit_result.test(subject).historyLt(trt,floor(6000*scale_factor):floor(7500*scale_factor)));   
%         
%        rLook= sum( xsit_result.test(subject).historyRt(trt,floor(500*scale_factor):floor(2000*scale_factor)));%250:2000
%        rLook=rLook+sum( xsit_result.test(subject).historyRt(trt,floor(2000*scale_factor):floor(3500*scale_factor)));
%        rLook=rLook+sum( xsit_result.test(subject).historyRt(trt,floor(4500*scale_factor):floor(6000*scale_factor)));
%        rLook=rLook+sum( xsit_result.test(subject).historyRt(trt,floor(6000*scale_factor):floor(7500*scale_factor)));
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

figure(1)% Plot total looking time during test trial
blockNames{1}=[blockNames{1}; parVal];
sts{1} = [sts{1}; mean(mean(targLookTime+dstrLookTime))/(scale_factor*slice_wordOnset4)];
errY{1}=[errY{1}; std(mean(targLookTime+dstrLookTime,2))/(scale_factor*slice_wordOnset4)];
b = barwitherr(errY{1}, sts{1});% Plot with errorbars
set(gca,'xticklabel',blockNames{1})
ylim([0 8000]);
title ('Avg Looking Time during a test trial');
ylabel('Looking time (ms)');
xlabel('atn_{sa} -> atn_c');

figure(2)% Plot Target vs Distractor looking time during test trial
blockNames{2}=[blockNames{2}; parVal];
sts{2} = [sts{2}; mean(mean(targLookTime))/(scale_factor*slice_wordOnset4)    mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4)];
errY{2} =[errY{2}; std(mean(targLookTime,2))/(scale_factor*slice_wordOnset4)    std(mean(dstrLookTime,2))/(scale_factor*slice_wordOnset4)];
barwitherr(errY{2}, sts{2});% Plot with errorbars
set(gca,'xticklabel',blockNames{2})
legend('Target','Distractor');
xlabel('atn_{sa} -> atn_c');
ylabel('Avg looking time per test trial');


figure(3)%Plot Targ vs Dist looking time for Strong learners
blockNames{3}=[blockNames{3}; parVal];
sts{3} = [sts{3}; mean(mean(targLookTime(goodLearners()==1,:)))/scale_factor    mean(mean(dstrLookTime(goodLearners()==1,:)))/scale_factor];
errY{3}=[errY{3}; std(mean(targLookTime(goodLearners()==1,:),2))/scale_factor     std(mean(dstrLookTime(goodLearners()==1,:),2))/scale_factor];
%b=bar(sts)
b = barwitherr(errY{3}, sts{3});% Plot with errorbars
set(gca,'xticklabel',blockNames{3});
legend('Target ','Distractor ');
xlabel('atn_{sa} -> atn_c');
title ('Strong Learners');
ylabel('Looking time (ms)');

figure(4)%Plot Targ vs Dist looking time for Weak learners
blockNames{4}=[blockNames{4}; parVal];
sts{4} = [sts{4}; mean(mean(targLookTime(goodLearners()==0,:)))/scale_factor   mean(mean(dstrLookTime(goodLearners()==0,:)))/scale_factor];
errY{4}=[errY{4}; std(mean(targLookTime(goodLearners()==0,:),2))/scale_factor     std(mean(dstrLookTime(goodLearners()==0,:),2))/scale_factor];
b = barwitherr(errY{4}, sts{4});% Plot with errorbars
set(gca,'xticklabel',blockNames{4})
title ('Weak Learners');
legend('Target','Distractor');
xlabel('atn_{sa} -> atn_c');
ylabel('Looking time (ms)');

% figure(5);%Plot Strong vs Weak looking time during over TEST trials
% plot(mean(targLookTime,1)/scale_factor,plotStyle{param});%
% hold on
% plot(mean(dstrLookTime,1)/scale_factor,plotStyle{param+compStyle});%
% legend(legendInfo);
% xlabel('test trial');
% ylabel('total looking time Target vs Distractor');

figure(6)%Plot Proportion of Strong/ Weak learners
blockNames{6}=[blockNames{6}; parVal];
sts{6} = [sts{6}; sum(goodLearners)/numSubjects   (numSubjects-sum(goodLearners))/numSubjects];
b = bar(sts{6});% Plot with errorbars
set(gca,'xticklabel',blockNames{6})
title ('Proportion of Learners');
legend('Strong','Weak');
xlabel('atn_{sa} -> atn_c');

figure(7)% Avg # of Words Learnt
blockNames{7}=[blockNames{7}; parVal];
sts{7} = [sts{7}; mean(sum(LearntWords(goodLearners()==1,:),2))     mean(sum(LearntWords(goodLearners()==0,:),2))];
errY{7}=[errY{7}; std(sum(LearntWords(goodLearners()==1,:),2))     std(sum(LearntWords(goodLearners()==0,:),2))];
b = barwitherr(errY{7}, sts{7});% Plot with errorbars
set(gca,'xticklabel',blockNames{7})
ylim([0 6]);
title ('Avg # of Words Learnt');
legend('Strong','Weak');
xlabel('atn_{sa} -> atn_c');
grid on

xx=['Avg Looking time per test trial is ',num2str(mean(mean(targLookTime+dstrLookTime))/(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Looking time to Target per test trial  is ',num2str(mean(mean(targLookTime))/(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Looking time to Distractor per test trial is ',num2str(mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Proportion of time looking correctly (Target/Total) is ',num2str(mean(sum (correct_proportion,2)))]; disp(xx);
xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);%
xx=['Avg # of Words Learnt by Strong Learners is ',num2str(mean(sum(LearntWords(goodLearners()==1,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt by Weak Learners is ',num2str(mean(sum(LearntWords(goodLearners()==0,:),2)))]; disp(xx);%  
%
%% Training trials analysis NEW
tLookAtLearntWords=zeros(numSubjects,3);%3 blocks
tLookAtNonLrntWords=zeros(numSubjects,3);%3 blocks
for subject=1:numSubjects
    savestate_historyL = xsit_result.train(subject).historyL;
    savestate_historyR = xsit_result.train(subject).historyR;    
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

strongLtime=mean(tLookAtLearntWords((goodLearners()==1),:))/(scale_factor);
SLerr=std(tLookAtLearntWords((goodLearners()==1),:))/(scale_factor);

strongNLtime=mean(tLookAtNonLrntWords((goodLearners()==1),:))/(scale_factor);
SNLerr=std(tLookAtNonLrntWords((goodLearners()==1),:))/(scale_factor);

weakLtime=mean(tLookAtLearntWords((goodLearners()==0),:))/(scale_factor);
WLerr=std(tLookAtLearntWords((goodLearners()==0),:))/(scale_factor);

weakNLtime=mean(tLookAtNonLrntWords((goodLearners()==0),:))/(scale_factor);
WNLerr=std(tLookAtNonLrntWords((goodLearners()==0),:))/(scale_factor);


% % % % % % % % % % figure(8)% Plot looking time at leant vs non-learnt words for Strong learners
% % % % % % % % % % blockNames{8}=[blockNames{8}; parVal];%'1to10','11to20','21to30'
% % % % % % % % % % sts{8} = [sts{8}; strongLtime(1) strongNLtime(1); strongLtime(2) strongNLtime(2);strongLtime(3) strongNLtime(3)];
% % % % % % % % % % errY{8} =[errY{8}; SLerr(1) SNLerr(1); SLerr(2) SNLerr(2); SLerr(3) SNLerr(3)];
% % % % % % % % % % barwitherr(errY{8}, sts{8});% Plot with errorbars
% % % % % % % % % % set(gca,'xticklabel',blockNames{8});
% % % % % % % % % % title ('Strong Learners');
% % % % % % % % % % ylabel('Looking time (ms)');
% % % % % % % % % % legend('Avg Learned Words','Avg NonLearned Words');
% % % % % % % % % % xlabel('atn_{sa} -> atn_c');
% % % % % % % % % % 
% % % % % % % % % % figure(9)% Plot looking time at leant vs non-learnt words for Weak learners
% % % % % % % % % % blockNames{9}=[blockNames{9}; parVal];%'1to10','11to20','21to30'
% % % % % % % % % % sts{9} = [sts{9}; weakLtime(1) weakNLtime(1); weakLtime(2) weakNLtime(2);weakLtime(3) weakNLtime(3)];
% % % % % % % % % % errY{9} =[errY{9}; WLerr(1) WNLerr(1); WLerr(2) WNLerr(2); WLerr(3) WNLerr(3)];
% % % % % % % % % % barwitherr(errY{9}, sts{9});% Plot with errorbars
% % % % % % % % % % set(gca,'xticklabel',blockNames{9});
% % % % % % % % % % title ('Weak Learners');
% % % % % % % % % % ylabel('Looking time (ms)');
% % % % % % % % % % legend('Avg Learnt Words','Avg Non-Learnt Words');
% % % % % % % % % % xlabel('atn_{sa} -> atn_c');


targLookTimeTraining=zeros(numSubjects,size(xsit_result.train(1).historyL,1));%number of training trials
dstrLookTimeTraining=zeros(numSubjects,size(xsit_result.train(1).historyL,1));
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
    savestate_historyL = xsit_result.train(subject).historyL;
    savestate_historyR = xsit_result.train(subject).historyR;
    %%%%% looking during training
    for tr=1:size(savestate_historyL,1)
        s1=char(xsit_result.train(subject).training_pair(tr,2*nFeatures+3));
        word1_ON=500;
        word1_OFF=1500;%1500
        word2_ON=2000;
        word2_OFF=3000;%3000
        if (strcmp(s1,'P'))% on left
           targLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyR(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking to target object right side
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyL(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking wrong way
        elseif (strcmp(s1,'X'))
           targLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyL(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking to target object right side
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyR(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking wrong way
            
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
            if (round(ldata(tr,time-1)) == 1)
                nlooks(side,tr) = nlooks(side,tr)+1;
            end
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
plot(mean(TotalLookTimeS)/scale_factor,plotStyle{param});%
hold on
plot(mean(TotalLookTimeW)/scale_factor,plotStyle{param+compStyle});%
legend(legendInfo);
xlabel('training trial');
ylabel('total looking time Strong vs Weak learners');
%ylabel('total looking time Strong learners');
title({'- Strong; -. Weak';'atn_{sa} -> atn_c'})
%hold off 

figure (12);% Plot number of fixations/looks over training trials
plot(mean(totnlooksS),plotStyle{param});% number of fixations/looks over training trials
hold on
plot(mean(totnlooksW),plotStyle{param+compStyle});% number of fixations/looks over training trials
xlabel('training trial');
ylabel('number of fixations/looks Strong vs Weak learners');
%ylabel('number of fixations/looks Strong learners');
legend(legendInfo);
title({'- Strong; -. Weak';'atn_{sa} -> atn_c'})
%hold off

figure (13);%Plot mean look duration of each fixation
%errorbar(mean(meanlookdurS)/scale_factor,std(meanlookdurS)/scale_factor, plotStyle{param});% mean look duration % multiped by timing scale factor
plot(mean(meanlookdurS)/scale_factor,plotStyle{param});%
hold on
%errorbar(mean(meanlookdurW)/scale_factor,std(meanlookdurW)/scale_factor, plotStyle{param+compStyle});% mean look duration % multiped by timing scale factor
plot(mean(meanlookdurW)/scale_factor,plotStyle{param+compStyle});%
xlabel('training trial');
ylabel('mean look duration Strong vs Weak learners');
%ylabel('mean look duration Strong learners');
legend(legendInfo);
title({'- Strong; -. Weak';'atn_{sa} -> atn_c'})
%hold off

figure (14);%Plot duration of longest look per trial
%errorbar(mean(totlonglookdurS)/scale_factor,std(totlonglookdurS)/scale_factor,plotStyle{param}); %duration of longest look per trial
plot(mean(totlonglookdurS)/scale_factor,plotStyle{param}); %duration of longest look per trial
hold on
%errorbar(mean(totlonglookdurW)/scale_factor,std(totlonglookdurW)/scale_factor, plotStyle{param+compStyle}); %duration of longest look per trial
plot(mean(totlonglookdurW)/scale_factor,plotStyle{param+compStyle}); %duration
xlabel('training trial');
ylabel('duration of longest look Strong vs Weak learners');
legend(legendInfo);
title({'- Strong; -. Weak';'atn_{sa} -> atn_c'})
%hold off

figure (15);%Plot Mean propertion of looking Varying vs repeated C
plot(mean(tLookVaryingS./(tLookVaryingS+tLookRepeatedS)),plotStyle{param}); %duration of longest look per trial
hold on
plot(mean(tLookRepeatedS./(tLookVaryingS+tLookRepeatedS)), plotStyle{param+compStyle}); %duration of longest look per trial
xlabel('Block Number');
ylabel('Proportion Looking Time Strong Learners');
legend(legendInfo);
title({'- Varying; -. Repeated';'atn_{sa} -> atn_c'});
xlim([0.5 6.5]);

figure (16);%Plot Mean propertion of looking Varying vs repeated WEAK
plot(mean(tLookVaryingW./(tLookVaryingW+tLookRepeatedW)),plotStyle{param}); %duration of longest look per trial
hold on
plot(mean(tLookRepeatedW./(tLookVaryingW+tLookRepeatedW)), plotStyle{param+compStyle}); %duration of longest look per trial
ylabel('Proportion Looking Time Weak Learners');
xlabel('Block Number');
legend(legendInfo);
title({'- Varying; -. Repeated';'atn_{sa} -> atn_c'});
xlim([0.5 6.5]);

summaF=mean(mean([TotalLookTimeS; TotalLookTimeW])) /scale_factor;
xx=['Avg looking time per training trial is ',num2str(summaF)]; disp(xx);
var1= mean(mean(TotalLookTimeS))/scale_factor;
xx=['Avg looking time per training trial per Strong learner is ',num2str(var1)]; disp(xx);
var2= mean(mean(TotalLookTimeW))/scale_factor;
xx=['Avg looking time per training trial per Weak learner is ',num2str(var2)]; disp(xx);


%%%%%%%% Association hwf trace analysis
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
        for jj=1:size(AsocMat,1)
           temp(jj)=AsocMat(jj,cell2mat(xsit_result.train(subject).Words(kk)));         
        end
        maxAsocnVal(subject,kk) = max(temp);
        NxtmaxAsocnVal(subject,kk) = max(temp(temp<max(temp)));
        ratioMax(subject,kk)= maxAsocnVal(subject,kk)./NxtmaxAsocnVal(subject,kk);
        prodtMR(subject,kk)=ratioMax(subject,kk).*maxAsocnVal(subject,kk);
        
        
        [maxIn(kk) indIn(kk)] = max(inputMapping(:,kk));
        [maxAs(kk) indAs(kk)] = max(AsocMat(:,kk));       
        if (abs(indIn(kk)-indAs(kk)) <= 3)%if association is correct i..e same as input?
           corrAsocn(subject, kk)=1; % wrongAssocn = 6-corrAsocn
        end
    end 
end

% SLer=[];WLer=[];SNon=[];WNon=[];SLer2=[];WLer2=[];SNon2=[];WNon2=[];
% for subject=1:numSubjects 
%     if(goodLearners(subject)==1)
%         SLer=[SLer maxAsocnVal(subject,LearntWords(subject,:)==1)];
%         SNon=[SNon maxAsocnVal(subject,LearntWords(subject,:)==0)];
%         
%         SLer2=[SLer2 ratioMax(subject,LearntWords(subject,:)==1)];
%         SNon2=[SNon2 ratioMax(subject,LearntWords(subject,:)==0)];
%      elseif (goodLearners(subject)==0)
%          WLer=[WLer maxAsocnVal(subject,LearntWords(subject,:)==1)];
%          WNon=[WNon maxAsocnVal(subject,LearntWords(subject,:)==0)];
%          
%          WLer2=[WLer2 ratioMax(subject,LearntWords(subject,:)==1)];
%          WNon2=[WNon2 ratioMax(subject,LearntWords(subject,:)==0)];
%     end
% end
% if ((size(SLer,2)+size(WLer,2)+size(SNon,2)+size(WNon,2))./numSubjects ~= nObjects), disp('ERROR ERROR ERROR ERROR'), end
% 
% figure(17)%Plot learnt vs nonleart association strength for Strong learners
% blockNames{17}=[blockNames{17}; parVal];
% sts{17} = [sts{17}; mean(SLer)    mean(SNon)];
% errY{17}=[errY{17}; std(SLer)     std(SNon)];
% b = barwitherr(errY{17}, sts{17});% Plot with errorbars
% set(gca,'xticklabel',blockNames{17});
% legend('Learnt','NonLearnt ');
% xlabel('atn_sa -> atn_c');
% title ('Avg association strength: Strong Learners');
% 
% 
% figure(18)%Plot learnt vs nonleart association strength for Weak learners
% blockNames{18}=[blockNames{18}; parVal];
% sts{18} = [sts{18}; mean(WLer)    mean(WNon)];
% errY{18}=[errY{18}; std(WLer)     std(WNon)];
% b = barwitherr(errY{18}, sts{18});% Plot with errorbars
% set(gca,'xticklabel',blockNames{18});
% legend('Learnt','NonLearnt ');
% xlabel('atn_sa -> atn_c');
% title ('Avg association strength: Weak Learners');
% 
% 
% figure(19)%Plot learnt vs nonleart Avg Ratio of 2 Maximums for Strong learners
% blockNames{19}=[blockNames{19}; parVal];
% sts{19} = [sts{19}; mean(SLer2)    mean(SNon2)];
% errY{19}=[errY{19}; std(SLer2)     std(SNon2)];
% b = barwitherr(errY{19}, sts{19});% Plot with errorbars
% set(gca,'xticklabel',blockNames{19});
% legend('Learnt','NonLearnt ');
% xlabel('atn_sa -> atn_c');
% title ('Avg Ratio of 2 Maximums: Strong Learners');
% ylim([0 5]);
% 
% figure(20)%Plot learnt vs nonleart Avg Ratio of 2 Maximums for Strong learners
% blockNames{20}=[blockNames{20}; parVal];
% sts{20} = [sts{20}; mean(WLer2)    mean(WNon2)];
% errY{20}=[errY{20}; std(WLer2)     std(WNon2)];
% b = barwitherr(errY{20}, sts{20});% Plot with errorbars
% set(gca,'xticklabel',blockNames{20});
% legend('Learnt','NonLearnt ');
% xlabel('atn_sa -> atn_c');
% title ('Avg Ratio of 2 Maximums: Weak Learners');
% ylim([0 5]);
% 
% figure(21)%Plot learnt vs nonleart association strength for Weak learners
% blockNames{21}=[blockNames{21}; parVal];
% sts{21} = [sts{21}; mean(sum(corrAsocn,2))/nObjects];
% errY{21}=[errY{21}; std(sum(corrAsocn,2))/nObjects];
% b = barwitherr(errY{21}, sts{21});% Plot with errorbars
% set(gca,'xticklabel',blockNames{21});
% xlabel('atn_sa -> atn_c');
% title ('Avg proportion of correctly associated words in Memory');
% 
% 
%  varb=mean(SLer); xx=['Avg association strength for Learnt words in Strong learners ',num2str(varb)]; disp(xx);
%  varb=mean(WLer); xx=['Avg association strength for Learnt words in Weak learners ',num2str(varb)]; disp(xx);
%  varb=mean(SNon); xx=['Avg association strength for NonLearnt words in Strong learners ',num2str(varb)]; disp(xx);
%  varb=mean(WNon); xx=['Avg association strength for NonLearnt words in Weak learners ',num2str(varb)]; disp(xx);
%   varb=mean(SLer2); xx=['Avg Ratio of 2Maximums for Learnt words in Strong learners ',num2str(varb)]; disp(xx);
%  varb=mean(WLer2); xx=['Avg Ratio of 2Maximums for Learnt words in Weak learners ',num2str(varb)]; disp(xx);
%  varb=mean(SNon2); xx=['Avg Ratio of 2Maximums for NonLearnt words in Strong learners ',num2str(varb)]; disp(xx);
%  varb=mean(WNon2); xx=['Avg Ratio of 2Maximums for NonLearnt words in Weak learners ',num2str(varb)]; disp(xx);
% % %wordsAssoc=sum(corrAsocn,2);
%  varb=mean(sum(corrAsocn,2))/nObjects; xx=['Avg proportion of correctly associated words in Memory ',num2str(varb)]; disp(xx);
% varb=mean(sum(corrAsocn((goodLearners()==1),:),2)); xx=['Avg # of correctly associated words for Strong in Memory ',num2str(varb)]; disp(xx);
% varb=mean(sum(corrAsocn((goodLearners()==0),:),2)); xx=['Avg # of correctly associated words for Weak in Memory ',num2str(varb)]; disp(xx);
% %Learnt non learnt by strong weak 
% varb=mean(corrWordsStrong); xx=['Avg # of correctly associated words in Memory for Strong counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% %std(corrWordsStrong);
% varb=mean(corrWordsWeak); xx=['Avg # of correctly associated words in Memory for Weak counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% %std(corrWordsWeak);
% varb=sum((corrAsocn(:,1:6).*LearntWords(:,1:6)));
% xx=['# of subjects with correctly associated word in Memory counted as Learnt thru looking for ',num2str(numSubjects),' subjects is ' num2str(varb)]; disp(xx);

end
% for subject=1:numSubjects
%     inputMapping=zeros(100,10);
%     for kk=1:nObjects
%         inputMapping(cell2mat(xsit_result.train(subject).Feature1(kk)),cell2mat(xsit_result.train(subject).Words(kk)))=1;
%     end
%     figure(21)
%     
%     lsp=subplot(1,2,1);
%     surface(inputMapping)
%     hold on
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

    
