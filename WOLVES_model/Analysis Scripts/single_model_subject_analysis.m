%% Data Analysis 
%%% load('abc_results.mat')
%sim.saveSettings('bbc.json')
clear all; close all;
nObjects=6;nFeatures=2; scale_factor=1/8;  %check with auto file
plotStyle = {'k-','b-','g-','c-','r-','m-','y-','k-.','b-.','g-.','c-.','r-.','m-.','y-.','b:','w.'};
compStyle=7;
blockNames{10}=[];sts{10}=[];errY{10}=[];
% xsit_result.Names={[1] [2] [3] [4] [5] [6]};
%xsit_result.train =load('5atn_13wf_xsitNF_surprise3_1_train.mat');
%xsit_result.test =load('5atn_13wf_xsitNF_surprise3_1_test.mat');

simName = 'fixed_orders_3_Hab_Smith_Yu_2008_'%'bbc_H50F30Asc50_Tzero_Smith_Yu_2013_'; %test4_Smith_Yu_2008_
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
for subject=1:numSubjects %%SET SUBJECT HERE
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
%% %PLOTS
%%% PLOTS

tLookAtLearntWords=zeros(numSubjects,3);%3 blocks
tLookAtNonLrntWords=zeros(numSubjects,3);%3 blocks
for subject=1:numSubjects %SET SUBJECT HERE
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
for subject=1:numSubjects %SET SUBJECT HERE
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
        %%%%%%%%%%%%%% for larissa
        TotalLookTimeF(subject,:)=TotalLookTimeS(indS,:);
        totnlooksF(subject,:) = totnlooksS(indS,:);
        meanlookdurF(subject,:)=meanlookdurS(indS,:);
        totlonglookdurF(subject,:)=totlonglookdurS(indS,:); 
        %%%%%%%%%%%%%%
        
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
        %%%%%%% for larissa
         TotalLookTimeF(subject,:)=TotalLookTimeW(indW,:);
        totnlooksF(subject,:) = totnlooksW(indW,:);
        meanlookdurF(subject,:)=meanlookdurW(indW,:);
        totlonglookdurF(subject,:)=totlonglookdurW(indW,:); 
        %%%%%
        indW=indW+1;
     else
         disp('ERROR goodLearners bad data');
     end
     

lookR = savestate_historyR';
lookL = savestate_historyL';
vecL = lookL(:);
vecR = lookR(:);
c(subject,:,:)= [vecL(3000:9000), vecR(3000:9000)];
end

a1=squeeze(c(1,:,:))
a2=squeeze(c(2,:,:))
a3=squeeze(c(3,:,:))
size(a1)


subplot(2,1,1);
area(a1)
title('Stimuli Order 1')

subplot(2,1,2);
area(a2)
title('Stimuli Order 2')

subplot(3,1,3);
area(a3)
title('Stimuli Order 3')
legend('Left Look', 'Right Look');

%%multi file code
figure(11);%Plot Strong vs Weak looking time during a training trial
plot((TotalLookTimeF(1,:))/scale_factor);%
hold on
plot((TotalLookTimeF(2,:))/scale_factor);%
hold on
plot((TotalLookTimeF(3,:))/scale_factor);%
legend('Order 1', 'Order 2','Order 3');
xlabel('per training trial');
ylabel('Total looking time of a learner');
%ylabel('total looking time Strong learners');
%hold off 


figure (19);% Plot number of fixations/looks over training trials
plot(totnlooksF(1,:));% number of fixations/looks over training trials
hold on
plot(totnlooksF(2,:));% number of fixations/looks over training trials
hold on
plot(totnlooksF(3,:));
xlabel('per training trial');
ylabel('number of fixations/looks of a learner');
%ylabel('number of fixations/looks Strong learners');
legend('Order 1', 'Order 2','Order 3');
%ylim([1.5 3.5])
%hold off


summaF=mean(mean(totnlooksS));
xx=['number of fixations/looks Strong learners',num2str(summaF)]; disp(xx);
summaF=mean(mean(totnlooksW));
xx=['number of fixations/looks Weak learners ',num2str(summaF)]; disp(xx);

figure (13);%Plot mean look duration of each fixation
plot((meanlookdurF(1,:))/scale_factor);% mean look duration % multiped by timing scale factor
hold on
plot((meanlookdurF(2,:))/scale_factor);% mean look duration % multiped by timing scale factor
hold on
plot((meanlookdurF(3,:))/scale_factor);
xlabel('per training trial');
ylabel('look duration of a learner');
%ylabel('mean look duration Strong learners');
legend('Order 1', 'Order 2','Order 3');
%ylim([1000 2000])
%hold off



summaF=mean(mean(meanlookdurS)/scale_factor)/1000;
xx=['mean look duration Strong learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanlookdurW)/scale_factor)/1000;
xx=['mean look duration Strong weak learners ',num2str(summaF)]; disp(xx);



figure (14);%Plot duration of longest look per trial
plot((totlonglookdurF(1,:))/scale_factor);% mean look duration % multiped by timing scale factor
hold on
plot((totlonglookdurF(2,:))/scale_factor);% mean look duration % multiped by timing scale factor
hold on
plot((totlonglookdurF(3,:))/scale_factor);
xlabel('per training trial');
ylabel('duration of the longest look of a learner');
legend('Order 1', 'Order 2','Order 3');

%hold off


