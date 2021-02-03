%% Data Analysis 
clear all; close all;
nObjects=12;nFeatures=2; scale_factor=1/8;  %check with auto file
plotStyle = {'k-','b-','g-','c-','r-','m-','y-','k-.','b-.','g-.','c-.','r-.','m-.','y-.','b:','w.'};
compStyle=7;
blockNames{10}=[];sts{10}=[];errY{10}=[];
% xsit_result.Names={[1] [2] [3] [4] [5] [6]};
%xsit_result.train =load('5atn_13wf_xsitNF_surprise3_1_train.mat');
%xsit_result.test =load('5atn_13wf_xsitNF_surprise3_1_test.mat');
%base_tauM_stim525_wmc9_hwmf15_Vlach_Johnson_2013_
simName = 'wf_wmFront_wf_testTask_'          % base2011_Vlach_Johnson_2013_    %wmc2_sim_Vlach_Johnson_2013_
OutName = [simName,'results.mat'];
xsit_result = load (OutName);
numSubjects=size(xsit_result.test,1);
%% Test trials
correct_proportion=zeros(numSubjects,2);
targLookTime=zeros(numSubjects,33);
dstrLookTime=zeros(numSubjects,33);
goodLearners=999*ones(numSubjects,1);
LearntWords= 999*ones(numSubjects,nObjects);
%word_On = floor([500 2000 4500 6000]*scale_factor);%[0 450 875 1300 ];    
%word_Off = floor([1500 3000 5500 7000]*scale_factor);%[375 650 1075 1500 ];
slice_wordOnset4=1;%6/8;
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
MassedTargLookTime=zeros(numSubjects,1);
MassedDistLookTime=zeros(numSubjects,1);
IleavedTargLookTime=zeros(numSubjects,1);
IleavedDistLookTime=zeros(numSubjects,1);
for subject=1:numSubjects
    lcorrect=0;rcorrect=0;targWord=zeros(nObjects,1);dstrWord=zeros(nObjects,1);
    for trt=1:33
          lLook= sum( xsit_result.train(subject).historyL(trt,:));%250:2000 
          rLook= sum( xsit_result.train(subject).historyR(trt,:));%250:2000
       
       for kk=1:nObjects
          if (xsit_result.train(subject).Words{kk}== xsit_result.train(subject).training_pair(trt,2*nFeatures+1))%word index
               s1= char(xsit_result.train(subject).training_pair(trt,2*nFeatures+3));
               
               if ( strcmp(s1,'L')) 
                       targWord(kk)=targWord(kk)+lLook;
                       dstrWord(kk)=dstrWord(kk)+rLook;
               elseif ( strcmp(s1,'R'))
                       targWord(kk)=targWord(kk)+rLook;
                       dstrWord(kk)=dstrWord(kk)+lLook;
               else
                       disp('ERROR reading wfTrain_pair_char');
               end
          end
       end
%%%%%
       s1= char(xsit_result.train(subject).training_pair(trt,2*nFeatures+3));
       if ( strcmp(s1,'L')) 
           targLookTime(subject,trt)=lLook;
           dstrLookTime(subject,trt)=rLook;
           lcorrect=lcorrect+lLook/(lLook+rLook);
       elseif ( strcmp(s1,'R'))
           targLookTime(subject,trt)=rLook;
           dstrLookTime(subject,trt)=lLook;
           rcorrect=rcorrect+rLook/(lLook+rLook);
       else
           disp('ERROR reading wfTrain_pair char');
       end
    
    end
    for kk=1:nObjects
        if (targWord(kk)>dstrWord(kk))
            LearntWords(subject,kk)=1;
        else
            LearntWords(subject,kk)=0;
        end
        
    end
    lcorrect = lcorrect/(nObjects);     rcorrect = rcorrect/(nObjects);
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
perform_CSL = sum(LearntWords,2);
freq_CSL=zeros(1,nObjects);
for subject=1:numSubjects
    for trt=1:nObjects
        if trt==perform_CSL(subject)
            freq_CSL(trt)= freq_CSL(trt)+1;
        end
    end
end

figure(1)% Plot total looking time during test trial
sts = [mean(mean(targLookTime+dstrLookTime))/(scale_factor*slice_wordOnset4)];
errY=[std(mean(targLookTime+dstrLookTime,2))/(scale_factor*slice_wordOnset4)];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel','Train')
title ('Avg Looking Time');
ylabel('Looking time (ms)');


figure(2)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime))/(scale_factor*slice_wordOnset4)    mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4)];
errY =[std(mean(targLookTime,2))/(scale_factor*slice_wordOnset4)    std(mean(dstrLookTime,2))/(scale_factor*slice_wordOnset4)];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames)
%legend('Target','Distractor');
ylabel('Avg looking time per train trial');

figure(5);%Plot  looking time during over TEST trials
errorbar(mean(targLookTime,1)/scale_factor,std(targLookTime,1)/scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTime,1)/scale_factor,std(dstrLookTime,1)/scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('train trial');
ylabel('total looking time Target vs Distractor');

xx=['Avg Looking time per test trial is ',num2str(mean(mean(targLookTime+dstrLookTime))/(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Looking time to Target per trial per subject is ',num2str(mean(mean(targLookTime))/(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Looking time to Distractor per trial per subject is ',num2str(mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4))]; disp(xx);
xx=['Proportion of time looking correctly (Target/Total) is ',num2str(mean(sum (correct_proportion,2)))]; disp(xx);
xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);% 
xx=['Avg # of Words Learnt by Strong Learners is ',num2str(mean(sum(LearntWords(goodLearners()==1,:),2)))]; disp(xx);%  
xx=['Avg # of Words Learnt by Weak Learners is ',num2str(mean(sum(LearntWords(goodLearners()==0,:),2)))]; disp(xx);%  
% stanard deviation std(sum(LearntWords(goodLearners()==1,:),2))
%% Training trials analysis NEW

savestate_historyL = xsit_result.train(subject).historyL;
savestate_historyR = xsit_result.train(subject).historyR; 
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
        word2_ON=2500;
        word2_OFF=3500;%3000
        if (strcmp(s1,'L'))% on left
           targLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyL(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking to target object correct side
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyR(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking wrong way
        elseif (strcmp(s1,'R'))
           targLookTimeTraining(subject,tr) = sum(savestate_historyR(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking to target object
            + sum(savestate_historyR(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking to target object 
            
            dstrLookTimeTraining(subject,tr) = sum(savestate_historyL(tr,floor(word1_ON*scale_factor):floor(word1_OFF*scale_factor))) ... % first audio presentation, Looking wrong way
           + sum(savestate_historyL(tr,floor(word2_ON*scale_factor):floor(word2_OFF*scale_factor)));% 2nd audio presentation, Looking wrong way
            
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
errorbar(mean(TotalLookTimeS)/scale_factor,std(TotalLookTimeS)/scale_factor,plotStyle{1});%
hold on
errorbar(mean(TotalLookTimeW)/scale_factor,std(TotalLookTimeW)/scale_factor,plotStyle{2+compStyle});%
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
%hold off

figure (13);%Plot mean look duration of each fixation
errorbar(mean(meanlookdurS)/scale_factor,std(meanlookdurS)/scale_factor, plotStyle{1});% mean look duration % multiped by timing scale factor
hold on
errorbar(mean(meanlookdurW)/scale_factor,std(meanlookdurW)/scale_factor, plotStyle{2+compStyle});% mean look duration % multiped by timing scale factor
xlabel('per training trial');
ylabel('mean look duration Strong vs Weak learners');
%ylabel('mean look duration Strong learners');
legend('Strong','Weak');
%hold off

figure (14);%Plot duration of longest look per trial
errorbar(mean(totlonglookdurS)/scale_factor,std(totlonglookdurS)/scale_factor,plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(totlonglookdurW)/scale_factor,std(totlonglookdurW)/scale_factor, plotStyle{2+compStyle}); %duration of longest look per trial
xlabel('per training trial');
ylabel('duration of longest look Strong vs Weak learners');
%ylabel('duration of longest look Strong learners');
legend('Strong','Weak');
%hold off

figure(15);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(mean(targLookTimeTraining(:,:))/scale_factor,std(targLookTimeTraining(:,:))/scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTimeTraining(:,:))/scale_factor,std(dstrLookTimeTraining(:,:))/scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Training Trial');
ylabel('Strong Learners: Looking Time when words are ON');

figure(16);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(mean(targLookTimeTraining((goodLearners()==0),:))/scale_factor,std(targLookTimeTraining((goodLearners()==0),:))/scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTimeTraining((goodLearners()==0),:))/scale_factor,std(dstrLookTimeTraining((goodLearners()==0),:))/scale_factor,plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('Training Trial');
ylabel('Weak Learners: Looking Time when words are ON');


summaF=mean(mean([TotalLookTimeS; TotalLookTimeW])) /scale_factor;
xx=['Avg looking time per training trial is ',num2str(summaF)]; disp(xx);

var1= mean(mean(TotalLookTimeS))/scale_factor;
xx=['Avg looking time per training trial per Strong learner is ',num2str(var1)]; disp(xx);
var2= mean(mean(TotalLookTimeW))/scale_factor;
xx=['Avg looking time per training trial per Weak learner is ',num2str(var2)]; disp(xx);


%%%%%%%% Association hwf trace analysis
corrAsocn=zeros(numSubjects,nObjects);
cS=1;cW=1;
% for subject=1:numSubjects  
%     inputMapping=zeros(100,10);
%     for kk=1:nObjects
%         inputMapping(cell2mat(xsit_result.train(subject).Feature1(kk)),cell2mat(xsit_result.train(subject).Words(kk)))=1;
%     end 
%     
%     AsocMat=squeeze(xsit_result.train(subject).hwf(1,:,:));
%     
%     for kk=1:nObjects        
%         temp=[];
%         for jj=1:nObjects%(features=size(Feature1))
%            temp(jj)=AsocMat(cell2mat(xsit_result.train(subject).Feature1(jj)),kk);         
%         end
%         maxAsocnVal(subject,kk) = max(temp);
%         NxtmaxAsocnVal(subject,kk) = max(temp(temp<max(temp)));
%         ratioMax(subject,kk)= maxAsocnVal(subject,kk)./NxtmaxAsocnVal(subject,kk);
%         prodtMR(subject,kk)=ratioMax(subject,kk).*maxAsocnVal(subject,kk);
%         
%         
%         [maxIn(kk) indIn(kk)] = max(inputMapping(:,kk));
%         [maxAs(kk) indAs(kk)] = max(AsocMat(:,kk));       
%         if (abs(indIn(kk)-indAs(kk)) <= 3)%if association is correct i..e same as input?
%            corrAsocn(subject, kk)=1; % wrongAssocn = 6-corrAsocn
%         end
%     end 
% end

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
% % varb=mean(sum(corrAsocn((goodLearners()==1),:),2)); xx=['Avg # of correctly associated words for Strong in Memory ',num2str(varb)]; disp(xx);
% varb=mean(sum(corrAsocn((goodLearners()==0),:),2)); xx=['Avg # of correctly associated words for Weak in Memory ',num2str(varb)]; disp(xx);
% %Learnt non learnt by strong weak 
% varb=mean(corrWordsStrong); xx=['Avg # of correctly associated words in Memory for Strong counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% %std(corrWordsStrong);
% varb=mean(corrWordsWeak); xx=['Avg # of correctly associated words in Memory for Weak counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% %std(corrWordsWeak);
% varb=sum((corrAsocn(:,1:6).*LearntWords(:,1:6)));
% xx=['# of subjects with correctly associated word in Memory counted as Learnt thru looking for ',num2str(numSubjects),' subjects is ' num2str(varb)]; disp(xx);

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

    
