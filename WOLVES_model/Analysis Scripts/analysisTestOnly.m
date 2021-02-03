%% Data Analysis 
clear all; close all;
nObjects=6;nFeatures=2; scale_factor=3000/8000;  %check with auto file
plotStyle = {'k-','b-','g-','c-','r-','m-','y-','k-.','b-.','g-.','c-.','r-.','m-.','y-.','b:','w.'};
compStyle=7;
% xsit_result.Names={[1] [2] [3] [4] [5] [6]};
%xsit_result.train =load('5atn_13wf_xsitNF_surprise3_1_train.mat');
%xsit_result.test =load('5atn_13wf_xsitNF_surprise3_1_test.mat');

simName = 'Bhwfwf2_wfatnf1.3_xsitNF_sub_';
OutName = [simName,'results.mat']
xsit_result = load (OutName);
numSubjects=size(xsit_result.test,1);
%numSubjects=28;
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
          if (xsit_result.Names{kk}== xsit_result.test(subject).test_pair(trt,2*nFeatures+1))%word index
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
sts = [mean(mean(targLookTime+dstrLookTime))/(scale_factor*slice_wordOnset4)];
errY=[std(mean(targLookTime+dstrLookTime,2))/(scale_factor*slice_wordOnset4)];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel','Test')
title ('Avg Looking Time');
ylabel('Looking time (ms)');


figure(2)% Plot Target vs Distractor looking time during test trial
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime))/(scale_factor*slice_wordOnset4)    mean(mean(dstrLookTime))/(scale_factor*slice_wordOnset4)];
errY =[std(mean(targLookTime,2))/(scale_factor*slice_wordOnset4)    std(mean(dstrLookTime,2))/(scale_factor*slice_wordOnset4)];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames)
%legend('Target','Distractor');
ylabel('Avg looking time per test trial');


figure(3)%Plot Targ vs Dist looking time for Strong learners
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime(goodLearners()==1,:)))/scale_factor    mean(mean(dstrLookTime(goodLearners()==1,:)))/scale_factor];
errY=[ std(mean(targLookTime(goodLearners()==1,:),2))/scale_factor     std(mean(dstrLookTime(goodLearners()==1,:),2))/scale_factor];
%b=bar(sts)
b = barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
%legend('Target ','Distractor ');
title ('Strong Learners');
ylabel('Looking time (ms)');

figure(4)%Plot Targ vs Dist looking time for Weak learners
blockNames={'Target'; 'Distractor'};
sts = [ mean(mean(targLookTime(goodLearners()==0,:)))/scale_factor   mean(mean(dstrLookTime(goodLearners()==0,:)))/scale_factor];
errY=[std(mean(targLookTime(goodLearners()==0,:),2))/scale_factor     std(mean(dstrLookTime(goodLearners()==0,:),2))/scale_factor];
b = barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames)
title ('Weak Learners');
%legend('Target','Distractor');
%xlabel('Fixed order # ');
ylabel('Looking time (ms)');

figure(5);%Plot Strong vs Weak looking time during over TEST trials
errorbar(mean(targLookTime,1)/scale_factor,std(targLookTime,1)/scale_factor,plotStyle{1});%
hold on
errorbar(mean(dstrLookTime,1)/scale_factor,std(dstrLookTime,1)/scale_factor,plotStyle{2+compStyle});%
legend('Strong','Weak');
xlabel('test trial');
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
tLookAtLearntWords=zeros(numSubjects,3);%3 blocks
tLookAtNonLrntWords=zeros(numSubjects,3);%3 blocks
for subject=1:numSubjects
    savestate_historyL = xsit_result.test(subject).historyLt;
    savestate_historyR = xsit_result.test(subject).historyRt;
    nHasLWord=zeros(3,1); tLookAtLearnt=zeros(3,1);
    nHasUWord=zeros(3,1); tLookAtNonLrnt=zeros(3,1);
    for kk=1:nObjects
        if (LearntWords(subject,kk)==1)%if the subject has learnt the kk word
           for block=1:3
                for tr=1:4%each trial may be 10 only
                    tt=4*(block-1)+tr;
                    if(xsit_result.Names{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial name1
                        nHasLWord(block)=nHasLWord(block)+1;
                        s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
                        if ( strcmp(s1,'L'))% on left
                            tLookAtLearnt(block)=tLookAtLearnt(block)+sum(savestate_historyL(tt,:));%add time historyL
                        elseif (strcmp(s1,'R'))
                            tLookAtLearnt(block)=tLookAtLearnt(block)+sum(savestate_historyR(tt,:));% add time historyR
                        end
                    end
                end
           end
        elseif (LearntWords(subject,kk)==0)%if word is NOT learnt
           for block=1:3
                for tr=1:4%each trial may be 10 only
                    tt=4*(block-1)+tr;
                    if(xsit_result.Names{kk}==xsit_result.test(subject).test_pair(tt,2*nFeatures+1))%object exists in trial
                        nHasUWord(block)=nHasUWord(block)+1;
                        s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
                        if ( strcmp(s1,'L'))% on left
                            tLookAtNonLrnt(block)=tLookAtNonLrnt(block)+sum(savestate_historyL(tt,:));%add time historyL
                        elseif (strcmp(s1,'R'))
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


figure(8)% Plot looking time at leant vs non-learnt words for Strong learners
blockNames={'1to4','5to8','9to12'};%
sts = [ strongLtime(1) strongNLtime(1); strongLtime(2) strongNLtime(2);strongLtime(3) strongNLtime(3)];
errY =[ SLerr(1) SNLerr(1); SLerr(2) SNLerr(2); SLerr(3) SNLerr(3)];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Strong Learners');
ylabel('Looking time (ms)');
legend('Avg Learned Words','Avg NonLearned Words');

figure(9)% Plot looking time at leant vs non-learnt words for Weak learners
blockNames={'1to4','5to8','9to12'};%
sts = [ weakLtime(1) weakNLtime(1); weakLtime(2) weakNLtime(2);weakLtime(3) weakNLtime(3)];
errY =[ WLerr(1) WNLerr(1); WLerr(2) WNLerr(2); WLerr(3) WNLerr(3)];
barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
title ('Weak Learners');
ylabel('Looking time (ms)');
legend('Avg Learnt Words','Avg Non-Learnt Words');
%
%% old
totnlooksS=zeros(sum(goodLearners),size(savestate_historyL,1));
meanlookdurS =zeros(sum(goodLearners),size(savestate_historyL,1));
totlonglookdurS = zeros(sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeS=zeros(sum(goodLearners),size(savestate_historyL,1));

totnlooksW=zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
meanlookdurW =zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
totlonglookdurW = zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeW=zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
tLookRepeatedS= zeros(sum(goodLearners),12);
tLookVaryingS=zeros(sum(goodLearners),12);
tLookRepeatedW=zeros(numSubjects-sum(goodLearners),12);
tLookVaryingW=zeros(numSubjects-sum(goodLearners),12);
indS=1;indW=1; summa=0;
for subject=1:numSubjects
    savestate_historyL = xsit_result.test(subject).historyLt;
    savestate_historyR = xsit_result.test(subject).historyRt;
    nlooks=zeros(2,size(savestate_historyL,1)); %L/R %%SAVE US!! 
    longlookdur=zeros(2,size(savestate_historyL,1));

    for side=1:2
        if side == 1
            ldata = savestate_historyL;
        else
            ldata = savestate_historyR;
        end
        for tr=1:size(ldata,1)
            look=0;
            templonglookdur=0;
            for time=1:size(ldata,2)
                if (round(ldata(tr,time)) == 1)
                    if look == 0
                        nlooks(side,tr) = nlooks(side,tr)+1;
                        if templonglookdur > longlookdur(side,tr)
                            longlookdur(side,tr) = templonglookdur;
                        end
                        templonglookdur=0;
                    end
                    look = 1;
                    templonglookdur = templonglookdur+1;
                else
                    look = 0;
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
     if(goodLearners(subject)==1)
         for blockz=1:12
              TinA=blockz;
              
            tLookRepeatedS(indS,blockz)=sum(sum(savestate_historyL(TinA,:)));
            tLookVaryingS(indS,blockz)=sum(sum(savestate_historyR(TinA,:)));
          end
        totnlooksS(indS,:)=sum(nlooks,1);
        meanlookdurS(indS,:)=(sum(savestate_historyL')+sum(savestate_historyR'))./totnlooksS(indS,:);
        totlonglookdurS(indS,:)=max(longlookdur,[],1);
        TotalLookTimeS(indS,:)=sum(savestate_historyL')+sum(savestate_historyR');
        indS=indS+1;
     elseif (goodLearners(subject)==0)
         for blockz=1:12
              TinA=blockz;
             
            tLookRepeatedW(indW,blockz)=sum(sum(savestate_historyL(TinA,:)));
            tLookVaryingW(indW,blockz)=sum(sum(savestate_historyR(TinA,:)));
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
xlabel('per test trial');
ylabel('total looking time Strong vs Weak learners');
%ylabel('total looking time Strong learners');
%hold off 

figure (12);% Plot number of fixations/looks over training trials
errorbar(mean(totnlooksS),std(totnlooksS),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(mean(totnlooksW),std(totnlooksW),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per test trial');
ylabel('number of fixations/looks Strong vs Weak learners');
%ylabel('number of fixations/looks Strong learners');
legend('Strong','Weak');
%hold off

figure (13);%Plot mean look duration of each fixation
errorbar(mean(meanlookdurS)/scale_factor,std(meanlookdurS)/scale_factor, plotStyle{1});% mean look duration % multiped by timing scale factor
hold on
errorbar(mean(meanlookdurW)/scale_factor,std(meanlookdurW)/scale_factor, plotStyle{2+compStyle});% mean look duration % multiped by timing scale factor
xlabel('per test trial');
ylabel('mean look duration Strong vs Weak learners');
%ylabel('mean look duration Strong learners');
legend('Strong','Weak');
%hold off

figure (14);%Plot duration of longest look per trial
errorbar(mean(totlonglookdurS)/scale_factor,std(totlonglookdurS)/scale_factor,plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(totlonglookdurW)/scale_factor,std(totlonglookdurW)/scale_factor, plotStyle{2+compStyle}); %duration of longest look per trial
xlabel('per test trial');
ylabel('duration of longest look Strong vs Weak learners');
%ylabel('duration of longest look Strong learners');
legend('Strong','Weak');
%hold off

figure (15);%Plot Mean propertion of looking Varying vs repeated C
errorbar(mean(tLookVaryingS./(tLookVaryingS+tLookRepeatedS)),std(tLookVaryingS./(tLookVaryingS+tLookRepeatedS)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeatedS./(tLookVaryingS+tLookRepeatedS)),std(tLookRepeatedS./(tLookVaryingS+tLookRepeatedS)), plotStyle{2+compStyle}); %duration of longest look per trial
xlabel('per test trial');
xlim([0.1 12.9]);
ylabel('Proportion Looking Time Strong Learners');
%ylabel('duration of longest look Strong learners');
legend('Varying','Repeated');


figure (16);%Plot Mean propertion of looking Varying vs repeated WEAK
errorbar(mean(tLookVaryingW./(tLookVaryingW+tLookRepeatedW)),std(tLookVaryingW./(tLookVaryingW+tLookRepeatedW)),plotStyle{1}); %duration of longest look per trial
hold on
errorbar(mean(tLookRepeatedW./(tLookVaryingW+tLookRepeatedW)),std(tLookRepeatedW./(tLookVaryingW+tLookRepeatedW)), plotStyle{2+compStyle}); %duration of longest look per trial
xlabel('per test trial');
xlim([0.1 12.9]);
ylabel('Proportion Looking Time Weak Learners');
%ylabel('duration of longest look Strong learners');
legend('Varying','Repeated');


summaF=mean(mean([TotalLookTimeS; TotalLookTimeW])) /scale_factor;
xx=['Avg looking time per test trial is ',num2str(summaF)]; disp(xx);
var1= mean(mean(TotalLookTimeS))/scale_factor;
xx=['Avg looking time per test trial per Strong learner is ',num2str(var1)]; disp(xx);
var2= mean(mean(TotalLookTimeW))/scale_factor;
xx=['Avg looking time per test trial per Weak learner is ',num2str(var2)]; disp(xx);


%%%%%%%% Association hwf trace analysis
corrAsocn=zeros(numSubjects,nObjects);
cS=1;cW=1;
for subject=1:numSubjects  
    inputMapping=zeros(100,10);
    for kk=1:nObjects
        inputMapping(cell2mat(xsit_result.test(subject).Feature1(kk)),cell2mat(xsit_result.test(subject).Words(kk)))=1;
    end 
    
    AsocMat=squeeze(xsit_result.test(subject).hwft(1,:,:));
    
    for kk=1:nObjects        
        temp=[];
        for jj=1:nObjects%(features=size(Feature1))
           temp(jj)=AsocMat(cell2mat(xsit_result.test(subject).Feature1(jj)),kk);         
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
% varb=mean(sum(corrAsocn((goodLearners()==1),:),2)); xx=['Avg # of correctly associated words for Strong in Memory ',num2str(varb)]; disp(xx);
% varb=mean(sum(corrAsocn((goodLearners()==0),:),2)); xx=['Avg # of correctly associated words for Weak in Memory ',num2str(varb)]; disp(xx);
% %Learnt non learnt by strong weak 
% varb=mean(corrWordsStrong); xx=['Avg # of correctly associated words in Memory for Strong counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% %std(corrWordsStrong);
% varb=mean(corrWordsWeak); xx=['Avg # of correctly associated words in Memory for Weak counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% %std(corrWordsWeak);
% varb=sum((corrAsocn(:,1:6).*LearntWords(:,1:6)));
% xx=['# of subjects with correctly associated word in Memory counted as Learnt thru looking for ',num2str(numSubjects),' subjects is ' num2str(varb)]; disp(xx);

for subject=1:numSubjects
    inputMapping=zeros(100,10);
    for kk=1:nObjects
        inputMapping(cell2mat(xsit_result.test(subject).Feature1(kk)),cell2mat(xsit_result.test(subject).Words(kk)))=1;
    end
    figure(21)
    
    %lsp=subplot(1,2,1);
    surface(inputMapping)
    hold on
    [mA, iA] = max(squeeze(xsit_result.test(subject).hwft(1,:,:)));
    surface(squeeze(xsit_result.test(subject).hwft(1,:,:)));
    title([num2str(mA)]);
    shading flat; 
    
    %rsp= subplot(1,2,2);
    %hold on
    
    %title('Input');
    
%     surface(squeeze(xsit_result.test(subject).hwft(1,:,:)));
%     title('Test');
    %shading flat;
    
    if (goodLearners(subject)==1), suptitle(['Strong # ' num2str(subject)]),
    elseif (goodLearners(subject)==0), suptitle(['Weak # ' num2str(subject)]), end
    pause(7);
    clf;
end