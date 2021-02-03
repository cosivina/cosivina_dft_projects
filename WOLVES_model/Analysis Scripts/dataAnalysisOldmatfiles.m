%% Data Analysis 
clear all; close all;
%simName = 'B_xsitNew_16_sub_';
numSubjects=28;nObjects=6; scale_factor=3000/8000;  %check with auto file
plotStyle = {'k-','b-','g-','c-','r-','m-','y-','k-.','b-.','g-.','c-.','r-.','m-.','y-.','b:','w.'};
for param=1:1
    parVal=2+2*(param-1);
    legendInfo{2*param -1}= [plotStyle{param} num2str(parVal)];
    legendInfo{2*param}= [plotStyle{param+7} num2str(parVal)];
simName = ['1hwf_wf_' num2str(parVal) '_xsit_sub71_']; %''
OutName = [simName,'results.mat'];
xsit_result = load (OutName);
%% Test trials
correct_proportion=zeros(numSubjects,2);
targLookTime=zeros(numSubjects,1);
dstrLookTime=zeros(numSubjects,1);
goodLearners=999*ones(numSubjects,1);
slice_wordOnset4=6/8;
%word_On = floor([500 2000 4500 6000]*scale_factor);%[0 450 875 1300 ];    
%word_Off = floor([1500 3000 5500 7000]*scale_factor);%[375 650 1075 1500 ];
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
for subject=1:numSubjects
    lcorrect=0;rcorrect=0;
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
       
       s1= xsit_result.test(subject).object_side(trt);
       if ( strcmp(s1,'L')) 
           targLookTime(subject)=targLookTime(subject)+lLook;
           dstrLookTime(subject)=dstrLookTime(subject)+rLook;
           lcorrect=lcorrect+lLook/(lLook+rLook);
       elseif ( strcmp(s1,'R'))
           targLookTime(subject)=targLookTime(subject)+rLook;
           dstrLookTime(subject)=dstrLookTime(subject)+lLook;
           rcorrect=rcorrect+rLook/(lLook+rLook);
       else
           disp('ERROR reading object_side');
       end
    end
    lcorrect = lcorrect/(2*nObjects);     rcorrect = rcorrect/(2*nObjects);
    correct_proportion(subject,1)=lcorrect;
    correct_proportion(subject,2)=rcorrect;
    targLookTime(subject) =targLookTime(subject)./(2*nObjects);
    dstrLookTime(subject) =dstrLookTime(subject)./(2*nObjects);
    if (targLookTime(subject)>dstrLookTime(subject))
        targTimeS=targTimeS+targLookTime(subject);
        distTimeS=distTimeS+dstrLookTime(subject);
        goodLearners(subject)=1;
    else
        targTimeW=targTimeW+targLookTime(subject);
        distTimeW=distTimeW+dstrLookTime(subject);
        goodLearners(subject)=0;
    end
end
targTimeS=targTimeS/(scale_factor*(sum(goodLearners)));
distTimeS=distTimeS/(scale_factor*(sum(goodLearners)));
targTimeW=targTimeW/(scale_factor*(numSubjects-sum(goodLearners)));
distTimeW=distTimeW/(scale_factor*(numSubjects-sum(goodLearners)));

figure(4)%strong targ vd dist
blockNames={'Target','Distractor'};
sts = [targTimeS; distTimeS];
b=bar(sts)
set(gca,'xticklabel',blockNames)
title ('Strong');
ylabel('Looking time (ms)');

figure(5)
blockNames={'Target','Distractor'};
sts = [targTimeW; distTimeW];
b=bar(sts)
set(gca,'xticklabel',blockNames)
title ('Weak');
ylabel('Looking time (ms)');


% figure
% plot(targLookTime/(scale_factor*slice_wordOnset4),'b');
% hold on;
% plot(dstrLookTime/(scale_factor*slice_wordOnset4),'r');
% hold off;
%%%plot(goodLearners,'g');
%%%MULtifile code 
targ(param)=mean(targLookTime)/(scale_factor*slice_wordOnset4);
dista(param)=mean(dstrLookTime)/(scale_factor*slice_wordOnset4);
parma(param)= parVal;
% xlabel('subject');
% ylabel('Total time looking at Target(blue) vs Distractor(red)');

 xx=['Total looking time to target per test trial per subject is ',num2str(mean(targLookTime)/(scale_factor*slice_wordOnset4))]; disp(xx);
 xx=['Total looking time to distractor per trial per subject is ',num2str(mean(dstrLookTime)/(scale_factor*slice_wordOnset4))]; disp(xx);
 xx=['Mean proportion of time looking correctly is ',num2str(mean(sum (correct_proportion,2)))]; disp(xx);
xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);% 

%% Training trials analysis
savestate_historyL = xsit_result.train(subject).historyL;
totnlooksS=zeros(sum(goodLearners),size(savestate_historyL,1));
meanlookdurS =zeros(sum(goodLearners),size(savestate_historyL,1));
totlonglookdurS = zeros(sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeLeftS=zeros(sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeRightS=zeros(sum(goodLearners),size(savestate_historyL,1));

totnlooksW=zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
meanlookdurW =zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
totlonglookdurW = zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeLeftW=zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
TotalLookTimeRightW=zeros(numSubjects-sum(goodLearners),size(savestate_historyL,1));
indS=1;indW=1; summa=0;
for subject=1:numSubjects
    savestate_historyL = xsit_result.train(subject).historyL;
    savestate_historyR = xsit_result.train(subject).historyR;
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
        totnlooksS(indS,:)=sum(nlooks,1);
        meanlookdurS(indS,:)=(sum(savestate_historyL')+sum(savestate_historyR'))./totnlooksS(indS,:);
        totlonglookdurS(indS,:)=max(longlookdur,[],1);
        TotalLookTimeLeftS(indS,:)=sum(savestate_historyL');
        TotalLookTimeRightS(indS,:)=sum(savestate_historyR');
        indS=indS+1;
        summa=summa+sum(sum(savestate_historyL'));
        summa=summa+sum(sum(savestate_historyR'));
     elseif (goodLearners(subject)==0)
        totnlooksW(indW,:)=sum(nlooks,1);
        meanlookdurW(indW,:)=(sum(savestate_historyL')+sum(savestate_historyR'))./totnlooksW(indW,:);
        totlonglookdurW(indW,:)=max(longlookdur,[],1);
        TotalLookTimeLeftW(indW,:)=sum(savestate_historyL');
        TotalLookTimeRightW(indW,:)=sum(savestate_historyR');
        indW=indW+1;
        summa=summa+sum(sum(savestate_historyL'));
        summa=summa+sum(sum(savestate_historyR'));
     else
         disp('ERROR goodLearners bad data');
     end
end
%%multi file code
sumlook(param)=(summa /scale_factor) /(numSubjects*size(savestate_historyL,1));
sumlookStr(param)=mean((mean(TotalLookTimeLeftS+TotalLookTimeRightS,2)))/scale_factor;
sumlookWek(param)=mean(mean(TotalLookTimeLeftW+TotalLookTimeRightW,2))/scale_factor;

 summaF=(summa /scale_factor) /(numSubjects*size(savestate_historyL,1));
 xx=['total looking time per training trial per subject is ',num2str(summaF)]; disp(xx);
 var1= mean((mean(TotalLookTimeLeftS+TotalLookTimeRightS,2)))/scale_factor;
 xx=['total looking time per training trial per Strong learner is ',num2str(var1)]; disp(xx);
 var2= mean(mean(TotalLookTimeLeftW+TotalLookTimeRightW,2))/scale_factor;
 xx=['total looking time per training trial per Weak learner is ',num2str(var2)]; disp(xx);


figure(1);
plot(mean(TotalLookTimeLeftS+TotalLookTimeRightS)/scale_factor,plotStyle{param});%
hold on
plot(mean(TotalLookTimeLeftW+TotalLookTimeRightW)/scale_factor,plotStyle{param+7});%
legend(legendInfo)
xlabel('per training trial');
ylabel('total looking time Strong vs Weak learners');
%ylabel('total looking time Strong learners');
%hold off 

figure (2);
plot(mean(totnlooksS),plotStyle{param});% number of fixations/looks over training trials
hold on
plot(mean(totnlooksW),plotStyle{param+7});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('number of fixations/looks Strong vs Weak learners');
%ylabel('number of fixations/looks Strong learners');
legend(legendInfo)
%hold off

mean(mean(totnlooksW))
mean(mean(meanlookdurW))/scale_factor
mean(mean(totlonglookdurW))/scale_factor

figure (3);
plot(mean(meanlookdurS)/scale_factor,plotStyle{param});% mean look duration % multiped by timing scale factor
hold on
plot(mean(meanlookdurW)/scale_factor,plotStyle{param+7});% mean look duration % multiped by timing scale factor
xlabel('per training trial');
ylabel('mean look duration Strong vs Weak learners');
%ylabel('mean look duration Strong learners');
legend(legendInfo)
%hold off

figure (4);
plot(mean(totlonglookdurS),plotStyle{param}); %duration of longest look per trial
hold on
plot(mean(totlonglookdurW),plotStyle{param+7}); %duration of longest look per trial
xlabel('per training trial');
ylabel('duration of longest look Strong vs Weak learners');
%ylabel('duration of longest look Strong learners');
legend(legendInfo)
%hold off
end

figure
plot(parma,sumlook,'r')
hold on
plot(parma,sumlookStr,'g')
hold on
plot(parma,sumlookWek,'b')
legend('r-total, g-strong, b-weak');
xlabel('amplitude wf n -> atn\_c n');
ylabel('Total time looking  per training trial')

figure
plot(parma,targ+dista,'r')
hold on
plot(parma,targ,'g')
hold on
plot(parma,dista,'b')
legend('r-total, g-target, b-distractor');
xlabel('amplitude wf n -> atn\_c n');
ylabel('Total time looking per test trial')