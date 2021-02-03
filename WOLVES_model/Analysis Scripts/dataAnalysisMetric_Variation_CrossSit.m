%% Data Analysis File
clear all; 
close all;
nObjects=12;nFeatures=2;nTrainTrials=30; nTestTrials=12; scale_factor=8; MIN_LOOK_DURATION=200/scale_factor; %check with auto file 160
TRAIN_DUR=4000; TEST_DUR=1000;
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;%blockNames{2}=[];sts{2}=[];errY{2}=[];
%TASK= 'NEAR';% 'FAR';% 

%% test vars
numSubjects=294;
correct_proportion=zeros(2,numSubjects,12);
targLookTime=zeros(2,numSubjects,nTestTrials);
dstrLookTime=zeros(2,numSubjects,nTestTrials);
LearntWords= NaN(2,numSubjects,nObjects);

%% training vars
corrLookTimeTraining=zeros(2,numSubjects,nTrainTrials);%number of training trials
incorrLookTimeTraining=zeros(2,numSubjects,nTrainTrials);
totnlooks=zeros(2,numSubjects,nTrainTrials);
meanlookdur =zeros(2,numSubjects,nTrainTrials);
TotalLookTime=zeros(2,numSubjects,nTrainTrials);
totlonglookdur = zeros(2,numSubjects,nTrainTrials);
mLookCorrect= zeros(2,numSubjects,nObjects);mLookIncorrect=zeros(2,numSubjects,nObjects);

for iterSim=1:2 
    if iterSim ==1
        label='NEAR';
        simName = 'wPPR_1k15k_1secTest_NEAR_Metric_Variation_2019_results'
        xsit_result = load (simName);
    else
        label='FAR';
        simName = 'wPPR_1k15k_1secTest_FAR_Metric_Variation_2019_results'
        xsit_result = load (simName);
    end
    legendInfo{iterSim}= label;
    
    %xsit_result.sim.saveSettings('tesw7.json'); %visdiff ('tesw2.json', 'tesw7.json')
    numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 

% DATA SMOOTHENING
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
    
    word_On=1; word_Off= floor(TEST_DUR/scale_factor);
    vis_On = 1;vis_Off = floor(TEST_DUR/scale_factor);

    targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
    for subject=1:numSubjects
        lcorrect=0;rcorrect=0;targWord=zeros(nObjects,1);dstrWord=zeros(nObjects,1);
        twoSuccess=zeros(nObjects,1);
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
                targLookTime(iterSim,subject,trt)=lLook;
                dstrLookTime(iterSim,subject,trt)=rLook;
            elseif ( strcmp(s1,'R'))
                targLookTime(iterSim,subject,trt)=rLook;
                dstrLookTime(iterSim,subject,trt)=lLook;
            else
                  disp('ERROR reading test_pair char');
            end
        end%% trials loop

        for kk=1:nObjects
            if (targWord(kk)>dstrWord(kk))
                LearntWords(iterSim,subject,kk)=1;
            else
                LearntWords(iterSim,subject,kk)=0;
            end
        end   
    end

    

% disp('t-test statistics between Target and Distractor Looking');
% [h,p,ci,stats] = ttest(mean(targLookTime,2),mean(dstrLookTime,2),'Tail','right')
% disp('t-test statistics for words learnt');
% [h,p,ci,stats] = ttest(nansum(LearntWords,2),3,'Tail','right')

% xx=['Avg Looking time per test trial is ',num2str(mean(mean(targLookTime+dstrLookTime))*(scale_factor/1000))]; disp(xx);
% xx=['Looking time to Target per test trial  is ',num2str(mean(mean(targLookTime))*(scale_factor/1000))]; disp(xx);
% xx=['Looking time to Distractor per test trial per is ',num2str(mean(mean(dstrLookTime))*(scale_factor/1000))]; disp(xx);
% xx=['Proportion of time looking correctly (Target/Total) is ',num2str(nanmean(nanmean(correct_proportion,2)))]; disp(xx);
% xx=['Number of Strong Learners is ',num2str(sum(goodLearners))]; disp(xx);% Number models classified as strong vs weak
% xx=['Number of Weak Learners is ',num2str(numSubjects-sum(goodLearners))]; disp(xx);% 
% xx=['Looking time to Target per by strong test trial  is ',num2str(mean(mean(targLookTime(goodLearners()==1,:)))*(scale_factor/1000))]; disp(xx);
% xx=['Looking time to Distractor by strong per test trial per is ',num2str(mean(mean(dstrLookTime(goodLearners()==1,:)))*(scale_factor/1000))]; disp(xx);
% xx=['Looking time to Target per by Weak test trial  is ',num2str(mean(mean(targLookTime(goodLearners()==0,:)))*(scale_factor/1000))]; disp(xx);
% xx=['Looking time to Distractor by Weak per test trial per is ',num2str(mean(mean(dstrLookTime(goodLearners()==0,:)))*(scale_factor/1000))]; disp(xx);
% 
% xx=['Avg # of Words Learnt by Strong Learners is ',num2str(mean(sum(LearntWords(goodLearners()==1,:),2)))]; disp(xx);%  
% xx=['Avg # of Words Learnt by Weak Learners is ',num2str(mean(sum(LearntWords(goodLearners()==0,:),2)))]; disp(xx);%  
% xx=['Avg # of Words Learnt is ',num2str(mean(sum(LearntWords(),2)))]; disp(xx);%
% % stanard deviation std(sum(LearntWords(goodLearners()==1,:),2))
% xx=['Avg Looking time per test trial by Strong is ',num2str(mean(mean(targLookTime(goodLearners()==1,:)+dstrLookTime(goodLearners()==1,:)))*(scale_factor/1000))]; disp(xx);
% xx=['Avg Looking time per test trial by Weak is ',num2str(mean(mean(targLookTime(goodLearners()==0,:)+dstrLookTime(goodLearners()==0,:)))*(scale_factor/1000))]; disp(xx);


%% TRAINING CONDITION ANALYSIS trials analysis
    word_On = floor([500 2000]/scale_factor);%XSIT            
    word_Off = floor([1500 3000 ]/scale_factor);word_Len=floor(1000/scale_factor);
    C_word_On = word_On+ floor(500/scale_factor); C_word_Off = floor([2000 4000]/scale_factor);
    vis_On=1;vis_Off=(TRAIN_DUR/scale_factor); nFix_limit=10;

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
               corrLookTimeTraining(iterSim,subject,tr) = sum(savestate_historyL(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking to target object
                + sum(savestate_historyR(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking to target object right side

                incorrLookTimeTraining(iterSim,subject,tr) = sum(savestate_historyR(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking wrong way
               + sum(savestate_historyL(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking wrong way
            elseif (strcmp(s1,'X'))
               corrLookTimeTraining(iterSim,subject,tr) = sum(savestate_historyR(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking to target object
                + sum(savestate_historyL(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking to target object right side

                incorrLookTimeTraining(iterSim,subject,tr) = sum(savestate_historyL(tr,C_word_On(1):C_word_Off(1))) ... % first audio presentation, Looking wrong way
               + sum(savestate_historyR(tr,C_word_On(2):C_word_Off(2)));% 2nd audio presentation, Looking wrong way            
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

    %     for blockz=1:nObjects
    %         TinA=(nObjects-1)*(blockz-1)+1;
    %         TinB=(nObjects-1)*(blockz);
    %         tLookRepeated(subject,blockz)=sum(sum(savestate_historyL(TinA:TinB,:)));
    %         tLookVarying(subject,blockz)=sum(sum(savestate_historyR(TinA:TinB,:)));
    % 
    %         mLookCorrect(subject,blockz)= mean(corrLookTimeTraining(subject,TinA:TinB));
    %         mLookIncorrect(subject,blockz)= mean(incorrLookTimeTraining(subject,TinA:TinB));
    %     end

        totnlooks(iterSim,subject,:)=sum(nlooks,1);
        meanLukhadur(iterSim,subject,:)=nanmean(nanmean(all_look_dur,3),1);
        totlonglookdur(iterSim,subject,:)=max(longlookdur,[],1);    
        TotalLookTime(iterSim,subject,:)=sum(savestate_historyL')+sum(savestate_historyR');    
        meanlookdur(iterSim,subject,:)= TotalLookTime(iterSim,subject,:)./totnlooks(iterSim,subject,:);
         %% calculate entropy in looking on very trial
         total_trialLook_duration=nansum(nansum(all_look_dur,3),1);%1x30 from 2x30x10
         mean_trialLook_duration=nanmean(nanmean(all_look_dur,3),1);%1x30 from 2x30x10
         pdf=NaN(size(all_look_dur));%2x30x10
         variancA=NaN(size(all_look_dur));%2x30x10
         EntropySub(iterSim,subject,:)= 0;
         VarianceSub(iterSim,subject,:)=0;

         for trial=1:size(nlooks,2) 
            for side=1:size(nlooks,1)
                 pdf_side(side,:)=abs(all_look_dur(side,trial,:))./total_trialLook_duration(trial);
                 variance_side(side,:)= (all_look_dur(side,trial,:)./total_trialLook_duration(trial)) .*((all_look_dur(side,trial,:)-mean_trialLook_duration(trial)).^2) ;
            end
            pdf=[pdf_side(1,:) pdf_side(2,:)];
            EntropySub(iterSim,subject,trial)= -1* nansum( pdf(:).*log2(pdf(:)) );  
            variancA=[variance_side(1,:) variance_side(2,:)];
            VarianceSub(iterSim,subject,trial)= nansum(variancA);
         end     
    end

    %% TRACE ANALYSIS
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
    %         a_cv=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cv=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
    %         C_inTr= C_inTr+ nanmean([nanmean(a_cv(a_cv>0.001))  nanmean(b_cv(b_cv>0.001))]);
            a_cv=inputMapping1(xx1(kk),yy(kk));b_cv=inputMapping2(xx2(kk),yy(kk));
            C_inTr= C_inTr+ nanmean([a_cv b_cv]);
            inputMapping1(xx1(kk),yy(kk))=0;
            inputMapping2(xx2(kk),yy(kk))=0;
            for jj=1:6
                inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
                inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
            end
            a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
            W_inTr = W_inTr + nanmean([nanmean(a_in(a_in>0.001)) nanmean(b_in(b_in>0.001))]);
        end
        Correct_inTrace(iterSim,subject)=C_inTr/nObjects;
        Wrong_inTrace(iterSim, subject)=W_inTr/nObjects;
        InCorr_assocs(iterSim,subject)=nanmean([as_count1-1 as_count2-1]);
        EntropyTrace(iterSim,subject)= nanmean( [entropy(inputMapping1) entropy(inputMapping2)] ); 
    end


end

correct_proportion=targLookTime./(targLookTime+dstrLookTime);
disp('t-test statistics between Correct proportions ');
[h,p,ci,stats] = ttest(nanmean(correct_proportion(1,1:293,:),3),nanmean(correct_proportion(2,1:293,:),3))

 disp('t-test statistics for words learnt');
 [h,p,ci,stats] = ttest(nansum(LearntWords(1,:,:),3),nansum(LearntWords(2,:,:),3))

figure(1)% Plot total looking time during test trial
%ylim([0 1]);
%x = [1 2];
%bar(x,squeeze(nanmean(nanmean(targLookTime,3),2))*scale_factor/1000)
barwitherr([nanstd(nanmean(targLookTime(1,:,:),3))*scale_factor/1000; nanstd(nanmean(targLookTime(2,:,:),3))*scale_factor/1000], squeeze(nanmean(nanmean(targLookTime,3),2))*scale_factor/1000);
hold on
%bar(x,sts)
title('At Test')
set(gca,'xticklabel',legendInfo,'fontsize',16);
ylabel('Proportion time looking at target');
grid on

figure(2)% Plot total looking time during test trial
ylim([0 1]);
x = [1 2];
barwitherr([nanstd(nanmean(correct_proportion(1,:,:),3)); nanstd(nanmean(correct_proportion(2,:,:),3))], squeeze(nanmean(nanmean(correct_proportion,3),2)));
%bar(x,squeeze(nanmean(nanmean(correct_proportion,3),2)))
hold on
title('At Test')
set(gca,'xticklabel',legendInfo,'fontsize',16);
ylabel('Proportion correct looking');
grid on

figure(3)% Plot total looking time during test trial
%ylim([0 1]);
%x = [1 2];
barwitherr([nanstd(sum(LearntWords(1,:,:),3)); nanstd(sum(LearntWords(2,:,:),3))], squeeze(nanmean(sum(LearntWords,3),2)));
%bar(x,squeeze(nanmean(sum(LearntWords,3),2))) 
hold on
%bar(x,sts)
title('At Test')
set(gca,'xticklabel',legendInfo,'fontsize',16);
ylabel('Words Learnt');



figure(5);%Plot looking time during over TEST trials
errorbar(squeeze(nanmean(targLookTime(1,:,:),2))*scale_factor/1000,(squeeze(nanstd(targLookTime(1,:,:)))*scale_factor/1000)./sqrt(length(targLookTime(1,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(nanmean(dstrLookTime(1,:,:),2))*scale_factor/1000,(squeeze(nanstd(dstrLookTime(1,:,:)))*scale_factor/1000)./sqrt(length(dstrLookTime(1,:,:))),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('test trial');
ylabel('NEAR total looking time Target vs Distractor');
ylim([0 1]);
title('NEAR')

figure(6);%Plot looking time during over TEST trials
errorbar(squeeze(nanmean(targLookTime(2,:,:),2))*scale_factor/1000,(squeeze(nanstd(targLookTime(2,:,:)))*scale_factor/1000)./sqrt(length(targLookTime(2,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(nanmean(dstrLookTime(2,:,:),2))*scale_factor/1000,(squeeze(nanstd(dstrLookTime(2,:,:)))*scale_factor/1000)./sqrt(length(dstrLookTime(2,:,:))),plotStyle{2+compStyle});%
legend('Target','Distractor');
xlabel('test trial');
ylabel('FAR total looking time Target vs Distractor');
ylim([0 1]);
title('FAR')





%% TRAINING CONDITION ANALYSIS trials analysis

% word_On = floor([500 2000]/scale_factor);%XSIT            
% word_Off = floor([1500 3000 ]/scale_factor);word_Len=floor(1000/scale_factor);
% C_word_On = word_On+ floor(500/scale_factor); C_word_Off = floor([2000 4000]/scale_factor);
% vis_On=1;vis_Off=(TRAIN_DUR/scale_factor); nFix_limit=10;

% tLookAtLearntWords=NaN(numSubjects,nTrainTrials);
% tLookAtNonLrntWords=NaN(numSubjects,nTrainTrials);
% for subject=1:numSubjects
%     savestate_historyL = xsit_result.train(subject).historyL(:,vis_On:vis_Off);
%     savestate_historyR = xsit_result.train(subject).historyR(:,vis_On:vis_Off);    
%     nHasLWord=zeros(nTrainTrials,1); tLookAtLearnt=zeros(nTrainTrials,1);
%     nHasUWord=zeros(nTrainTrials,1); tLookAtNonLrnt=zeros(nTrainTrials,1);
%     for kk=1:nObjects
%         if (LearntWords(iterSim,subject,kk)==1)%if the subject has learnt the kk word      
%             for tt=1:nTrainTrials
%                 if(xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+1))%object exists in trial name1
%                     nHasLWord(tt)=1;
%                     s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
%                     if (strcmp(s1,'P'))% on left
%                         tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
%                     elseif (strcmp(s1,'X'))
%                         tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
%                     end
%                 elseif (xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+2))%object exists in trial
%                     nHasLWord(tt)=1;
%                     s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
%                     if (strcmp(s1,'X'))% on left
%                          tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
%                     elseif (strcmp(s1,'P'))
%                         tLookAtLearnt(tt)=tLookAtLearnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
%                     end
%                 end
%             end
%            
%         elseif (LearntWords(iterSim,subject,kk)==0)%if word is NOT learnt     
%             for tt=1:nTrainTrials
%                 if(xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+1))%object exists in trial
%                     nHasUWord(tt)=1;
%                     s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
%                     if (strcmp(s1,'P'))% on left
%                         tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
%                     elseif (strcmp(s1,'X'))
%                         tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
%                     end
%                 elseif (xsit_result.train(subject).Words{kk}==xsit_result.train(subject).training_pair(tt,2*nFeatures+2))%object exists in trial
%                     nHasUWord(tt)=1;
%                     s1=char(xsit_result.train(subject).training_pair(tt,2*nFeatures+3));
%                     if (strcmp(s1,'X'))% on left
%                          tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyL(tt,:));%add time historyL
%                     elseif (strcmp(s1,'P'))
%                         tLookAtNonLrnt(tt)=tLookAtNonLrnt(tt)+sum(savestate_historyR(tt,:));% add time historyR
%                     end
%                 end
%             end
%         else
%             disp('Error in LearntWords array: supurious data ');
%         end
%      end
%  
%     for tt=1:nTrainTrials
%         if nHasLWord(tt)==0; tLookAtLearnt(tt)=NaN;end % take only those trials wherein the word existed
%         if nHasUWord(tt)==0; tLookAtNonLrnt(tt)=NaN;end
%         tLookAtLearntWords(subject,tt)=tLookAtLearnt(tt);
%         tLookAtNonLrntWords(subject,tt)=tLookAtNonLrnt(tt);
%     end
% end



%  strongLtime=nanmean(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000;
%  SLerr=(nanstd(tLookAtLearntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==1,:)));
%  %
%  strongNLtime=nanmean(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000;
%  SNLerr=(nanstd(tLookAtNonLrntWords(goodLearners()==1,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==1,:)));
% %
%  weakLtime=nanmean(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000;
%  WLerr=(nanstd(tLookAtLearntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtLearntWords(goodLearners()==0,:)));
% % 
%  weakNLtime=nanmean(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000;
%  WNLerr=(nanstd(tLookAtNonLrntWords(goodLearners()==0,:))*scale_factor/1000)./sqrt(length(tLookAtNonLrntWords(goodLearners()==0,:)));
% %%% above nan lengths are NOT CORRECT... length(x(~isnan(x)))
% 
% figure(8);%Plot looking to learnt vs unlearnt 
% errorbar(strongLtime,SLerr,plotStyle{1});%
% hold on
% errorbar(strongNLtime,SNLerr,plotStyle{2+compStyle});%
% xlabel('per training trial');
% ylabel('Looking time (s)');
% legend('Learned Words','Non-Learned Words');
% title('Strong Learners')
% 
% figure(801)% Plot looking time at leant vs non-learnt words for Strong learners
% blockNames={'1to10','11to20','21to30'};%
% sts = [ mean(strongLtime(1:10)) mean(strongNLtime(1:10)); mean(strongLtime(11:20)) mean(strongNLtime(11:20)); mean(strongLtime(21:30)) mean(strongNLtime(21:30))];
% errY =[ mean(SLerr(1:10)) mean(SNLerr(1:10)); mean(SLerr(11:20)) mean(SNLerr(11:20)); mean(SLerr(21:30)) mean(SNLerr(21:30))];
% barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% title ('Strong Learners');
% ylabel('Looking time (ms)');
% legend('Avg Learned Words','Avg NonLearned Words');
% 
% figure(9);%Plot looking to learnt vs unlearnt 
% errorbar(weakLtime,WLerr,plotStyle{1});%
% hold on
% errorbar(weakNLtime,WNLerr,plotStyle{2+compStyle});%
% xlabel('per training trial');
% ylabel('Looking time (s)');
% legend('Learned Words','Non-Learned Words');
% title('Weak Learners')
% 
% figure(901)% Plot looking time at leant vs non-learnt words for Weak learners
% blockNames={'1to10','11to20','21to30'};%
% sts = [ mean(weakLtime(1:10)) mean(weakNLtime(1:10)); mean(weakLtime(11:20)) mean(weakNLtime(11:20)); mean(weakLtime(21:30)) mean(weakNLtime(21:30))];
% errY =[ mean(WLerr(1:10)) mean(WNLerr(1:10)); mean(WLerr(11:20)) mean(WNLerr(11:20)); mean(WLerr(21:30)) mean(WNLerr(21:30))];
% barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% title ('Weak Learners');
% ylabel('Looking time (ms)');
% legend('Avg Learnt Words','Avg Non-Learnt Words');


%%



figure (10);% Plot entropy in looking fixation durations
errorbar(squeeze(mean(VarianceSub(1,:,:)))*scale_factor/1000,(squeeze(std(VarianceSub(1,:,:))*scale_factor/1000))./sqrt( length( VarianceSub(1,:,:) )),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(squeeze(mean(VarianceSub(2,:,:)))*scale_factor/1000,(squeeze(std(VarianceSub(2,:,:))*scale_factor/1000))./sqrt( length( VarianceSub(2,:,:) )),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('Variance Near vs Far Condition');
%ylabel('number of fixations/looks Strong learners');
legend('NEAR','FAR');
%ylim([0 2.5]);
%ylim([1.5 3.5])
%hold off
summaF=mean(mean(VarianceSub(1,:,:)));
xx=['Variance NEAR',num2str(summaF)]; disp(xx);
summaF=mean(mean(VarianceSub(2,:,:)));
xx=['Variance FAR',num2str(summaF)]; disp(xx);

figure (1001);% Plot entropy in looking fixation durations
errorbar(squeeze(mean(EntropySub(1,:,:))),squeeze(std(EntropySub(1,:,:)))./sqrt( length( EntropySub(1,:,:) )),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(squeeze(mean(EntropySub(2,:,:))),squeeze(std(EntropySub(2,:,:)))./sqrt( length( EntropySub(2,:,:) )),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('Entropy Near vs far learners');
%ylabel('number of fixations/looks Strong learners');
legend('NEAR','FAR');
summaF=mean(mean(EntropySub(1,:,:)));
xx=['Entropy NEAR',num2str(summaF)]; disp(xx);
summaF=mean(mean(EntropySub(2,:,:)));
xx=['Entropy FAR',num2str(summaF)]; disp(xx);

figure(11);%Plot Strong vs Weak looking time during a training trial
errorbar(squeeze(mean(TotalLookTime(1,:,:)))*scale_factor/1000, (squeeze(std(TotalLookTime(1,:,:)))*scale_factor/1000)./sqrt(length(TotalLookTime(1,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(TotalLookTime(2,:,:)))*scale_factor/1000,(squeeze(std(TotalLookTime(2,:,:)))*scale_factor/1000)./sqrt(length(TotalLookTime(2,:,:))),plotStyle{2+compStyle});%
legend('NEAR','FAR');
xlabel('per training trial');
ylabel('total looking time NEAR vs FAR');


figure (12);% Plot number of fixations/looks over training trials
errorbar(squeeze(mean(totnlooks(1,:,:))),(squeeze(std(totnlooks(1,:,:))))./sqrt(length(totnlooks(1,:,:))),plotStyle{1});% number of fixations/looks over training trials
hold on
errorbar(squeeze(mean(totnlooks(2,:,:))),(squeeze(std(totnlooks(2,:,:))))./sqrt(length(totnlooks(2,:,:))),plotStyle{2+compStyle});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('number of fixations/looks NEAR vs FAR condition');
%ylabel('number of fixations/looks Strong learners');
legend('NEAR','FAR');
summaF=mean(mean(totnlooks(1,:,:)));
xx=['number of fixations/looks NEAR ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totnlooks(2,:,:)));
xx=['number of fixations/looks FAR  ',num2str(summaF)]; disp(xx);


figure (13);%Plot mean look duration of each fixation
errorbar(squeeze(mean(meanlookdur(1,:,:)))*scale_factor/1000, (squeeze(std(meanlookdur(1,:,:)))*scale_factor/1000)./sqrt(length(meanlookdur(1,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(meanlookdur(2,:,:)))*scale_factor/1000,(squeeze(std(meanlookdur(2,:,:)))*scale_factor/1000)./sqrt(length(meanlookdur(2,:,:))),plotStyle{2+compStyle});%
xlabel('per training trial');
ylabel('mean look duration Near Vs Far Condition');
%ylabel('mean look duration Strong learners');
legend('NEAR','FAR');
%hold off
summaF=mean(mean(meanlookdur(1,:,:))*scale_factor)/1000;
xx=['mean look duration NEAR learners ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanlookdur(2,:,:))*scale_factor)/1000;
xx=['mean look duration FAR weak learners ',num2str(summaF)]; disp(xx);


figure (1301);%Plot mean look duration of each fixation
errorbar(squeeze(mean(meanLukhadur(1,:,:)))*scale_factor/1000, (squeeze(std(meanLukhadur(1,:,:)))*scale_factor/1000)./sqrt(length(meanLukhadur(1,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(meanLukhadur(2,:,:)))*scale_factor/1000,(squeeze(std(meanLukhadur(2,:,:)))*scale_factor/1000)./sqrt(length(meanLukhadur(2,:,:))),plotStyle{2+compStyle});%
xlabel('per training trial (from durations)');
ylabel('mean look duration Near vs Far Condition');
%ylabel('mean look duration Strong learners');
legend('NEAR','FAR');
%hold off
summaF=mean(mean(meanLukhadur(1,:,:))*scale_factor)/1000;
xx=['mean look duration NEAR  (indiv calcs) ',num2str(summaF)]; disp(xx);
summaF=mean(mean(meanLukhadur(2,:,:))*scale_factor)/1000;
xx=['mean look duration FAR (indiv calcs) ',num2str(summaF)]; disp(xx);


figure (14);%Plot duration of longest look per trial
errorbar(squeeze(mean(totlonglookdur(1,:,:)))*scale_factor/1000, (squeeze(std(totlonglookdur(1,:,:)))*scale_factor/1000)./sqrt(length(totlonglookdur(1,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(totlonglookdur(2,:,:)))*scale_factor/1000,(squeeze(std(totlonglookdur(2,:,:)))*scale_factor/1000)./sqrt(length(totlonglookdur(2,:,:))),plotStyle{2+compStyle});%
xlabel('per training trial');
ylabel('duration of longest look Near vs Far learners');

%ylabel('duration of longest look Strong learners');
legend('NEAR','FAR');
%hold off
summaF=mean(mean(totlonglookdur(1,:,:))*scale_factor)/1000;
xx=['Longest look duration NEAR ',num2str(summaF)]; disp(xx);
summaF=mean(mean(totlonglookdur(2,:,:))*scale_factor)/1000;
xx=['Longest look duration FAR ',num2str(summaF)]; disp(xx);


figure(15);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(squeeze(mean(corrLookTimeTraining(1,:,:)))*scale_factor/1000, (squeeze(std(corrLookTimeTraining(1,:,:)))*scale_factor/1000)./sqrt(length(corrLookTimeTraining(1,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(incorrLookTimeTraining(1,:,:)))*scale_factor/1000,(squeeze(std(incorrLookTimeTraining(1,:,:)))*scale_factor/1000)./sqrt(length(corrLookTimeTraining(1,:,:))),plotStyle{2+compStyle});%
legend('Correct','Incorrect');
xlabel('Training Trial');
ylabel('NEAR Condition: Looking Time when words are ON');
ylim([0.2 1.75]);

figure(16);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(squeeze(mean(corrLookTimeTraining(2,:,:)))*scale_factor/1000, (squeeze(std(corrLookTimeTraining(2,:,:)))*scale_factor/1000)./sqrt(length(corrLookTimeTraining(2,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(incorrLookTimeTraining(2,:,:)))*scale_factor/1000,(squeeze(std(incorrLookTimeTraining(2,:,:)))*scale_factor/1000)./sqrt(length(corrLookTimeTraining(2,:,:))),plotStyle{2+compStyle});%
legend('Correct','Incorrect');
xlabel('Training Trial');
ylabel('FAR Condition: Looking Time when words are ON');
ylim([0.2 1.75]);

figure(1591);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(squeeze(mean(mLookCorrect(1,:,:)))*scale_factor/1000, (squeeze(std(mLookCorrect(1,:,:)))*scale_factor/1000)./sqrt(length(mLookCorrect(1,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(mLookIncorrect(1,:,:)))*scale_factor/1000,(squeeze(std(mLookIncorrect(1,:,:)))*scale_factor/1000)./sqrt(length(mLookIncorrect(1,:,:))),plotStyle{2+compStyle});%
legend('Correct','Incorrect');
xlabel('Training Trial');
ylabel('NEAR Condition: Looking Time when words are ON');
ylim([0.2 1.75]);

figure(1591);%Plot Target vs Distractor looking time during a training trial when words are ON
errorbar(squeeze(mean(mLookCorrect(2,:,:)))*scale_factor/1000, (squeeze(std(mLookCorrect(2,:,:)))*scale_factor/1000)./sqrt(length(mLookCorrect(2,:,:))),plotStyle{1});%
hold on
errorbar(squeeze(mean(mLookIncorrect(2,:,:)))*scale_factor/1000,(squeeze(std(mLookIncorrect(2,:,:)))*scale_factor/1000)./sqrt(length(mLookIncorrect(2,:,:))),plotStyle{2+compStyle});%
legend('Correct','Incorrect');
xlabel('Training Trial');
ylabel('FAR Condition: Looking Time when words are ON');
ylim([0.2 1.75]);


summaF=mean(mean(TotalLookTime(1,:,:))) *(scale_factor/1000);
xx=['NEAR Condition: Avg looking time per training trial is ',num2str(summaF)]; disp(xx);

summaF=mean(mean(TotalLookTime(2,:,:))) *(scale_factor/1000);
xx=['FAR Condition: Avg looking time per training trial is ',num2str(summaF)]; disp(xx);




%% TRACE ANALYSIS    
figure(21);%Entropy
blockNames={'NEAR';'FAR'};
sts = [mean(EntropyTrace(1,:)); mean(EntropyTrace(2,:))  ];
errY =[std(EntropyTrace(1,:))/sqrt(length(EntropyTrace(1,:))); std(EntropyTrace(2,:))/sqrt(length(EntropyTrace(2,:)))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',16);
title ('Entropy in the traces');
%ylabel('Looking time per test trial');
%ylim([0 4]);

figure(22);%My own Entropy: No of incorrect traces
blockNames={'NEAR';'FAR'};
sts = [mean(InCorr_assocs(1,:)); mean(InCorr_assocs(2,:))  ];
errY =[std(InCorr_assocs(1,:))/sqrt(length(InCorr_assocs(1,:))); std(InCorr_assocs(2,:))/sqrt(length(InCorr_assocs(2,:)))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',16);
title ('Proportion of incorrect assocs in the traces');
ylabel('# of incorrect assocs per word');

figure(221);%Strength of correct assocs in the traces
blockNames={'NEAR';'FAR'};
sts = [nanmean(Correct_inTrace(1,:)); nanmean(Correct_inTrace(2,:))  ];
errY =[nanstd(Correct_inTrace(1,:))/sqrt(length(Correct_inTrace(1,:))); nanstd(Correct_inTrace(2,:))/sqrt(length(Correct_inTrace(2,:)))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',16);
title ('Strength of correct assocs in the traces');

figure(222);%Strength of correct assocs in the traces
blockNames={'NEAR';'FAR'};
sts = [nanmean(Wrong_inTrace(1,:)); nanmean(Wrong_inTrace(2,:))  ];
errY =[nanstd(Wrong_inTrace(1,:))/sqrt(length(Wrong_inTrace(1,:))); nanstd(Wrong_inTrace(2,:))/sqrt(length(Wrong_inTrace(2,:)))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',16);
title ('Strength of INCORRECT assocs in the traces');




%%%%%%% Association hwf trace analysis
% corrAsocn=zeros(numSubjects,nObjects);
% cS=1;cW=1;
% for subject=1:numSubjects  
%     
%     AsocMat=squeeze(xsit_result.train(subject).hwf(1,:,:));
%     inputMapping=zeros(size(AsocMat));
%      for kk=1:nObjects
%          inputMapping(cell2mat(xsit_result.train(subject).Feature1(kk)),cell2mat(xsit_result.train(subject).Words(kk)))=1;
%      end 
%     
%     for kk=1:nObjects        
%         temp=[];
%         temp=AsocMat(:,cell2mat(xsit_result.train(subject).Words(kk)));  
%         maxAsocnVal(subject,kk) = max(temp);
%         [temp2 in2] = max(temp); 
%         NxtmaxAsocnVal(subject,kk)= max(max (temp(1:max(1,in2-5))), max (temp(min(size(temp),in2+5):size(temp))));
%         %NxtmaxAsocnVal(subject,kk) = max(temp(temp<max(temp)));
%         ratioMax(subject,kk)= maxAsocnVal(subject,kk)./NxtmaxAsocnVal(subject,kk);
%         prodtMR(subject,kk)=ratioMax(subject,kk).*maxAsocnVal(subject,kk);
%         
%         
%         [maxIn(kk), indIn(kk)] = max(inputMapping(:,cell2mat(xsit_result.train(subject).Words(kk))));
%         [maxAs(kk) indAs(kk)] = max(AsocMat(:,cell2mat(xsit_result.train(subject).Words(kk))));       
%         if (abs(indIn(kk)-indAs(kk)) <= 2)%if association is correct i..e same as input?
%            corrAsocn(subject, kk)=1; % wrongAssocn = 6-corrAsocn
%         end
%     end 
% end
% 
% % 
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
% varb=mean(sum(corrAsocn((goodLearners()==1),:),2)); xx=['Avg # of correctly associated words for Strong in Memory ',num2str(varb)]; disp(xx);
% varb=mean(sum(corrAsocn((goodLearners()==0),:),2)); xx=['Avg # of correctly associated words for Weak in Memory ',num2str(varb)]; disp(xx);
% % %Learnt non learnt by strong weak 
% varb=mean(corrWordsStrong); xx=['Avg # of correctly associated words in Memory for Strong counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% std(corrWordsStrong);
% varb=mean(corrWordsWeak); xx=['Avg # of correctly associated words in Memory for Weak counted as Learnt thru looking ',num2str(varb)]; disp(xx);
% std(corrWordsWeak);
% varb=sum((corrAsocn(:,1:6).*LearntWords(:,1:6)));
% xx=['# of subjects with correctly associated word in Memory counted as Learnt thru looking for ',num2str(numSubjects),' subjects is ' num2str(varb)]; disp(xx);


% for subject=1:numSubjects
%     inputMapping1=squeeze(xsit_result.train(subject).hwm_c(1,:,:));
%     inputMapping2=squeeze(xsit_result.train(subject).hwm_c(2,:,:));
%     Repeated_side(subject) = mean([mean(sum(inputMapping1(:,1:50),2))   mean(sum(inputMapping2(:,1:50),2))  ]);
%     Varying_side(subject) = mean([mean(sum(inputMapping1(:,51:100),2))  mean(sum(inputMapping2(:,51:100),2)) ]);
%     
% end
%     
% figure(20132)% Plot Target vs Distractor looking time during test trial
% blockNames={'Repeated'; 'Varying'};
% sts = [  nanmean(Repeated_side((goodLearners()==1))) nanmean(Repeated_side((goodLearners()==0)));   nanmean(Varying_side((goodLearners()==1))) nanmean(Varying_side((goodLearners()==0)));];
% errY =[  nanstd(Repeated_side((goodLearners()==1))) nanstd(Repeated_side((goodLearners()==0)));   nanstd(Varying_side((goodLearners()==1))) nanstd(Varying_side((goodLearners()==0)));];
% b=barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames);
% legend('Learners', 'Non-Learners');
% title ('Strength of scene memory trace');
% %ylabel('Proportion');
% %ylim([0 0.6]);
% set(gca,'xticklabel',blockNames,'fontsize',16);
% grid on
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
% %     surface(squeeze(xsit_result.test(subject).hwft(1,:,:)));
% %     title('Test');
% %     shading flat;
%     
% %     if (goodLearners(subject)==1), suptitle(['Strong # ' num2str(subject)]),
% %     elseif (goodLearners(subject)==0), suptitle(['Weak # ' num2str(subject)]), end
%     pause(7);
%     clf;
% end

%% TIME COURSE LOOKING IMAGEMAP LEFT RIGHT  
% subplot(2,1,1);
% lookL = xsit_result.train(1).historyL(28:30,vis_On:vis_Off)';
% lookR = xsit_result.train(1).historyR(28:30,vis_On:vis_Off)';
% vecL = lookL(:);
% vecR = lookR(:);
% d = [vecL, vecR];%size(d)
% area(d);
% %title('Run 1')
% set(gca,'fontsize',12);
% %edit change color
% subplot(2,1,2);
% lookL = xsit_result.train(9).historyL(28:30,vis_On:vis_Off)';
% lookR = xsit_result.train(9).historyR(28:30,vis_On:vis_Off)';
% vecL = lookL(:);
% vecR = lookR(:);
% c = [vecL, vecR];%size(c)
% area(c);
% %title('Run 2')
% set(gca,'fontsize',12);




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

% % %%% fixations against other measures
% fix_count_All=mean(totnlooks,2);
% TargDis_All=mean((targLookTime./(targLookTime+dstrLookTime)),2);
% Words_All=mean(LearntWords,2);
% T1=fix_count_All(:);
% S1=InCorr_assocs(:); S2=TargDis_All(:); S3=Words_All(:); S4= Correct_inTrace(:)+Wrong_inTrace(:);
% figure (31)
% scatter(T1,S1);
% yb = scatstat1(T1,S1,0,@mean);
% plot(T1,yb,'bo')
% xlabel('fixation count')
% ylabel('# incorrect associations');
% % %[fitobject,gof]=fit(T1,yb,'poly1')
% % %scatter(mean(sync_time,2),mean(targLookTime./(targLookTime+dstrLookTime),2))
% % 
% % figure (33)
% % scatter(T1,S2);
% % %hold on
% % yb = scatstat1(T1,S2,0,@mean);
% % plot(T1,yb,'bo')
% % xlabel('fixation count')
% % ylabel('Prop time looking to target');
% % 
% % figure (33)
% % scatter(T1,S3);
% % %hold on
% % yb = scatstat1(T1,S3,0,@mean);
% % plot(T1,yb,'bo')
% % xlabel('fixation count')
% % ylabel('words learnt');
% % 
% % 
% % figure (34)
% % scatter(T1,S4);
% % yb = scatstat1(T1,S4,0,@nanmean);
% % plot(T1,yb,'bo')
% % xlabel('fixation count')
% % ylabel('association strength');
% % [fitobject2,gof2]=fit(T1,yb,'poly1');
% % %[h,p,ci,stats] = ttest2(fix_count_All((goodLearners()==1)),fix_count_All((goodLearners()==0)))

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

% ((2.75*12) + (3.82*6))/18
% sqrt(((a-b)^2 + (a-b)^2 + (a-b)^2)/3)
% mean ([abs((man-mac)/man)*100  abs((man-mac)/man)*100 ])