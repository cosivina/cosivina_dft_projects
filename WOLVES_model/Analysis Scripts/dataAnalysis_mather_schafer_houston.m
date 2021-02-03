%% Data Analysis File
clear all; 
close all;
TASK = 'Exp2';
nFeatures=2;scale_factor=8;legendInfo=[];MIN_LOOK_DURATION=200/scale_factor;nFix_limit=10;rmseErr=0;
if strcmp (TASK ,'Exp1')
    nObjects=31; nTestTrials=30;  TEST_DUR=6000; pre_fam_trials=1;blockSize=6;nPhases=6; 
    words_On = floor([600 3200]/scale_factor);word_Len=floor((500)/scale_factor);
        mather_dataTB(1,:)     =   [56.3   59.5    59.3    57.2    60.7];%silent
        mather_errorTB(1,:)    =   [1.3    1.5     1.6     1.6     1.6 ];%silent
        mather_dataTB(2,:)     =   [52.5   54.5    56.2    60.9    58.2];%labelling
        mather_errorTB(2,:)    =   [1.6    1.6     1.7     1.4     1.4 ];%labelling
else
    nObjects=28; nTestTrials=24;  TEST_DUR=5000; pre_fam_trials=4;blockSize=5;nPhases=5;
    words_On = floor([1500 3000]/scale_factor);word_Len=floor((700)/scale_factor);
        mather_dataTB(1,:)     =   [56.3   59.5    59.3    57.2    ];%silent
        mather_errorTB(1,:)    =   [1.3    1.5     1.6     1.6     ];%silent
        mather_dataTB(2,:)     =   [52.5   54.5    56.2    60.9    ];%labelling
        mather_errorTB(2,:)    =   [1.6    1.6     1.7     1.4     ];%labelling
end
vis_On=1; vis_Off = floor(TEST_DUR/scale_factor);trial_Blocks = floor(nTestTrials/blockSize);
plotStyle = {'k:o','b-+','g-*','r-s','c-x','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
for i=1:2 
    if i==1
        label='Silent';
        simName = 'wPPR_noisewmf2_iors13_12h80h_Silent_Mather_Schafer_Houston_Price_Exp2_results'
        xsit_result = load (simName);
        %     % add another simulation
%              simName = 'wPPR_noisewmf2_iors13_12h50h_Silent_Mather_Schafer_Houston_Price_results';
%              xsit_result1 = load (simName);
%              xsit_result.test = [xsit_result.test ; xsit_result1.test]; 
    else%if i==2
        label='Labelling';
        simName = 'wPPR_noisewmf2_iors13_12h80h_Labelling_Mather_Schafer_Houston_Price_Exp2_results'
        xsit_result = load (simName);
        %     % add another simulation
%              simName = 'wPPR_noisewmf2_iors13_12h50h_Labelling_Mather_Schafer_Houston_Price_results';
%              xsit_result1 = load (simName);
%              xsit_result.test = [xsit_result.test ; xsit_result1.test]; 
% %     else
% %         label='Labelling Novel';
% %         simName = 'wPPR_fix50_noisewmf2_iors13_12h60h_LabellingNovel_Mather_Schafer_Houston_Price_results'
% %         xsit_result = load (simName);
% % 
     end
    legendInfo{i}= label;
    
 %if i==2    xsit_result.sim.saveSettings('tes10.json'); end %visdiff ('wolvesPaperPR.json', 'tes10.json')
    numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 
    fam_Look=zeros(numSubjects,nTestTrials,vis_Off);
    novel_Look=zeros(numSubjects,nTestTrials,vis_Off);
    off_Look=zeros(numSubjects,nTestTrials,vis_Off);
    nov_pref_TrialPhase = zeros(numSubjects,nTestTrials-pre_fam_trials,nPhases);
    for subject=1:numSubjects
        for trt=1:nTestTrials
            lLook(:)= xsit_result.test(subject).historyLt(trt,vis_On:vis_Off);%full trial
            rLook(:)= xsit_result.test(subject).historyRt(trt,vis_On:vis_Off);%
            for j=1:length(lLook)
               if  (round(lLook(j)) + round(rLook(j))) > 0.5
                   oLook(j)=0;  else; oLook(j)=1; end
            end
            off_Look (subject,trt,:) = oLook;
            fam_Location = char (xsit_result.test(subject).fam_word(trt,1));
            if ( strcmp(fam_Location,'L'))
                fam_Look (subject,trt,:) = lLook;
                novel_Look(subject,trt,:) = rLook;            
            elseif (strcmp(fam_Location,'R'))
                fam_Look (subject,trt,:) = rLook;
                novel_Look(subject,trt,:) = lLook;
            end
        end
    end
    nov_pref = (sum(novel_Look,3)./ (sum(fam_Look,3) + sum(novel_Look,3)))*100;
    nov_pref_Phase = (novel_Look./ (fam_Look + novel_Look))*100;
    phaseSize = vis_Off/nPhases;
    for phase=1:nPhases
        phaseOn = ((phase-1)*phaseSize)+1;  phaseOff = (phase*phaseSize);
        nov_pref_TrialPhase(:,:,phase) =  mean(nov_pref_Phase(:,pre_fam_trials+1:nTestTrials,phaseOn:phaseOff),3);
    end
    

   if strcmp (TASK ,'Exp1')
       for subject=1:numSubjects
            for tb = 1:trial_Blocks         
                if tb == 1
                    nov_pref_tb(subject,tb) = mean( nov_pref(subject,(tb-1)*blockSize+1+pre_fam_trials:(tb)*blockSize), 2);
                else
                    nov_pref_tb(subject,tb) = mean( nov_pref(subject,(tb-1)*blockSize+1:(tb)*blockSize), 2);
                end
            end
       end
   else
       for subject=1:numSubjects
           for tb = 1:trial_Blocks         
                  nov_pref_tb(subject,tb) = mean( nov_pref(subject,(tb-1)*blockSize+1+pre_fam_trials:(tb)*blockSize+pre_fam_trials), 2);
           end
       end
   end
    

    figure(1);%Plot trial wise preferntial looking percentage to novel object during testing 
    hold on
    errorbar(mean(nov_pref),(std(nov_pref)./sqrt(length(nov_pref))), plotStyle{i});%
    set(gca,'fontsize',18);
    legend(legendInfo);
    xlabel('Test Trial');
    ylabel('% novelty preference');
    ylim([0 100]);
    hline(50);
    grid on  
    summaF=mean(mean(nov_pref(:,pre_fam_trials+1:nTestTrials)));xx=['mean % novelty preference in ',label,' condition is ',num2str(summaF)]; disp(xx);
    
    
    figure(2);%%Plot trial-block wise preferntial looking percentage to novel object during testing 
    errorbar(squeeze(mean(nov_pref_tb)),(squeeze(std(nov_pref_tb))./squeeze(sqrt(length(nov_pref_tb)))), plotStyle{i},'LineWidth',3);%
    set(gca,'fontsize',18);
    legend(legendInfo);
    xlabel('Trial Block');
    ylabel('% novelty preference');
    if strcmp (TASK ,'Exp1') xlim([0.5 5.5]); else xlim([0.5 4.5]); end;
    ylim([45 65]);
    hold on
    hline(50);
    grid on  
    
%     figure(201);%% Infant Data Plot trial-block wise preferntial looking percentage to novel object during testing 
%     errorbar(mather_dataTB(i,:),mather_errorTB(i,:), plotStyle{i},'LineWidth',3);%
%     hold on
%     set(gca,'fontsize',18);
%     legend(legendInfo);
%     xlabel('Trial Block');
%     ylabel('% novelty preference');
%     if strcmp (TASK ,'Exp1') xlim([0.5 5.5]); else xlim([0.5 4.5]); end;
%     ylim([45 65]);  
%     hline(50);
%     grid on 

    
    figure(3);%%Plotwithin trial phase wise preferntial looking percentage to novel object during testing 
    hold on
    errorbar(squeeze(mean(mean(nov_pref_TrialPhase,2))),(squeeze(std(mean(nov_pref_TrialPhase,2)))./squeeze(sqrt(length(mean(nov_pref_TrialPhase,2))))), plotStyle{i},'LineWidth',3);%
    set(gca,'fontsize',18);
    legend(legendInfo);
    xlabel('Trial Phase');
    ylabel('% novelty preference');
    if strcmp (TASK ,'Exp1') xlim([0.5 6.5]);ylim([40 70]); else xlim([0.5 5.5]);ylim([40 60]); end;
    hold on
    hline(50);
    grid on
    
    figure (4);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse
    hold on
    rectangle('Position',[words_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    rectangle('Position',[words_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    plot(squeeze(nanmean(nanmean(nov_pref_Phase,2),1)),plotStyle{i});
    hold on
    legend(legendInfo);
    xlabel('time');
    ylabel('proportion looking');
    set (gca, 'FontSize',18);
    hline(50,'k--')
    grid on
    xlim([vis_On vis_Off ]);
    ylim([25 75]);
    if strcmp (TASK ,'Exp1') xticklabels({'1.6','3.2','4.8'});    %xticks([0 200 400 600 ])
    else xticklabels({'0.8','1.6','2.4','3.2','4.0','4.8'}); end;   
    
    figure (41);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse
    hold on
    rectangle('Position',[words_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    rectangle('Position',[words_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    plot(squeeze(nanmean(nanmean(nov_pref_Phase(:,2:7,:),2),1)),plotStyle{i});
    hold on
    legend(legendInfo);
    xlabel('time');
    ylabel('proportion looking - start 6 trials');
    set (gca, 'FontSize',18);
    hline(50,'k--')
    grid on
    xlim([vis_On vis_Off ]);
    ylim([25 75]);
    if strcmp (TASK ,'Exp1') xticklabels({'1.6','3.2','4.8'});    %xticks([0 200 400 600 ])
    else xticklabels({'0.8','1.6','2.4','3.2','4.0','4.8'}); end;   
   
    figure (42);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse
    hold on
    rectangle('Position',[words_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    rectangle('Position',[words_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    plot(squeeze(nanmean(nanmean(nov_pref_Phase(:,8:16,:),2),1)),plotStyle{i});
    hold on
    legend(legendInfo);
    xlabel('time');
    ylabel('proportion looking - mid 8-16 trials');
    set (gca, 'FontSize',18);
    hline(50,'k--')
    grid on
    xlim([vis_On vis_Off ]);
    ylim([25 75]);
    if strcmp (TASK ,'Exp1') xticklabels({'1.6','3.2','4.8'});    %xticks([0 200 400 600 ])
    else xticklabels({'0.8','1.6','2.4','3.2','4.0','4.8'}); end;  
    
    figure (43);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse
    hold on
    rectangle('Position',[words_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    rectangle('Position',[words_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    plot(squeeze(nanmean(nanmean(nov_pref_Phase(:,nTestTrials-5:nTestTrials,:),2),1)),plotStyle{i});
    hold on
    legend(legendInfo);
    xlabel('time');
    ylabel('proportion looking - end 6 trials');
    set (gca, 'FontSize',18);
    hline(50,'k--')
    grid on
    xlim([vis_On vis_Off ]);
    ylim([25 75]);
    if strcmp (TASK ,'Exp1') xticklabels({'1.6','3.2','4.8'});    %xticks([0 200 400 600 ])
    else xticklabels({'0.8','1.6','2.4','3.2','4.0','4.8'}); end;
    

    
    %%%%%%%% Fixation data anlaysis
    corrLookTimeTraining=zeros(numSubjects,nTestTrials);%number of training trials
    incorrLookTimeTraining=zeros(numSubjects,nTestTrials);
    totnlooks=zeros(numSubjects,nTestTrials);
    meanlookdur =zeros(numSubjects,nTestTrials);
    TotalLookTime=zeros(numSubjects,nTestTrials);
    totlonglookdur = zeros(numSubjects,nTestTrials);
    for subject=1:numSubjects
        savestate_historyLt = xsit_result.test(subject).historyLt(:,vis_On:vis_Off);
        savestate_historyRt = xsit_result.test(subject).historyRt(:,vis_On:vis_Off);    
        % create the off-looking history Vector
        for trt=1:nTestTrials
            for j=1:TEST_DUR/scale_factor
               if  (round(savestate_historyLt(trt,j)) + round(savestate_historyRt(trt,j))) > 0; savestate_historyOt(trt,j)=0;               
               else savestate_historyOt(trt,j)=1; end
            end
            fam_Location = char (xsit_result.test(subject).fam_word(trt,1));
            if ( strcmp(fam_Location,'L'))
                analyse_fam_Look (trt,:) = savestate_historyLt(trt,:);
                analyse_novel_Look(trt,:) = savestate_historyRt(trt,:);            
            elseif (strcmp(fam_Location,'R'))
                analyse_fam_Look (trt,:) = savestate_historyRt(trt,:);
                analyse_novel_Look(trt,:) = savestate_historyLt(trt,:);
            end
        end

        %% no# of looks and fixation-durations calculation
        nlooks=zeros(2,nTestTrials); %L/R
        longlookdur=zeros(2,nTestTrials);
        all_look_dur=NaN(2,nTestTrials,nFix_limit);

        for side=1:2     
            if side == 1
                ldata = analyse_fam_Look;
            else
                ldata = analyse_novel_Look;
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



        totnlooks(subject,:)=sum(nlooks,1);
        meanLukhadur(subject,:)=nanmean(nanmean(all_look_dur,3),1);
        totlonglookdur(subject,:)=max(longlookdur,[],1);    
        TotalLookTime(subject,:)=sum(savestate_historyLt')+sum(savestate_historyRt');    
        meanlookdur(subject,:)= TotalLookTime(subject,:)./totnlooks(subject,:);
        
        totnlooksFam(subject,:)=nlooks(1,:);
        totnlooksNov(subject,:)=nlooks(2,:);
        meanLukhadurFam(subject,:)=nanmean(all_look_dur(1,:,:),3);
        meanLukhadurNov(subject,:)=nanmean(all_look_dur(2,:,:),3);
        totlonglookdurFam(subject,:)=max(longlookdur(1,:),[],1);
        totlonglookdurNov(subject,:)=max(longlookdur(2,:),[],1);
        TotalLookTime(subject,:)=sum(savestate_historyLt')+sum(savestate_historyRt');    
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

    figure (10);% Plot variance in looking fixation durations
    hold on
    errorbar(mean(VarianceSub)*scale_factor/1000,(std(VarianceSub)*scale_factor/1000)./sqrt( length(VarianceSub)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('Variance in looking');
    legend(legendInfo);
    ylim([20 120]);
    grid on
    summaF=mean(mean(VarianceSub));xx=['Variance ',num2str(summaF)]; disp(xx);

    figure (11);% Plot entropy in looking fixation durations
    hold on
    errorbar(mean(EntropySub),(std(EntropySub))./sqrt( length(EntropySub)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('Entropy in looking');
    legend(legendInfo);
    ylim([0.5 2.5]);
    grid on
    summaF=mean(mean(EntropySub));xx=['Entropy ',num2str(summaF)]; disp(xx);

    figure (12);% Plot total looking time
    hold on
    errorbar(mean(TotalLookTime)*scale_factor/1000,(std(TotalLookTime)*scale_factor/1000)./sqrt( length(TotalLookTime)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('Total looking time');
    legend(legendInfo);
    ylim([3 6]);
    grid on
    summaF=mean(mean(TotalLookTime))*scale_factor/1000;xx=['TotalLookTime ',num2str(summaF)]; disp(xx);

    figure (13);% Plot entropy in looking fixation durations
    hold on
    errorbar(mean(totnlooks),(std(totnlooks))./sqrt( length(totnlooks)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('number of fixations/looks');
    legend(legendInfo);
    ylim([3 5]);
    grid on
    summaF=mean(mean(totnlooks));xx=['number of fixations/looks ',num2str(summaF)]; disp(xx);

    figure (131);% Plot entropy in looking fixation durations
    hold on
    errorbar(mean(totnlooksFam),(std(totnlooksFam))./sqrt( length(totnlooksFam)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('number of fixations/looks to Familiar');
    legend(legendInfo);
    ylim([3 5]);
    grid on
    summaF=mean(mean(totnlooksFam));xx=['number of fixations/looks to Familiar ',num2str(summaF)]; disp(xx);
    
    figure (132);% Plot entropy in looking fixation durations
    hold on
    errorbar(mean(totnlooksNov),(std(totnlooksNov))./sqrt( length(totnlooksNov)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('number of fixations/looks to Novel');
    legend(legendInfo);
    ylim([3 5]);
    grid on
    summaF=mean(mean(totnlooksNov));xx=['number of fixations/looks to Novel ',num2str(summaF)]; disp(xx);
    
    figure (14);% Plot looking  durations
    hold on
    errorbar(mean(meanlookdur)*scale_factor/1000,(std(meanlookdur)*scale_factor/1000)./sqrt( length(meanlookdur)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('mean look duration');
    legend(legendInfo);
    ylim([0.75 2.5]);
    grid on
    summaF=mean(mean(meanlookdur))*scale_factor/1000;xx=['mean look duration ',num2str(summaF)]; disp(xx);

    figure (15);% Plot looking  durations
    hold on
    errorbar(mean(meanLukhadur)*scale_factor/1000,(std(meanLukhadur)*scale_factor/1000)./sqrt( length(meanLukhadur)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('mean look duration (indic calcs)');
    legend(legendInfo);
    ylim([0.75 2.5]);
    grid on
    summaF=mean(mean(meanLukhadur))*scale_factor/1000;xx=['mean look duration (indic calcs) ',num2str(summaF)]; disp(xx);

    figure (151);% Plot looking  durations
    hold on
    errorbar(nanmean(meanLukhadurFam)*scale_factor/1000,(nanstd(meanLukhadurFam)*scale_factor/1000)./sqrt( length(meanLukhadurFam)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('mean look duration to Familiar (indic calcs)');
    legend(legendInfo);
    ylim([0.75 2.5]);
    grid on
    summaF=nanmean(nanmean(meanLukhadurFam))*scale_factor/1000;xx=['mean look duration to Familiar (indic calcs) ',num2str(summaF)]; disp(xx);
    
    
    figure (152);% Plot looking  durations
    hold on
    errorbar(nanmean(meanLukhadurNov)*scale_factor/1000,(nanstd(meanLukhadurNov)*scale_factor/1000)./sqrt( length(meanLukhadurNov)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('mean look duration to Novel (indic calcs)');
    legend(legendInfo);
    ylim([0.75 2.5]);
    grid on
    summaF=nanmean(nanmean(meanLukhadurNov))*scale_factor/1000;xx=['mean look duration to Novel (indic calcs) ',num2str(summaF)]; disp(xx);
    
    
    figure (16);% Plot duration of longest look
    hold on
    errorbar(mean(totlonglookdur)*scale_factor/1000,(std(totlonglookdur)*scale_factor/1000)./sqrt( length(totlonglookdur)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('duration of longest look');
    legend(legendInfo);
    ylim([0.75 2.75]);
    grid on
    summaF=mean(mean(totlonglookdur))*scale_factor/1000;xx=['mean duration of longest look ',num2str(summaF)]; disp(xx);
    
    
    figure (161);% Plot duration of longest look
    hold on
    errorbar(nanmean(totlonglookdurFam)*scale_factor/1000,(nanstd(totlonglookdurFam)*scale_factor/1000)./sqrt( length(totlonglookdurFam)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('duration of longest look to Familiar');
    legend(legendInfo);
    ylim([0.75 2.75]);
    grid on
    summaF=mean(mean(totlonglookdurFam))*scale_factor/1000;xx=['mean duration of longest look to Familiar ',num2str(summaF)]; disp(xx);
    
    figure (162);% Plot duration of longest look
    hold on
    errorbar(nanmean(totlonglookdurNov)*scale_factor/1000,(nanstd(totlonglookdurNov)*scale_factor/1000)./sqrt( length(totlonglookdurNov)),plotStyle{i});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('duration of longest look to Novel');
    legend(legendInfo);
    ylim([0.75 2.75]);
    grid on
    summaF=mean(mean(totlonglookdurNov))*scale_factor/1000;xx=['mean duration of longest look to Novel ',num2str(summaF)]; disp(xx);
    
    
    rmseErr = rmseErr+sum((squeeze(mean(nov_pref_tb)) - mather_dataTB(i,:)).^2)/trial_Blocks;
    peErr(i)= mean(abs(squeeze(mean(nov_pref_tb)) - mather_dataTB(i,:) )./mather_dataTB(i,:))*100;
end
 rmse = sqrt(rmseErr)
 pe = mean (peErr)