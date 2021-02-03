%% Yu Zhong & Fricker (2012)
clear all;close all;
    %datafile: wfChanges_conwmf0_1k15k_hwmc1_fix48_wrdlen2k_new_Yu_Zhong_Fricker_2012_results
    simName = 'wPPR_1k15k_wrdlen2k_pretrain2k_Yu_Zhong_Fricker_2012_results'
    xsit_result = load (simName);
    %OutName1 = 'wfChanges_conwmf0_1k15k_hwmc1_fix63_wrdlen2k_Yu_Zhong_Fricker_2012_results'
    %xsit_result1 = load (OutName1);
    %xsit_result.train = [xsit_result.train ; xsit_result1.train];  
    %xsit_result.test = [xsit_result.test ; xsit_result1.test]; 
    plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
    nObjects = 18; scale_factor=8;%%total number of objects   
    vis_Off = 1000/scale_factor;
    nTrials=nObjects;nFeatures=2;  %check with auto file
    numSubjects=size(xsit_result.test,1);
    inA=0;inB=0;inC=0;inD=0; corrMapped=zeros(numSubjects,1);
    corrMappedWrd=zeros(numSubjects,2); Look=[];
    for subject=1:numSubjects
        Testset = xsit_result.test(subject).Testset; %size
        %Testset, has 3 distractors, 1 target/paired-target and its location
        Words=cell2mat(xsit_result.train(subject).Words);
        for trt=1:nTrials %%
            Look(1)= sum( xsit_result.test(subject).historyLt(trt,1:vis_Off));%left extreme
            Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,1:vis_Off));%left middle
            Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,1:vis_Off));%right middle
            Look(4)= sum( xsit_result.test(subject).historyRt(trt,1:vis_Off));%right extreme
            targLoc = Testset(trt,5);%the spatial location of target
            %next, add looking time to all 3 distracrors
            dstrLook=0;for i=1:4;if i~=targLoc; dstrLook= dstrLook + Look(i); end; end
            if (Look(targLoc) > dstrLook) %WARNING divide by 3
                corrMapped(subject)= corrMapped(subject)+1;
                if (Testset(trt,targLoc))== 1 || (Testset(trt,targLoc))==2 || (Testset(trt,targLoc))==3
                    corrMappedWrd(subject,1)=corrMappedWrd(subject,1)+1;
                else
                    corrMappedWrd(subject,2)=corrMappedWrd(subject,2)+1;
                end
            end   
        end
        if corrMapped(subject) > 13
            learnerType(subject) = 1;
        elseif corrMapped(subject) < 8
            learnerType(subject) = -1;
        else
            learnerType(subject) = 0;
        end
    end
    %(mean(corrMapped)/nObjects)
    %sum(learnerType==-1)

    figure(1);% proprtion correct response
    blockNames={'Yu Zhang & Fricker (2012)'; 'WOLVES Model'};
    sts = [0.58 mean(corrMapped)/nObjects];
    errY =[0.05 (std(corrMapped)/nObjects)/sqrt(length(corrMapped))];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);

    figure(2);% proprtion correct response
    blockNames={'Pre-trained words';'All other words'};
    sts = [0.9187 mean(corrMappedWrd(:,1))/3; 0.5812 mean(corrMappedWrd(:,2))/15 ];
    errY =[0.05 (std(corrMappedWrd(:,1))/3)/sqrt(length(corrMappedWrd(:,1)));0.05 (std(corrMappedWrd(:,2))/15)/sqrt(length(corrMappedWrd(:,2))) ];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'FontSize',14)
    ylabel('Proportion Correct');
    ylim ([0 1]);
    legend('Yu Zhang & Fricker (2012)', 'WOLVES Model');
    hold on
    a=hline(0.25,'k:');
    set(a,'LineWidth',2.0);
    grid on

    rmse = sqrt(((0.9187 - mean(corrMappedWrd(:,1))/3)^2 + (0.5812-mean(corrMappedWrd(:,2))/15)^2)./3);
    pe = mean ([abs(0.9187 - mean(corrMappedWrd(:,1))/3)/0.9187*100 abs(0.5812-mean(corrMappedWrd(:,2))/15)/0.5812*100 ] );
    disp (['RMSE in words learnt = ',num2str(rmse)]); 
    disp (['MAPE in words learnt = ',num2str(pe)]);
    
    % training analysis
    nTrials=27; vis_Off=11250/scale_factor;
%   %1900 word measurement as specified in paper
    word_On  = [2250 4500 6750 9000];    %word_Off = [4500 6750 9000 11250];%
    off_setter = 350;
    word_On = word_On + off_setter;
    word_Off = word_On + 1900;%1900
    word_On = floor(word_On/scale_factor); word_Off = floor(word_Off/scale_factor);
    
    propLook=NaN(numSubjects,nObjects,6);% 6 occurances of each object
    for subject=1:numSubjects
        occurance=zeros(nObjects,1);
        training_pair = xsit_result.train(subject).training_pair; % the objects/words presneted in the trial (actually from the sequence 1-18)
        training_order = xsit_result.train(subject).training_order; %the order in which objects, word were presneted 
        for trt=1:nTrials 
            obj_at_loc = training_order(trt, 1:4);
            word_seq = training_order(trt, 5:8);         
            for w=1:4 % lets take word w  %it was presneted as w temporally.. therefore time on will be word_On(w):word_Off(w)
                look(1) = sum( xsit_result.train(subject).historyL(trt,word_On(w):word_Off(w)));
                look(2) = sum( xsit_result.train(subject).historyLWB(trt,word_On(w):word_Off(w)));
                look(3) = sum( xsit_result.train(subject).historyRWB(trt,word_On(w):word_Off(w)));
                look(4) = sum( xsit_result.train(subject).historyR(trt,word_On(w):word_Off(w)));
                currword = training_pair(trt,w);
                occurance(currword) = occurance(currword)+1;
                for locs =1:4
                    if obj_at_loc(locs) == word_seq(w) % if this is the object of this currently played word
                         propLook(subject,currword,occurance(currword)) = look(locs)/((mean(look(look~=look(locs)))+look(locs)));  %looking to target vs mean looking to distractors
                    end
                end
            end
        end
    end
    kak1= 1.0:0.05:1.25;kak2= 1.0:0.025:1.125; kak3= 1.0:0.0125:1.0625;%%
    figure(34);%Plot Mean proportion correct looking at every learning instance
    errorbar(squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1)),squeeze(std(nanmean(propLook(learnerType==1,:,:),2),1))./sqrt(size(propLook(learnerType==1),2)),plotStyle{1});%./sqrt(length(mean(tarLook))
    hold on
    errorbar(squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1)),squeeze(std(nanmean(propLook(learnerType==0,:,:),2),1))./sqrt(size(propLook(learnerType==0),2)),plotStyle{2});%./sqrt(length(mean(tarLook))
    hold on
    errorbar(squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1)),squeeze(std(nanmean(propLook(learnerType==-1,:,:),2),1))./sqrt(size(propLook(learnerType==-1),2)),plotStyle{3});%./sqrt(length(mean(tarLook))
    hold on
    xlim([0.5 6.5]);
    ylim([0 1]);
    ylabel('Proportion of Time on Target');
    legend('Strong','Average','Weak');
    set (gca, 'FontSize',14);
    grid on
    
    lt_1_emp = [0.41; 0.38; 0.5; 0.59; 0.65; 0.71];
    lt_0_emp = [0.35; 0.31; 0.36; 0.41; 0.45; 0.51];
    lt_11_emp = [0.38; 0.30; 0.32; 0.37; 0.35; 0.41];
    rmse = sqrt((sum((lt_1_emp - squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1))).^2) + sum((lt_0_emp - squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1))).^2) + sum((lt_11_emp - squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1))).^2))./18);
    pe = mean ([mean(abs(lt_1_emp - squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1)))./lt_1_emp)*100  mean(abs(lt_0_emp - squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1)))./lt_1_emp)*100 mean(abs(lt_11_emp - squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1)))./lt_1_emp)*100 ] );
    disp (['RMSE in prop time on target = ',num2str(rmse)]); 
    disp (['MAPE in prop time on target = ',num2str(pe)]);
    

    nTrainTrials=27;TRAIN_DUR=11250;;scale_factor=8;MIN_LOOK_DURATION=160/scale_factor;vis_On=1;vis_Off=floor(TRAIN_DUR/scale_factor);nFix_limit=10;
    totnlooks=zeros(numSubjects,nTrainTrials);
    meanlookdur =zeros(numSubjects,nTrainTrials);
    TotalLookTime=zeros(numSubjects,nTrainTrials);
    totlonglookdur = zeros(numSubjects,nTrainTrials);
    for subject=1:numSubjects
        savestate_history(1,:,:) = xsit_result.train(subject).historyL(:,vis_On:vis_Off);
        savestate_history(2,:,:) = xsit_result.train(subject).historyLWB(:,vis_On:vis_Off);
        savestate_history(3,:,:) = xsit_result.train(subject).historyRWB(:,vis_On:vis_Off);
        savestate_history(4,:,:) = xsit_result.train(subject).historyR(:,vis_On:vis_Off);    
        %% no# of looks and fixation-durations calculation
        nlooks=zeros(4,nTrainTrials); %L/LWB/RWB/R
        longlookdur=zeros(4,nTrainTrials);
        all_look_dur=NaN(4,nTrainTrials,nFix_limit);

        for side=1:4   
            ldata = squeeze(savestate_history(side,:,:));
            
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
        TotalLookTime(subject,:)=squeeze(sum(sum(savestate_history,3)));    
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
errorbar(mean(VarianceSub)*scale_factor/1000,(std(VarianceSub)*scale_factor/1000)./sqrt( length( VarianceSub )),plotStyle{1});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('Variance ');
summaF=mean(mean(VarianceSub));xx=['Variance Strong learners ',num2str(summaF)]; disp(xx);

figure (1001);% Plot entropy in looking fixation durations
errorbar(mean(EntropySub)*scale_factor/1000,(std(EntropySub)*scale_factor/1000)./sqrt( length( EntropySub )),plotStyle{1});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('Entropy ');
summaF=mean(mean(EntropySub));xx=['Entropy Strong learners ',num2str(summaF)]; disp(xx);


figure(11);%Plot Strong vs Weak looking time during a training trial
errorbar(mean(TotalLookTime)*scale_factor/1000,(std(TotalLookTime)*scale_factor/1000)./sqrt( length( TotalLookTime )),plotStyle{1});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('total looking time ');


figure (12);% Plot number of fixations/looks over training trials
errorbar(mean(totnlooks)*scale_factor/1000,(std(totnlooks)*scale_factor/1000)./sqrt( length( totnlooks )),plotStyle{1});% number of fixations/looks over training trials
xlabel('per training trial');
ylabel('number of fixations/looks');
summaF=mean(mean(totnlooks));xx=['number of fixations/looks learners ',num2str(summaF)]; disp(xx);
sim_fix_cont(1)=mean(mean(totnlooks(learnerType()==1,:)));xx=['number of fixations/looks Strong learners ',num2str(sim_fix_cont(1))]; disp(xx);
sim_fix_cont(2) =mean(mean(totnlooks(learnerType()==0,:)));xx=['number of fixations/looks Average learners ',num2str(sim_fix_cont(2))]; disp(xx);
sim_fix_cont(3) =mean(mean(totnlooks(learnerType()==-1,:)));xx=['number of fixations/looks Weak learners ',num2str(sim_fix_cont(3))]; disp(xx);
emp_fix_cont = [2.05 2.19 2.15]*4;
rmse = sqrt(mean((emp_fix_cont - sim_fix_cont).^2));
pe = mean(abs((emp_fix_cont - sim_fix_cont)/emp_fix_cont))*100;
disp (['RMSE in fixation count = ',num2str(rmse)]); 
disp (['MAPE in fixation count = ',num2str(pe)]);



figure (13);%Plot mean look duration of each fixation
errorbar(mean(meanlookdur)*scale_factor/1000,(std(meanlookdur)*scale_factor/1000)./sqrt( length( meanlookdur )),plotStyle{1});% 
xlabel('per training trial');
ylabel('mean look duration learners');
%hold off
summaF=mean(mean(meanlookdur)*scale_factor)/1000; xx=['mean look duration Strong learners ',num2str(summaF)]; disp(xx);

figure (1301);%Plot mean look duration of each fixation
errorbar(mean(meanLukhadur)*scale_factor/1000,(std(meanLukhadur)*scale_factor/1000)./sqrt( length( meanLukhadur )),plotStyle{1});% 
xlabel('per training trial (from durations)');
ylabel('mean look duration ');
summaF=mean(mean(meanLukhadur)*scale_factor)/1000;xx=['mean look duration Strong learners (indiv calcs) ',num2str(summaF)]; disp(xx);

figure (14);%Plot duration of longest look per trial
errorbar(mean(totlonglookdur)*scale_factor/1000,(std(totlonglookdur)*scale_factor/1000)./sqrt( length( totlonglookdur )),plotStyle{1});% 
xlabel('per training trial');
ylabel('duration of longest look ');
summaF=mean(mean(totlonglookdur)*scale_factor)/1000; xx=['Longest look duration Strong learners ',num2str(summaF)]; disp(xx);
