%clear all; close all;
Exp=1;
if Exp==1
    simName = 'wPPR_noisewmf2_iors13_vis20_known30_hab25_12h80h_Mather_Plunkett_2012_results'
    %xsit_result = load (simName);
%     % add another simulation
%     simName = 'wPPR_noisewmf2_iors13_known10_12h80h_hab25_Mather_Plunkett_2012_results';
%     xsit_result1 = load (simName);
%     xsit_result.train = [xsit_result.train ; xsit_result1.train];  
%     xsit_result.test = [xsit_result.test ; xsit_result1.test]; 
    plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
    nObjects = 2; scale_factor=8;%%total number of objects 
    nBlocks=2;nTrials=3;nFeatures=2;repeats=nTrials/nObjects;  %check with auto file
    vis_Off = floor(8000/scale_factor); %t_max = 1900/scale_factor;% scale with Experiment training trial Duration
    word_On =  floor([3633 5633]/scale_factor); word_Off = floor((600+[3633 5633])/scale_factor);word_Len=floor((600)/scale_factor);
    numSubjects=size(xsit_result.test,1)

    for subject=1:10
        test_pair = xsit_result.test(subject).test_pair;
        % test_pair(trt,4) is location of known,
        % test_pair(trt,5) is loc of familiar
        % test_pair(trt,6) is loc of novel
        % 
        looksToKnown=zeros(nTrials,vis_Off);looksToFamiliar=zeros(nTrials,vis_Off);looksToNovel=zeros(nTrials,vis_Off);looksOff=zeros(nTrials,vis_Off);
        for ib=1:nBlocks
            for trt=1:nTrials
                look(1,:) = xsit_result.test(subject).historyLt((ib-1)*nTrials+trt,1:vis_Off);
                look(2,:) = xsit_result.test(subject).historyCt((ib-1)*nTrials+trt,1:vis_Off);
                look(3,:) = xsit_result.test(subject).historyRt((ib-1)*nTrials+trt,1:vis_Off);
                for j=1:size(look(1,:))
                   if  (round(look(1,j)) + round(look(2,j)) + round(look(3,j))) > 0.5
                       look(4,j)=0; 
                   else; look(4,j)=1; end
                end

                looksToKnown(trt,:)    = looksToKnown(trt,:) + look(test_pair((ib-1)*nTrials+trt,4),:);
                looksToFamiliar(trt,:) = looksToFamiliar(trt,:) + look(test_pair((ib-1)*nTrials+trt,5),:);
                looksToNovel(trt,:)    = looksToNovel(trt,:) + look(test_pair((ib-1)*nTrials+trt,6),:);
                looksOff(trt,:)        = looksOff(trt,:) + look(4,:);
            end
        end
        % both blocks put together
        propKnown(subject,:,:)= looksToKnown(:,:)./(looksToKnown(:,:) + looksToFamiliar(:,:) + looksToNovel(:,:) );%+ looksOff(:,:);
        propFamiliar(subject,:,:)= looksToFamiliar(:,:)./(looksToKnown(:,:) + looksToFamiliar(:,:) + looksToNovel(:,:) );%+ looksOff(:,:);
        propNovel(subject,:,:)= looksToNovel(:,:)./(looksToKnown(:,:) + looksToFamiliar(:,:) + looksToNovel(:,:) );%+ looksOff(:,:);
        TotLooking(subject,:,:)=looksToKnown(:,:)+looksToFamiliar(:,:)+looksToNovel(:,:);
    end


    figure (1);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
    rectangle('Position',[word_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    rectangle('Position',[word_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    plot(squeeze(nanmean(nanmean(propNovel,2),1))*100,'LineWidth',3,'Color','b')
    hold on
    plot(squeeze(nanmean(nanmean(propFamiliar,2),1))*100,'LineWidth',3,'Color','g')
    hold on
    plot(squeeze(nanmean(nanmean(propKnown,2),1))*100,'LineWidth',3,'Color','r')
    hold on
    vline(floor((4000)/scale_factor),'k--')
    hline(33,'k--')
    legend('Novel','Prexposed','Known');
    xlabel('time');
    ylabel('proportion of subjects looking');
    set (gca, 'FontSize',16);
    ylim([0 60]);
    grid on
    xticklabels({'0','1600','3200','4800','6400','8000'})

    
    figure (2);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
    rectangle('Position',[word_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    rectangle('Position',[word_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
    plot(squeeze(nanmean(propNovel(:,1,:),1))*100,'LineWidth',3,'Color','r')
    hold on
    plot(squeeze(nanmean(propNovel(:,2,:),1))*100,'LineWidth',3,'Color','g')
    hold on
    plot(squeeze(nanmean(propNovel(:,3,:),1))*100,'LineWidth',3,'Color','b')
    hold on
    vline(floor((4000)/scale_factor),'k--')
    hline(33,'k--')
    legend('1st trial','2nd trial','3rd trial');
    xlabel('time');
    ylabel('proportion of subjects looking');
    set (gca, 'FontSize',16);
    ylim([20 60]);
    grid on
    xticklabels({'0','1600','3200','4800','6400','8000'})
    
    
    
    for trial=1:3
        figure (11+trial);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
        rectangle('Position',[word_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
        rectangle('Position',[word_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
        plot(squeeze(nanmean(propNovel(:,trial,:)))*100,'LineWidth',3,'Color','b')
        hold on
        plot(squeeze(nanmean(propFamiliar(:,trial,:)))*100,'LineWidth',3,'Color','g')
        hold on
        plot(squeeze(nanmean(propKnown(:,trial,:)))*100,'LineWidth',3,'Color','r')
        hold on
        vline(floor((4000)/scale_factor),'k--')
        hline(33,'k--')
        legend('Novel','Prexposed','Known');
        xlabel('time');
        ylabel('proportion of subjects looking');
        set (gca, 'FontSize',16);
        ylim([0 60]);
        title([' trial ', num2str(trial)])
        grid on
        xticklabels({'0','1600','3200','4800','6400','8000'})  
    end
    
    %% trace analysis for Mather_Plunkett
    strength_Assoc= zeros(numSubjects,6); incorr_Assoc= zeros(numSubjects,6);
    for subject=1:numSubjects    
        inputMapping1=squeeze(xsit_result.test(subject).hwf(1,:,:));
        inputMapping2=squeeze(xsit_result.test(subject).hwf(2,:,:));
        for kk=1:6
            xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
            xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
            yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
        end
        for kk=1:6
            a_cv=inputMapping1(xx1(kk),yy(kk));b_cv=inputMapping2(xx2(kk),yy(kk));
            strength_Assoc(subject, kk) = nanmean([a_cv b_cv]);
            for rkk=1:6
                if rkk ~= kk %nanmean(a_in)
                    a_cv=inputMapping1(xx1(rkk),yy(kk));
                    b_cv=inputMapping2(xx2(rkk),yy(kk)); 
                    incorr_Assoc(subject, kk) = incorr_Assoc(subject, kk) + nanmean([(a_cv>0.001) (b_cv>0.001)]);
                end
            end
        end
    end
    %%mean(strength_Assoc)
    %%mean(incorr_Assoc)
    %propNovel(subject,:,:)  OffLooking(subject,:,:) (:,:,500:1000)
    totLok=mean(mean(TotLooking,3),2)*scale_factor/nBlocks;
    propNov=mean(mean(propNovel,3),2); %(:,:,500:1000)
    T1=totLok(:);
    S1=propNov(:);

    figure (33)
    scatter(T1,S1);
    %hold on
    yb = scatstat1(T1,S1,0,@mean);
    plot(T1,yb,'bo')
    xlabel('Total Looking time')
    ylabel('Novety pref');
    set (gca, 'FontSize',16);


    figure(35);%Plot preferntial looking percentage to novel object during testing against strength of known word after 4000 msec
    %%calcuated from... mean(mean(mean(propNovel(:,:,500:1000),3),1)); std(mean(mean(propNovel(:,:,500:1000),3),1))
    Strength_KnownWord=[0.10      0.15       0.20        0.25    0.30];%       0.35        0.40        0.45];
    means_propNov =[0.3729      0.3786      0.3796      0.3851  0.3991];%       0.3895      0.3883      0.3937];
    stds_propNov  =[0.0045/20  0.0054/20  0.0118/20  0.0045/20  0.0058/20 ];%   0.0054/20   0.0030/15   0.0118/15];
    hold on
    errorbar(Strength_KnownWord, means_propNov,stds_propNov, 'LineWidth',3,'Color','b');%
    set(gca,'fontsize',18);
    xlabel('known word memory strength');
    ylabel('% novelty preference');
    xlim([0 0.35]);
    %ylim([0 100]);
    %hline(50);
    grid on 
    
    figure(36);%Plot preferntial looking percentage to novel object during testing against habituation parameter after 4000 msec
    %calcuated from... mean(mean(mean(propNovel(:,:,500:1000),3),1)); std(mean(mean(propNovel(:,:,500:1000),3),1));
    hab_param =     [1              1.5             2.5             3.5];
    means_propNov_Hab =[0.3664      0.3772          0.3796      0.3855];
    stds_propNov_Hab  =[0.1024/20   0.1000/20     0.0979/20    0.0943/20];
    hold on
    errorbar(hab_param, means_propNov_Hab,stds_propNov_Hab, 'LineWidth',3,'Color','b');%
    set(gca,'fontsize',18);
    xlabel('hwm\_c -> wm\_c');
    ylabel('% novelty preference');
    xlim([0.5 4]);
    %ylim([0 100]);
    %hline(50);
    grid on 

    figure(37);%  Plot trialwise preferntial looking percentage to novel object during testing against habituation parameter after 4000 msec
    %calcuated from...   mean(mean(propNovel(:,:,500:1000),3),1)    std(mean(propNovel(:,:,500:1000),3),1)
    hab_param =     [1              1.5             2.5             3.5];
    means_propNov_Hab1 =[0.3649      0.3774          0.3799      0.3840];
    stds_propNov_Hab1  =[0.1796/20   0.1619/20     0.1556/20    0.1561/20];
    
    means_propNov_Hab2 =[0.3631      0.3846          0.3699      0.3766];
    stds_propNov_Hab2  =[0.1768/20   0.1692/20     0.1585/20    0.1564/20];
    
    means_propNov_Hab3 =[0.3713      0.3696          0.3888      0.3958];
    stds_propNov_Hab3  =[0.1780/20   0.1577/20     0.1668/20    0.1634/20];
    errorbar(hab_param, means_propNov_Hab1,stds_propNov_Hab1, 'LineWidth',3,'Color','r');%
    hold on
    errorbar(hab_param, means_propNov_Hab2,stds_propNov_Hab2, 'LineWidth',3,'Color','g');%
    hold on
    errorbar(hab_param, means_propNov_Hab3,stds_propNov_Hab3, 'LineWidth',3,'Color','b');%
    set(gca,'fontsize',18);
    xlabel('hwm\_c -> wm\_c');
    legend('Trial 1','Trial 2','Trial 3');
    ylabel('% novelty preference');
    xlim([0.5 4]);
    %ylim([0 100]);
    %hline(50);
    grid on 
    
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
    
    
else

%% Experiment 2
%%%%#######################
%%%%##########################
clear all; close all;
simName = 'wPPR_noisewmf2_iors13_12h80h_hab35_Mather_Plunkett_2012_Exp2_results'

xsit_result = load (simName);
% add another simulation
% simName = 'wPPR_noisewmf2_iors13_12h80h_Mather_Plunkett_2012_Exp2_results';
% xsit_result1 = load (simName);
% xsit_result.train = [xsit_result.train ; xsit_result1.train];  
% xsit_result.test = [xsit_result.test ; xsit_result1.test]; 
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
nObjects = 2; scale_factor=8;%%total number of objects 
nBlocks=2;nTrials=3;nFeatures=2;repeats=nTrials/nObjects;  %check with auto file
vis_Off = floor(8000/scale_factor); %t_max = 1900/scale_factor;% scale with Experiment training trial Duration
word_On =  floor([3633 5633]/scale_factor); word_Off = floor((600+[3633 5633])/scale_factor);word_Len=floor((600)/scale_factor);
numSubjects=size(xsit_result.test,1)

for subject=1:numSubjects
    test_pair = xsit_result.test(subject).test_pair;
    % test_pair(trt,3) is loc of familiar
    % test_pair(trt,4) is loc of novel
    % 
looksToFamiliar=zeros(nTrials,vis_Off);looksToNovel=zeros(nTrials,vis_Off);looksOff=zeros(nTrials,vis_Off);
    for ib=1:nBlocks
        for trt=1:nTrials
            look(1,:) = xsit_result.test(subject).historyLt((ib-1)*nTrials+trt,1:vis_Off);
            look(2,:) = xsit_result.test(subject).historyRt((ib-1)*nTrials+trt,1:vis_Off);
            for j=1:size(look(1,:))
               if  (round(look(1,j)) + round(look(2,j))) > 0.5
                   look(3,j)=0; 
               else; look(3,j)=1; end
            end

            looksToFamiliar(trt,:) = looksToFamiliar(trt,:) + look(test_pair((ib-1)*nTrials+trt,3),:);
            looksToNovel(trt,:)    = looksToNovel(trt,:) + look(test_pair((ib-1)*nTrials+trt,4),:);
            looksOff(trt,:)       = looksOff(trt,:) + look(3,:);
        end
    end
    
    propFamiliar(subject,:,:)= looksToFamiliar(:,:)./( looksToFamiliar(:,:) + looksToNovel(:,:) );%+ looksOff;
    propNovel(subject,:,:)= looksToNovel(:,:)./( looksToFamiliar(:,:) + looksToNovel(:,:) );%+ looksOff;
    TotLooking(subject,:,:)=looksToFamiliar(:,:)+looksToNovel(:,:);
end


figure (2);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot(squeeze(nanmean(nanmean(propNovel,2),1))*100,'LineWidth',3,'Color','b')
hold on
plot(squeeze(nanmean(nanmean(propFamiliar,2),1))*100,'LineWidth',3,'Color','g')
hold on
vline(floor((4000)/scale_factor),'k--')
hline(50,'k--')
legend('Novel','Prexposed');
xlabel('time');
ylabel('proportion of subjects looking');
set (gca, 'FontSize',16);
%ylim([0 60]);
grid on
xticklabels({'0','1600','3200','4800','6400','8000'})

novelProp=squeeze(nanmean(propNovel,1))*100;
figure (3);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot(novelProp(1,:),'LineWidth',3,'Color','r')
hold on
plot(novelProp(2,:),'LineWidth',3,'Color','g')
hold on
plot(novelProp(3,:),'LineWidth',3,'Color','b')
hold on
vline(floor((4000)/scale_factor),'k--')
hline(50,'k--')
legend('1st trial','2nd trial','3rd trial');
xlabel('time');
ylabel('proportion of subjects looking at novel');
set (gca, 'FontSize',16);
ylim([20 90]);
grid on
xticklabels({'0','1600','3200','4800','6400','8000'})


totLok=mean(mean(TotLooking,3),2)*scale_factor/nBlocks;
propNov=mean(mean(propFamiliar,3),2); %%(:,:,500:1000)

T1=totLok(:);
S1=propNov(:);

figure (332)
scatter(T1,S1);
%hold on
yb = scatstat1(T1,S1,0,@mean);
plot(T1,yb,'bo')
xlabel('Total Looking time')
ylabel('Novety pref');
set (gca, 'FontSize',16);

    figure(362);%Plot preferntial looking percentage to novel object during testing against habituation parameter after 4000 msec
    %calcuated from... mean(mean(mean(propNovel(:,:,500:1000),3),1))     std(mean(mean(propNovel(:,:,500:1000),3),1))
    hab_param =     [1              1.5             2.5             3.5];
    means_propNov_Hab =[0.5449      0.5487          0.5767      0.6330];
    stds_propNov_Hab  =[0.1030/20   0.0107/20     0.0133/20    0.0586/20];
    hold on
    errorbar(hab_param, means_propNov_Hab,stds_propNov_Hab, 'LineWidth',3,'Color','b');%
    set(gca,'fontsize',18);
    xlabel('hwm\_c -> wm\_c');
    ylabel('% novelty preference');
    xlim([0.5 4]);
    %ylim([0 100]);
    %hline(50);
    grid on 

    figure(372);%  Plot trialwise preferntial looking percentage to novel object during testing against habituation parameter after 4000 msec
    %calcuated from...   mean(mean(propNovel(:,:,500:1000),3),1)    std(mean(propNovel(:,:,500:1000),3),1)
    hab_param =     [1              1.5             2.5             3.5];
    means_propNov_Hab1 =[0.5479      0.5487          0.5790      0.6537];
    stds_propNov_Hab1  =[0.1088/20   0.1042/20     0.1003/20    0.2511/20];
    
    means_propNov_Hab2 =[0.5563      0.5593          0.5886      0.5668];
    stds_propNov_Hab2  =[0.1246/20   0.1114/20     0.0992/20    0.2931/20];
    
    means_propNov_Hab3 =[0.5307      0.5380          0.5623      0.6784];
    stds_propNov_Hab3  =[0.1117/20   0.0996/20     0.0826/20    0.2695/20];
    errorbar(hab_param, means_propNov_Hab1,stds_propNov_Hab1, 'LineWidth',3,'Color','r');%
    hold on
    errorbar(hab_param, means_propNov_Hab2,stds_propNov_Hab2, 'LineWidth',3,'Color','g');%
    hold on
    errorbar(hab_param, means_propNov_Hab3,stds_propNov_Hab3, 'LineWidth',3,'Color','b');%
    set(gca,'fontsize',18);
    xlabel('hwm\_c -> wm\_c');
    legend('Trial 1','Trial 2','Trial 3');
    ylabel('% novelty preference');
    xlim([0.5 4]);
    %ylim([0 100]);
    %hline(50);
    grid on 
end
