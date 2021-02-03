
close all;
%paramValues=paramValues/10;%log10(paramValues)%paramValues/1000;
if strcmp(TASK, 'XSIT')
    word_On = floor([500 2000 4500 6000]/scale_factor);%XSIT    
    word_Off = floor([1500 3000 5500 7000]/scale_factor);word_Len=floor(1000/scale_factor);
elseif strcmp(TASK, 'XTRAP')
    word_On = floor([8 1800 3500 5200 6900]/scale_factor);  %% XTRAP    
    word_Off = floor((745+[8 1800 3500 5200 6900])/scale_factor);word_Len=floor(745/scale_factor);
end
vis_On = 1;vis_Off = floor(TEST_DUR/scale_factor);
word_Len=floor(1000/scale_factor);
compStyle=7;
parameter='scene trace hab';%

%test_TotalLookingTime_Mean= [6.7422    6.2439    5.9922    5.8249    5.6964];
%test_TotalLookingTime_Error=[0.0079    0.0079    0.0049    0.0033    0.0032];
%paramValues=[1 2 3 4 5];

figure(1);%Plot Target vs Distractor looking time during testing
errorbar((paramValues),test_TotalLookingTime_Mean,test_TotalLookingTime_Error,'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
%set(gca,'xticklabel',paramValues,'fontsize',12)
title('At Test');
ylabel('Total looking time (secs)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
%xlim([0.75 6.25]);
ylim([4 8])
xlabel(parameter);
grid on

figure(2);%Plot Target vs Distractor looking time during testing
errorbar((paramValues),test_TargetLookingTime_Mean,test_TargetLookingTime_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),test_DistractorLookingTime_Mean,test_DistractorLookingTime_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('Target','Distractor');
title('Test Trial');
ylabel('Looking time (s)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
%xlim([1.75 4.5]);
ylim([0 6]);
xlabel(parameter);
grid on
Q=test_TargetLookingTime_Mean+test_DistractorLookingTime_Mean;
QE=sqrt(test_TargetLookingTime_Error.^2 + test_DistractorLookingTime_Error.^2);
R=test_TargetLookingTime_Mean./Q;
delR= R.* sqrt(test_TargetLookingTime_Error./abs(test_TargetLookingTime_Mean) + QE./abs(Q)); 
S=test_DistractorLookingTime_Mean./Q;
delS= S.* sqrt(test_DistractorLookingTime_Error./abs(test_DistractorLookingTime_Mean) + QE./abs(Q)); 

figure(201);%Plot Target vs Distractor looking time during testing
errorbar((paramValues),test_TargetLookingTime_Mean./(test_TargetLookingTime_Mean+test_DistractorLookingTime_Mean),test_TargetLookingTime_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),test_DistractorLookingTime_Mean./(test_TargetLookingTime_Mean+test_DistractorLookingTime_Mean),test_DistractorLookingTime_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('Target','Distractor');
title('Test Trial');
ylabel('Proportion looking time');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
%xlim([1.75 4.5]);
ylim([0.2 0.8]);
xlabel(parameter);
grid on

figure(31);%Plot Target vs Distractor looking time during testing
for i=1:length(paramValues)
plot(test_STargetLookingTimePerTrial_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
hold on
%plot(test_SDistractorLookingTimePerTrial_Mean(i,:),plotStyle{i+compStyle},'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
end
xlabel('Test Trial');
ylabel('Strong Learners: Looking Time at Test');
ylim([2 6]);
legend(legendInfo);
title({'- Target; :Distractor';parameter})

figure(32);%Plot Target vs Distractor looking time during testing
for i=1:length(paramValues)
plot(test_WTargetLookingTimePerTrial_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
hold on
%plot(test_WDistractorLookingTimePerTrial_Mean(i,:),plotStyle{i+compStyle},'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
end
xlabel('Test Trial');
ylabel('Weak Learners: Looking Time at Test');
ylim([2 6]);
legend(legendInfo);
title({'- Target; :Distractor';parameter})

figure(4);%Plot Proportion of Strong/ Weak learners
errorbar((paramValues),test_ProportionStrong_Mean,test_ProportionStrong_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),test_ProportionWeak_Mean,test_ProportionWeak_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('Strong', 'Weak');
ylabel('Proportion of Learners');
ylim([0 1]);
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
xlabel(parameter);
grid on

%test_WordsLearnt_Mean(1)=test_WordsLearnt_Mean(1)+0.2;
%test_WordsLearnt_Mean(2)=test_WordsLearnt_Mean(2)-0.2;
figure(5)% Avg # of Words Learnt
hold on
errorbar((paramValues),test_WordsLearnt_Mean/6,test_WordsLearnt_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
title('Test Trial');
xlabel(parameter)
ylabel('Proportion of words learnt');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
ylim([0 1]);
grid on

figure(51)% Scatter plot
s=scatter(test_sync_time_Mean,test_correct_proportion);%
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
set (gca, 'FontSize',12);
title('Test Trial');
xlabel('sync time');
ylabel('looking to target');
grid on

figure(510)% 
errorbar((paramValues),test_sync_time_Mean,test_sync_time_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
title('Test Trial');
xlabel(parameter);
ylabel('sync time');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
grid on

% figure (1008); % plot sync time vs learning
% legend('sync mean');
% T1=targLookTime(:); S1=sync_time(:);
% scatter(S1,T1);
% ylabel('looking to target');
% xlabel('sync time');

figure(1);%Temporal analysis Proportion of looking
hold on
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
for i=1:length(paramValues)
plot((test_PropLookAtTarg_Mean(i,:)),plotStyle{i},'MarkerSize',1,'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
end
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
title('Looking to the target object at test');
legend(legendInfo);
xlabel('time (secs)');
ylabel('Proportion of looking');
ylim([0 1]);
xticklabels({'0','1.6','3.2','4.8','6.4','8.0'})

figure(611);%Temporal analysis Proportion of looking
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
for i=1:length(paramValues)
plot((test_SPropLookAtObjs_Mean(i,:)),plotStyle{i});%
hold on
% plot(test_SPropLookAtTarg_Mean(i,:),plotStyle{i},'MarkerSize',1,'LineWidth',1.5);%
% hold on
set (gca, 'FontSize',12);
end
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
title('Both Objects');
legend(legendInfo);
xlabel('time');
ylabel('Strong Learners: Proportion of looking');
ylim([0 1]);


figure(61);%Temporal analysis Proportion of looking
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
for i=1:length(paramValues)
%plot(test_SPropLookAtObjs_Mean(i,:),plotStyle{i});%
%hold on
plot((test_SPropLookAtTarg_Mean(i,:)),plotStyle{i},'MarkerSize',1,'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
end
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
title('Looking to the target object at test');
legend(legendInfo);
xlabel('time (secs)');
ylabel('Strong Learners: Proportion of looking');
ylim([0 1]);
xticklabels({'0','1.6','3.2','4.8','6.4','8.0'})

figure(612);%Temporal analysis Proportion of looking
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
for i=1:length(paramValues)
%plot(test_SPropLookAtObjs_Mean(i,:),plotStyle{i});%
%hold on
plot((test_SPropLookAtDist_Mean(i,:)),plotStyle{i},'MarkerSize',1,'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
end
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
title('Distractor Object');
legend(legendInfo);
xlabel('time');
ylabel('Strong Learners: Proportion of looking');
ylim([0 1]);


figure(613);%%Temporal analysis Proportion of looking
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
for i=1:length(paramValues)
plot((test_WPropLookAtObjs_Mean(i,:)),plotStyle{i});%
hold on
% plot(test_WPropLookAtTarg_Mean(i,:),plotStyle{i},'MarkerSize',1,'LineWidth',1.5);%
% hold on
set (gca, 'FontSize',12);
end
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
title('Both Objects');
legend(legendInfo);
xlabel('time');
ylabel('Weak Learners: Proportion of looking');
ylim([0 1]);

figure(614);%%Temporal analysis Proportion of looking
rectangle('Position',[word_On(1),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(3),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(4),0,word_Len,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
for i=1:length(paramValues)
%plot(test_WPropLookAtObjs_Mean(i,:),plotStyle{i});%
%hold on
plot((test_WPropLookAtTarg_Mean(i,:)),plotStyle{i},'MarkerSize',1,'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
end
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
title('Target Object');
legend(legendInfo);
xlabel('time');
ylabel('Weak Learners: Proportion of looking');
ylim([0 1]);



figure(17);%Plot Target vs Distractor looking time during testing
errorbar((paramValues),train_TotalLookingTime_Mean,train_TotalLookingTime_Error,'LineWidth',1.5);%
hold on
set (gca, 'FontSize',12);
%set(gca,'xticklabel',paramValues,'fontsize',12)
ylabel('Looking time (seconds) at training');
xlabel(parameter);
ylim([0 4]);
xlim([paramValues(1)-0.5 paramValues(length(paramValues))+0.5]);
grid on

figure(7);%%Plot Total looking time over training trials
for i=1:length(paramValues)
    if isnan(train_WTotalLookTime_Mean(i,:))
        plot(train_STotalLookTime_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_STotalLookTime_Mean(i,:); train_WTotalLookTime_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_WTotalLookTime_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
title(parameter);
legend(legendInfo);
xlabel('training trial');
ylabel('total looking time');
ylim([2.75 3.25]);
grid on

figure(8);%%Plot number of fixations/looks over training trials
for i=1:length(paramValues)
    if isnan(train_Wtotnlooks_Mean(i,:))
        plot(train_Stotnlooks_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_Stotnlooks_Mean(i,:); train_Wtotnlooks_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_Wtotnlooks_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
title(parameter);
legend(legendInfo);
xlabel('training trial');
ylabel('number of fixations/looks');
ylim([0 5]);
grid on



figure(9);%%Plot mean look duration of each fixation
for i=1:length(paramValues)
    if isnan(train_Wmeanlookdur_Mean(i,:))
        plot(train_Smeanlookdur_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_Smeanlookdur_Mean(i,:); train_Wmeanlookdur_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_Wmeanlookdur_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
%title({'- Strong';'-. Weak'});
legend(legendInfo);
xlabel('training trial');
ylabel('mean look duration');
ylim([0 4]);
grid on


figure(10);%%Plot duration of longest look per trial
for i=1:length(paramValues)
    if isnan(train_Wtotlonglookdur_Mean(i,:))
        plot(train_Stotlonglookdur_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_Stotlonglookdur_Mean(i,:); train_Wtotlonglookdur_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_Wtotlonglookdur_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
%title({'- Strong';'-. Weak'});
legend(legendInfo);
xlabel('training trial');
ylabel('duration of longest look');
ylim([0 4]);
grid on

figure(1088);%Plot Target vs Distractor looking time during testing
errorbar((paramValues),[1.2  1.1595    1.0765    1.0172    0.8775    0.8265],[0.015 0.01 0.005 0.008 0.011 0.013],'LineWidth',1.5);%
hold on
set (gca, 'FontSize',16);
%set(gca,'xticklabel',paramValues,'fontsize',12)
title('At Training');
ylabel('mean fixation length (secs)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
%xlim([0.75 6.25]);
ylim([0.5 2])
xlabel(parameter);
grid on




figure(10051);%%Looking time to learnt words
for i=1:length(paramValues)
    if isnan(train_Wlook_to_Learnt_words_Mean(i,:))
        plot(train_Slook_to_Learnt_words_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_Slook_to_Learnt_words_Mean(i,:); train_Wlook_to_Learnt_words_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_Wtotlonglookdur_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
title('Learnt');
legend(legendInfo);
xlabel('training trial');
ylabel('Look time to Learnt words');
ylim([0 4]);
grid on


figure(10052);%%Looking time to Non-Learnt words
for i=1:length(paramValues)
    if isnan(train_Wlook_to_NonLearnt_words_Mean(i,:))
        plot(train_Slook_to_NonLearnt_words_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_Slook_to_NonLearnt_words_Mean(i,:); train_Wlook_to_NonLearnt_words_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_Wtotlonglookdur_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
title('Non-Learnt');
legend(legendInfo);
xlabel('training trial');
ylabel('Look time to Non-Learnt words');
ylim([0 4]);
grid on

figure(13);%Plot 
errorbar((paramValues),totML_Mean,totML_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),totMN_Mean,totMN_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('Learnt','Non-Learnt');
title('After Training');
ylabel('word-feature association trace strength (max)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
%xlim([1.75 4.5]);
ylim([0 1]);
xlabel(parameter);
grid on


figure(14);%Plot max trace 
errorbar((paramValues),totMS_Mean,totMS_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),totMW_Mean,totMW_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('Strong','Weak');
title('After Training');
ylabel('word-feature association trace strength (max)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
%xlim([1.75 4.5]);
ylim([0 1]);
xlabel(parameter);
grid on

figure(141);%Plot # of Incorrect associations ratio in trace
errorbar((paramValues),train_InCorr_assocs_Mean,train_InCorr_assocs_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),train_SInCorr_assocs_Mean,train_SInCorr_assocs_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),train_WInCorr_assocs_Mean,train_WInCorr_assocs_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('All','Strong','Weak');
title('After Training');
ylabel('Proportion of InCorrect Associations (trace)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
ylim([0 1]);
xlabel(parameter);
grid on

figure(15);%Plot Entropy trace
errorbar((paramValues),train_TraceEntropy_Mean,train_TraceEntropy_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),train_STraceEntropy_Mean,train_STraceEntropy_Error,'LineWidth',1.5);%
hold on
errorbar((paramValues),train_WTraceEntropy_Mean,train_WTraceEntropy_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('All','Strong','Weak');
title('After Training');
ylabel('Entropy (trace)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
ylim([0 1]);
xlabel(parameter);
grid on

figure(151);%Plot trace mean strengths
errorbar((paramValues),train_Correct_inTrace_Strength_Mean,train_Correct_inTrace_Strength_Error,'LineWidth',1.5);%; 
hold on
errorbar((paramValues),train_Wrong_inTrace_Strength_Mean,train_Wrong_inTrace_Strength_Error,'LineWidth',1.5);%;  ;
hold on
errorbar((paramValues),(train_Correct_inTrace_Strength_Mean-train_Wrong_inTrace_Strength_Mean)./train_Correct_inTrace_Strength_Mean,train_Correct_inTrace_Strength_Error,'LineWidth',1.5);%
set (gca, 'FontSize',12);
legend('Correct Associations', 'Wrong Associations','Difference');
title('After Training');
ylabel('Mean Strength of word-feature Association (trace)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
ylim([0 1]);
xlabel(parameter);
grid on

% figure(16);%Plot trace mean strengths
% bar(train_Wrong_inTrace_Mean,'LineWidth',1.5);%
% hold on
% %errorbar((paramValues),train_Wrong_inTrace_Mean,train_Wrong_inTrace_Error,'LineWidth',1.5);%
% set (gca, 'FontSize',12);
% %legend('Correct','Wrong');
% title('After Training');
% ylabel('Association strength (trace)');
% xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
% ylim([0 0.7]);
% xlabel('old wf params                  new wf params ');
% grid on

figure(1711);%%Plot mean look duration of each fixation
for i=1:length(paramValues)
    if isnan(train_WLookEntropy_Mean(i,:))
        plot(train_SLookEntropy_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_SLookEntropy_Mean(i,:); train_WLookEntropy_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_Wmeanlookdur_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
legend(legendInfo);
grid on
title('At Training');
ylabel('Entropy (looking)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
xlim([1 30]);
xlabel('trial ');


figure(1712);%%Plot mean look duration of each fixation
for i=1:length(paramValues)
    if isnan(train_WLookVariance_Mean(i,:))
        plot(train_SLookVariance_Mean(i,:),plotStyle{i},'LineWidth',1.5);%
    else
        plot(mean([train_SLookVariance_Mean(i,:); train_WLookVariance_Mean(i,:)]),plotStyle{i},'LineWidth',1.5);%
    end
hold on
%plot(train_Wmeanlookdur_Mean(i,:),plotStyle{i+compStyle});%
%hold on
set (gca, 'FontSize',12);
end
legend(legendInfo);
grid on
title('At Training');
ylabel('Variance (looking)');
xlim([paramValues(1)-0.1 paramValues(length(paramValues))+0.1]);
xlim([1 30]);
xlabel('trial ');




test_TotalLookingTime_Mean,test_TotalLookingTime_Error
test_correct_proportion
test_ProportionStrong_Mean
test_WordsLearnt_Mean

mean([train_Smeanlookdur_Mean train_Wmeanlookdur_Mean],2)
mean([train_Stotlonglookdur_Mean train_Stotlonglookdur_Mean],2)
mean([train_Stotnlooks_Mean train_Wtotnlooks_Mean],2)
train_STotalLookTime_Mean(:,1)-train_STotalLookTime_Mean(:,30)
nanmean([train_STraceEntropy_Mean train_WTraceEntropy_Mean],2)

% % % figure(100);%Plots 
% % % x1=[ 1.8641    3.0427    3.6315    3.8942    4.0641];
% % % y1=[ 3.9167    4.4571    4.3111    3.2833    2.5833]./6;
% % % %yEr5=[0.0054    0.0053    0.0053    0.0060    0.0067    0.0079    0.0081];
% % % %errorbar(x1,y5,yEr5,'LineWidth',1.5);%
% % % plot(x1,y1,'LineWidth',1.5);%
% % % hold on 
% % % ylim([0.35 1]);
% % % set (gca, 'FontSize',12);
% % % title('Varying Spatial Attention for Object Consolidation');
% % % ylabel('Proportion (at test)');
% % % xlim([1.75 4.25]);
% % % xlabel('Number of fixations per trial (at training)');
% % %  legend('Preferntial Looking to Target','Words Learnt'); 
% % % grid on
% % % %legend('hwf -> wf','wf -> atn\_f','word -> wf','tau\_Build');
% % % %legend('atn\_sa -> atn_c','hwm\_c -> wm\_c','hwm\_f -> wm\_f', 'hwf -> wf','wf -> atn\_f','word -> wf','tau\_Build', 'tau\_Decay','noise');
% % % legend('atn\_sa -> atn_c','hwm\_c -> wm\_c','hwm\_f -> wm\_f')
% % % %legend('tau\_Build','tau\_Decay');

barColorMap(1,:) = [.2 .71 .3];	% Green Color for segment 1.
barColorMap(2,:) = [.25 .55 .79];	% Blue Color for segment 2.
barColorMap(3,:) = [.9 .1 .14];	% Red Color for segment 3.
barColorMap(4,:) = [.9 .9 .14];	% Yellow Color for segment 4.
barColorMap(5,:) = [.25 .3 .3];	% brownish Color for segment 5.   
% figure(1)% Plot total looking time during test trial
% blockNam{1}=[blockNam{1}; parVal];
% stsA{1} = [stsA{1};(mean(mean(targLookTime+dstrLookTime))*(scale_factor*slice_wordOnset4))/1000];
% errorY{1}=[errorY{1};(std(mean(targLookTime+dstrLookTime,2))*(scale_factor*slice_wordOnset4))/1000];
% %b = barwitherr(errorY{1}, stsA{1});% Plot with errorbars
% title ('Looking time at test');
% ylabel('Looking time (ms)');
% ylim([0 8]);
% xlabel('parameter');
% grid on
% hold on;
% barFontSize = 11;
% for b = 1 : length(stsA{1})
% 	% Plot one single bar as a separate bar series.
% 	handleToThisBarSeries(b) = bar(b, stsA{1}(b), 'BarWidth', 0.6);
%     set(gca,'xticklabel','');
% 	% Apply the color to this bar series.
% 	set(handleToThisBarSeries(b), 'FaceColor', barColorMap(b,:));
% 	% Place text atop the bar
% 	barTopper = sprintf('%.3f', stsA{1}(b));
%     barLabel = sprintf('%.1f', blockNam{1}(b));
% 	text(b-0.1, stsA{1}(b)+0.3, barTopper, 'FontSize', barFontSize);
%     text(b, -0.3, barLabel, 'FontSize', barFontSize);
% 	hold on;
% end

figure(22)% Plot
blockNames={'Correct Associations'; 'Wrong Associations';'Difference'};
sts = [ train_Correct_inTrace_Strength_Mean;  train_Wrong_inTrace_Strength_Mean;(train_Correct_inTrace_Strength_Mean-train_Wrong_inTrace_Strength_Mean)./train_Correct_inTrace_Strength_Mean];
errY =[ train_Correct_inTrace_Strength_Error;  train_Wrong_inTrace_Strength_Error;train_Correct_inTrace_Strength_Error-train_Wrong_inTrace_Strength_Error];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',12);
%legend('Old Model params', 'New wf params');
title('After Training (word->wf= 3.1)');
ylabel('Mean association strength (trace)');
ylim([0 1]);
grid on

figure(23);%Plot Entropy trace
blockNames={'All'; 'Strong'; 'Weak'};
sts = [ train_InCorr_assocs_Mean;   train_SInCorr_assocs_Mean;  train_WInCorr_assocs_Mean];
errY =[ train_InCorr_assocs_Error;  train_SInCorr_assocs_Error;  train_WInCorr_assocs_Error];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',12);
%legend('Old Model params', 'New wf params');
title('After Training');
ylabel('Proportion of Incorrect Associations (trace)');
ylim([0 1]);
grid on

figure(24);%Plot Entropy trace
blockNames={'All'; 'Strong'; 'Weak'};
sts = [ train_TraceEntropy_Mean;   train_STraceEntropy_Mean;  train_WTraceEntropy_Mean];
errY =[ train_TraceEntropy_Error;  train_STraceEntropy_Error;  train_WTraceEntropy_Error];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',12);
%legend('Old Model params', 'New wf params');
title('After Training');
ylabel('Entropy (trace)');
ylim([0 1]);
grid on

figure(25);%Plot
blockNames={'Target'; 'Distractor'};
sts = [ test_TargetLookingTime_Mean;  test_DistractorLookingTime_Mean];
errY =[ test_TargetLookingTime_Error;  test_DistractorLookingTime_Error];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'fontsize',12);
%legend('Old Model params', 'New wf params');
title('After Training');
ylabel('Looking at test (ms)');
ylim([0 4]);
grid on


x = 2*pi*[0 1 .1:.2:.9];
y = cos(x);
cs = csapi(x,y);
fnplt(cs,2);
axis([-1 7 -1.2 1.2])
hold on
plot(x,y,'o')
hold off

figure(200)% Plot Target vs Distractor looking time during test trial
blockNames={'Looking to target'; 'Incorrect associations'; 'Trace strength'};
sts = [  0.58 0.62;  0.69 0.48; 0.40 0.31];
errY =[  0.03  0.03; 0.006  0.007; 0.004 0.003];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
legend('WOLVES [Multi-peak]', 'WOLVES [Winner Take All]');
%title ('Performance at test');
ylabel('Proportion');
ylim([0 0.8]);
set(gca,'xticklabel',blockNames,'fontsize',12);
grid on 

%Hwf trace strength
X= [0  2  4 5 6];
y= [4.7     4.5      4.1 3.6 3.1];
yE = [0.01       0.01             0.01 0.01      0.01];

%word influence
X1= [1 2 3 4 5 6 7];
y1=[3.0000    2.6689    2.6944    4.4571    5.1824    5.1883    5.0976]/6;
yE1= [0.01        0.01          0.01        0.01      0.01         0.01         0.01] ;


figure(1)% Avg # of Words Learnt

errorbar(X,y,yE,'LineWidth',1.5);%
hold on
%errorbar(X1,y1,yE1,'LineWidth',1.5);%
set (gca, 'FontSize',16);
title('Training');
xlabel('# incorrect associations')
ylabel('Number of words learnt');
ylim([0 6]);
grid on

% Assoc
% 3.68/ (3.68 + 2.58) 
% 4.45
% 
% WTA
% 3.88 / (3.88 + 2.37)
% 
% 5.36 


figure(2013)% Plot Target vs Distractor looking time during test trial
blockNames={'Looking to Target'; 'Words Learnt'; 'Learners'};
sts = [  0.54 0.53;  4.0/6 3.83/6; 0.74 0.70];
errY =[  0.02  0.02; 0.01  0.01; 0.01 0.01];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
legend('WOLVES [Yu & Smith (2011)]', 'WOLVES [Smith  & Yu (2013)]');
%title ('Performance at test');
ylabel('Proportion');
ylim([0 1]);
set(gca,'xticklabel',blockNames,'fontsize',12);
grid on 


figure(20132)% Plot Target vs Distractor looking time during test trial
blockNames={'Incorrect traces (Prop.)'; 'Correct traces (Strength)'};
sts = [  2.6524/6 2.7612/6;   0.2592 0.2494];
errY =[  0.05/6  0.1/6;  0.0029 0.0062];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames);
legend('Learners', 'Non-Learners');
%title ('Performance at test');
%ylabel('Proportion');
ylim([0.2 0.5]);
set(gca,'xticklabel',blockNames,'fontsize',16);
grid on




x = [0.5           1.0          1.5     2.0];
y1 =[0.5186    0.5308    0.5446    0.5524];
y1E=[0.0173    0.0138    0.0137   0.0161]
y2=[0.5918    0.6488    0.7031   0.75];

figure (20133);%Plot Mean propertion of looking Varying vs repeated C
errorbar(y1,y1E,plotStyle{1}); %duration of longest look per trial
hold on
%errorbar(mean(tLookRepeated((goodLearners()==1),:)),std(tLookRepeated((goodLearners()==1),:)./(tLookVarying((goodLearners()==1),:)+tLookRepeated((goodLearners()==1),:)))./sqrt(length(tLookVarying((goodLearners()==1),:))), plotStyle{2+compStyle} ); %duration of longest look per trial
xlabel('trial block');
xlim([0.5 6.5]);
%ylim([0.3 0.65]);
ylabel('Proportion looking time Strong Learners');
legend('Varying Objects','Repeated Objects');
set (gca, 'FontSize',12);
title('Habituation over training');