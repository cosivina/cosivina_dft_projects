%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis for Trueswell, Medina, Hafri & Gleitman (2013)
clear all; close all; % This script generates an output file
%Global Variables
scale_factor = 8; nFeatures = 2;
%Plotting Variables
plotStyle = {'k-o','b--+','g-*','c-x','r-o','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts=[];errY=[];
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);
%%  
%% Raw Data File Name
simName='WPPR_wlen20h_spchange4_Trueswell_Medina_Hafri_Gleitman_2013_results'%
%simName ='wfChanges_conwmf0_1k15k_hmwc1_fix48_locChange_Trueswell_Medina_Hafri_Gleitman_2013_results'% PR sim
%simName = 'WPPR_wlen1250_Trueswell_Medina_Hafri_Gleitman_2013_results'% 
%simName = 'WPPR_vislen5k_w3_len1500_spchange_pred_Trueswell_Medina_Hafri_Gleitman_2013_results' %prediction

% loading the file
xsit_result = load( [simName,'.mat']);
numSubjects=size(xsit_result.train,1)
%xsit_result.sim.saveSettings('simfile.json'); %visdiff ('simfile.json', 'wolvesPaperPR.json')

%Experiment (Task) Variables
new_sim_formatting = false; % new format sims with task_variables saved at simulation
if (new_sim_formatting == true)  
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrials = xsit_result.train(1).maxTR; repeats = nTrials/nObjects;  %check with auto file
    %word_On = xsit_result.train(1).word1_On; word_Off = xsit_result.train(1).word1_Off; 
    One_look_test_len =  floor(1000/scale_factor);% mean(word_Off - word_On);
    vis_Off = xsit_result.train(1).visuals_Off;
else %old sims, need manually set task_variables below
    nObjects = 12;%%total number of objects 
    nTrials = 60; repeats = nTrials/nObjects;  %check with auto file
    One_look_test_len = floor(1000/scale_factor);
    vis_Off = floor(2400/scale_factor);
end


%% EXP TRIAL ANALYSIS
empirical_data_A = [0.18 0.22 0.27 0.31 0.33]; empirical_error_A = [0.08 0.1  0.1 0.15 0.2];
empirical_data_B = [0.205 0.47]; empirical_error_B= [0.03 0.13];
corrResp=zeros(numSubjects,repeats,nObjects);
LI_Incorrect =NaN(numSubjects,repeats,nObjects); LI_Correct =NaN(numSubjects,repeats,nObjects);
Incorr_lookingTarg=NaN(numSubjects,repeats,nObjects,vis_Off);Incorr_lookingDist=NaN(numSubjects,repeats,nObjects,vis_Off);
Corr_lookingTarg=NaN(numSubjects,repeats,nObjects,vis_Off);Corr_lookingDist=NaN(numSubjects,repeats,nObjects,vis_Off);TimeCourse=[];Look=[];
for subject = 1:numSubjects
    training_pair =  xsit_result.train(subject).training_pair;% 
    %training_pair, has 4 distractors, 1 target/paired-target, its LOCATION, and word/target sequence no     
    for rep_appear=1:repeats
        for wordy=1:nObjects
            trt=(rep_appear-1)*nObjects+wordy;
            TimeCourse(1,:) = xsit_result.train(subject).historyL(trt,1:vis_Off);%left extreme
            TimeCourse(2,:) = xsit_result.train(subject).historyLWB(trt,1:vis_Off);%left middle
            TimeCourse(3,:) = xsit_result.train(subject).historyC(trt,1:vis_Off);
            TimeCourse(4,:) = xsit_result.train(subject).historyRWB(trt,1:vis_Off);%right middle
            TimeCourse(5,:) = xsit_result.train(subject).historyR(trt,1:vis_Off);%right extreme

            Look(1)= sum( xsit_result.train(subject).historyL(trt,1:One_look_test_len));%left extreme
            Look(2)= sum( xsit_result.train(subject).historyLWB(trt,1:One_look_test_len));%left middle
            Look(3)= sum( xsit_result.train(subject).historyC(trt,1:One_look_test_len));
            Look(4)= sum( xsit_result.train(subject).historyRWB(trt,1:One_look_test_len));%right middle
            Look(5)= sum( xsit_result.train(subject).historyR(trt,1:One_look_test_len));%right extreme
            targLoc = training_pair(trt,6);%the spatial location of target
            %next, add looking time to all 4 distracrors
            dstrLook=0;for i=1:5;if i~=targLoc; dstrLook= dstrLook + Look(i); end; end
            if (Look(targLoc) > dstrLook) 
                corrResp(subject,rep_appear,wordy) = 1; % trials looked correctly
            end
            if rep_appear > 1  %previous instance analysis              
                dis=randperm(5,1); while dis==targLoc; dis=randperm(5,1);end
                if corrResp(subject,rep_appear-1,wordy) == 0 % previous instance looked incorrectly
                    LI_Incorrect(subject,rep_appear,wordy) = corrResp(subject,rep_appear,wordy);
                    Incorr_lookingTarg(subject,rep_appear,wordy,:) = TimeCourse(targLoc,:);
                    Incorr_lookingDist(subject,rep_appear,wordy,:) = TimeCourse(dis,:);
                elseif corrResp(subject,rep_appear-1,wordy) == 1 % previous instance looked correctly
                    LI_Correct(subject,rep_appear,wordy) = corrResp(subject,rep_appear,wordy);
                    Corr_lookingTarg(subject,rep_appear,wordy,:) = TimeCourse(targLoc,:);
                    Corr_lookingDist(subject,rep_appear,wordy,:) = TimeCourse(dis,:);
                end
            end
        end
    end
end

%% Data for Output File Saving 
measurement_i = 'Mean proportion correct at a learning instance';
empirical_mean_i = empirical_data_A;
mean_i = mean(mean(corrResp,3));
SE_i =   SE(mean(corrResp,3));
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


figure (131);%Plot Mean proportion correct looking at every learning instance
errorbar(mean(mean(corrResp,3)),1.96*SE(mean(corrResp,3)),plotStyle{5},'LineWidth',2); %duration of longest look per trial
hold on %Error bars show 95% confidence intervals
errorbar(empirical_data_A, empirical_error_A,plotStyle{2}, 'LineWidth',2); %d
xlabel('Learning Instance');
xlim([0.5 5.5]);
ylim([0 1]);
ylabel('Proportion Correct');
legend('WOLVES Model','Trueswell et al (2013)');
set (gca, 'FontSize',14);
grid on
movegui('east');

%% Data for Output File Saving 
measurement_i = 'Mean proportion correct at current versus previous learning instance';
empirical_mean_i = empirical_data_B;
mean_i = [ nanmean(nanmean(nanmean(LI_Incorrect,3),2)) nanmean(nanmean(nanmean(LI_Correct,3),2))];
SE_i =   [SE(nanmean(nanmean(LI_Incorrect,3),2)) SE(nanmean(nanmean(LI_Correct,3),2))];
RMSE_i = RMSE(empirical_mean_i, mean_i); MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


figure(132);%Previous Learning Instance
blockNames={'Incorrect'; 'Correct'};
sts = [empirical_data_B(1) nanmean(nanmean(nanmean(LI_Incorrect,3),2)) ; empirical_data_B(2) nanmean(nanmean(nanmean(LI_Correct,3),2)) ];
errY =[empirical_error_B(1) 1.96*SE(nanmean(nanmean(LI_Incorrect,3),2));empirical_error_B(2)  1.96*SE(nanmean(nanmean(LI_Correct,3),2)) ];
b=barwitherr(errY, sts);% %Error bars show 95% confidence intervals
set(gca,'xticklabel',blockNames,'FontSize',14)
ylabel('Proportion Correct');
ylim ([0 1]);
ylabel('Mean proportion correct');
legend('Trueswell et al (2013)','WOLVES Model');
xlabel('Previous Learning Instance');
hold on
a=hline(0.2,'k:');
set(a,'LineWidth',2.0);
grid on
movegui('west');

%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])


disp('t-test statistics for previously incorrect response');
[h,p,ci,stats] = ttest((nanmean(nanmean(LI_Incorrect,3),2)),0.2,'Tail','right');% against chance =0.2
disp(['h = ', num2str(h),'         p-value = ', num2str(p)] );


figure(15);%Previous Learning Instances Seperated
blockNames={'Incorrect'; 'Correct'};
sts = [nanmean(nanmean(LI_Incorrect,3),1) ; nanmean(nanmean(LI_Correct,3),1) ];
errY =[1.96*nanstd(nanmean(LI_Incorrect,3),1)./sqrt(length(nanmean(nanmean(LI_Incorrect,3),2)));1.96*nanstd( nanmean(LI_Correct,3),1)./sqrt(length(nanmean(nanmean(LI_Correct,3),2))) ];
b=barwitherr(errY, sts);% %Error bars show 95% confidence intervals
set(gca,'xticklabel',blockNames,'FontSize',14)
ylabel('Proportion Correct');
ylim ([0 1]);
ylabel('Mean proportion correct');
legend('Instance 1','Instance 2','Instance 3','Instance 4','Instance 5');
xlabel('Previous Learning Instance');
hold on
a=hline(0.2,'k:');
set(a,'LineWidth',2.0);
grid on
movegui('north');

figure (143);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
%rectangle('Position',[1,0,1500/8,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot(squeeze(nanmean(nanmean(nanmean(Incorr_lookingTarg,3),2),1)),plotStyle{1})
hold on
plot(squeeze(nanmean(nanmean(nanmean(Incorr_lookingDist,3),2),1)),plotStyle{2})
hold on
plot(squeeze(nanmean(nanmean(nanmean(Corr_lookingTarg,3),2),1)), plotStyle{3})
hold on
plot(squeeze(nanmean(nanmean(nanmean(Corr_lookingDist,3),2),1)),plotStyle{4})
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Prev. Incorrect Target','Prev. Incorrect Competitor','Prev. Correct Target','Prev. Correct Competitor');
xlabel('time');
ylabel('Proportion of Looks');
set (gca, 'FontSize',14);
grid on
ylim([0 1]);
movegui('south');

figure (144);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
%rectangle('Position',[1,0,1500/8,1],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot(squeeze(nanmean(nanmean(nanmean(Incorr_lookingTarg,3),2),1))-squeeze(nanmean(nanmean(nanmean(Incorr_lookingDist,3),2),1)),plotStyle{1})
hold on
%plot(squeeze(nanmean(nanmean(nanmean(Incorr_lookingDist,3),2),1)),plotStyle{2})
%hold on
plot(squeeze(nanmean(nanmean(nanmean(Corr_lookingTarg,3),2),1))-squeeze(nanmean(nanmean(nanmean(Corr_lookingDist,3),2),1)), plotStyle{3})
hold on
%plot(squeeze(nanmean(nanmean(nanmean(Corr_lookingDist,3),2),1)),plotStyle{4})
%vline(word_On(1),'r-.','On');vline(word_Off(1),'g:','Off');vline(word_On(2),'r-.','On');vline(word_Off(2),'g:','Off');vline(word_On(3),'r-.','On');vline(word_Off(3),'g:','Off');vline(word_On(4),'r-.','On');vline(word_Off(4),'g:','Off');
legend('Prev. Incorrect Target Advantage','Prev. Correct Target Advantage');
xlabel('time');
ylabel('Target Minus Competitor Looks');
set (gca, 'FontSize',14);
grid on
%ylim([0 1]);
