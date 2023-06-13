%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn
% This script generates an output file
%% Data Analysis File
clear all; close all; 
%Global Variables 
nFeatures=2;scale_factor=8;legendInfo=[];MIN_LOOK_DURATION=200/scale_factor;nFix_limit=10;rmseErr=0;
%Plotting Variables
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_Young = [];  MAPE_Young = []; RMSE_Old = [];  MAPE_Old = [];
T = table (Measure, Mean, Standard_Error, RMSE_Young, MAPE_Young, RMSE_Old, MAPE_Old);
            
%Experiment (Task) Variables
TASK = 'Exp1'; SIM_TYPE = {'Silent', 'Labelling'};

if strcmp (TASK ,'Exp1')
    nObjects=31; nTestTrials=30;  TEST_DUR=6000; pre_fam_trials=1;blockSize=6;nPhases=6; 
    words_On = floor([600 3200]/scale_factor);word_Len=floor((500)/scale_factor);
    
    %EMPIRICAL DATA
    mather_dataTB(1,:)     =   [56.2   52.5    53.3    58.0    56.0];%silent young
    mather_dataTB(2,:)     =   [50.4   53.5    57.9    51.5    57.0];%labelling young
    mather_errorTB(1,:)    =   [1.3    1.5     1.6     1.6     1.6 ];%silent errorbar
    mather_errorTB(2,:)    =   [1.6    1.6     1.7     1.4     1.4 ];%labelling errorbar
    empirical_mean_young= [mather_dataTB(1,:)  mather_dataTB(2,:)];
    mather_dataTB(1,:)     =   [56.2   59.5    59.3    57.2    60.3];%silent  old
    mather_dataTB(2,:)     =   [52.5   54.6    56.2    60.9    58.2];%labelling  old
    mather_errorTB(1,:)    =   [1.3    1.5     1.6     1.6     1.6 ];%silent errorbar
    mather_errorTB(2,:)    =   [1.6    1.6     1.7     1.4     1.4 ];%labelling errorbar
    empirical_mean_old= [mather_dataTB(1,:)  mather_dataTB(2,:)];
else
    nObjects=28; nTestTrials=24;  TEST_DUR=5000; pre_fam_trials=4;blockSize=5;nPhases=5;
    words_On = floor([1500 3000]/scale_factor);word_Len=floor((700)/scale_factor);
    %EMPIRICAL DATA
    mather_dataTB(1,:)     =   [48.0   54.4    55.0    53.1    ];%silent young
    mather_dataTB(2,:)     =   [52.2   49.4    54.2    58.3    ];%labelling young
    mather_errorTB(1,:)    =   [1.3    1.5     1.6     1.6     ];%silent errorbar
    mather_errorTB(2,:)    =   [1.6    1.6     1.7     1.4     ];%labelling errorbar
    empirical_mean_young= [mather_dataTB(1,:)  mather_dataTB(2,:)];

    mather_dataTB(1,:)     =   [53.0   55.6    56.0    56.4    ];%silent old
    mather_dataTB(2,:)     =   [51.5   53.3    55.2    48.7    ];%labelling old
    mather_errorTB(1,:)    =   [1.3    1.5     1.6     1.6     ];%silent errorbar
    mather_errorTB(2,:)    =   [1.6    1.6     1.7     1.4     ];%labelling errorbar
    empirical_mean_old= [mather_dataTB(1,:)  mather_dataTB(2,:)];
end
vis_On=1; vis_Off = floor(TEST_DUR/scale_factor);trial_Blocks = floor(nTestTrials/blockSize);
novelty_pref_sim_mean=[]; novelty_pref_sim_SE =[]; empirical_mean_i =[];               
for i=1:2
    %% Raw Data File Name 
    if i==1
        label=SIM_TYPE{1};
        simName = 'WPR_Silent_MSHP_Exp1_2011_results';
        xsit_result = load (simName);
        % CHECK IF FILES ARE THERE % if not( exist(simName, 'file') == 2), disp(simName);  end                       
    else%if i==2
        label=SIM_TYPE{2};
        simName = 'WPR_Labelling_MSHP_Exp1_2011_results';
        xsit_result = load (simName); 
    end
    legendInfo{i}= label;

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

    novelty_pref_sim_mean = [novelty_pref_sim_mean squeeze(mean(nov_pref_tb))];
    novelty_pref_sim_SE   = [novelty_pref_sim_SE   squeeze(SE(nov_pref_tb))];



    figure(2);%%Plot trial-block wise preferntial looking percentage to novel object during testing 
    errorbar(squeeze(mean(nov_pref_tb)),(squeeze(std(nov_pref_tb))./squeeze(sqrt(length(nov_pref_tb)))), plotStyle{i},'LineWidth',3);%
    set(gca,'fontsize',18);
    legend(legendInfo);
    xlabel('Trial Block');
    ylabel('% novelty preference');
    if strcmp (TASK ,'Exp1'), xlim([0.5 5.5]); else, xlim([0.5 4.5]); end
    ylim([45 65]);
    hold on
    a=hline(50,'k-');
    set(a,'LineWidth',2.0);
    grid on

    figure(3);%%Plotwithin trial phase wise preferntial looking percentage to novel object during testing 
    hold on
    errorbar(squeeze(mean(mean(nov_pref_TrialPhase,2))),(squeeze(std(mean(nov_pref_TrialPhase,2)))./squeeze(sqrt(length(mean(nov_pref_TrialPhase,2))))), plotStyle{i},'LineWidth',3);%
    set(gca,'fontsize',18);
    legend(legendInfo);
    xlabel('Trial Phase');
    ylabel('% novelty preference');
    if strcmp (TASK ,'Exp1'), xlim([0.5 6.5]);ylim([40 70]); else, xlim([0.5 5.5]);ylim([40 60]); end
    hold on
    a=hline(50,'k-');
    set(a,'LineWidth',2.0);
    grid on 

end

                
%% statistics for Output File Saving 
mean_i =  novelty_pref_sim_mean;
SE_i = novelty_pref_sim_SE;
measurement_i = 'novelty preference';
RMSE_young = RMSE(empirical_mean_young, mean_i);MAPE_young = MAPE(empirical_mean_young, mean_i);
xx=['RMSE_young = ', num2str(RMSE_young),' and ', 'MAPE_young = ', num2str(MAPE_young)]; disp(xx);
RMSE_old = RMSE(empirical_mean_old, mean_i);MAPE_old = MAPE(empirical_mean_old, mean_i);
xx=['RMSE_old = ', num2str(RMSE_old),' and ', 'MAPE_old = ', num2str(MAPE_old)]; disp(xx);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_young, MAPE_young,RMSE_old,MAPE_old}; T = [T; row_i];


%% write table T to output csv file% Name your output file
writetable(T,[simName  '_Analysis.csv'])




            


