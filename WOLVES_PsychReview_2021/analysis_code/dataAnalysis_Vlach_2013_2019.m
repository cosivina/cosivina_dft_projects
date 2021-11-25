%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis 
clear all; close all; % This script generates an output file
%Global Variables
scale_factor=8;MIN_LOOK_DURATION=200/scale_factor;nFeatures=2; 

%Plotting Variables
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);

%% Raw Data File Name             
simName = 'WPPR_50month_Vlach_CSWL_2013_17_19_results'
% loading the file
xsit_result = load ([simName '.mat']); 
numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 


%% Experiment Variables
new_sim_formatting = true; % new format sims with task_variables saved at simulation
if (new_sim_formatting == true)
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrainTrials = xsit_result.train(1).maxTR; 
    nTestTrials = xsit_result.train(1).maxTRt;
    TRAIN_DUR = xsit_result.train(1).t_max - floor(1000/scale_factor);
    TEST_DUR = xsit_result.train(1).t_maxt - floor(1000/scale_factor);
    word_On = xsit_result.test(1).word_On; 
    word_Off = xsit_result.test(1).word_Off; 
    vis_On = xsit_result.test(1).visuals_On;
    vis_Off = xsit_result.test(1).visuals_Off;
   else
    nObjects=12;nTrainTrials=36; nTestTrials=12;  TRAIN_DUR=4000; TEST_DUR=1000;
    %word_On = floor([500 2000 4500 6000]/scale_factor);    
    %word_Off = floor([1500 3000 5500 7000]/scale_factor);word_Len=floor(1000/scale_factor);
    vis_On = 1;vis_Off = floor(TEST_DUR/scale_factor);
end
Look_Smoothening_SY_2AFC;% DATA SMOOTHENING if necessary 
        
%% TEST CONDITION ANALYSIS
targLookTime=zeros(numSubjects,nTestTrials);
dstrLookTime=zeros(numSubjects,nTestTrials);
goodLearners=NaN(numSubjects,1);
LearntWords= NaN(numSubjects,nObjects);
targTimeS=0;distTimeS=0;targTimeW=0;distTimeW=0;
MassedTargLookTime=zeros(numSubjects,1);
MassedDistLookTime=zeros(numSubjects,1);
IleavedTargLookTime=zeros(numSubjects,1);
IleavedDistLookTime=zeros(numSubjects,1);
for subject=1:numSubjects
    lcorrect=0;rcorrect=0;targWord=zeros(nObjects,1);dstrWord=zeros(nObjects,1);
    for trt=1:nTestTrials
        lLook= sum( xsit_result.test(subject).historyLt(trt,vis_On:vis_Off));%full trial
        rLook= sum( xsit_result.test(subject).historyRt(trt,vis_On:vis_Off));%
        
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
        
        for ss= 1:nObjects/2
           if (xsit_result.train(subject).Words{ss}== xsit_result.test(subject).test_pair(trt,2*nFeatures+1))
               s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
               if ( strcmp(s1,'L')) 
                   MassedTargLookTime(subject)=MassedTargLookTime(subject)+lLook;
                   MassedDistLookTime(subject)= MassedDistLookTime(subject)+rLook;
               elseif ( strcmp(s1,'R'))
                   MassedTargLookTime(subject)= MassedTargLookTime(subject) + rLook;
                   MassedDistLookTime(subject)= MassedDistLookTime(subject)+ lLook;
               else
                   disp('ERROR reading test_pair char');
               end
           end
        end

        for ss= (nObjects/2)+1:nObjects
           if (xsit_result.train(subject).Words{ss}== xsit_result.test(subject).test_pair(trt,2*nFeatures+1))
               s1= char(xsit_result.test(subject).test_pair(trt,2*nFeatures+2));
               if ( strcmp(s1,'L')) 
                   IleavedTargLookTime(subject)= IleavedTargLookTime(subject) + lLook;
                   IleavedDistLookTime(subject)= IleavedDistLookTime(subject) + rLook;
                   lcorrect=lcorrect+lLook/(lLook+rLook);
               elseif ( strcmp(s1,'R'))
                   IleavedTargLookTime(subject)= IleavedTargLookTime(subject) + rLook;
                   IleavedDistLookTime(subject)= IleavedDistLookTime(subject) + lLook;
                   rcorrect=rcorrect+rLook/(lLook+rLook);
               else
                   disp('ERROR reading test_pair char');
               end
           end
        end%%massed/interleaved
    end%% trials loop
    
    for kk=1:nObjects
        if (targWord(kk)>dstrWord(kk))
            LearntWords(subject,kk)=1;
        else
            LearntWords(subject,kk)=0;
        end
    end  
    if (mean(targLookTime(subject,:)) > mean(dstrLookTime(subject,:)))
        goodLearners(subject)=1;
    else
        goodLearners(subject)=0;
    end    
end
correct_proportion=targLookTime./(targLookTime+dstrLookTime);

% perform_CSL = sum(LearntWords,2);
% freq_CSL=zeros(1,nObjects);
% for subject=1:numSubjects
%     for trt=1:nObjects
%         if trt==perform_CSL(subject)
%             freq_CSL(trt)= freq_CSL(trt)+1;
%         end
%     end
% end
% mean(perform_CSL)
% std(perform_CSL)./sqrt(length(perform_CSL))

%% Data for Output File Saving 
%% The simulation file above can correspond to only one of the three age groups below

%% Proportion looking to target at test vs Vlach & Johnson 2013 (16 month olds)
measurement_i = 'Proportion looking to target at test vs Vlach & Johnson 2013 (16 m)';
empirical_mean_i = [0.55 0.48];
mean_i = [mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
SE_i = [SE(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) SE(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

figure(271)% Plot Massed vs Interleaved looking time during test trial
blockNames={'Massed'; 'Interleaved'};
sts = [0.55 mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) ;  0.48 mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
errY =[0.04 SE(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime));     0.04 SE(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames, 'fontsize',18)
ylabel('Prop. Looking Time to Target');
legend('Vlach & Johnson 2013 (16 m)', 'WOLVES Model');
hold on
hline(0.05,'k:');
grid on
ylim([0 0.6]);

%% Proportion looking to target at test vs Vlach & Johnson 2013 (20 month olds)
measurement_i = 'Proportion looking to target at test vs Vlach & Johnson 2013 (20 m)';
empirical_mean_i = [0.54 0.55];
mean_i = [mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
SE_i = [SE(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) SE(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

figure(272)% Plot Massed vs Interleaved looking time during test trial
blockNames={'Massed'; 'Interleaved'};
sts = [0.54 mean(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime)) ;  0.55 mean(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
errY =[0.04 SE(MassedTargLookTime ./(MassedTargLookTime+MassedDistLookTime));    0.04  SE(IleavedTargLookTime ./(IleavedTargLookTime+IleavedDistLookTime))];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames, 'fontsize',18)
ylabel('Mean prop of looking time to target');
legend('Vlach & Johnson 2013 (20 m)', 'WOLVES Model');
grid on
%ylim([0 0.6]);

%% Proportion looking to target at test vs Vlach & DeBrock 2019 (47-58 month olds)
measurement_i = 'Proportion looking to target at test vs Vlach & Johnson 2019 (47-58 m)';
empirical_mean_i = [3.556 4.111]./6;
mean_i = [mean(sum(LearntWords(:,1:6),2)) mean(sum(LearntWords(:,7:12),2))]./6;
SE_i = [SE(sum(LearntWords(:,1:6),2)) SE(sum(LearntWords(:,7:12),2))]./6;
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

figure(273)% Plot Massed vs Interleaved looking time during test trial
blockNames={'Massed'; 'Interleaved'};
sts = [3.556 mean(sum(LearntWords(:,1:6),2)) ;  4.111 mean(sum(LearntWords(:,7:12),2)) ];
errY =[1.042./sqrt(18)  SE(sum(LearntWords(:,1:6),2)) ;  0.900./sqrt(18)  SE(sum(LearntWords(:,7:12),2)) ];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames, 'fontsize',16)
ylabel('Mean prop of looking time to target');
legend('Vlach & Johnson 2019 (47-58 m)', 'WOLVES Model');
grid on
hold on
hline(3)
%ylim([0 0.6]);


%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])