%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis %% To run the code below, you need:
%% n simulation files of name: i_Vlach_CSWL_2013_17_19_results.mat name_format where i= 1, 2,.. n and
%% n simulations of name: i_Vlach_DeBrockWOB_2017_results format
%% corresponding to n different tau values n = 6 in paper
% tau_decay values we used are [800, 1200, 2000, 3000, 3500, 5000]

clear all; close all; % This script generates an output file
n = 6; % number of tau values on which simulations have been performed
%Plotting Variables
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);

for i=1:n
    %Experiment Variables
    nObjects=12;nFeatures=2;nTrainTrials=36; nTestTrials=12;  TRAIN_DUR=4000; TEST_DUR=1000; scale_factor=8;
    compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
    %% Raw Data File Name iTH       
    simName = [num2str(i),'_Vlach_CSWL_2013_17_19_results','.mat']
    xsit_result.train=[];xsit_result.test=[];
    xsit_result = load (simName);
    numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);%
   
    vis_On = 1;vis_Off = floor(TEST_DUR/scale_factor);
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

    perform_CSL = sum(LearntWords,2);
    freq_CSL=zeros(1,nObjects);
    for subject=1:numSubjects
        for trt=1:nObjects
            if trt==perform_CSL(subject)
                freq_CSL(trt)= freq_CSL(trt)+1;
            end
        end
    end
    
    y_CSWL(i)=  mean(perform_CSL);
    yErr(i)=   SE(perform_CSL);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%% WOB File Analysis %%%%%%%%%%%
    xsit_result.train=[];    xsit_result.test=[];
    nTrials=12;nFeatures=2; scale_factor=1/8;
    OutName = [num2str(i),'_Vlach_DeBrockWOB_2017_results','.mat']
    xsit_result = load (OutName);
    numSubjects=size(xsit_result.test,1);
    freq_WOB=zeros(1,nTrials);
    perform_WOB=[];
    for subject=1:numSubjects
        corResp=0;
        for trt=1:nTrials
             lLook= sum( xsit_result.test(subject).historyLt(trt,:));%250:2000 
             rLook= sum( xsit_result.test(subject).historyRt(trt,:));%
            wordLocation = char (xsit_result.test(subject).test_word(trt,1));
            if ( strcmp(wordLocation,'L'))
                if (lLook > rLook)
                    corResp=corResp+1;
                end
            elseif (strcmp(wordLocation,'R'))
                if (rLook > lLook)
                    corResp=corResp+1;
                end
            end
        end
        perform_WOB(subject)= corResp;
        for trt=1:nTrials
            if trt==perform_WOB(subject)
                freq_WOB(trt)= freq_WOB(trt)+1;
            end
        end
    end
    x_WOB(i) =   mean(perform_WOB);
    xErr(i) =   SE(perform_WOB);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Data for Output File Saving 

% Vlach DeBrock Data empirical data fit line
x_emp_WOB = ([6 7 8 9 10 11]); % based on emprical range on the plot in Vlach 2017 
y_emp_CSWL = 4 + ((10-7)/(12-6)) * x_emp_WOB; % measured manually. intercept (4) and slope ((10-7)/(12-6))
x_emp_WOB = x_emp_WOB/nObjects;y_emp_CSWL = y_emp_CSWL/nObjects;% scale to 0-1 proportion correct

% WOLVES Model data points
% xTau_nKs = [0.8      1.2      2       3.0      3.5        5];%%parameters. tau_decay.k
% x_WOB = ([6.5000   7.000    8.2500  9.6915   10.490      11.3])/12;%% WOB 
% y_CSWL =([6.1858   6.3529   6.9875  7.9139   8.55      9.2813])/12;%;%CSWL 1k test duration
% yErr =   [0.1245   0.1439   0.1457  0.1240   0.1344    0.0886]/12;%


%% Relationship between Proprotion Corrects performance in CSWL and a word-object binding task via a regression fit
measurement_i = 'Relationship between Proprotion Corrects performance in CSWL and a word-object binding task via a regression fit';
empirical_mean_i = y_emp_CSWL;% =pred_CSWL, prediction of CSWL based on empirical fit*WOB from model sims
mean_i = y_CSWL/nObjects;%model_CSWL
SE_i = yErr/nObjects;% cswl_error
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);


figure(271)% Relationship between Proprotion Corrects performance in CSWL and a word-object binding task via a regression fit
plot(x_emp_WOB,y_emp_CSWL,'r','LineWidth',3);% vlach DeBrock data fit
hold on
errorbar(x_WOB/nObjects,y_CSWL/nObjects,yErr/nObjects,':ob','LineWidth',3,'MarkerEdgeColor','b',...
        'MarkerFaceColor','b','MarkerSize',6)% model data
xlim([0.45 1.05]);
ylim([0.3 1]);
xlabel('Word-Object Binding Memory');
ylabel('CSWL performance');
legend('Vlach & DeBrock 2017 (fit)','WOLVES Model');
set (gca, 'FontSize',16);
grid on


%% DATA FROM KACHERGIS Model
% kach_x_WOB = [0.5201517 0.5207173 0.5212066 0.5217403 0.5222435 0.5227589];
% kach_y_CSWL = [0.5238417 0.5249040 0.5264033 0.5277431 0.5288035 0.5300902];
% kach_yErr = [0.001   0.001   0.001  0.001  0.001    0.001];
% 
% errorbar(kach_x_WOB,kach_y_CSWL,kach_yErr,':ob','LineWidth',1.6,'MarkerEdgeColor','b',...
%         'MarkerFaceColor','b','MarkerSize',6)
% xlim([0.45 1.05]);
% ylim([0.3 1]);
% xlabel('Word-Object Binding Memory');
% ylabel('CSWL performance');
% legend('Vlach & DeBrock 2017 (fit)','Kachergis Model');
% set (gca, 'FontSize',14);
% grid on
% 
% %% DATA FROM STEVEN'S PURSUIT ON remember parameter 
% hold on
% pursuit_xWOB = [0.6002500 0.6737500 0.7470417 0.8285833 0.8985417 0.9737917]
% pursuit_y_CSWL = [0.5990833 0.6725833 0.7528333 0.8242917 0.9003333 0.9756667]
% pursuit_yErr = [ 0.001 0.001   0.001  0.001  0.001    0.001];
% 
% errorbar(pursuit_xWOB,pursuit_y_CSWL,pursuit_yErr,':ob','LineWidth',1.6,'MarkerEdgeColor','c',...
%         'MarkerFaceColor','c','MarkerSize',6)
% xlim([0.45 1.05]);
% ylim([0.3 1.05]);
% xlabel('Word-Object Binding Memory');
% ylabel('CSWL performance');
% legend('Vlach & DeBrock 2017 (fit)','Kachergis Model', 'Pursuit Model');
% set (gca, 'FontSize',16);
% grid on

%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])