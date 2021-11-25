%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis for %% Fitneva & Christiansen (2015)
clear all; close all; % This script generates an output file
%Global Variables
%Experiment Variables
nObjects = 10; %%total number of objects
nTrials=5;nFeatures=2;
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);
        


%% TEST CONDITION ANALYSIS
for s=1:3 %  loop over 3 age group simulation files
    %% Raw Data File Name 
    if s==1
        simName = 'wfChanges_conwmf0_12h3k_hwmc1_fix48_Fitneva_Christiansen_2015_results';% 4-year old's sim
    elseif s==2
        simName = 'wfChanges_conwmf0_12h55h_hwmc1_fix48_Fitneva_Christiansen_2015_results';% 10-year old's sim
    elseif s==3
        simName = 'wfChanges_conwmf0_1k15k_hwmc1_fix48_Fitneva_Christiansen_2015_results';% Adult sim
    end
    % loading the file
    xsit_result = load ([simName '.mat']);
    numSubjects=size(xsit_result.test,1)
    
    %% TEST CONDITION ANALYSIS   
    lCA=0;hCA=0;lCI=0;hCI=0;ini_Accu=zeros(numSubjects,1);ini_AccuTested=zeros(numSubjects,1);
    highIA=[];highII=[];lowIA=[];lowII=[];
    ini_NoAccu=zeros(numSubjects,1);ini_NoAccuTested=zeros(numSubjects,1);
    for subject=1:numSubjects
        Words =     cell2mat(xsit_result.train(subject).Words);
        Words_Fam = cell2mat(xsit_result.train(subject).Words_Fam);
        test_vector= cell2mat(xsit_result.test(subject).test_vector);
        for trt=1:nTrials
            lLook= sum( xsit_result.test(subject).historyLt(trt,:));%left
            rLook= sum( xsit_result.test(subject).historyRt(trt,:));%right      
            s1= char(test_vector(trt,4));
            win=0;
            if ( strcmp(s1,'L')) %%target on Left
                if (lLook > rLook); win=1; else; win=0; end
            elseif ( strcmp(s1,'R'))
                if (rLook > lLook); win=1; else; win=0; end
            end

            if ( Words(test_vector(trt,3)) == Words_Fam(test_vector(trt,3)) ) % if IA word
                ini_Accu(subject)= ini_Accu(subject)+win;
                ini_AccuTested(subject)= ini_AccuTested(subject)+1;
            else %II word
                ini_NoAccu(subject)= ini_NoAccu(subject)+win;
                ini_NoAccuTested(subject)= ini_NoAccuTested(subject)+1;
            end
        end
        subj=xsit_result.train(subject).subject;
        if mod(subj,2)== 0  %% Condition HIGH IA: only 4 pairings reshufffled
            if (ini_AccuTested(subject)>0)
                hCA=hCA+1;  highIA(hCA)= ini_Accu(subject)/ini_AccuTested(subject);
            end
            if (ini_NoAccuTested(subject)>0)
                hCI=hCI+1;  highII(hCI)= ini_NoAccu(subject)/ini_NoAccuTested(subject);
            end
        elseif mod(subj,2)== 1 %% Condition LOW IA
            if (ini_AccuTested(subject)>0)
                lCA=lCA+1; lowIA(lCA)= ini_Accu(subject)/ini_AccuTested(subject);
            end
            if (ini_NoAccuTested(subject)>0)
                lCI=lCI+1; lowII(lCI)= ini_NoAccu(subject)/ini_NoAccuTested(subject);
            end
        end
    end

    %% 3 in 1 figure
    hFig = figure(41);set(hFig, 'Position', [100 500 1800 300]);
    
    %subpanel
    subplot(1,3,s);
    blockNames={'Unchanged Pairs'; 'Changed Pairs'};
    sts = [ mean(lowIA)  mean(highIA) ; mean(lowII)  mean(highII)];
    errY =[ std(lowIA)./sqrt(length(lowIA)) std(highIA)./sqrt(length(highIA)); std(lowII)./sqrt(length(lowII)) std(highII)./sqrt(length(highII))];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'FontSize',16)
    ylabel('Mean proportion correct');
    legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
    ylim ([0 1]);
    grid on
    
    measurement_i = 'Mean proprotion correct looking at test: ';
    if s==1
        groupname = '4-year olds'; title(groupname);
        empirical_mean_i = [0.52 0.65 0.55 0.65];
        %rmse = sqrt(((mean(lowIA) - 0.52)^2 + (mean(highIA) - 0.65)^2 + (mean(lowII) - 0.55)^2 + (mean(highII) - 0.65)^2)./4)
        %pe = mean ([abs(mean(lowIA) - 0.52)/0.52*100 abs(mean(highIA) - 0.65)/0.65*100  abs(mean(lowII) - 0.55)/0.55*100 (mean(highII) - 0.65)/0.65*100])     
    elseif s==2
        groupname = '10-year olds'; title(groupname);
        empirical_mean_i = [0.75 0.75 0.63 0.57];
        %rmse = sqrt(((mean(lowIA) - 0.75)^2 + (mean(highIA) - 0.75)^2 + (mean(lowII) - 0.63)^2 + (mean(highII) - 0.57)^2)./4)
        %pe = mean ([abs(mean(lowIA) - 0.75)/0.75*100 abs(mean(highIA) - 0.75)/0.75*100  abs(mean(lowII) - 0.63)/0.63*100 (mean(highII) - 0.57)/0.57*100])
    elseif s==3
        groupname = 'Adults'; title(groupname);
        empirical_mean_i = [0.90 0.76 0.83 0.73];
        %rmse = sqrt(((mean(lowIA) - 0.90)^2 + (mean(highIA) - 0.76)^2 + (mean(lowII) - 0.83)^2 + (mean(highII) - 0.73)^2)./4)
        %pe = mean ([abs(mean(lowIA) - 0.90)/0.90*100 abs(mean(highIA) - 0.76)/0.76*100  abs(mean(lowII) - 0.83)/0.83*100 (mean(highII) - 0.73)/0.73*100])
    end
    
    %% Data for Output File Saving 
    mean_i = [mean(lowIA) mean(highIA)  mean(lowII) mean(highII)];
    SE_i =   [SE(lowIA) SE(highIA)  SE(lowII) SE(highII)];
    RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
    row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
    xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
    xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);
end

%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])


figure(4)% Plot Mean proportion correct looking standalone
blockNames={'Unchanged Pairs'; 'Changed Pairs'};
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'FontSize',16)
ylabel('Mean proportion correct');
legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
ylim ([0 1]);
grid on

figure(1)% Plot Mean proportion correct looking at words against 4YEAR OLDS
blockNames={'Low IA'; 'High IA'; 'Low II'; 'High II'};
sts = [ mean(lowIA) 0.52; mean(highIA) 0.65;  mean(lowII) 0.55; mean(highII) 0.65];
errY =[ std(lowIA)./sqrt(length(lowIA)) 0.1;  std(highIA)./sqrt(length(highIA)) 0.07; std(lowII)./sqrt(length(lowII)) 0.08; std(highII)./sqrt(length(highII)) 0.09];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'FontSize',12)
ylabel('Mean proportion correct');
ylim ([0 1]);
legend('WOLVES Model','Fitneva & Christiansen 2015 (4 yr old)');
grid on
%sqrt(((mean(lowIA) - 0.52)^2 + (mean(lowII) - 0.55)^2 + (mean(highIA) - 0.65)^2 + (mean(highII) - 0.65)^2)./4)

figure(2)% Plot Mean proportion correct looking at words against 10-yr olds
blockNames={'Low IA'; 'High IA'; 'Low II'; 'High II'};
sts = [ mean(lowIA) 0.75;  mean(highIA) 0.75; mean(lowII) 0.63; mean(highII) 0.57];
errY =[ std(lowIA)./sqrt(length(lowIA)) 0.1;  std(highIA)./sqrt(length(highIA)) 0.07; std(lowII)./sqrt(length(lowII)) 0.08; std(highII)./sqrt(length(highII)) 0.09];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'FontSize',14)
ylabel('Mean proportion correct');
ylim ([0 1]);
ylabel('Mean proportion correct');
legend('WOLVES Model','Fitneva & Christiansen 2015 (10 yr olds)');
grid on
%sqrt(((mean(lowIA) - 0.75)^2 + (mean(lowII) - 0.63)^2 + (mean(highIA) - 0.75)^2 + (mean(highII) - 0.57)^2)./4)


figure(3)% Plot Mean proportion correct looking at words against ADULTS
blockNames={'Low IA';  'High IA'; 'Low II'; 'High II'};
sts = [ mean(lowIA) 0.90;  mean(highIA) 0.76; mean(lowII) 0.83; mean(highII) 0.73];
errY =[ std(lowIA)./sqrt(length(lowIA)) 0.1;  std(highIA)./sqrt(length(highIA)) 0.07; std(lowII)./sqrt(length(lowII)) 0.08; std(highII)./sqrt(length(highII)) 0.09];
b=barwitherr(errY, sts);% Plot with errorbars
set(gca,'xticklabel',blockNames,'FontSize',14)
ylabel('Mean proportion correct');
ylim ([0 1]);
ylabel('Mean proportion correct');
legend('WOLVES Model','Fitneva & Christiansen 2015 (Adults)');
grid on
%sqrt(((mean(lowIA) - 0.90)^2 + (mean(lowII) - 0.83)^2 + (mean(highIA) - 0.76)^2 + (mean(highII) - 0.73)^2)./4)



%% EMPIRICAL DATA PLOT 3 in 1
% hFig = figure(51);%set(hFig, 'Position', [100 500 1800 300]);
% %4-yr olds
% subplot(1,3,1);
% blockNames={'Unchanged Pairs'; 'Changed Pairs'};
% sts = [ 0.52  0.65 ; 0.55  0.65];
% errY =[ 0.1   0.07;  0.08  0.09];
% b=barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames,'FontSize',16)
% ylabel('Mean proportion correct');
% legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
% ylim ([0 1]);
% title('4-year olds');
% grid on
% 
% %10-yr olds
% subplot(1,3,2);
% blockNames={'Unchanged Pairs'; 'Changed Pairs'};
% sts = [ 0.75  0.75 ; 0.63  0.57];
% errY =[ 0.1   0.07;  0.08  0.09];
% b=barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames,'FontSize',16)
% %ylabel('Mean proportion correct');
% %legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
% ylim ([0 1]);
% title('10-year olds');
% grid on
% 
% %Adults
% subplot(1,3,3);
% blockNames={'Unchanged Pairs'; 'Changed Pairs'};
% sts = [ 0.90  0.76 ; 0.83  0.73];
% errY =[ 0.1   0.07;  0.08  0.09];
% b=barwitherr(errY, sts);% Plot with errorbars
% set(gca,'xticklabel',blockNames,'FontSize',16)
% %ylabel('Mean proportion correct');
% %legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
% ylim ([0 1]);
% title('Adults');
% grid on