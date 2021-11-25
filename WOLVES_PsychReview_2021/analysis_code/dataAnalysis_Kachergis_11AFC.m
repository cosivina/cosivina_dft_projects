%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis for Kachergis, Yu, & Shiffrin (2012)
clear all; close all; % This script generates an output file
%Global Variables

%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, Standard_Error, RMSE_val, MAPE_val);
  

%% Raw Data File Name 
simName = 'WPPR_Kachergis_Yu_Shiffrin_2012_11AFC_results'  %'batch4_11AFC_Kachergis_Yu_Shiffrin_2012_results';
% loading the file
xsit_result = load ([simName '.mat']);
numSubjects = size(xsit_result.test,1)

%Experiment (Task) Variables
new_sim_formatting = true; % new format sims with task_variables saved at simulation
if (new_sim_formatting == true)
    nObjects = xsit_result.train(1).nObjects;%%total number of objects 
    nTrials = 2*xsit_result.train(1).maxTRt; 
    AFC = xsit_result.train(1).AFC;
    LCond=3;ECond=4; %number of Late & Early Conditions
    pairtypes = 4; % w1-o1,w7-o7,w1-o7,w7-o1
else
    nObjects = 12; %%total number of objects
    nTrials=24;
    AFC = 11; 
    LCond=3;ECond=4; %number of Late & Early Conditions
    pairtypes = 4; % w1-o1,w7-o7,w1-o7,w7-o1
end


%% TEST CONDITION ANALYSIS
hFig = figure(1);set(hFig, 'Position', [0 0 1200 300]);
rmse_Within=[]; pe_Within=[]; rmse_Cross=[]; pe_Cross=[];m_0=[];
for lateC=1:LCond
    inA=0;inB=0;inC=0;inD=0; 
    corrMapped=zeros(numSubjects,pairtypes); % 
    Early0=[];Early3=[];Early6=[];Early9=[];
    for subject=1:numSubjects
        Words =  cell2mat(xsit_result.train(subject).Words);
        %Sequence=  xsit_result.train(subject).Sequence;% first 6 indices for early ,lst 6 late (1-7 pairing)
        Testset = xsit_result.test(subject).Testset; %size wl be 24 i.e 2 times for 12 objects
        %Testset, has 10 distractors, 1 target/paired-target(11), 
        %its location (12), pairtype(13), and the word (14) 
        for trt=1:nTrials
            for ikch = 1:AFC
                Look(ikch)= sum(xsit_result.test(subject).historyK(trt,ikch,:));
            end    
            targLoc = Testset(trt,12);%the spatial location of target   
            %next, add looking time to all 10 distracrors
            dstrLook=0;for i=1:AFC;if i~=targLoc; dstrLook= dstrLook + Look(i); end; end
            %%% find which pairtype w1-o1,w7-o7,w1-o7,w7-o1
            pairType = Testset(trt,13) - 64; %subtracting ascii 64 base and avoid any ( strcmp(s1,'A'))       
            if (Look(targLoc) > dstrLook) %WARNING divide by 10
                corrMapped(subject,pairType)= corrMapped(subject,pairType)+1; 
            end     
        end
        
        %pause(3);
         %% LATE TRIALS       
        if lateC==1, lateSubj=[0,4,7,10]; M_late_empirical_within=0.71;M_late_empirical_cross=0.15;
        elseif lateC==2, lateSubj=[2,5,8,11]; M_late_empirical_within=0.74;M_late_empirical_cross=0.28;
        elseif lateC==3, lateSubj=[3,6,9,1]; M_late_empirical_within=0.84;M_late_empirical_cross=0.34;
        end
        iteration=mod(xsit_result.test(subject).subject,LCond*ECond);
        %%% 3 time late pair condition is [1,4,7,10]) %% 6 latee [2,5,8,11]) %% 9 late [3,6,9,0]
        if ismember(iteration, lateSubj)%%[1,4,7,10])%%[3,6,9,0]) %%
             if  ismember(iteration, [0,5,9]) %%  early condition  0 times
                 inA=inA+1;
                 Early0(inA, :)= corrMapped(subject, :)./(nTrials/pairtypes);
             elseif ismember(iteration, [2,6,10]) %%  early condition;  3 times
                 inB=inB+1;
                 Early3(inB, :)= corrMapped(subject, :)./(nTrials/pairtypes);
             elseif ismember(iteration, [3,7,11]) %%  early condition;  6 times
                 inC=inC+1;
                 Early6(inC, :)= corrMapped(subject, :)./(nTrials/pairtypes);
             elseif  ismember(iteration, [4,8,1]) %% early condition; wherein object appeared 9 times
                 inD=inD+1;
                 Early9(inD, :)= corrMapped(subject, :)./(nTrials/pairtypes);
             end
        end
    end
    subplot(1,3,lateC);
    if lateC==1 
        ylabel('Mean proportion correct');
    elseif lateC==2 
        xlabel('Early repetitions');
    end
    %plot (mean(Early0),'-or'); plot (mean(Early3),'-og'); plot (mean(Early6),'-ob'); plot (mean(Early9),'-oc');
    E0=mean(Early0); Err0=std(Early0)./sqrt(length(Early0));
    E3=mean(Early3); Err3=std(Early3)./sqrt(length(Early3));
    E6=mean(Early6); Err6=std(Early6)./sqrt(length(Early6));
    E9=mean(Early9); Err9=std(Early9)./sqrt(length(Early9));
    hold all

    x=[3 6 9];
    y=[E3(1),E6(1),E9(1)];
    yErr=[Err3(1),Err6(1),Err9(1)];
    %plot([E0(1),E3(1),E6(1),E9(1)],':sr','LineWidth',1);plot([E0(3),E3(3),E6(3),E9(3)],'-.db','LineWidth',1);plot([E0(4),E3(4),E6(4),E9(4)],'--^k','LineWidth',1);plot([E0(2),E3(2),E6(2),E9(2)],'-vg','LineWidth',1);
    errorbar(x,y,yErr,':sr','LineWidth',1,'MarkerEdgeColor','r',...
        'MarkerFaceColor','r','MarkerSize',10);
    y=[E3(3),E6(3),E9(3)];
    yErr=[Err3(3),Err6(3),Err9(3)];
    errorbar(x,y,yErr,'-.db','LineWidth',1,'MarkerEdgeColor','b',...
        'MarkerFaceColor','b','MarkerSize',10);
    y=[E3(4),E6(4),E9(4)];
    yErr=[Err3(4),Err6(4),Err9(4)];
    errorbar(x,y,yErr,'--^c','LineWidth',1,'MarkerEdgeColor','c',...
        'MarkerFaceColor','c','MarkerSize',10);
    y=[E3(2),E6(2),E9(2)];
    yErr=[Err3(2),Err6(2),Err9(2)];
    errorbar(x,y,yErr,'-vg','LineWidth',1,'MarkerEdgeColor','g',...
        'MarkerFaceColor','g','MarkerSize',10);
    title([num2str(lateC*3),' Late'])
    legend('w1-o1','w1-o7','w7-o1','w7-o7')
    xticks([0 3 6 9])
    xticklabels({0 3 6 9})
    x=0;
    y=mean(E0); yErr= mean(Err0);
    errorbar(x,y,yErr,'om','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',10);
    xlim([-1 10]) 
    ylim([0 1])
    set(gca,'fontsize',12);
    grid on
    h = hline(1/11);
    set(h,'LineWidth',2)
    
    
    %% Data for Output File Saving 
    measurement_i = ['Mean proportion correct at test: Within-stage (' num2str(3*lateC) '-Late)'];
    empirical_mean_i = M_late_empirical_within;
    mean_i =  mean([E3(1) E6(1) E9(1) E3(2) E6(2) E9(2)]);
    SE_i   =  SE([E3(1) E6(1) E9(1) E3(2) E6(2) E9(2)]);
    RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
    row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
    xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
    xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);
    
    rmse_Within = [rmse_Within RMSE_i];    pe_Within = [pe_Within MAPE_i]; 
    
    %% Data for Output File (Save) 
    measurement_i = ['Mean proportion correct at test: Cross-stage (' num2str(3*lateC) '-Late)'];
    empirical_mean_i = M_late_empirical_cross;
    mean_i =  mean([E3(3) E6(3) E9(3) E3(4) E6(4) E9(4)]);
    SE_i   =  SE([E3(3) E6(3) E9(3) E3(4) E6(4) E9(4)]);
    RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
    row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
    xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
    xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);
    
    rmse_Cross = [rmse_Cross RMSE_i];    pe_Cross = [pe_Cross MAPE_i];
    
    m_0 = [m_0 mean(E0)];
    
    %M_late_within = mean([E3(1) E6(1) E9(1) E3(2) E6(2) E9(2)]);
    %M_late_cross =  mean([E3(3) E6(3) E9(3) E3(4) E6(4) E9(4)]);
    %rmse_Within = sqrt(((M_late_within - M_late_empirical_within)^2 )./1)
    %rmse_Cross = sqrt((M_late_cross - M_late_empirical_cross)^2)
    %pe_Within = mean ([abs(M_late_within - M_late_empirical_within)/M_late_empirical_within*100  ] )%
    %pe_Cross = mean ([abs(M_late_cross - M_late_empirical_cross)/M_late_empirical_cross*100 ] )
end
%rmse_0=sqrt((m_0/3 - 0.4766)^2)
%pe_0 = mean ([abs(m_0/3 - 0.4766)/0.4766*100 ] )

%% Data for Output File Saving 
measurement_i = ['Mean proportion correct at test: No Early-stage'];
empirical_mean_i = 0.4766;
mean_i = mean(m_0);
SE_i =   SE(m_0);
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=[measurement_i,' = ', num2str(mean_i)]; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);
rmse_0 = RMSE_i; pe_0 = MAPE_i;

%% Data for Output File Saving 
measurement_i = ['Overall Accuracy'];
empirical_mean_i = NaN; mean_i = NaN;SE_i =   NaN;
RMSE_i = mean([rmse_Within rmse_Cross rmse_0]);MAPE_i = mean([pe_Within pe_Cross pe_0]);
row_i = {measurement_i, num2str(mean_i), num2str(SE_i), RMSE_i, MAPE_i}; T = [T; row_i];
xx=measurement_i; disp(xx);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);

%disp (['overall RMSE is ', num2str(mean([rmse_Within rmse_Cross rmse_0]))])
%disp (['overall MAPE is ',num2str(mean([pe_Within pe_Cross pe_0]))])


%% write table T to output csv file
writetable(T,[simName 'Analysis.csv'])
