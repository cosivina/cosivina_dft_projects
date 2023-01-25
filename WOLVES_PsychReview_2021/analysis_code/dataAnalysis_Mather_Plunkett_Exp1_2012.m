%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% Data Analysis 
clear all; close all; % This script generates an output file
%Global Variables
scale_factor=8;nFeatures=2;MIN_LOOK_DURATION=200/scale_factor; 
%Plotting Variables
plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;
%Output File Save Variables
Measure = {}; Mean = {}; Standard_Error = {}; RMSE_val = [];  MAPE_val = [];
T = table (Measure, Mean, RMSE_val, MAPE_val);
%Experiment (Task) Variables
simName = 'WPR_Mather_Plunkett_Exp1_2012_results'
xsit_result = load (simName);
nObjects = 2;%%total number of objects 
nBlocks=2;nTrials=3;nFeatures=2;repeats=nTrials/nObjects;  %check with auto file
vis_Off = floor(8000/scale_factor); %t_max = 1900/scale_factor;% scale with Experiment training trial Duration
word_On =  floor([3633 5633]/scale_factor); word_Off = floor((600+[3633 5633])/scale_factor);word_Len=floor((600)/scale_factor);
numSubjects=size(xsit_result.test,1);xx=['Number of Subjects is ',num2str(numSubjects)]; disp(xx);% 

for subject=1:numSubjects
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


figure (100);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
rectangle('Position',[word_On(1),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
rectangle('Position',[word_On(2),0,word_Len,100],'FaceColor',[.5 .5 .5 0.4],'EdgeColor',[1 1 1 0.25],'LineWidth',0.1);hold on
plot(squeeze(nanmean(squeeze(nanmean(propNovel,1)),1))*100,'LineWidth',3,'Color','b')
hold on
plot(squeeze(nanmean(squeeze(nanmean(propFamiliar,1)),1))*100,'LineWidth',3,'Color','g')
hold on
plot(squeeze(nanmean(squeeze(nanmean(propKnown,1)),1))*100,'LineWidth',3,'Color','r')
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

%mean pref to novel against pre-exposed: [pre post]
emp_prefN2P = [55.6 64.5];
%mean pref to novel against known: [pre post]
emp_prefN2K = [49.4 62.5];    %49.4% (SD = 12.3) to 62.5% (SD = 18.6),
%mean pref to known against pre-exposed: [pre post not provided]
emp_prefK2P = [55.7 ]; %55.7%, SD only pre given
empirical_mean_i = [emp_prefN2P emp_prefN2K emp_prefK2P];

sim_prefN2P_pre = squeeze(nanmean(squeeze(nanmean(propNovel(:,:,1:500),1)),1)/(squeeze(nanmean(squeeze(nanmean(propNovel(:,:,1:500),1)),1)) + squeeze(nanmean(squeeze(nanmean(propFamiliar(:,:,1:500),1)),1))))*100;
sim_prefN2P_post = squeeze(nanmean(squeeze(nanmean(propNovel(:,:,501:end),1)),1)/(squeeze(nanmean(squeeze(nanmean(propNovel(:,:,501:end),1)),1)) + squeeze(nanmean(squeeze(nanmean(propFamiliar(:,:,501:end),1)),1))))*100;

sim_prefN2K_pre = squeeze(nanmean(squeeze(nanmean(propNovel(:,:,1:500),1)),1)/(squeeze(nanmean(squeeze(nanmean(propNovel(:,:,1:500),1)),1)) + squeeze(nanmean(squeeze(nanmean(propKnown(:,:,1:500),1)),1))))*100;
sim_prefN2K_post = squeeze(nanmean(squeeze(nanmean(propNovel(:,:,501:end),1)),1)/(squeeze(nanmean(squeeze(nanmean(propNovel(:,:,501:end),1)),1)) + squeeze(nanmean(squeeze(nanmean(propKnown(:,:,501:end),1)),1))))*100;

sim_prefK2P_pre = squeeze(nanmean(squeeze(nanmean(propKnown(:,:,1:500),1)),1)/(squeeze(nanmean(squeeze(nanmean(propFamiliar(:,:,1:500),1)),1)) + squeeze(nanmean(squeeze(nanmean(propKnown(:,:,1:500),1)),1))))*100;

mean_i = [sim_prefN2P_pre sim_prefN2P_post sim_prefN2K_pre sim_prefN2K_post sim_prefK2P_pre];

measurement_i = 'proportion of subjects looking';
RMSE_i = RMSE(empirical_mean_i, mean_i);MAPE_i = MAPE(empirical_mean_i, mean_i);
xx=['RMSE = ', num2str(RMSE_i),' and ', 'MAPE = ', num2str(MAPE_i)]; disp(xx);
row_i = {measurement_i, num2str(mean_i),  RMSE_i, MAPE_i}; T = [T; row_i];

%% write table T to output csv file
writetable(T,[simName '_Analysis.csv'])


