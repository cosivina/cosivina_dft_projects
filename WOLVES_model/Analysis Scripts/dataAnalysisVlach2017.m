%% Data Analysis 
%clear all; close all;
nTrials=9;nFeatures=2; scale_factor=1/8;  %check with auto file
plotStyle = {'k-','b-','g-','c-','r-','m-','y-','k-.','b-.','g-.','c-.','r-.','m-.','y-.','b:','w.'};
compStyle=7;
blockNames{10}=[];sts{10}=[];errY{10}=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Vlach_DeBrockWR_2017%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simName = 'base2011_Vlach_DeBrockWR_2017_';
OutName = [simName,'results.mat'];
xsit_result = load (OutName);
numSubjects=size(xsit_result.test,1);
freq_WR=zeros(1,nTrials);
perform_WR=[];
for subject=1:numSubjects
    corResp=0;
    for trt=1:nTrials
        %activity = sum(sum(xsit_result.test(subject).historyW(trt,:,:),2),3) % if 
        wordTested = xsit_result.test(subject).test_word(trt,1);
        wordStatus = char (xsit_result.test(subject).test_word(trt,2));
        activity = sum(xsit_result.test(subject).historyW(trt,:,wordTested))/(scale_factor*1000);
         if ( strcmp(wordStatus,'F')) % familiar word
             if activity > 0.5
                 corResp=corResp+1;
             end
         elseif ( strcmp(wordStatus,'N')) % Novel word
             if activity < 0.5
                 corResp=corResp+1;
             end
         end          
    end
    perform_WR(subject)= corResp;
    for trt=1:nTrials
        if trt==perform_WR(subject)
            freq_WR(trt)= freq_WR(trt)+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Vlach_DeBrockOR_2017%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all; close all;
nTrials=9;nFeatures=2; scale_factor=1/8; 
simName = 'base2011_Vlach_DeBrockOR_2017_';
OutName = [simName,'results.mat'];
xsit_result = load (OutName);
numSubjects=size(xsit_result.test,1);
freq_OR=zeros(1,nTrials);
perform_OR=[];
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
    perform_OR(subject)= corResp;
    for trt=1:nTrials
        if trt==perform_OR(subject)
            freq_OR(trt)= freq_OR(trt)+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Vlach_DeBrockWOB_2017%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all; close all;
xsit_result.train=[];
xsit_result.test=[];
nTrials=12;nFeatures=2; scale_factor=1/8;
simName = 'wfChanges_conwmf0_hwmc1_fix48_xxx_12h1200_Vlach_DeBrockWOB_2017_'%'base2011_Vlach_DeBrockWOB_2017_';
OutName = [simName,'results.mat'];
xsit_result1 = load (OutName);
xsit_result.train = [xsit_result.train ; xsit_result1.train];  
xsit_result.test = [xsit_result.test ; xsit_result1.test];
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
mean(perform_WOB)
std(perform_WOB)./sqrt(length(perform_WOB))


figure(1)
sz = 100;
%scatter(freq_OR,freq_CSL(1:9),sz, 'filled')
scatter(perform_OR,perform_CSL,sz)
grid on
hold on
plotregression(perform_OR,perform_CSL,'Regression')
xlabel('Object Recognition Memory');
ylabel('CSWL performance');


figure(2)
sz = 100;
%scatter(freq_WR,freq_CSL(1:9),sz, 'filled')
scatter(perform_WR,perform_CSL,sz)
grid on
hold on
plotregression(perform_WR,perform_CSL,'Regression')
xlabel('Word (audio) Recognition Memory');
ylabel('CSWL performance');

figure(3)
sz = 133;
%scatter(freq_WOB,freq_CSL(1:9),sz, 'filled')
scatter(perform_WOB(1:sz),perform_CSL(1:sz),sz)
grid on
hold on
plotregression(perform_WOB(1:sz),perform_CSL(1:sz),'Regression')
xlabel('Word-Object Binding Memory');
ylabel('CSWL performance');


figure(4)
x1 = [6 12];%%exp data fit
y1 = [7 10];
plot(x1,y1,'-r','LineWidth',2);
hold on


%y = 4 + ((10-7)/(12-6)) * x;

%yVD_age_fit = 4 + ((8.7-6.5)/(73-27)) * xP;

% 2,3,4,6,(7,8)
% x = [6.5106 7.2766 8.6327 10.9091 11.4931 12];%WOB
% y = [6.0044 6.1955 7.0807 8.2018 ? ?];%CSL 8.8
% yErr =[0.1112 0.1511 0.1152 0.1163 0.15 0.15];%CSL_Error

xP =   [0.8      1.2      2       3.0      3.5        5]*3;%          5.5       6        7               8];%%parameters. tau_decay.k
x =    [6.5000   7.000    8.2500  9.6915   10.490      11.3];%   11.700    11.8    11.900];%         12];%% WOB 3x
y1000= [6.1858   6.3529   6.9875  7.9139   8.55      9.2813];%   9.5000    9.7083   10.166];%     10.208];%CSL 1k test duration
yErr = [0.1245   0.1439   0.1457  0.1240   0.1344    0.0886];%   0.1033    0.1299   0.0787];%     0.0796];
x = [6 7 8 9 10 11 12]
yVD = 4 + ((10-7)/(12-6)) * x;
plot(x,yVD,'r','LineWidth',2);% vlach DeBrock data fit
hold on
% errorbar(xP,y1000,yErr,':og','LineWidth',1.6,'MarkerEdgeColor','g',...
%         'MarkerFaceColor','g','MarkerSize',6)
% hold on
errorbar(x,y1000,yErr,':ob','LineWidth',1.6,'MarkerEdgeColor','b',...
        'MarkerFaceColor','b','MarkerSize',6)
%xlim([5.5 12.5]);
ylim([4 12]);
xlabel('Word-Object Binding Memory');
ylabel('CSWL performance');
legend('Vlach & DeBrock 2017 (fit)','WOLVES Model');
set (gca, 'FontSize',14);
grid on


rmse = sqrt((sum((y1000 - [4 + ((10-7)/(12-6)) * x]).^2))./6)
pe = mean ([mean(abs(y1000 - [4 + ((10-7)/(12-6)) * x])./[4 + ((10-7)/(12-6)) * x])*100  ] )


%param =[0   12h,12h 12h,2k   12h,35k   12h,5k   12h,55h   12h,6k   12h,7k   12h,8k];% lets run for 60, 65, 70 k and see how far we need to go
x2w =  [x             xx                7.6944    8.8056   9.0722     xx];%WOB 2x [desired range = 6-12]
y =    [x             6.5936             6.5233    6.6129   xx         xx];%CSL [designed range = 7-10]

%xP =  [0.8      1.2      2       3.0      3.5        5       5.5       6        7          8];%%parameters.. DROP 8k
x =    [6.5000   7.000    8.2500  9.6915   10.490      11.3];%   11.700    11.8    11.900];%         12];%% WOB 3x
y1000= [6.1858   6.3529   6.9875  7.9139   8.55      9.2813];%   9.5000    9.7083   10.166];%     10.208];%CSL 1k test duration
yErr = [0.1245   0.1439   0.1457  0.1240   0.1344    0.0886];%   0.1033    0.1299   0.0787];%     0.0796];



%% DATA FROM KACHERGIS ON gereral_params(0.20, 0.88, c) parameter set 2 where c = [0.96; 0.78; 0.6; 0.42; 0.24; 0.06]
v_wob = [11.423920 11.087354 10.161895  9.266802  8.312669];%kachergis in wob task
v_cswl = [11.701650 11.134953  9.595189  8.774666  8.363047]% kachergis in cswl task
plot (v_wob, v_cswl)
rmse = sqrt((sum((v_cswl - [4 + ((10-7)/(12-6)) * v_wob]).^2))./6)
mape=mean ([mean(abs(v_cswl - [4 + ((10-7)/(12-6)) * v_wob])./[4 + ((10-7)/(12-6)) * v_wob])*100  ] )

%% DATA FROM KACHERGIS ON gereral_params(0.347,18.31, c) parameter set 1 where c = [0.995 onwards; 5 points each subtracted by 0.18]
v_wob = [11.671103 11.502932 10.880329  9.732049  8.808695];
v_cswl = [11.594934 10.785198  9.845292  8.704619  8.855543];
%% DATA FROM KACHERGIS ON gereral_params(0.001, 15, c) parameter set 3 where c = [1 onwards; 5 points each subtracted by 0.18]
hold on
v_wob = [6.285714 5.170167 4.450466 3.966793 3.816368]
v_cswl = [6.382753 4.988066 4.975463 5.397850 5.529094]
plot(v_wob,v_cswl,'g','LineWidth',2);

%% DATA FROM STEVEN'S PURSUIT ON remember parameter where  rm=(0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 0.95, 1.0)
hold on
v_wob = [6.524  6.285  7.220  8.071  8.998  9.914 10.810 11.686]
v_cswl = [6.579  6.295  7.250  8.139  8.981  9.829 10.802 11.710]
plot(v_wob,v_cswl,'k','LineWidth',2);

figure (5)
%param =[12h8h   12h10h  12k12h  12k15h  12k 2k  12h30h   12h,35h  12h,5k   12h,55h   12h,6k   12h,7k   12h,8k];
mas_W = [3.1858  3.233   3.2118  3.3470  3.5125  3.7990   3.8881   4.4271   4.4340     4.6424   4.8958];%    4.8472]
int_W = [3.0000  3.0631  3.1412  3.3196  3.475   4.1148   4.3659   4.8542   5.0660     5.0660   5.2708];%    5.3611]
xP =    [0.8     1.0     1.2     1.5     2       3.0      3.5       5        5.5       6        7];%          8];%%may be DROP 8k
plot (xP,mas_W);
hold on
plot (xP,int_W);
legend('Massed','Interleaved');


figure (6)
age =      [12    14    16     20     55      73.6      120      300 ];
param =    [700   800   1000   1500   3000    4000      5500     15000 ];
plot (age,param);
xlabel('Age');
ylabel('tau\_Decay');
set (gca, 'FontSize',16);
grid on

