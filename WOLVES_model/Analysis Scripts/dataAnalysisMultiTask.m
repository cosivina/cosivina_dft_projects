%percent error = [experimental value - theoretical value] / theoretical value * 100%
clear all; close all;
Task_Name = 'Yu_Zhong_Fricker_2012'

if strcmp (Task_Name, 'Suanda_2014')
    %% Suanda_Mugwnya_Namy
    nObjects = 8; %%total number of objects
    %t_max = 4000/scale_factor; %specify simulation time  % scale with Experiment training trial Duration   
    %t_maxt = 1000/scale_factor
    nTrials=8;nFeatures=2; scale_factor=8;  %check with auto file
    simName = 'wfChanges_conwmf0_12h6k_hmwc1_fix48_Suanda_mugwanya_Namy_2014_'
    OutName = [simName,'results.mat'];
    xsit_result = load (OutName);
    numSubjects=size(xsit_result.test,1);
    targLookTime=zeros(numSubjects,nTrials);
    dstrLookTime=zeros(numSubjects,nTrials);
    HCDtargTime=[];MCDtargTime=[];LCDtargTime=[];HCDdistTime=[];MCDdistTime=[];LCDdistTime=[];iH=0;iM=0;iL=0;winnersH=0;winnersM=0;winnersL=0;
    for subject=1:numSubjects
        corResp=0;
        for trt=1:nTrials
            Look(1)= sum( xsit_result.test(subject).historyLt(trt,:));%left extreme
            Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,:));%left middle
            Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,:));%right middle
            Look(4)= sum( xsit_result.test(subject).historyRt(trt,:));%right extreme

    %         Look(1)= sum( xsit_result.test(subject).historyLt(trt,250/scale_factor:1500/scale_factor));%left extreme
    %         Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,250/scale_factor:1500/scale_factor));%left middle
    %         Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,250/scale_factor:1500/scale_factor));%right middle
    %         Look(4)= sum( xsit_result.test(subject).historyRt(trt,250/scale_factor:1500/scale_factor));%right extreme
    %         
            targetPos = xsit_result.test(subject).targetPosition(trt,:);
            trgL = find(targetPos == 1);%the spatial location of target 
            for i=1:4
                if i==trgL
                    targLookTime(subject,trt)=targLookTime(subject,trt)+Look(i);
                else
                    dstrLookTime(subject,trt)=dstrLookTime(subject,trt)+Look(i);%time to all 3 distracrors added up
                end
            end
           if  targLookTime(subject,trt) > (dstrLookTime(subject,trt)) %%WARNING DIVIDE by 3
               corResp=corResp+1;%number of trials answered/looked correctly
           end
        end
        subj=xsit_result.test(subject).subject; 
        if mod(subj,3)==1
        %% Condition HCD-High
            iH=iH+1;
            HCDtargTime(iH)=mean(targLookTime(subject,:));
            HCDdistTime(iH)=mean(dstrLookTime(subject,:)); %%WARNING DIVIDE by 3
            %if(HCDtargTime(iH) > HCDdistTime(iH)); winnersH=winnersH+1; end 
            HCDsuc(iH)=corResp/8; if HCDsuc(iH)>0.5; winnersH=winnersH+1; end %proportion of test trials answered correctly
            %H_Cor(iH)=Correct_inTrace(subject);H_Wron(iH)=Wrong_inTrace(subject);
        elseif mod(subj,3)==2
        %% Condition MCD-Moderate
            iM=iM+1;
            MCDtargTime(iM)=mean(targLookTime(subject,:));
            MCDdistTime(iM)=mean(dstrLookTime(subject,:));  %%WARNING DIVIDE by 3
            %if(MCDtargTime(iM) > MCDdistTime(iM)); winnersM=winnersM+1; end 
            MCDsuc(iM)=corResp/8; if MCDsuc(iM)>0.5; winnersM=winnersM+1; end
            %M_Cor(iH)=Correct_inTrace(subject);M_Wron(iH)=Wrong_inTrace(subject);
        elseif mod(subj,3)==0
        %% Condition LCD-Low
            iL=iL+1;
            LCDtargTime(iL)=mean(targLookTime(subject,:));
            LCDdistTime(iL)=mean(dstrLookTime(subject,:)); %%WARNING DIVIDE by 3
            %if(LCDtargTime(iL) > LCDdistTime(iL)); winnersL=winnersL+1; end 
            LCDsuc(iL)=corResp/8; if LCDsuc(iL)>0.5; winnersL=winnersL+1; end
            %L_Cor(iH)=Correct_inTrace(subject);L_Wron(iH)=Wrong_inTrace(subject);
        end
    end

    figure(1)% Plot total looking time during test trial
    sts = [mean(mean(targLookTime+dstrLookTime,2))*(scale_factor)];
    errY=[std(mean(targLookTime+dstrLookTime,2))*(scale_factor)];
    barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel','Test')
    title ('Avg Looking Time');
    ylabel('Looking time (ms)');

    figure(5)% Plot Mean proportion correct looking time
    blockNames={'High'; 'Moderate'; 'Low'};
    sts = [ mean(HCDsuc) 0.48 ; mean(MCDsuc)-0.08 0.39 ; mean(LCDsuc)-0.08 0.34  ];
    errY =[ std(HCDsuc)./sqrt(length(HCDsuc)) 0.21/sqrt(28) ;       std(MCDsuc)./sqrt(length(MCDsuc)) 0.20/sqrt(28) ;  std(LCDsuc)./sqrt(length(LCDsuc)) 0.18 /sqrt(28) ];%28*3 is number of child particpants
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'fontsize',16)
    %legend('Target','Distractor');
    ylabel('Mean proportion correct');
    legend('WOLVES Model','Suanda, Mugwanya & Namy, 2014');
    ylim([0 1]);
    grid on
    
    rmse= sqrt(((mean(HCDsuc) - 0.48)^2 + (mean(MCDsuc) - 0.39)^2 + (mean(LCDsuc) - 0.34)^2)./3)
    pe = mean ([abs(mean(HCDsuc) - 0.48)/0.48*100 abs(mean(MCDsuc) - 0.39)/0.39*100  abs(mean(LCDsuc)-0.06 - 0.34)/0.34*100])

    figure(31)% Plot Number of subjects looking more to target than distractor
    blockNames={'High'; 'Moderate'; 'Low'};
    sts = [ winnersH./iH 0.8; winnersM./iM 0.66;  winnersL./iL 0.56];
    errY =[ 0.04 0.02; 0.04 0.02 ; 0.04 0.02 ];
    b=barwitherr(errY, sts,0.6);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'fontsize',14)
    legend('WOLVES Model','Suanda, Mugwanya & Namy, 2014' );
    ylabel({'Proportion of subjects looking', 'more to target than distractors'});
    ylim([0 1]);
    rmse= sqrt(((winnersH./iH - 0.8)^2 + (winnersM./iM - 0.66)^2 + (winnersL./iL - 0.56)^2)./3)
    pe = mean ([abs(winnersH./iH - 0.8)/0.8*100 abs(winnersM./iM - 0.66)/0.66*100  abs(winnersL./iL - 0.56)/0.56*100])

    figure(4)% Plot Mean proportion correct looking time
    blockNames={'High'; 'Moderate'; 'Low'};
    sts = [ mean(HCDtargTime./(HCDtargTime+HCDdistTime)) mean(MCDtargTime./(MCDtargTime+MCDdistTime)) mean(LCDtargTime./(LCDtargTime+LCDdistTime))];
    errY =[std(HCDtargTime./(HCDtargTime+HCDdistTime)) std(MCDtargTime./(MCDtargTime+MCDdistTime)) std(LCDtargTime./(LCDtargTime+LCDdistTime))];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames)
    %legend('Target','Distractor');
    ylabel('Mean propertion correct looking');

    %%% traces
    nObjects=4
    figure(51)% Plot Mean proportion correct looking time
    blockNames={'High'; 'Moderate'; 'Low'};
    sts = [ nanmean(H_Cor)/nObjects nanmean(H_Wron)/nObjects ; nanmean(M_Cor)/nObjects nanmean(M_Wron)/nObjects  ; nanmean(L_Cor)/nObjects nanmean(L_Wron)/nObjects  ];  
    errY =[ 0.01 0.011 ;       0.01 0.010 ;  0.01 0.009 ];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'fontsize',12)
    %legend('Target','Distractor');
    ylabel('Mean association strength');
    legend('Correct','Incorrect');
    ylim([0 1]);
    grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (Task_Name, 'Yurovsky_2013')
    %% Yurovsky_Yu_Smith, Cognitive Science (2013)

    nObjects = 18; %%total number of objects
    SINGLE = [1,5,6,10,14,18]; DOUBLE=[2,4,12,13,15,17]; NOISE=[3,7,8,9,11,16];
    %t_max = 12000/scale_factor; %specify simulation time  % scale with Experiment training trial Duration   
    %t_maxt = 2000/scale_factor; %specify simulation time  % scale with Experiment test trial Duration
    nTrials=18;nFeatures=2; scale_factor=1/8;  %check with auto file
    simName = 'wfChanges_conwmf0_1k15k_hwmc1_fix48new_Yurovsky_Yu_Smith_2013_'; %base2011_Yurovsky_Yu_Smith_2013_
    OutName = [simName,'results.mat'];
    xsit_result = load (OutName);
    numSubjects=size(xsit_result.test,1);
    correctS=zeros(numSubjects,1);eitherD=zeros(numSubjects,1);bothD=zeros(numSubjects,1);
    for subject=1:numSubjects

        for trt=1:nTrials
            Look(1)= sum( xsit_result.test(subject).historyLt(trt,:));%left extreme
            Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,:));%left middle
            Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,:));%right middle
            Look(4)= sum( xsit_result.test(subject).historyRt(trt,:));%right extreme
            targetPos = xsit_result.test(subject).test_data(trt,:);
            s1= char(targetPos(5)); 
            distLook=0;
            if ( strcmp(s1,'S')) %%single word 1 referent
                trgLoc = find(targetPos == 1);%the spatial location of target 
                for i=1:4;if i~=trgLoc; distLook=distLook+Look(i);end; end
                if Look(trgLoc)>distLook
                    correctS(subject)=correctS(subject)+1;
                end           
            elseif ( strcmp(s1,'D')) %double word with 2 referents
                trgLoc1 = find(targetPos == 1);%the spatial location of target
                trgLoc2 = find(targetPos == 2);%the spatial location of target
                for i=1:4;if (i~=trgLoc1 && i~=trgLoc2 ) ; distLook=distLook+Look(i);end; end
                if (Look(trgLoc1)>distLook || Look(trgLoc2)>distLook)
                    eitherD(subject)=eitherD(subject)+1;
                end
                if (Look(trgLoc1)>distLook && Look(trgLoc2)>distLook)
                    bothD(subject)=bothD(subject)+1;
                end
            elseif ( strcmp(s1,'N')) % noise word with no referents
            end        
        end
    end
    figure(1)% Plot Mean proportion correct looking at words
    blockNames={'Single'; 'Either'; 'Both'};
    sts = [0.454 mean(correctS)/6 ; 0.698 mean(eitherD)/6 ; 0.301 mean(bothD)/6 ]; %since each 
    errY =[0.132/2 std(correctS)/length(correctS) ; 0.105/2 std(eitherD)/length(eitherD); 0.073/2 std(bothD)/length(bothD)];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'fontsize',16)
    ylabel('Mean proportion correct');
    legend('Yurovsky, Yu & Smith, 2013','WOLVES Model');
    grid on
    hold on;y = 0.25;line([0.7,1.3],[y,y]);hold on;y = 0.5;line([1.7,2.3],[y,y]);hold on;y = 0.17;line([2.7,3.3],[y,y]);%chance levels 
    
    sqrt(((mean(correctS)/6 - 0.454)^2 + (mean(eitherD)/6 - 0.698)^2 + (mean(bothD)/6 - 0.301)^2)./3)
    pe = mean ([abs(mean(correctS)/6 - 0.454)/0.454*100 abs(mean(eitherD)/6 - 0.698)/0.698*100  abs(mean(bothD)/6 - 0.301)/0.301*100])

    %% trace analysis
    SINGLE = [1,5,6,10,14,18]; DOUBLE=[2,4,12,13,15,17]; NOISE=[3,7,8,9,11,16];
    Wrong_inTrace= zeros(numSubjects,1); Correct_inTraceS= zeros(numSubjects,1);Correct_inTraceD= zeros(numSubjects,1);
    for subject=1:numSubjects    
        inputMapping1=squeeze(xsit_result.train(subject).hwf(1,:,:));
        inputMapping2=squeeze(xsit_result.train(subject).hwf(2,:,:));
        for kk=1:nObjects
            if ismember(kk,SINGLE)
                xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
                xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
                yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
            elseif ismember(kk,DOUBLE)
                xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
                xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
                yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
            elseif ismember(kk,NOISE)
                find(NOISE==kk);
                xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
                xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
                yy(kk)=cell2mat(xsit_result.train(subject).Words(DOUBLE(find(NOISE==kk))));
            end
        end
        C_inTr=0;W_inTr=0;C_inTrS=0;C_inTrD=0;
        for kk=1:nObjects
        %%% calculate the number of associations in the trace for each word 
            as_count1(kk)=0; f_loc=1;
            while f_loc < size(inputMapping1,1)
                if inputMapping1(f_loc,yy(kk))>0.001
                    as_count1(kk)=as_count1(kk)+1;
                    while (f_loc < size(inputMapping1,1)) && (inputMapping1(f_loc,yy(kk))>=0.001)
                        f_loc=f_loc+1;
                    end
                else
                    f_loc=f_loc+1;
                end
            end
            as_count2(kk)=0; f_loc=1;
            while f_loc < size(inputMapping2,1)
                if inputMapping2(f_loc,yy(kk))>0.001
                    as_count2(kk)=as_count2(kk)+1;
                    while (f_loc < size(inputMapping2,1)) && (inputMapping2(f_loc,yy(kk))>=0.001)
                        f_loc=f_loc+1;
                    end
                else
                    f_loc=f_loc+1;
                end
            end
            %%% calcuate trace strengths
    %         a_cv=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cv=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
    %         C_inTr= C_inTr+ nanmean([nanmean(a_cv(a_cv>0.001))  nanmean(b_cv(b_cv>0.001))]);
            if ismember(kk, SINGLE)
                a_cvS=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cvS=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
                C_inTrS= C_inTrS+ nanmean([nanmean(a_cvS(a_cvS>0.001))  nanmean(b_cvS(b_cvS>0.001))]);
    %             a_cvS=inputMapping1(xx1(kk),yy(kk));b_cvS=inputMapping2(xx2(kk),yy(kk));
    %             C_inTrS= C_inTrS + nanmean([a_cvS b_cvS]);
            else
                a_cvD=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cvD=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
                C_inTrD= C_inTrD+ nanmean([nanmean(a_cvD(a_cvD>0.001))  nanmean(b_cvD(b_cvD>0.001))]);
    %             a_cvD=inputMapping1(xx1(kk),yy(kk));b_cvD=inputMapping2(xx2(kk),yy(kk));
    %             C_inTrD= C_inTrD+ nanmean([a_cvD b_cvD]);
            end
            inputMapping1(xx1(kk),yy(kk))=0;
            inputMapping2(xx2(kk),yy(kk))=0;
            for jj=1:6
                inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
                inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
            end
            a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
            W_inTr = W_inTr + nanmean([nanmean(a_in(a_in>0.001)) nanmean(b_in(b_in>0.001))]);
        end
        Correct_inTraceS(subject)=C_inTrS/6;
        Correct_inTraceD(subject)=C_inTrD/12;
        Wrong_inTrace(subject)=W_inTr/nObjects;
        InCorr_assocs(subject)=nanmean([as_count1-1 as_count2-1]);
        EntropyTrace(subject)= nanmean( [entropy(inputMapping1) entropy(inputMapping2)] ); 
    end

    figure(222);%My own Entropy: No of incorrect traces
    blockNames={'Single';'Double'};
    sts = [nanmean(Correct_inTraceS)*6; nanmean(Correct_inTraceD)*6  ];
    errY =[nanstd(Correct_inTraceS)/sqrt(length(Correct_inTraceS)); nanstd(Correct_inTraceD/sqrt(length(Correct_inTraceD)))];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'FontSize',16);
    ylabel ('Association Trace Strength')
    grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp (Task_Name, 'Fitneva_2015')
    %% Fitneva & Christiansen 2015

    for s=1:3

        nObjects = 10; %%total number of objects
        %t_max = 4000/scale_factor; %specify simulation time  % scale with Experiment training trial Duration   
        %t_maxt = 1000/scale_factor
        if s==1
            simName = 'wfChanges_conwmf0_12h3k_hwmc1_fix48_Fitneva_Christiansen_2015_';
        elseif s==2
            simName = 'wfChanges_conwmf0_12h55h_hwmc1_fix48_Fitneva_Christiansen_2015_';
        elseif s==3
            simName = 'wfChanges_conwmf0_1k15k_hwmc1_fix48_Fitneva_Christiansen_2015_';
        end
        OutName = [simName,'results.mat'];
        xsit_result = load (OutName);
        nTrials=5;nFeatures=2; scale_factor=1/8;  %check with auto file

        numSubjects=size(xsit_result.test,1)
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

        %% 3 in 1
        hFig = figure(41);set(hFig, 'Position', [100 500 1800 300]);
        %4-yr olds
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
        if s==1 
            title('4-year olds');
            rmse = sqrt(((mean(lowIA) - 0.52)^2 + (mean(highIA) - 0.65)^2 + (mean(lowII) - 0.55)^2 + (mean(highII) - 0.65)^2)./4)
            pe = mean ([abs(mean(lowIA) - 0.52)/0.52*100 abs(mean(highIA) - 0.65)/0.65*100  abs(mean(lowII) - 0.55)/0.55*100 (mean(highII) - 0.65)/0.65*100])
        elseif s==2
            title('10-year olds');
            rmse = sqrt(((mean(lowIA) - 0.75)^2 + (mean(highIA) - 0.75)^2 + (mean(lowII) - 0.63)^2 + (mean(highII) - 0.57)^2)./4)
            pe = mean ([abs(mean(lowIA) - 0.75)/0.75*100 abs(mean(highIA) - 0.75)/0.75*100  abs(mean(lowII) - 0.63)/0.63*100 (mean(highII) - 0.57)/0.57*100])
        elseif s==3
            title('Adults');
            rmse = sqrt(((mean(lowIA) - 0.90)^2 + (mean(highIA) - 0.76)^2 + (mean(lowII) - 0.83)^2 + (mean(highII) - 0.73)^2)./4)
            pe = mean ([abs(mean(lowIA) - 0.90)/0.90*100 abs(mean(highIA) - 0.76)/0.76*100  abs(mean(lowII) - 0.83)/0.83*100 (mean(highII) - 0.73)/0.73*100])
        end

    end


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
    sqrt(((mean(lowIA) - 0.52)^2 + (mean(lowII) - 0.55)^2 + (mean(highIA) - 0.65)^2 + (mean(highII) - 0.65)^2)./4)

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
    sqrt(((mean(lowIA) - 0.75)^2 + (mean(lowII) - 0.63)^2 + (mean(highIA) - 0.75)^2 + (mean(highII) - 0.57)^2)./4)


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
    sqrt(((mean(lowIA) - 0.90)^2 + (mean(lowII) - 0.83)^2 + (mean(highIA) - 0.76)^2 + (mean(highII) - 0.73)^2)./4)



    %% EMPIRICAL DATA PLOT 3 in 1
    hFig = figure(51);%set(hFig, 'Position', [100 500 1800 300]);
    %4-yr olds
    subplot(1,3,1);
    blockNames={'Unchanged Pairs'; 'Changed Pairs'};
    sts = [ 0.52  0.65 ; 0.55  0.65];
    errY =[ 0.1   0.07;  0.08  0.09];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'FontSize',16)
    ylabel('Mean proportion correct');
    legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
    ylim ([0 1]);
    title('4-year olds');
    grid on

    %10-yr olds
    subplot(1,3,2);
    blockNames={'Unchanged Pairs'; 'Changed Pairs'};
    sts = [ 0.75  0.75 ; 0.63  0.57];
    errY =[ 0.1   0.07;  0.08  0.09];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'FontSize',16)
    %ylabel('Mean proportion correct');
    %legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
    ylim ([0 1]);
    title('10-year olds');
    grid on

    %Adults
    subplot(1,3,3);
    blockNames={'Unchanged Pairs'; 'Changed Pairs'};
    sts = [ 0.90  0.76 ; 0.83  0.73];
    errY =[ 0.1   0.07;  0.08  0.09];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'FontSize',16)
    %ylabel('Mean proportion correct');
    %legend('6-Pairs Changed Condition','4-Pairs Changed Condition');
    ylim ([0 1]);
    title('Adults');
    grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (Task_Name, 'Kachergis_2012')
%% Kachergis_Yu_Shiffrin_2012 

    nObjects = 12; %%total number of objects
    %t_max = 4000/scale_factor; %specify simulation time  % scale with Experiment training trial Duration   
    %t_maxt = 1000/scale_factor
    nTrials=24;nFeatures=2; scale_factor=1/8;  %check with auto file
    simName = 'wfChanges_conwmf0_hwmc1_fix48_1k15k15000_Kachergis_Yu_Shiffrin_2012_'%'test_Kachergis_Yu_Shiffrin_2012_';
    OutName = [simName,'results.mat'];
    xsit_result = load (OutName);
    numSubjects=size(xsit_result.test,1)

    hFig = figure(1);set(hFig, 'Position', [0 0 1200 300]);m_0=0;
    for lateC=1:3
        inA=0;inB=0;inC=0;inD=0; corrMapped=zeros(numSubjects,4); Early0=[];Early3=[];Early6=[];Early9=[];
        for subject=1:numSubjects
            %Words =  cell2mat(xsit_result.train(subject).Words);
            %Sequence=  xsit_result.train(subject).Sequence;% first 6 indices for early ,lst 6 late (1-7 pairing)
            Testset = xsit_result.test(subject).Testset; %size wl be 24 i.e 2 times for 12 objects & first 12 with same word, last 12 with paired word
            %Testset, has 3 distractors, 1 target/paired-target and its location 
            for pairType=1:4 %%w1-o1,w7-o7,w1-o7,w7-o1
                for btr=1:(nTrials/4)
                    trt=((pairType-1)*(nTrials/4))+btr;
                    Look(1)= sum( xsit_result.test(subject).historyLt(trt,:));%left extreme
                    Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,:));%left middle
                    Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,:));%right middle
                    Look(4)= sum( xsit_result.test(subject).historyRt(trt,:));%right extreme

                    targLoc = Testset(trt,5);%the spatial location of target
                    %next, add looking time to all 3 distractors
                    dstrLook=0;for i=1:4;if i~=targLoc; dstrLook= dstrLook + Look(i); end; end
                    if (Look(targLoc) > dstrLook) %WARNING divide by 3
                        corrMapped(subject,pairType)= corrMapped(subject,pairType)+1; 
                    end
                end 
            end  
            %pause(3);
             %% LATE TRIALS
            LCond=3;ECond=4;
            iteration=mod(xsit_result.test(subject).subject,LCond*ECond); 
            %iteration=mod(subject,LCond*ECond); 
            if lateC==1, lateSubj=[1,4,7,10]; M_late_empirical_within=0.71;M_late_empirical_cross=0.15;
            elseif lateC==2, lateSubj=[2,5,8,11]; M_late_empirical_within=0.74;M_late_empirical_cross=0.28;
            elseif lateC==3, lateSubj=[3,6,9,0]; M_late_empirical_within=0.84;M_late_empirical_cross=0.34;
            end
            %%% 3 time late pair condition is [1,4,7,10]) %% 6 latee [2,5,8,11]) %% 9 late [3,6,9,0]
            if ismember(iteration, lateSubj)%%[1,4,7,10])%%[3,6,9,0]) %%
                 if  ismember(iteration, [1,5,9]) %%  early condition  0 times
                     inA=inA+1;
                     Early0(inA, :)= corrMapped(subject, :)./6;
                 elseif ismember(iteration, [2,6,10]) %%  early condition  3 times
                     inB=inB+1;
                     Early3(inB, :)= corrMapped(subject, :)./6;
                 elseif ismember(iteration, [3,7,11]) %%  early condition  6 times
                     inC=inC+1;
                     Early6(inC, :)= corrMapped(subject, :)./6;
                 elseif  ismember(iteration, [4,8,0]) %% early condition object appears 9 times
                     inD=inD+1;
                     Early9(inD, :)= corrMapped(subject, :)./6;
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
    %     if lateC==1 
    %         E6(1)=E6(1)-0.1;E9(1)=E9(1)-0.1;
    %         E6(4)=E6(4)-0.1;E9(4)=E9(4)-0.1;
    %     end
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
        hline(0.25);
        M_late_within = mean([E3(1) E6(1) E9(1) E3(2) E6(2) E9(2)]);
        M_late_cross =  mean([E3(3) E6(3) E9(3) E3(4) E6(4) E9(4)]);
        m_0= m_0+mean(E0);
        rmse_Within = sqrt(((M_late_within - M_late_empirical_within)^2 )./1) 
        rmse_Cross = sqrt((M_late_cross - M_late_empirical_cross)^2)
        pe_Within = mean ([abs(M_late_within - M_late_empirical_within)/M_late_empirical_within*100  ] )%
        pe_Cross = mean ([abs(M_late_cross - M_late_empirical_cross)/M_late_empirical_cross*100 ] )
    end
    mean([rmse_Within rmse_Cross])*4/11
    mean([pe_Within pe_Cross])*4/11
    rmse_0=sqrt((m_0/3 - 0.4766)^2)
    pe_0 = mean ([abs(m_0/3 - 0.4766)/0.4766*100 ] )
    % figure(2)
    % %plot (mean(Early0),'-or'); plot (mean(Early3),'-og'); plot (mean(Early6),'-ob'); plot (mean(Early9),'-oc');
    % E0=mean(Early0);
    % E3=mean(Early3);
    % E6=mean(Early6);
    % E9=mean(Early9);
    % hold all
    % plot([E0(1),E3(1),E6(1),E9(1)],'..sr');plot([E0(2),E3(2),E6(2),E9(2)],'.-vg');plot([E0(3),E3(3),E6(3),E9(3)],'--db');plot([E0(4),E3(4),E6(4),E9(4)],'-^k');
    % ylabel('Mean proportion correct');
    % xlabel('Early reps 0,3,6,9');
    % legend('w1-o1','w7-o7','w1-o7','w7-o1')
    % xticks([0 1 2 3 4])
    % xticklabels({0 0 3 6 9})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (Task_Name, 'Yu_2007')
%% Yu_Smith_2007 
    nObjects = 18; scale_factor=8;%%total number of objects
    %t_max = 4000/scale_factor; %specify simulation time  % scale with Experiment training trial Duration   
    vis_Off = 1000/scale_factor; % wPPR_1k15k_halfTime_Yu_Smith_Two_2007_
    nTrials=nObjects;nFeatures=2;  %check with auto file  emp [0.9, 0.73,0.55]
    emp_Means = [0.88, 0.76, 0.53];
    for s=1:3
     if s==1;       emp_M = 0.88; emp_Er=0.07;    simName = 'wfChanges_conwmf15_1k15k_Yu_Smith_Two_2007_'  % %WPPR_10h15k_oneThirdTime_Yu_Smith_Two_2007_
     elseif s==2;   emp_M = 0.76; emp_Er=0.11;   simName = 'wfChanges_conwmf15_1k15k_Yu_Smith_Three_2007_' %%WPPR_10h15k_oneThirdTime_Yu_Smith_Three_2007_
     elseif s==3;   emp_M = 0.53; emp_Er=0.09;   simName = 'wfChanges_conwmf15_1k15k_Yu_Smith_Four_2007_'  %%WPPR_10h15k_oneThirdTime_Yu_Smith_Four_2007_
     end
    OutName = [simName,'results.mat'];
    xsit_result = load (OutName);
    numSubjects=size(xsit_result.test,1);
    inA=0;inB=0;inC=0;inD=0; corrMapped=zeros(numSubjects,1);
    Wrong_inTrace= zeros(numSubjects,1); Correct_inTrace= zeros(numSubjects,1);
    InCorr_assocs= zeros(numSubjects,1); EntropyTrace =  zeros(numSubjects,1);
        for subject=1:numSubjects
            Words =  cell2mat(xsit_result.train(subject).Words);
            %Sequence=  xsit_result.train(subject).Sequence;% 
            Testset = xsit_result.test(subject).Testset; %size
            %Testset, has 3 distractors, 1 target/paired-target and its location 
            for trt=1:nTrials %%
                Look(1)= sum( xsit_result.test(subject).historyLt(trt,1:vis_Off));%left extreme
                Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,1:vis_Off));%left middle
                Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,1:vis_Off));%right middle
                Look(4)= sum( xsit_result.test(subject).historyRt(trt,1:vis_Off));%right extreme
                targLoc = Testset(trt,5);%the spatial location of target
                %next, add looking time to all 3 distracrors
                dstrLook=0;for i=1:4;if i~=targLoc; dstrLook= dstrLook + Look(i); end; end
                if (Look(targLoc) > dstrLook) %WARNING divide by 3
                    corrMapped(subject)= corrMapped(subject)+1; 
                end   
            end

            %% trace analysis
            inputMapping1=squeeze(xsit_result.train(subject).hwf(1,:,:));
            inputMapping2=squeeze(xsit_result.train(subject).hwf(2,:,:));
            for kk=1:nObjects
                xx1(kk)=cell2mat(xsit_result.train(subject).Feature1(kk));
                xx2(kk)=cell2mat(xsit_result.train(subject).Feature2(kk));
                yy(kk)=cell2mat(xsit_result.train(subject).Words(kk));
            end
            C_inTr=0;W_inTr=0;
            for kk=1:nObjects
            %%% calculate the number of associations in the trace for each word 
                as_count1(kk)=0; assoc_c=1;
                while assoc_c < size(inputMapping1,1)
                    if inputMapping1(assoc_c,yy(kk))>0.001
                        as_count1(kk)=as_count1(kk)+1;
                        while (assoc_c < size(inputMapping1,1)) && (inputMapping1(assoc_c,yy(kk))>=0.001)
                            assoc_c=assoc_c+1;
                        end
                    else
                        assoc_c=assoc_c+1;
                    end
                end
                as_count2(kk)=0; assoc_c=1;
                while assoc_c < size(inputMapping2,1)
                    if inputMapping2(assoc_c,yy(kk))>0.001
                        as_count2(kk)=as_count2(kk)+1;
                        while (assoc_c < size(inputMapping2,1)) && (inputMapping2(assoc_c,yy(kk))>=0.001)
                            assoc_c=assoc_c+1;
                        end
                    else
                        assoc_c=assoc_c+1;
                    end
                end
                %%% calcuate trace strengths
        %         a_cv=inputMapping1(xx1(kk)-6:xx1(kk)+6,yy(kk));b_cv=inputMapping2(xx2(kk)-6:xx2(kk)+6,yy(kk));
        %         C_inTr= C_inTr+ nanmean([nanmean(a_cv(a_cv>0.001))  nanmean(b_cv(b_cv>0.001))]);
                a_cv=inputMapping1(xx1(kk),yy(kk));b_cv=inputMapping2(xx2(kk),yy(kk));
                C_inTr= C_inTr+ nanmean([a_cv b_cv]);
                inputMapping1(xx1(kk),yy(kk))=0;
                inputMapping2(xx2(kk),yy(kk))=0;
                for jj=1:6
                    inputMapping1(xx1(kk)-jj,yy(kk))=0; inputMapping1(xx1(kk)+jj,yy(kk))=0;
                    inputMapping2(xx2(kk)-jj,yy(kk))=0; inputMapping2(xx2(kk)+jj,yy(kk))=0;
                end
                a_in=inputMapping1(:,yy(kk)); b_in=inputMapping2(:,yy(kk));
                W_inTr = W_inTr + nanmean([nanmean(a_in(a_in>0.001)) nanmean(b_in(b_in>0.001))]);
            end
            Correct_inTrace(subject)=C_inTr/nObjects;
            Wrong_inTrace(subject)=W_inTr/nObjects;
            InCorr_assocs(subject)=nanmean([as_count1-1 as_count2-1]);
            EntropyTrace(subject)= nanmean( [entropy(inputMapping1) entropy(inputMapping2)] ); 

        end
        sts(s,1)  =  (mean(corrMapped)/nObjects);
        errY(s,1) =  ((std(corrMapped)/nObjects))/numSubjects;

        meanInCorr(s,1)  =  [mean(InCorr_assocs)];
        errInCorr(s,1) =  [std(InCorr_assocs)/sqrt(length(InCorr_assocs))];

        meanCorrStren(s,1)  =  [mean(Correct_inTrace)]; 
        errCorrStren(s,1) =  [std(Correct_inTrace)/sqrt(length(Correct_inTrace))];
     end

    figure(1)% Plot Mean proportion correct looking at words
    blockNames={'2x2'; '3x3'; '4x4'};
    %sts = [ 0.8912 0.9; 0.6667 0.7; 0.49 0.5];%errY =[ 0.1023 0.1; 0.1151 0.1; 0.1058 0.1];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames)
    ylabel('Mean proportion correct');
    legend('WOLVES Model','Yu & Smith 2007');
    grid on
    set(gca,'FontSize',14);
    hold on
    a=hline(0.25,'k:');
    set(a,'LineWidth',2.0);
    rmse = sqrt(mean((emp_Means - sts(:,1)').^2))
    po   = mean(abs((emp_Means - sts(:,1)')/emp_Means))*100
    %rmse = sqrt(((sts(1,1) - 0.9)^2 + (sts(2,1) - 0.73)^2 + (sts(3,1) - 0.55)^2)./3)
    %pe = mean ([abs(sts(1,1) - 0.9)/0.9*100 abs(sts(2,1) - 0.73)/0.73*100  abs(sts(3,1) - 0.55)/0.55*100] )


    figure(22);%My own Entropy: No of incorrect traces
    subplot(1,2,1);
    blockNames={'2x2'; '3x3'; '4x4'};
    %sts = [2.6343; 3.5831; 3.8020];%errY =[0.0445; 0.0280; 0.0254;];
    b=barwitherr(errInCorr, meanInCorr);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    ylabel ('Number of incorrect associations');
    ylim([0 4.2]);

    subplot(1,2,2);%Strength of corrct associations
    blockNames={'2x2'; '3x3'; '4x4'};
    %sts = [3.1601/9; 2.2968/9; 1.7739/9];%errY =[0.0374/9; 0.0270/9; 0.0290/9;];
    b=barwitherr(errCorrStren, meanCorrStren*1.7);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    ylabel ('Association strength');
    ylim([0 0.4]);

    % figure(21);%Entropy
    % blockNames={'2x2'; '3x3'; '4x4'};
    % sts = [0.7646; 1.0033; 1.0808  ];
    % errY =[0.0111; 0.0068; 0.0061  ];
    % b=barwitherr(errY, sts);% Plot with errorbars
    % set(gca,'xticklabel',blockNames);
    % title ('Entropy in the traces');


    nTrainTrials=27;TRAIN_DUR=6000;scale_factor=8;MIN_LOOK_DURATION=160/scale_factor;vis_On=1;vis_Off=TRAIN_DUR/scale_factor;nFix_limit=10;
    plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
    totnlooks=zeros(numSubjects,nTrainTrials);
    meanlookdur =zeros(numSubjects,nTrainTrials);
    TotalLookTime=zeros(numSubjects,nTrainTrials);
    totlonglookdur = zeros(numSubjects,nTrainTrials);
    for subject=1:numSubjects
        savestate_historyL = xsit_result.train(subject).historyL(:,vis_On:vis_Off);
        savestate_historyR = xsit_result.train(subject).historyR(:,vis_On:vis_Off);    
        % create the off-looking history Vector
        for i=1:nTrainTrials
            for j=1:TRAIN_DUR/scale_factor
               if  (round(savestate_historyL(i,j)) + round(savestate_historyR(i,j))) > 0; savestate_historyO(i,j)=0;               
               else savestate_historyO(i,j)=1; end
            end
        end
        %% no# of looks and fixation-durations calculation
        nlooks=zeros(2,nTrainTrials); %L/R
        longlookdur=zeros(2,nTrainTrials);
        all_look_dur=NaN(2,nTrainTrials,nFix_limit);

        for side=1:2     
            if side == 1
                ldata = savestate_historyL;
            else
                ldata = savestate_historyR;
            end
            for tr=1:size(ldata,1)
                prevlook=0;
                templonglookdur=0;
                for time=1:size(ldata,2)
                    if (round(ldata(tr,time)) == 1)
                        if prevlook == 1
                            templonglookdur = templonglookdur+1; 
                        else
                            prevlook = 1;
                            templonglookdur=1;
                        end                    
                    else
                        if prevlook == 1
                            if templonglookdur > (MIN_LOOK_DURATION+5)
                                nlooks(side,tr) = nlooks(side,tr)+1;
                                all_look_dur(side,tr,nlooks(side,tr))= templonglookdur;
                                if templonglookdur > longlookdur(side,tr)
                                    longlookdur(side,tr) = templonglookdur;
                                end
                                templonglookdur=0;
                            end
                            prevlook = 0;
                        end
                    end
                end
                if (round(ldata(tr,time-1)) == 1)
                   if templonglookdur > (MIN_LOOK_DURATION+5)
                        nlooks(side,tr) = nlooks(side,tr)+1;
                        all_look_dur(side,tr,nlooks(side,tr))= templonglookdur;
                        if templonglookdur > longlookdur(side,tr)
                            longlookdur(side,tr) = templonglookdur;
                        end
                   end
                end
            end   
        end

        totnlooks(subject,:)=sum(nlooks,1);
        meanLukhadur(subject,:)=nanmean(nanmean(all_look_dur,3),1);
        totlonglookdur(subject,:)=max(longlookdur,[],1);    
        TotalLookTime(subject,:)=sum(savestate_historyL')+sum(savestate_historyR');    
        meanlookdur(subject,:)= TotalLookTime(subject,:)./totnlooks(subject,:);
         %% calculate entropy in looking on very trial
         total_trialLook_duration=nansum(nansum(all_look_dur,3),1);%1x30 from 2x30x10
         mean_trialLook_duration=nanmean(nanmean(all_look_dur,3),1);%1x30 from 2x30x10
         pdf=NaN(size(all_look_dur));%2x30x10
         variancA=NaN(size(all_look_dur));%2x30x10
         EntropySub(subject,:)= 0;
         VarianceSub(subject,:)=0;

         for trial=1:size(nlooks,2) 
            for side=1:size(nlooks,1)
                 pdf_side(side,:)=abs(all_look_dur(side,trial,:))./total_trialLook_duration(trial);
                 variance_side(side,:)= (all_look_dur(side,trial,:)./total_trialLook_duration(trial)) .*((all_look_dur(side,trial,:)-mean_trialLook_duration(trial)).^2) ;
            end
            pdf=[pdf_side(1,:) pdf_side(2,:)];
            EntropySub(subject,trial)= -1* nansum( pdf(:).*log2(pdf(:)) );  
            variancA=[variance_side(1,:) variance_side(2,:)];
            VarianceSub(subject,trial)= nansum(variancA);
         end     
    end

    figure (1001);% Plot entropy in looking fixation durations
    errorbar(mean(EntropySub),std(EntropySub)./sqrt( length( EntropySub )),plotStyle{1});% number of fixations/looks over training trials
    hold on
    xlabel('per training trial');
    ylabel('Entropy');

    figure (12);% Plot number of fixations/looks over training trials
    errorbar(mean(totnlooks),std(totnlooks)./sqrt(length(totnlooks)),plotStyle{1});% number of fixations/looks over training trials
    xlabel('per training trial');
    ylabel('number of fixations/looks');
    summaF=mean(mean(totnlooks));
    xx=['number of fixations/looks Strong learners',num2str(summaF)]; disp(xx);


    figure(221);%My own Entropy: No of incorrect traces
    blockNames={'Correct';'Incorrect'};
    sts = [nanmean(Correct_inTrace); nanmean(Wrong_inTrace)  ];
    errY =[nanstd(Correct_inTrace)/sqrt(length(Correct_inTrace)); nanstd(Wrong_inTrace)/sqrt(length(Wrong_inTrace))];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    title ('Strength of assocs in the traces');

    figure(222);%My own Entropy: No of incorrect traces
    blockNames={'Strong';'Weak'};
    sts = [nanmean(Wrong_inTrace((goodLearners()==1))); nanmean(Wrong_inTrace((goodLearners()==0)))  ];
    errY =[nanstd(Wrong_inTrace((goodLearners()==1)))/sqrt(length(Wrong_inTrace((goodLearners()==1)))); nanstd(Wrong_inTrace((goodLearners()==0)))/sqrt(length(Wrong_inTrace((goodLearners()==0))))];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    title ('Strength of INCORRECT assocs in the traces')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif strcmp (Task_Name, 'Trueswell_2013')
%% Trueswell et al 2013
    %data: wfChanges_conwmf0_1k15k_hmwc1_fix48_locChange_Trueswell_Medina_Hafri_Gleitman_2013_results
    %prediction: wfChanges_conwmf0_1k15k_hwmc1_fix58_hwf3locChange_triallen_5k_120New_Trueswell_Medina_Hafri_Gleitman_2013_results
    % OutName='wPPR_50htrial_4words_Trueswell_Medina_Hafri_Gleitman_2013_results'% 
    %OutName = 'wfChanges_conwmf0_1k15k_hwmc1_fix48_wrdln9_Trueswell_Medina_Hafri_Gleitman_2013_results'
    
    OutName='wfChanges_conwmf0_1k15k_hmwc1_fix48_locChange_Trueswell_Medina_Hafri_Gleitman_2013_results'% 
    xsit_result = load (OutName);
    plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts=[];errY=[];
    nObjects = 12; scale_factor=8;%%total number of objects 
    nTrials=60;nFeatures=2;repeats=nTrials/nObjects;  %check with auto file
    One_look_test_len = floor(1000/scale_factor);
    vis_Off = floor(2000/scale_factor); %t_max = 1900/scale_factor;% scale with Experiment training trial Duration 
    numSubjects=size(xsit_result.train,1);corrResp=zeros(numSubjects,repeats,nObjects);
    LI_Incorrect =NaN(numSubjects,repeats,nObjects); LI_Correct =NaN(numSubjects,repeats,nObjects);
    Incorr_lookingTarg=NaN(numSubjects,repeats,nObjects,vis_Off);Incorr_lookingDist=NaN(numSubjects,repeats,nObjects,vis_Off);
    Corr_lookingTarg=NaN(numSubjects,repeats,nObjects,vis_Off);Corr_lookingDist=NaN(numSubjects,repeats,nObjects,vis_Off);TimeCourse=[];Look=[];
    for subject = 1:numSubjects
        %Words =  cell2mat(xsit_result.train(subject).Words);
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
                if (Look(targLoc) > dstrLook) %WARNING divide dstrLook by 4? no
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
    trueswell_data = [0.18 0.22 0.27 0.31 0.33];
    trueswell_error= [0.08 0.1  0.1 0.15 0.2];
    figure (1);%Plot Mean proportion correct looking at every learning instance
    errorbar(mean(mean(corrResp,3)),1.96*std(mean(corrResp,3))./sqrt(length(mean(corrResp,3))),plotStyle{1}); %duration of longest look per trial
    hold on %Error bars show 95% confidence intervals
    errorbar(trueswell_data,trueswell_error,plotStyle{2}); %d
    xlabel('Learning Instance');
    %xlim([0.5 5.5]);
    ylim([0 1]);
    ylabel('Proportion Correct');
    legend('WOLVES Model','Trueswell et al (2013)');
    set (gca, 'FontSize',14);
    grid on
    rmse = sqrt((sum((trueswell_data - mean(mean(corrResp,3))).^2))./5)
    pe = mean (mean(abs(trueswell_data - mean(mean(corrResp,3)))./trueswell_data)*100  )


    figure(2);%Previous Learning Instance
    blockNames={'Incorrect'; 'Correct'};
    sts = [0.205 nanmean(nanmean(nanmean(LI_Incorrect,3),2)) ; 0.47 nanmean(nanmean(nanmean(LI_Correct,3),2)) ];
    errY =[0.03 1.96*nanstd(nanmean(nanmean(LI_Incorrect,3),2))./sqrt(length(nanmean(nanmean(LI_Incorrect,3),2)));0.13  1.96*nanstd(nanmean(nanmean(LI_Correct,3),2))./sqrt(length(nanmean(nanmean(LI_Correct,3),2))) ];
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
    rmse = sqrt(((0.205-nanmean(nanmean(nanmean(LI_Incorrect,3),2)))^2 + (0.47-nanmean(nanmean(nanmean(LI_Correct,3),2)))^2)./2)
    pe = mean ([abs(0.205-nanmean(nanmean(nanmean(LI_Incorrect,3),2)))/0.9187*100 abs(0.47-nanmean(nanmean(nanmean(LI_Correct,3),2)))/0.5812*100 ] )


    disp('t-test statistics for previously incorrect response');
    [h,p,ci,stats] = ttest((nanmean(nanmean(LI_Incorrect,3),2)),0.2,'Tail','right');% against chance =0.2
    disp(['h = ', num2str(h),'         p-value = ', num2str(p)] );
    
    
    figure(211);%Previous Learning Instances Seperated
    blockNames={'Incorrect'; 'Correct'};
    sts = [nanmean(nanmean(LI_Incorrect,3),1) ; nanmean(nanmean(LI_Correct,3),1) ];
    %  sts= [NaN    0.2021    0.2464    0.2722    0.2665;
    %        NaN    0.2914    0.3990    0.5886    0.6710]

    %  sts= [ NaN    0.1879    0.2186    0.1825    0.1758;
    %         NaN    0.3771    0.4227    0.6343    0.6492]
    errY =[1.96*nanstd(nanmean(LI_Incorrect,3),1)./sqrt(length(nanmean(nanmean(LI_Incorrect,3),2)));1.96*nanstd( nanmean(LI_Correct,3),1)./sqrt(length(nanmean(nanmean(LI_Correct,3),2))) ];
    b=barwitherr(errY, sts);% %Error bars show 95% confidence intervals
    set(gca,'xticklabel',blockNames,'FontSize',14)
    ylabel('Proportion Correct');
    ylim ([0 1]);
    ylabel('Mean proportion correct');
    %legend('Instance 1','Instance 2','Instance 3','Instance 4','Instance 5','Instance 6','Instance 7','Instance 8','Instance 9','Instance 10');
    legend('Instance 1','Instance 2','Instance 3','Instance 4','Instance 5');
    xlabel('Previous Learning Instance');
    hold on
    a=hline(0.2,'k:');
    set(a,'LineWidth',2.0);
    grid on


    figure (24);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
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

    figure (25);%Plot proportions of subjects looking to objects, target, distractor on test timrcourse 
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
    cross_correlations=zeros(12,12);
    training_pair=xsit_result.train(2).training_pair;
    for i=1:60
        ia=training_pair(i,7);   
        for j=1:4
            ja=training_pair(i,j);
        cross_correlations(ia,ja) = cross_correlations(ia,ja) +1;
        end
    end
    for i=1:12
        cross_correlations(i,i)=0;
    end
    cross_correlations
    
    
%% TRUESWELL TASK CODE PAIRINGS GENERATION
% T = readtable('Trueswell_Pairings_Chosen.xlsx'); % Loads the file into workspace
% T(1:5,:) % An overview of the table
% M=T;
% L=M;
% T.Number(10)
% counter=1;
% for i=1:5
%     for j=1:12
%         M(counter,:) = T(((j-1)*5)+i,:);
%         counter=counter+1;
%     end   
% end
% suff='.bmp'
% for k=1:60
%     %cell2mat(M.Image1(1))
%     L.Image1(1) = {string(M.Image1(1)) + string('.bmp')};
% end
% filename='testdata.xlsx';
% xlswrite(filename,L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp (Task_Name, 'Yu_Zhong_Fricker_2012')
%% Yu Zhong & Fricker (2012)
    %data: wfChanges_conwmf0_1k15k_hwmc1_fix63_wrdlen2k_new_Yu_Zhong_Fricker_2012_results
    simName = 'wPPR_1k50k_fix63_wrdlen2k_pretrain10k_Yu_Zhong_Fricker_2012_results'
    xsit_result = load (simName);
    plotStyle = {'k-o','b-+','g-*','c-x','r-s','m-d','y->','k:o','b:+','g:*','c:x','r:s','m:d','y:>','b:<','w.<'};compStyle=7;blockNames{10}=[];sts{10}=[];errY{10}=[];
    nObjects = 18; scale_factor=8;%%total number of objects   
    vis_Off = 1000/scale_factor;
    nTrials=nObjects;nFeatures=2;  %check with auto file
    numSubjects=size(xsit_result.test,1);
    inA=0;inB=0;inC=0;inD=0; corrMapped=zeros(numSubjects,1);
    corrMappedWrd=zeros(numSubjects,2); Look=[];
    for subject=1:numSubjects
        Testset = xsit_result.test(subject).Testset; %size
        %Testset, has 3 distractors, 1 target/paired-target and its location
        Words=cell2mat(xsit_result.train(subject).Words);
        for trt=1:nTrials %%
            Look(1)= sum( xsit_result.test(subject).historyLt(trt,1:vis_Off));%left extreme
            Look(2)= sum( xsit_result.test(subject).historyLWBt(trt,1:vis_Off));%left middle
            Look(3)= sum( xsit_result.test(subject).historyRWBt(trt,1:vis_Off));%right middle
            Look(4)= sum( xsit_result.test(subject).historyRt(trt,1:vis_Off));%right extreme
            targLoc = Testset(trt,5);%the spatial location of target
            %next, add looking time to all 3 distracrors
            dstrLook=0;for i=1:4;if i~=targLoc; dstrLook= dstrLook + Look(i); end; end
            if (Look(targLoc) > dstrLook) %WARNING divide by 3
                corrMapped(subject)= corrMapped(subject)+1;
                if (Testset(trt,targLoc))== 1 || (Testset(trt,targLoc))==2 || (Testset(trt,targLoc))==3
                   % Testset(trt,4)
                    corrMappedWrd(subject,1)=corrMappedWrd(subject,1)+1;
                else
                    corrMappedWrd(subject,2)=corrMappedWrd(subject,2)+1;
                end
            end   
        end
        if corrMapped(subject) > 13
            learnerType(subject) = 1;
        elseif corrMapped(subject) < 8
            learnerType(subject) = -1;
        else
            learnerType(subject) = 0;
        end
    end
    %(mean(corrMapped)/nObjects)
    %sum(learnerType==-1)

    figure(1);% proprtion correct response
    blockNames={'Yu Zhang & Fricker (2012)'; 'WOLVES Model'};
    sts = [0.58 mean(corrMapped)/nObjects];
    errY =[0.05 (std(corrMapped)/nObjects)/sqrt(length(corrMapped))];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);

    figure(11);% proprtion correct response
    blockNames={'Pre-trained words';'All other words'};
    sts = [0.9187 mean(corrMappedWrd(:,1))/3; 0.5812 mean(corrMappedWrd(:,2))/15 ];
    errY =[0.05 (std(corrMappedWrd(:,1))/3)/sqrt(length(corrMappedWrd(:,1)));0.05 (std(corrMappedWrd(:,2))/15)/sqrt(length(corrMappedWrd(:,2))) ];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames,'FontSize',14)
    ylabel('Proportion Correct');
    ylim ([0 1]);
    legend('Yu Zhang & Fricker (2012)', 'WOLVES Model');
    hold on
    a=hline(0.25,'k:');
    set(a,'LineWidth',2.0);
    grid on

    rmse = sqrt(((0.9187 - mean(corrMappedWrd(:,1))/3)^2 + (0.5812-mean(corrMappedWrd(:,2))/15)^2)./3);
    pe = mean ([abs(0.9187 - mean(corrMappedWrd(:,1))/3)/0.9187*100 abs(0.5812-mean(corrMappedWrd(:,2))/15)/0.5812*100 ] );
    disp (['RMSE = ',num2str(rmse)]); 
    disp (['MAPE = ',num2str(pe)]);
    % training analysis
    nTrials=27; vis_Off=11250/scale_factor;
%     %2k word
    word_On(1)  = floor((2250+350)/scale_factor);% 
    word_Off(1) = floor(4500/scale_factor);%
    word_On(2)  = floor((4500+350)/scale_factor);%
    word_Off(2) = floor(6750/scale_factor);%
    word_On(3)  = floor((6750+350)/scale_factor);%
    word_Off(3) = floor(9000/scale_factor);%
    word_On(4)  = floor((9000+350)/scale_factor);%
    word_Off(4) = floor(11250/scale_factor);%    
    
    tarLook=NaN(numSubjects,nObjects,6);
    distLook=NaN(numSubjects,nObjects,6);
    for subject=1:numSubjects
        occurance=zeros(nObjects,1);
        training_pair = xsit_result.train(subject).training_pair; % the objects/words presneted in the trial (actually from the sequence 1-18)
        training_order = xsit_result.train(subject).training_order; %the order in which objects, word were presneted 
        for r = 1:size(training_pair,1) 
              for c = 1:size(training_pair,2)
                  trial_object(r,c) = training_pair(r, training_order(r,c));
                  trial_word(r,c) = training_pair(r, training_order(r,c+4));
              end
        end
        for trt=1:nTrials 
            for w=1:4 % lets take word w  %it was presneted as w temporally.. therefore time on will be word_On(w):word_Off(w)
                currword=trial_word(trt,w); %and this is the currently played word    
                %and lets find where it is placed spatially
                if currword == trial_object(trt,1)
                    occurance(currword)=occurance(currword)+1;
                    tarLook(subject,currword,occurance(currword))  = sum( xsit_result.train(subject).historyL(trt,word_On(w):word_Off(w)));%left extreme
                    distLook(subject,currword,occurance(currword)) = sum( xsit_result.train(subject).historyLWB(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyRWB(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyR(trt,word_On(w):word_Off(w)));
                elseif currword == trial_object(trt,2)  
                    occurance(currword)=occurance(currword)+1;
                    tarLook(subject,currword,occurance(currword))  = sum( xsit_result.train(subject).historyLWB(trt,word_On(w):word_Off(w)));%left middle
                    distLook(subject,currword,occurance(currword)) = sum( xsit_result.train(subject).historyL(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyRWB(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyR(trt,word_On(w):word_Off(w)));
                elseif currword == trial_object(trt,3)
                    occurance(currword)=occurance(currword)+1;
                    tarLook(subject,currword,occurance(currword))  = sum( xsit_result.train(subject).historyRWB(trt,word_On(w):word_Off(w)));%right middle
                    distLook(subject,currword,occurance(currword)) = sum( xsit_result.train(subject).historyLWB(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyL(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyR(trt,word_On(w):word_Off(w)));
                elseif currword == trial_object(trt,4)
                    occurance(currword)=occurance(currword)+1;
                    tarLook(subject,currword,occurance(currword))  = sum( xsit_result.train(subject).historyR(trt,word_On(w):word_Off(w)));%right extreme
                    distLook(subject,currword,occurance(currword)) = sum( xsit_result.train(subject).historyLWB(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyRWB(trt,word_On(w):word_Off(w)))+ ... 
                                                                        sum( xsit_result.train(subject).historyL(trt,word_On(w):word_Off(w)));
                end
            end  
        end
    end
    distLook=distLook./3;
    propLook=tarLook./(tarLook+distLook);
    kak1= 1.0:0.05:1.25;kak2= 1.0:0.025:1.125; kak3= 1.0:0.0125:1.0625;%%*kak1'
    figure(2);%Plot Mean proportion correct looking at every learning instance
    errorbar(squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1)),squeeze(std(nanmean(propLook(learnerType==1,:,:),2),1))./sqrt(size(propLook(learnerType==1),2)),plotStyle{1});%./sqrt(length(mean(tarLook))
    hold on
    errorbar(squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1)),squeeze(std(nanmean(propLook(learnerType==0,:,:),2),1))./sqrt(size(propLook(learnerType==0),2)),plotStyle{2});%./sqrt(length(mean(tarLook))
    hold on
    errorbar(squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1)),squeeze(std(nanmean(propLook(learnerType==-1,:,:),2),1))./sqrt(size(propLook(learnerType==-1),2)),plotStyle{3});%./sqrt(length(mean(tarLook))
    hold on
    xlim([0.5 6.5]);
    ylim([0 1]);
    ylabel('Proportion of Time on Target');
    legend('Strong','Average','Weak');
    set (gca, 'FontSize',14);
    grid on
    lt_1_emp = [0.41; 0.38; 0.5; 0.59; 0.65; 0.71];
    lt_0_emp = [0.35; 0.31; 0.36; 0.41; 0.45; 0.51];
    lt_11_emp = [0.38; 0.30; 0.32; 0.37; 0.35; 0.41];

    %rmse = sqrt((sum((lt_1_emp - squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1))).^2) + sum((lt_0_emp - squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1))).^2) + sum((lt_11_emp - squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1))).^2))./18)
    %pe = mean ([mean(abs(lt_1_emp - squeeze(nanmean(nanmean(propLook(learnerType==1,:,:),2),1)))./lt_1_emp)*100  mean(abs(lt_0_emp - squeeze(nanmean(nanmean(propLook(learnerType==0,:,:),2),1)))./lt_1_emp)*100 mean(abs(lt_11_emp - squeeze(nanmean(nanmean(propLook(learnerType==-1,:,:),2),1)))./lt_1_emp)*100 ] )

    figure(21);%Entropy
    blockNames={'Strong';'Average';'Weak'};
    set (gca, 'FontSize',18);
    sts = [mean(EntropyTrace((learnerType()==1))); mean(EntropyTrace((learnerType()==0))); mean(EntropyTrace((learnerType()==-1)))  ];
    errY =[std(EntropyTrace((learnerType()==1)))/sqrt(length(EntropyTrace((learnerType()==1)))); std(EntropyTrace((learnerType()==0)))/sqrt(length(EntropyTrace((learnerType()==0)))); std(EntropyTrace((learnerType()==-1)))/sqrt(length(EntropyTrace((learnerType()==-1))));];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    title ('Entropy in the traces');
    %ylabel('Looking time per test trial');
    %ylim([0 4]);

    figure(22);%My own Entropy: No of incorrect traces
    blockNames={'Strong';'Average';'Weak'};
    set (gca, 'FontSize',14);
    sts = [mean(InCorr_assocs((learnerType()==1))); mean(InCorr_assocs((learnerType()==0))); mean(InCorr_assocs((learnerType()==-1)))  ];
    errY =[std(InCorr_assocs((learnerType()==1)))/sqrt(length(InCorr_assocs((learnerType()==1)))); std(InCorr_assocs((learnerType()==0)))/sqrt(length(InCorr_assocs((learnerType()==0)))); std(InCorr_assocs((learnerType()==-1)))/sqrt(length(InCorr_assocs((learnerType()==-1))));];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    title ('Proportion of incorrect assocs in the traces');
    ylabel('# of incorrect assocs/# of words (18)');

    figure(221);%correct traces
    blockNames={'Strong';'Average';'Weak'};
    set (gca, 'FontSize',14);
    sts = [mean(Correct_inTrace((learnerType()==1))); mean(Correct_inTrace((learnerType()==0))); mean(Correct_inTrace((learnerType()==-1)))  ];
    errY =[std(Correct_inTrace((learnerType()==1)))/sqrt(length(Correct_inTrace((learnerType()==1)))); std(Correct_inTrace((learnerType()==0)))/sqrt(length(Correct_inTrace((learnerType()==0)))); std(Correct_inTrace((learnerType()==-1)))/sqrt(length(Correct_inTrace((learnerType()==-1))));];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    title ('Strength of correct assocs in the traces');

    figure(222);%My own Entropy: No of incorrect traces
    blockNames={'Strong';'Average';'Weak'};
    set (gca, 'FontSize',14);
    sts = [mean(Wrong_inTrace((learnerType()==1))); mean(Wrong_inTrace((learnerType()==0))); mean(Wrong_inTrace((learnerType()==-1)))  ];
    errY =[std(Wrong_inTrace((learnerType()==1)))/sqrt(length(Wrong_inTrace((learnerType()==1)))); std(Wrong_inTrace((learnerType()==0)))/sqrt(length(Wrong_inTrace((learnerType()==0)))); std(Wrong_inTrace((learnerType()==-1)))/sqrt(length(Wrong_inTrace((learnerType()==-1))));];
    b=barwitherr(errY, sts);% Plot with errorbars
    set(gca,'xticklabel',blockNames);
    title ('Strength of INCORRECT assocs in the traces');


end


%%%%%%%%

%rmse=sqrt(((3.6-3.36)^2 + (3.25-3.36)^2 + (2.5-2.89)^2 +(2.67-2.89)^2  )./4)
%pe = mean ([abs(3.6-3.36)/3.36*100 abs(3.25-3.36)/3.25*100 abs(2.5-2.89)/2.5*100  abs(2.67-2.89)/2.67*100 ] ) %

