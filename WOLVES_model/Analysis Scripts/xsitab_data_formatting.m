%% Save Word-Object Mappings, Training and Test Trials -SY
% data=readtable('XSITSY_mappings.xlsx');%'Format','%s%s%u%f%f%s'
% %data1 = data(:,{'version','TRIAL_INDEX','condition','TRIAL_DURATION','Other','VideoXVD','Object_Left','Object_Right','First_Word','Second_Word','VideoAVI','SoundWAV1','SoundWAV1'});
% Order(1).Words =    {'manu'; 'bosa';    'gake'; 'regli'; 'kaki'; 'colat'};
% Order(1).Objects =  {'teal'; 'yellow';  'cyan'; 'green'; 'brown'; 'Rblue'};
% Order(1).training_pair=strings(30,5);
% Order(1).test_pair=strings(12,4);
% Order(2).Words =    {'manu';    'bosa';  'gake';    'regli';    'kaki'; 'colat'};
% Order(2).Objects =  {'Yellow';  'Brown'; 'Rblue';   'Cyan';     'Green'; 'Teal'};
% Order(2).training_pair=strings(30,5);
% Order(2).test_pair=strings(12,4);
% for j=1:2
%     for k=1:42
%         i=(j-1)*42+k;
%         if strcmpi(data.version(i),'Order1')
%             if k<=30
%                 Order(1).training_pair(k,1) = cell2mat(data.Object_Left(i));
%                 Order(1).training_pair(k,2) = cell2mat(data.Object_Right(i));
%                 Order(1).training_pair(k,3) = cell2mat(data.First_Word(i));
%                 Order(1).training_pair(k,4) = cell2mat(data.Second_Word(i));
%                 for m=1:6
%                     if strcmpi(cell2mat(Order(1).Words(m)), cell2mat(data.First_Word(i)))
%                         if strcmpi(cell2mat(Order(1).Objects(m)), cell2mat(data.Object_Left(i)))
%                             Order(1).training_pair(k,5) = 'P';
%                         else
%                             Order(1).training_pair(k,5) = 'X';
%                         end
%                         break;
%                     end
%                 end
%             elseif k>30
%                 Order(1).test_pair(k-30,1) = cell2mat(data.Object_Left(i));
%                 Order(1).test_pair(k-30,2) = cell2mat(data.Object_Right(i));
%                 Order(1).test_pair(k-30,3) = cell2mat(data.First_Word(i));
%                 for m=1:6
%                     if strcmpi(cell2mat(Order(1).Words(m)), cell2mat(data.First_Word(i)))
%                         if strcmpi(cell2mat(Order(1).Objects(m)), cell2mat(data.Object_Left(i)))
%                             Order(1).test_pair(k-30,4) = 'L';
%                         else
%                             Order(1).test_pair(k-30,4) = 'R';
%                         end
%                         break;
%                     end
%                 end
%             end
%         elseif strcmpi(data.version(i),'Order2')
%             if k<=30
%                 Order(2).training_pair(k,1) = cell2mat(data.Object_Left(i));
%                 Order(2).training_pair(k,2) = cell2mat(data.Object_Right(i));
%                 Order(2).training_pair(k,3) = cell2mat(data.First_Word(i));
%                 Order(2).training_pair(k,4) = cell2mat(data.Second_Word(i));
%                 for m=1:6
%                     if strcmpi(cell2mat(Order(2).Words(m)), cell2mat(data.First_Word(i)))
%                         if strcmpi(cell2mat(Order(2).Objects(m)), cell2mat(data.Object_Left(i)))
%                             Order(2).training_pair(k,5) = 'P';
%                         else
%                             Order(2).training_pair(k,5) = 'X';
%                         end
%                         break;
%                     end
%                 end
%             elseif k>30
%                 Order(2).test_pair(k-30,1) = cell2mat(data.Object_Left(i));
%                 Order(2).test_pair(k-30,2) = cell2mat(data.Object_Right(i));
%                 Order(2).test_pair(k-30,3) = cell2mat(data.First_Word(i));
%                 for m=1:6
%                     if strcmpi(cell2mat(Order(2).Words(m)), cell2mat(data.First_Word(i)))
%                         if strcmpi(cell2mat(Order(2).Objects(m)), cell2mat(data.Object_Left(i)))
%                             Order(2).test_pair(k-30,4) = 'L';
%                         else
%                             Order(2).test_pair(k-30,4) = 'R';
%                         end
%                         break;
%                     end
%                 end
%             end
%         end
%     end
% end

clear all
close all
diary output1.txt
load Order.mat
delete 'output.txt'
%historyL[trial_time] % For a subject in a trial; create a historyL array 
                     % that shows whether the kid was looking at left object during a trial 
%historyR[trial_time] % historyR for looking history to right side
%historyOFF[trial_time] % historyOFF for off-looking .. may be needed
% next to decide where the kid is looking
% find/get the screen size.. say 1024x768
% decide left side bounding box and right side bounding box-interest areas
% take a kid id.. take an order..take a condition.. take a trial..
% take a sample event and classify/value it based on the x,y bounding conditions of left, right and off 
% find a mapping between LR samples counts 8 sc or 4 sec trial time..
% save Order.condition.subject_id.trial.historyR[trialNumbr,LRSamleCount]
% save Order.test.subject_id.trial.test_pair(trialNumbr,[leftObj, Right Obj, Word])
% save Order.Train.subject_id.trial.training_pair(trialNumbr,[leftObj, Right Obj, Word_1, word_2])
% Save Order.mappings[Object,Word]
%ReadVariableNames true
%Delimiter \t  .. \b

%ds = spreadsheetDatastore('test2.xlsx');%%tabularTextDatastore
%first_word	second_word	left_obj	right_obj	TIMESTAMP	TRIAL_START_TIME	RIGHT_GAZE_X	RIGHT_GAZE_Y	RIGHT_FIX_INDEX	RIGHT_INTEREST_AREAS	RIGHT_INTEREST_AREA_LABEL
%ds.SelectedVariableNames = {'RECORDING_SESSION_LABEL','TRIAL_INDEX','version','condition','TRIAL_START_TIME','RIGHT_INTEREST_AREA_LABEL'};
%ds.SelectedFormats{strcmpi(ds.SelectedVariableNames,'version')} = '%s';
%ds.SelectedFormats{strcmpi(ds.SelectedVariableNames,'TRIAL_INDEX')} = '%d';
%T = readall(ds);

opts = detectImportOptions('rawdata3.txt'); %getvaropts(opts,{'TRIAL_INDEX'})
opts = setvartype(opts,{'TRIAL_INDEX','RIGHT_FIX_INDEX'},'int32');
opts.SelectedVariableNames = {'RECORDING_SESSION_LABEL','TRIAL_INDEX','version','condition','TIMESTAMP','TRIAL_START_TIME','RIGHT_FIX_INDEX','RIGHT_INTEREST_AREA_LABEL','SAMPLE_MESSAGE'};
T = readtable('rawdata3.txt',opts);
T(1:5,1:9)
% % ROW=1  %%% checking if i can use sound play timings too
% % prev_STAM= T.TIMESTAMP(ROW)
% % for ROW = 1:size(T,1)
% %     if (contains(T.SAMPLE_MESSAGE(ROW),'PLAY_SOU'))
% %         T.SAMPLE_MESSAGE(ROW)
% %         T.TIMESTAMP(ROW)-prev_STAM
% %         prev_STAM=T.TIMESTAMP(ROW);
% %     end
% % end

expName='XSITAB_';
LEARNING_TRIALS_COUNTER=30;TESTING_TRIALS_COUNTER=12; 
TRAIN_DUR=4000; TEST_DUR=8000; MIN_LOOK_DURATION=80;
ROW=1;subject=0;lCount=0;
while (ROW <= size(T,1)) % for each subject
    curr_Label= T.RECORDING_SESSION_LABEL(ROW);
    if (strcmpi(T.version(ROW), 'Order1'))
        curr_Order= 1;
    elseif (strcmpi(T.version(ROW), 'Order2'))
        curr_Order= 2;
    end
    subject=subject+1;
    lCount=0;
    tCount=0;
    while strcmpi(curr_Label,T.RECORDING_SESSION_LABEL(ROW))
        curr_Condition= T.condition(ROW);
        if (strcmpi(curr_Condition, 'LEARN'))
           %disp('Expected Learn OK');
           historyL=0;historyR=0;historyO=0;
           while (ROW <= size(T,1) && strcmpi(T.condition(ROW), 'LEARN'))
                curr_Trial = T.TRIAL_INDEX(ROW);
                timeIndex=0;
                lCount=lCount+1;
                time_Stam=T.TIMESTAMP(ROW);
                t_start_time=T.TRIAL_START_TIME(ROW);
                row_ind=ROW;
                while (ROW <= size(T,1) && eq(T.TRIAL_INDEX(ROW),curr_Trial)) %strcmpi
                    timeIndex=timeIndex+1; 
                    if  strcmpi(T.RIGHT_INTEREST_AREA_LABEL(ROW),'LEFT')%% 2d ehistoy ma not be needed
                        historyL(curr_Trial,timeIndex)=1; historyR(curr_Trial,timeIndex)=0; historyO(curr_Trial,timeIndex)=0;
                        historyL(curr_Trial,timeIndex+1)=1; historyR(curr_Trial,timeIndex+1)=0; historyO(curr_Trial,timeIndex+1)=0;
                    elseif  strcmpi(T.RIGHT_INTEREST_AREA_LABEL(ROW),'RIGHT')
                        historyL(curr_Trial,timeIndex)=0; historyR(curr_Trial,timeIndex)=1; historyO(curr_Trial,timeIndex)=0;
                        historyL(curr_Trial,timeIndex+1)=0; historyR(curr_Trial,timeIndex+1)=1; historyO(curr_Trial,timeIndex+1)=0;
                    else
                        historyL(curr_Trial,timeIndex)=0; historyR(curr_Trial,timeIndex)=0; historyO(curr_Trial,timeIndex)=1;
                        historyL(curr_Trial,timeIndex+1)=0; historyR(curr_Trial,timeIndex+1)=0; historyO(curr_Trial,timeIndex+1)=1;
                    end
                    ROW=ROW+1;
                    timeIndex=timeIndex+1; %%doing it twice bcoz the samling rate is twice the event rate
                end     
                  %output size error due to varying timeIndex
%                 historyL(curr_Trial,:)=round(squeeze(resample(ehistoryL(curr_Trial,:),TRAIN_DUR,timeIndex)));
%                 historyR(curr_Trial,:)=round(squeeze(resample(ehistoryR(curr_Trial,:),TRAIN_DUR,timeIndex)));
%                 historyO(curr_Trial,:)=round(squeeze(resample(ehistoryO(curr_Trial,:),TRAIN_DUR,timeIndex)));
                %T.TRIAL_INDEX(ROW)
                %T.TIMESTAMP(ROW)-time_Stam
                %T.TRIAL_START_TIME(ROW)-t_start_time
                %ROW-row_ind
                %% SMOOTHEN
                %%% Smoothing data..remove looks less than theshold
                for side=1:2 
                    nlooks =0;ldata=0;cdata=0;odata=0;
                    if side == 1
                        ldata = historyL(curr_Trial,:);
                        cdata = historyR(curr_Trial,:);
                        odata = historyO(curr_Trial,:);
                    else
                        ldata = historyR(curr_Trial,:);
                        cdata = historyL(curr_Trial,:);
                        odata = historyO(curr_Trial,:);
                    end                      
                    for tr=1:size(ldata,1) %% or size(rdata,1) 
                        looking=0;
                        templonglookdur=0;
                        lastLookEndTime=1;
                        for time=1:size(ldata,2) %% or size(rdata,2) %sum(round(odata(tr,lastLookEndTime:time))) > prev_look_threshold ||
                            if (round(ldata(tr,time)) == 1)%% if model looking to left
                                if looking == 0 %% if not previously looking left 
                                    if (sum(round(cdata(tr,lastLookEndTime:time))) > MIN_LOOK_DURATION) || (sum(round(odata(tr,lastLookEndTime:time))) > 2*MIN_LOOK_DURATION) || (nlooks == 0) %%%% if it was looking off or to some other object before this look or away or no looks so far
                                        nlooks = nlooks+1;
                                        %templonglookdur=0;
                                    else
                                        %templonglookdur= templonglookdur + (time-lastLookEndTime); %%% add the gap between consecutive looks to same side
                                        %fix the arrays
                                        ldata(tr,lastLookEndTime:time)=1;
                                        odata(tr,lastLookEndTime:time)=0;
                                    end
                                end
                                looking = 1;
                                %templonglookdur = templonglookdur+1;
                            else
                                if looking == 1; lastLookEndTime=time;end
                                looking = 0;
                            end
                        end
                    end  
                    if side == 1
                             historyL(curr_Trial,:)=ldata;
                             historyR(curr_Trial,:)=cdata;
                             historyO(curr_Trial,:)=odata;
                    else
                             historyR(curr_Trial,:)=ldata;
                             historyL(curr_Trial,:)=cdata;
                             historyO(curr_Trial,:)=odata;
                    end
                end                
            end
        elseif (strcmpi(curr_Condition, 'TEST'))
            historyLt=0;historyRt=0;historyOt=0;
            %disp('Expected Test OK');
            while (ROW <= size(T,1) && strcmpi(T.condition(ROW), 'TEST'))
                curr_Trial = T.TRIAL_INDEX(ROW);
                test_Index=curr_Trial-LEARNING_TRIALS_COUNTER;
                timeIndex=0;
                 tCount=tCount+1;
                while (ROW <= size(T,1) && eq(T.TRIAL_INDEX(ROW),curr_Trial) )    %strcmpi
                    timeIndex=timeIndex+1;
                    if strcmpi(T.RIGHT_INTEREST_AREA_LABEL(ROW),'LEFT')
                        historyLt(test_Index,timeIndex)=1; historyRt(test_Index,timeIndex)=0; historyOt(test_Index,timeIndex)=0;
                        historyLt(test_Index,timeIndex+1)=1; historyRt(test_Index,timeIndex+1)=0; historyOt(test_Index,timeIndex+1)=0;
                    elseif strcmpi(T.RIGHT_INTEREST_AREA_LABEL(ROW),'RIGHT')
                        historyLt(test_Index,timeIndex)=0; historyRt(test_Index,timeIndex)=1; historyOt(test_Index,timeIndex)=0;
                        historyLt(test_Index,timeIndex+1)=0; historyRt(test_Index,timeIndex+1)=1; historyOt(test_Index,timeIndex+1)=0;
                    else
                        historyLt(test_Index,timeIndex)=0; historyRt(test_Index,timeIndex)=0; historyOt(test_Index,timeIndex)=1;
                        historyLt(test_Index,timeIndex+1)=0; historyRt(test_Index,timeIndex+1)=0; historyOt(test_Index,timeIndex+1)=1;
                    end
                    ROW=ROW+1;
                    timeIndex=timeIndex+1;  %%doing it twice bcoz the samling rate is twice the event rate
                end
                %% SMOOTHEN
                %%% Smoothing data..remove looks less than theshold
                for side=1:2 
                    nlooks =0;ldata=0;cdata=0;odata=0;
                    if side == 1
                        ldata = historyLt(test_Index,:);
                        cdata = historyRt(test_Index,:);
                        odata = historyOt(test_Index,:);
                    else
                        ldata = historyRt(test_Index,:);
                        cdata = historyLt(test_Index,:);
                        odata = historyOt(test_Index,:);
                    end                      
                    for tr=1:size(ldata,1) %% or size(rdata,1) 
                        looking=0;
                        templonglookdur=0;
                        lastLookEndTime=1;
                        for time=1:size(ldata,2) %% or size(rdata,2) %sum(round(odata(tr,lastLookEndTime:time))) > prev_look_threshold ||
                            if (round(ldata(tr,time)) == 1)%% if model looking to left
                                if looking == 0 %% if not previously looking left 
                                    if (sum(round(cdata(tr,lastLookEndTime:time))) > MIN_LOOK_DURATION) || (sum(round(odata(tr,lastLookEndTime:time))) > 2*MIN_LOOK_DURATION) || (nlooks == 0) %%%% if it was looking off or to some other object before this look or away or no looks so far
                                        nlooks = nlooks+1;
                                        %templonglookdur=0;
                                    else
                                        %templonglookdur= templonglookdur + (time-lastLookEndTime); %%% add the gap between consecutive looks to same side
                                        %fix the arrays
                                        ldata(tr,lastLookEndTime:time)=1;
                                        odata(tr,lastLookEndTime:time)=0;
                                    end
                                end
                                looking = 1;
                                %templonglookdur = templonglookdur+1;
                            else
                                if looking == 1; lastLookEndTime=time;end
                                looking = 0;
                            end
                        end
                    end  
                    if side == 1
                             historyLt(test_Index,:)=ldata;
                             historyRt(test_Index,:)=cdata;
                             historyOt(test_Index,:)=odata;
                    else
                             historyRt(test_Index,:)=ldata;
                             historyLt(test_Index,:)=cdata;
                             historyOt(test_Index,:)=odata;
                    end
                end
            end
        else
            disp('Expected Learn or Test, found NA: ERROR at ROW!!!'); disp(ROW); ROW= ROW+1;
        end
        if(ROW > size(T,1))%SANITY CHECK
            disp('ROW indxx crossed bounds');
            break;
        end
    end    
        %% save subject data
        disp (['subject ',num2str(subject), ' processed with ',num2str(lCount),' learning trials and ',num2str(tCount),' test trials: moving to next']);
        OutName = [expName,num2str(subject),'_train.mat'];
        matTrain=matfile(OutName,'writable',true);
        matTrain.historyL = historyL;
        matTrain.historyR = historyR;
        matTrain.historyO = historyO;
        % if we know Order,we know these
        %matTrain.training_pair = training_pair; 
        %matTrain.Objects = Objects;
        %matTrain.Words = Words;
        matTrain.Order = curr_Order;
        matTrain.ID = curr_Label;
        %matTrain.subject = subject;
        
        OutName = [expName,num2str(subject),'_test.mat'];
        matTest=matfile(OutName,'writable',true);
        matTest.historyLt = historyLt;
        matTest.historyRt = historyRt;
        matTest.historyOt = historyOt;
        %matTest.test_pair = test_pair;
         
end

numSubjects=subject;
train=[];test=[];
%% Save Data          
for subject=1:numSubjects
    try
        OutName1 = [expName,num2str(subject),'_train.mat'];              
        OutName2 = [expName,num2str(subject),'_test.mat'];
        tempTrn=load(OutName1);
        tempTst=load(OutName2);
        %[tempTrn(:).subject] = subject;
        train = [train; tempTrn];
        test = [test; tempTst];
        delete(OutName1);
        delete(OutName2);
    catch
        disp('Error on concatenating subject number ');
        disp(subject);
        continue;
    end
end
OutName = [expName,'results.mat'];
save(OutName,'Order','train','test');% add orders and training and testing pairs here




