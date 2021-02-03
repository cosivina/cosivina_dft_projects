nObjects = 18; %%total number of objects
maxTR = 27; % total # of training trials
maxTRt = nObjects; % total # of test trials
t_max = floor(12000/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
t_maxt = floor(3000/scale_factor); %specify simulation time  % scale with Experiment test trial Duration
handle_historyL=sim.getElement('history lookL'); %reset lookingduration history variable                 
handle_historyL.timeSlots = t_max; 
handle_historyR=sim.getElement('history lookR');                       
handle_historyR.timeSlots = t_max; 
handle_historyLWB=sim.getElement('history lookLWB'); %reset lookingduration history variable                 
handle_historyLWB.timeSlots = t_max; 
handle_historyRWB=sim.getElement('history lookRWB');                       
handle_historyRWB.timeSlots = t_max; 
handle_historyC=sim.getElement('history lookC');                       
handle_historyC.timeSlots = t_max; 

handle_historyLt=sim.getElement('history lookLt');                       
handle_historyLt.timeSlots = t_maxt; 
handle_historyRt=sim.getElement('history lookRt');                       
handle_historyRt.timeSlots = t_maxt;
handle_historyLWBt=sim.getElement('history lookLWBt'); %reset lookingduration history variable                 
handle_historyLWBt.timeSlots = t_maxt; 
handle_historyRWBt=sim.getElement('history lookRWBt');                       
handle_historyRWBt.timeSlots = t_maxt; 
handle_historyCt=sim.getElement('history lookCt');                       
handle_historyCt.timeSlots = t_maxt;

if (mode == 1 )      
        gui.init();
        delete(gui.visualizations{1}.axesHandle);
        gui.visualizations{1}.axesProperties{2} = [-t_max, 10];
        gui.visualizations{1}.plotProperties{1}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{2}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{3}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{4}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{5}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.init(gui.figureHandle);
end
parfor subject = 1:numSubjects
                
    sim2 = 5;
    if (mode == 1 )             
        sim2 = sim;
    elseif (mode == 0)
        sim2 = sim;
    elseif (mode == 2)
        sim2 = sim.copy();
    end
    
    visuals_On = 0;
    visuals_Off = t_max;
    word1_On  = floor(2000/scale_factor);
    word1_Off = floor(3000/scale_factor);
    word2_On  = floor(4000/scale_factor);
    word2_Off = floor(5000/scale_factor);
    word3_On  = floor(6000/scale_factor);
    word3_Off = floor(7000/scale_factor);
    word4_On  = floor(8000/scale_factor);
    word4_Off = floor(9000/scale_factor);
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[]; 
    Words=[]; 
    for j=1:nFeatures
        chosen_indices = randperm(size(Property{j},2),nObjects);
        TempFe=[];
        for i=1:nObjects
            TempFe = [TempFe  Property{j}(chosen_indices(i))];
            %TempFe = [TempFe  Property{j}(i)];
        end
        Feature{j}=TempFe;
    end
    % Choose Words to associate on random
    chosen_indices = randperm(size(Names,2),nObjects);
    for i=1:nObjects
        Words = [Words  Names(chosen_indices(i))];
        %Words = [Words  Names(i)];
    end
    
%% DESIGN TRAINING TRIALS
    wordCombAll=[];
    SINGLE = [1,5,6,10,14,18]; DOUBLE=[2,4,12,13,15,17]; NOISE=[3,7,8,9,11,16]; Double2=[3,7,8,9,11,16]; 
    wordCombAll(1,:)= [9 2 4 8];      wordCombAll(2,:)= [2 3 5 1];      wordCombAll(3,:)= [3 4 6 10];
    wordCombAll(4,:)= [4 9 7 15];     wordCombAll(5,:)= [5 6 8 12];     wordCombAll(6,:)= [6 7 5 13];
    wordCombAll(7,:)= [15 8 10 14];    wordCombAll(8,:)= [8 9 4 15];    wordCombAll(9,:)= [9 10 12 1];
    wordCombAll(10,:)= [10 14 13 7]; wordCombAll(11,:)= [11 12 14 18]; wordCombAll(12,:)= [12 13 11 16];
    wordCombAll(13,:)= [4 11 16 2];  wordCombAll(14,:)= [14 15 18 3];  wordCombAll(15,:)= [15 16 17 11];
    wordCombAll(16,:)= [16 17 1 5];   wordCombAll(17,:)= [17 18 3 6];   wordCombAll(18,:)= [18 1 2 7];
    wordCombAll(19,:)= [1 10 6 14];   wordCombAll(20,:)= [2 11 7 15];   wordCombAll(21,:)= [3 12 8 17];
    wordCombAll(22,:)= [4 13 9 11];   wordCombAll(23,:)= [5 14 1 18];   wordCombAll(24,:)= [6 17 3 10];
    wordCombAll(25,:)= [7 16 2 17];   wordCombAll(26,:)= [8 16 13 12];   wordCombAll(27,:)= [9 18 5 13];
    objCombAll=wordCombAll;
    for tC=1:maxTR
         noisyInd = find (ismember(wordCombAll(tC,:), NOISE));
         doublyInd = find (ismember(wordCombAll(tC,:), DOUBLE));
         for xC=1:size(noisyInd,2)
             newInd = find(DOUBLE == wordCombAll(tC,doublyInd(xC)));
             objCombAll(tC,noisyInd(xC))=Double2(newInd);   
         end
    end

%% RUN TRAINING TRIALS   
     for i = 1 : nFeatures
        n=num2str(i);
        handle_hcon_f{subject}(i)=sim2.getElement(['hcon_f' n]);
        handle_hwm_f{subject}(i)=sim2.getElement(['hwm_f' n]);
        handle_hwf{subject}(i)=sim2.getElement(['hwf' n]);
        handle_hwm_c{subject}(i)=sim2.getElement(['hwm_c' n]);       
    end
    handle_hcon_s{subject}=sim2.getElement('hcon_s');
    handle_hwm_s{subject}=sim2.getElement('hwm_s');
    handle_hword{subject}=sim2.getElement('hword');
    
    savestate_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
    savestate_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
    savestate_hcon_s{subject}=zeros(fieldSize_spt);
    savestate_hwm_s{subject}=zeros(fieldSize_spt);
    savestate_hword{subject}=zeros(fieldSize_wd);
    savestate_historyL{subject}=zeros(maxTR,t_max);
    savestate_historyR{subject}=zeros(maxTR,t_max);
    savestate_historyLWB{subject}=zeros(maxTR,t_max);
    savestate_historyRWB{subject}=zeros(maxTR,t_max);
    
    for tr=1:maxTR
        %sim2.init();
        sim2.t =sim2.tZero;
        ObjPos=randperm(4);
        WordPos=randperm(4);
        if tr > 1
            for i = 1 : nFeatures
                n=num2str(i);
                handle_hcon_f{subject}(i).output = squeeze(savestate_hcon_f{subject}(i,:));
                handle_hwm_f{subject}(i).output = squeeze(savestate_hwm_f{subject}(i,:));
                handle_hwf{subject}(i).output = squeeze(savestate_hwf{subject}(i,:,:));
                handle_hwm_c{subject}(i).output = squeeze(savestate_hwm_c{subject}(i,:,:));
            end
            handle_hcon_s{subject}.output = savestate_hcon_s{subject};
            handle_hwm_s{subject}.output = savestate_hwm_s{subject};
            handle_hword{subject}.output = savestate_hword{subject};
        end
        
        while sim2.t <= t_max
            t = sim2.t;
            
            if t == visuals_On
                for f=1:nFeatures
                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(objCombAll(tr,ObjPos(1)))),cell2mat(Feature{f}(objCombAll(tr,ObjPos(2)))), ...
                        cell2mat(Feature{f}(objCombAll(tr,ObjPos(3)))),cell2mat(Feature{f}(objCombAll(tr,ObjPos(4))))});                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX'}, ...
                         {spt_L, spt_Lm, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
                         { vstrength, vstrength, vstrength, vstrength }); 
                    
                end 
            end
            
            if t== word1_On
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(1)))),wstrength});
            end
            if t== word1_Off
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(1)))),0});
            end
            
            if t== word2_On
                sim2.setElementParameters({'Word_Label_2'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(2)))),wstrength});
            end
            if t== word2_Off
                sim2.setElementParameters({'Word_Label_2'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(2)))),0});
            end
            
            if t== word3_On
                sim2.setElementParameters({'Word_Label_3'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(3)))),wstrength});
            end
            if t== word3_Off
                sim2.setElementParameters({'Word_Label_3'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(3)))),0});
            end
            
            if t== word4_On
                sim2.setElementParameters({'Word_Label_4'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(4)))),wstrength});
            end
            if t== word4_Off
                sim2.setElementParameters({'Word_Label_4'}, {'position','amplitude'},... 
                    {cell2mat(Words(wordCombAll(tr,WordPos(4)))),0});
            end
            
            if t == visuals_Off
                for i = 1 : nFeatures
                    n=num2str(i);
                    savestate_hcon_f{subject}(i,:) = sim2.getComponent(['hcon_f' n], 'output');
                    savestate_hwm_f{subject}(i,:) = sim2.getComponent(['hwm_f' n], 'output');
                    savestate_hwf{subject}(i,:,:) = sim2.getComponent(['hwf' n], 'output');
                    savestate_hwm_c{subject}(i,:,:) = sim2.getComponent(['hwm_c' n], 'output');
                end
                savestate_hcon_s{subject} = sim2.getComponent('hcon_s', 'output');
                savestate_hwm_s{subject} = sim2.getComponent('hwm_s', 'output');
                savestate_hword{subject} = sim2.getComponent('hword', 'output');
                savestate_historyL{subject}(tr,:) = flipud(sim2.getComponent('history lookL', 'output'));
                savestate_historyR{subject}(tr,:) = flipud(sim2.getComponent('history lookR', 'output'));
                savestate_historyLWB{subject}(tr,:) = flipud(sim2.getComponent('history lookLWB', 'output'));
                savestate_historyRWB{subject}(tr,:) = flipud(sim2.getComponent('history lookRWB', 'output'));
                
               for f=1:nFeatures
                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(objCombAll(tr,ObjPos(1)))),cell2mat(Feature{f}(objCombAll(tr,ObjPos(2)))), ...
                        cell2mat(Feature{f}(objCombAll(tr,ObjPos(3)))),cell2mat(Feature{f}(objCombAll(tr,ObjPos(4))))});                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX'}, ...
                         {spt_L, spt_Lm, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
                         { 0, 0, 0, 0 }); 
                    
                end
            end
           if (mode == 1) && (mod(t,gui_speed)==0)
                if ~gui.pauseSimulation
                    sim2.step();
                end
                if gui.quitSimulation
                    gui.close();
                    break;
                end             
                gui.step();
                %gui.checkAndUpdateControls();
                %gui.updateVisualizations();                  
           else
                sim2.step();
           end
         end
   OutName = [simName,num2str(subject),'_train.mat'];
    matTrain=matfile(OutName,'writable',true);

    matTrain.hcon_f = savestate_hcon_f{subject};
    matTrain.hwm_f = savestate_hwm_f{subject};
    matTrain.hwm_c = savestate_hwm_c{subject};
    matTrain.hwf = savestate_hwf{subject};
    matTrain.hcon_s = savestate_hcon_s{subject};
    matTrain.hwm_s = savestate_hwm_s{subject};
    matTrain.hword = savestate_hword{subject};
%STOPPED TRAIN HISTORY TOO BIG
%     matTrain.historyL = savestate_historyL{subject};
%     matTrain.historyR = savestate_historyR{subject};
%     matTrain.historyLWB = savestate_historyLWB{subject};
%     matTrain.historyRWB = savestate_historyRWB{subject};

    matTrain.Feature1 = Feature{1};
    matTrain.Feature2 = Feature{2};
    matTrain.Words = Words;     
    end
%% TEST TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if (mode == 1), tts('Test Trials start now','Microsoft Hazel Desktop - English (Great Britain)'),end
    visuals_On = 0;
    visuals_Off = t_maxt;
%     word1_On = 0;
%     word1_Off = t_maxt;
    word1_On = floor([0 1100 2200]/scale_factor);%[0 450 875 1300 ];    
    word1_Off = floor([800 1900 3000]/scale_factor);%[375 650 1075 1500 ]; 
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyLWBt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRWBt{subject}=zeros(maxTRt,t_maxt);
    test_data=[];
    testSet=[];
    for trt=1:maxTRt
        ObjPos=randperm(4);
        if ismember(trt, SINGLE) 
            testSet(trt,1)= trt;%referent
            inXx = find(SINGLE == trt);
            inX= randi(6);while(inXx==inX); inX= randi(6);end
            testSet(trt,2)= SINGLE(inX);    %single distrctor
            testSet(trt,3)= DOUBLE(inX);    %doublle distrctor
            inX= mod(inX,6)+1;
            testSet(trt,4)= Double2(inX);    %another double distrctor
            test_data(trt,1:4) =ObjPos;
            test_data(trt,5) ='S';
        elseif ismember(trt, DOUBLE)
            testSet(trt,1)= trt;            %referent 
            inXx = find(DOUBLE == trt);
            testSet(trt,2)= Double2(inXx);   %another referent 
            inX= randi(6);while(inXx==inX); inX= randi(6);end  
            testSet(trt,3)= DOUBLE(inX);    % double distractor
            testSet(trt,4)= SINGLE(inX);    % single distractor
            test_data(trt,1:4) =ObjPos;
            test_data(trt,5) ='D';
        elseif ismember(trt, NOISE)
            inX= randi(6); testSet(trt,1)= SINGLE(inX); 
            inX= randi(6); testSet(trt,2)= DOUBLE(inX);
            inX= mod(inX,6)+1; testSet(trt,3)= Double2(inX);
            inX= mod(inX,6)+1; testSet(trt,4)= DOUBLE(inX);
            test_data(trt,1:4) =ObjPos;
            test_data(trt,5) ='N';
        end        
        sim2.init();  
        
        for i = 1 : nFeatures
                n=num2str(i);
                handle_hcon_f{subject}(i).output = squeeze(savestate_hcon_f{subject}(i,:));
                handle_hwm_f{subject}(i).output = squeeze(savestate_hwm_f{subject}(i,:));
                handle_hwf{subject}(i).output = squeeze(savestate_hwf{subject}(i,:,:));
                handle_hwm_c{subject}(i).output = squeeze(savestate_hwm_c{subject}(i,:,:));
        end
        handle_hcon_s{subject}.output = savestate_hcon_s{subject};
        handle_hwm_s{subject}.output = savestate_hwm_s{subject};
        handle_hword{subject}.output = savestate_hword{subject};
        
        while sim2.t <= t_maxt

            t = sim2.t;
            
            if t == visuals_On
                for f=1:nFeatures
                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(testSet(trt,ObjPos(1)))),cell2mat(Feature{f}(testSet(trt,ObjPos(2)))), ...
                        cell2mat(Feature{f}(testSet(trt,ObjPos(3)))),cell2mat(Feature{f}(testSet(trt,ObjPos(4))))});                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX'}, ...
                         {spt_L, spt_Lm, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
                         { vstrength, vstrength, vstrength, vstrength }); 
                    
                end 
            end
            
            if isalmost(t, word1_On,tolerance) 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(trt)),wstrength});
            end
            if isalmost(t, word1_Off,tolerance) 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {0,0});
            end
            

            if t == visuals_Off
                for i = 1 : nFeatures
                    n=num2str(i);
                    savestate_hcon_f{subject}(i,:) = sim2.getComponent(['hcon_f' n], 'output');
                    savestate_hwm_f{subject}(i,:) = sim2.getComponent(['hwm_f' n], 'output');
                    savestate_hwf{subject}(i,:,:) = sim2.getComponent(['hwf' n], 'output');
                    savestate_hwm_c{subject}(i,:,:) = sim2.getComponent(['hwm_c' n], 'output');
                end
                savestate_hcon_s{subject} = sim2.getComponent('hcon_s', 'output');
                savestate_hwm_s{subject} = sim2.getComponent('hwm_s', 'output');
                savestate_hword{subject} = sim2.getComponent('hword', 'output');
                savestate_historyLt{subject}(trt,:) = flipud(sim2.getComponent('history lookLt', 'output'));
                savestate_historyRt{subject}(trt,:) = flipud(sim2.getComponent('history lookRt', 'output'));
                savestate_historyLWBt{subject}(trt,:) = flipud(sim2.getComponent('history lookLWBt', 'output'));
                savestate_historyRWBt{subject}(trt,:) = flipud(sim2.getComponent('history lookRWBt', 'output'));
                
               for f=1:nFeatures
                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(testSet(trt,ObjPos(1)))),cell2mat(Feature{f}(testSet(trt,ObjPos(2)))), ...
                        cell2mat(Feature{f}(testSet(trt,ObjPos(3)))),cell2mat(Feature{f}(testSet(trt,ObjPos(4))))});                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX'}, ...
                         {spt_L, spt_Lm, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
                         { 0, 0, 0, 0 }); 
                    
                end
            end

           if (mode == 1) && (mod(t,2)==0)
                if ~gui.pauseSimulation
                    sim2.step();
                end
                if gui.quitSimulation
                    gui.close();
                    break;
                end             
                gui.step();
                %gui.checkAndUpdateControls();
                %gui.updateVisualizations();                  
           else
                sim2.step();
           end
        end        
    end
    OutName = [simName,num2str(subject),'_test.mat'];
    matTest=matfile(OutName,'writable',true);

    matTrain.hcon_ft = savestate_hcon_f{subject};
    matTest.hwm_ft = savestate_hwm_f{subject};
    matTest.hwm_ct = savestate_hwm_c{subject};
    matTest.hwft = savestate_hwf{subject};
    matTrain.hcon_st = savestate_hcon_s{subject};
    matTrain.hwm_st = savestate_hwm_s{subject};
    matTest.hwordt = savestate_hword{subject};
    
    matTest.historyLt = savestate_historyLt{subject};
    matTest.historyRt = savestate_historyRt{subject};
    matTest.historyLWBt = savestate_historyLWBt{subject};
    matTest.historyRWBt = savestate_historyRWBt{subject};
    matTest.test_data = test_data;
end



