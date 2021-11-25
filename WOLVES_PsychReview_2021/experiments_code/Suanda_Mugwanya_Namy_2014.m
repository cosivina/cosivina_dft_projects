%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

nObjects = 8; %%total number of objects
maxTR = 16; % total # of training trials
maxTRt = nObjects; % total # of test trials
t_max = floor((4000+1000)/scale_factor); %specify simulation time = scaled with empirical Training trial duration in millisecs + 1000 ms break between trials   
t_maxt = floor((1000+1000)/scale_factor); %specify simulation time = scaled with empirical Test trial duration in millisecs + 1000 ms break between trials

History_Reset;
sim.init();
Gui_History_Reset;

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
    visuals_Off = floor((4000)/scale_factor);
    word1_On = floor([500, 2000]/scale_factor); word1_Off = floor([1500, 3000]/scale_factor); 
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
    if mod(subject,3)==1
    %% Condition HCD-High
        Object_pairs=[];
        Object_pairs(1,:)= [1 2]; Object_pairs(2,:)= [1 3]; Object_pairs(3,:)= [1 4]; Object_pairs(4,:)= [1 5];
        Object_pairs(5,:)= [2 3]; Object_pairs(6,:)= [2 4]; Object_pairs(7,:)= [2 6]; Object_pairs(8,:)= [3 4];
        Object_pairs(9,:)= [3 7]; Object_pairs(10,:)= [4 8]; Object_pairs(11,:)= [5 6]; Object_pairs(12,:)= [5 7];
        Object_pairs(13,:)= [5 8]; Object_pairs(14,:)= [6 7]; Object_pairs(15,:)= [6 8]; Object_pairs(16,:)= [7 8];
    elseif mod(subject,3)==2
    %% Condition MCD-Moderate
        Object_pairs=[];
        Object_pairs(1,:)= [1 2]; Object_pairs(2,:)= [1 2]; Object_pairs(3,:)= [1 3]; Object_pairs(4,:)= [1 4];
        Object_pairs(5,:)= [2 3]; Object_pairs(6,:)= [2 4]; Object_pairs(7,:)= [3 4]; Object_pairs(8,:)= [3 4];
        Object_pairs(9,:)= [5 6]; Object_pairs(10,:)= [5 6]; Object_pairs(11,:)= [5 7]; Object_pairs(12,:)= [5 8];
        Object_pairs(13,:)= [6 7]; Object_pairs(14,:)= [6 8]; Object_pairs(15,:)= [7 8]; Object_pairs(16,:)= [7 8];
    elseif mod(subject,3)==0
    %% Condition LCD-Low
        Object_pairs=[];
        Object_pairs(1,:)= [1 2]; Object_pairs(2,:)= [1 2]; Object_pairs(3,:)= [1 2]; Object_pairs(4,:)= [1 3];
        Object_pairs(5,:)= [2 4]; Object_pairs(6,:)= [3 4]; Object_pairs(7,:)= [3 4]; Object_pairs(8,:)= [3 4];
        Object_pairs(9,:)= [5 6]; Object_pairs(10,:)= [5 6]; Object_pairs(11,:)= [5 6]; Object_pairs(12,:)= [5 7];
        Object_pairs(13,:)= [6 8]; Object_pairs(14,:)= [7 8]; Object_pairs(15,:)= [7 8]; Object_pairs(16,:)= [7 8];
    end
    Sequence=randperm(size(Object_pairs,1));
    
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
    %training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%
    %% RUN TRAINING TRIALS
    for tr=1:maxTR

        sim2.t =sim2.tZero;
        pos =  Sequence(tr); % Finds postion of a random tuple
        inF=randperm(2);
        inW=randperm(2);
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
            
            % STEP 1 Add sesame attention 3 second for first 4 trials
            if isalmost (t, visuals_On, tolerance)
                for f=1:nFeatures
                    
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                        strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(Object_pairs(pos,inF(1)))),cell2mat(Feature{f}(Object_pairs(pos,inF(2))))});                

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set spatial dimension value
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'positionX', 'positionX'}, {spt_L,spt_R});

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set left right stimuli up
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'amplitude', 'amplitude'}, vstrength * ones(1, 2));                     
                end
            end

            for www=1:length(word1_On)
                if isalmost(t, word1_On(www),tolerance) 
                    sim2.setElementParameters({['Word_Label_',num2str(www)]}, {'position','amplitude'},... 
                        {cell2mat(Words(Object_pairs(pos,inW(www)))),wstrength});
                end
                if isalmost(t, word1_Off(www),tolerance)
                    sim2.setElementParameters({['Word_Label_',num2str(www)]}, {'position','amplitude'},... 
                        {cell2mat(Words(Object_pairs(pos,inW(www)))),0});
                end
            end

            if isalmost (t, visuals_Off, tolerance)
                for f=1:nFeatures
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... 
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'amplitude', 'amplitude'}, 0 * ones(1, 2));% stop left right stimuli

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... 
                        strcat('Feature_',num2str(f),'_Right')}, ...
                        {'positionY', 'positionY'},{0,0}); % Unset feature dimension value

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... 
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'positionX', 'positionX'}, {0,0});% Unset spatial dimension value
                end
            end
            
            if isalmost(t, t_max,tolerance)
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
                 savestate_historyR{subject}(tr,:) =flipud(sim2.getComponent('history lookR', 'output'));
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
           else
                sim2.step();
           end
        end % Time sim.t
    end % Trial maxTR loop
    try
    OutName = [simName,num2str(subject),'_train.mat'];
    matTrain=matfile(OutName,'writable',true);

    matTrain.hcon_f = savestate_hcon_f{subject};
    matTrain.hwm_f = savestate_hwm_f{subject};
    matTrain.hwm_c = savestate_hwm_c{subject};
    matTrain.hwf = savestate_hwf{subject};
    matTrain.hcon_s = savestate_hcon_s{subject};
    matTrain.hwm_s = savestate_hwm_s{subject};
    matTrain.hword = savestate_hword{subject};

    matTrain.historyL = savestate_historyL{subject};
    matTrain.historyR = savestate_historyR{subject};
    %matTrain.training_pair = training_pair{subject};
    matTrain.Feature1 = Feature{1};
    matTrain.Feature2 = Feature{2};
    matTrain.Words = Words;
    
    matTrain.nObjects = nObjects;
    matTrain.maxTR = maxTR;
    matTrain.maxTRt = maxTRt;
    matTrain.t_max = t_max;
    matTrain.t_maxt = t_maxt;
    matTrain.visuals_On = visuals_On;
    matTrain.visuals_Off = visuals_Off;
    matTrain.word1_On = word1_On;
    matTrain.word1_Off = word1_Off;
    catch
        disp(['error writing train file for ',num2str(subject)])
    end
    
    if (mode == 1); disp ('Training Phase Complete - Test phase begins now'); end
    %% TEST TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %if (mode == 1), tts('Test Trials start now','Microsoft Hazel Desktop - English (Great Britain)'),end
    visuals_On = 0;
    visuals_Off = floor(1000/scale_factor);
    word1_On = 1; %250  
    word1_Off = floor(1000/scale_factor);
    %% Design Test Trials %%%%%%%%%%%%%%%%%%
    testSet=[];
    testSet(1,:)=[1 6 7 8]; testSet(2,:)=[2 5 7 8];
    testSet(3,:)=[3 5 6 8]; testSet(4,:)=[4 5 6 7];
    testSet(5,:)=[5 2 3 4]; testSet(6,:)=[6 1 3 4];
    testSet(7,:)=[7 1 2 4]; testSet(8,:)=[8 1 2 3];
    trial_order=randperm(size(testSet,1));
    targetPosition=[];
    %% RUN TEST TRIALS 
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyLWBt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRWBt{subject}=zeros(maxTRt,t_maxt);
    for trt=1:maxTRt
        inF=randperm(4); %1==inF(x) will be  target at loc x .. find(inF==1) get index indicating correct spt_LRM
        targetPosition(trt,:)=inF;
        pos =trial_order(trt);
        %sim2.init();  
        sim2.t =sim2.tZero;
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
            
            if isalmost (t, visuals_On, tolerance)
                for f=1:nFeatures
                    n=num2str(f);                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(testSet(pos,inF(1)))),cell2mat(Feature{f}(testSet(pos,inF(2)))), ...
                        cell2mat(Feature{f}(testSet(pos,inF(3)))),cell2mat(Feature{f}(testSet(pos,inF(4))))});                

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
            
            if isalmost (t, word1_On, tolerance)
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos)),wstrength});
            end
            if isalmost (t, word1_Off, tolerance)
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos)),0});
            end
            

            if isalmost (t, visuals_Off, tolerance)             
               for f=1:nFeatures
                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(testSet(pos,inF(1)))),cell2mat(Feature{f}(testSet(pos,inF(2)))), ...
                        cell2mat(Feature{f}(testSet(pos,inF(3)))),cell2mat(Feature{f}(testSet(pos,inF(4))))});                

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
            
            if isalmost(t, t_maxt,tolerance)
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
           else
                sim2.step();
           end
        end        
    end
    try
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
    matTest.targetPosition = targetPosition;
    
    matTest.visuals_On = visuals_On;
    matTest.visuals_Off = visuals_Off;
    matTest.word_On = word_On;
    matTest.word_Off = word_Off;
    catch
        disp(['error writing test file for ',num2str(subject)])
    end
end
