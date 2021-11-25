%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

nObjects = 12; %%total number of objects
maxTR = 12; % 12 in exp actually total # of training trials
maxTRt = 12; % total # of test trials
t_max = floor((3000+1000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
t_maxt = floor((1000+1000)/scale_factor); %specify simulation time  % scale with Experiment test trial Duration

History_Reset;
sim.init();
Gui_History_Reset;

%% WOB MEMORY TASK (Word Object Binding Task)
%% TRAINING PHASE %skipped on a different set of objects/words

%% LEARNING PHASE
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
    visuals_Off = floor(3000/scale_factor);
    words_On = floor([1]/scale_factor);%[;
    words_Off = floor([3000]/scale_factor);
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[]; 
    Distra=[];
    for j=1:nFeatures
        chosen_indices = randperm(size(Property{1},2));%% Create all 18 objects
        TempFe=[];
        for i=1:nObjects
            TempFe = [TempFe  Property{j}(chosen_indices(i))];
            %TempFe = [TempFe  Property{j}(i)];
        end
        Feature{j}=TempFe;
        
        %adding distractor here only
        TempDe=[];
        for i=1:nObjects
            TempDe = [TempDe  Property{j}(chosen_indices(i+nObjects-6))];
            %TempDe = [TempDe  Property{j}(i+nObjects-6)];
        end
        Distra{j}=TempDe;
    end
    Words=[]; 
    chosen_indices = randperm(size(Names,2), nObjects);%% Create  12 words for learning
    for i=1:nObjects
        Words = [Words  Names(chosen_indices(i))];
        %Words = [Words  Names(i)];
    end
    
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
    savestate_historyC{subject}=zeros(maxTR,t_max);
    %% RUN TRAINING TRIALS
     for tr=1:maxTR

        sim2.t =sim2.tZero;
        pos = tr; % Finds postion of a random tuple
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

            if isalmost(t, visuals_On,tolerance)
                for f=1:nFeatures
                     sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                    {cell2mat(Feature{f}(pos)),spt_C, vstrength});                    
                end
            end
            if isalmost(t, words_On,tolerance)        
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos)),wstrength});               
            end
            if isalmost(t, words_Off,tolerance)           
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {0,0});               
            end
            if isalmost(t, visuals_Off,tolerance)   
                 for f=1:nFeatures
                    sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                    {0, 0, 0});    
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
                savestate_historyC{subject}(tr,:) = flipud(sim2.getComponent('history lookC', 'output'));
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

    matTrain.historyC = savestate_historyC{subject};
    matTrain.Feature1 = Feature{1};
    matTrain.Feature2 = Feature{2};
    matTrain.Words = Words;

%% TESTING PHASE
    visuals_On = 0;
    visuals_Off = floor(1000/scale_factor);
    words_On = 1;
    words_Off = floor(1000/scale_factor);%%NOT defined
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    test_word{subject}=zeros(maxTRt,1);

    for trt=1:maxTRt
        pos=trt;

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
            %
            if isalmost(t, visuals_On,tolerance) 
                if rand > 0.5
                    for f=1:nFeatures                   
                         sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(pos)),spt_L, vstrength});   
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Right')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Distra{f}(pos)),spt_R, vstrength}); 
                    end
                    test_word{subject}(trt,1)=  'L';
                else
                     for f=1:nFeatures
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(pos)),spt_R, vstrength});   
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Right')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Distra{f}(pos)),spt_L, vstrength});
                     end
                     test_word{subject}(trt,1)=  'R';
                 end
             end
            
            if isalmost(t, words_On,tolerance)                 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos)),wstrength});               
            end
            if isalmost(t, words_Off,tolerance)              
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos)),0}); 
            end    
            if isalmost(t, visuals_Off,tolerance) 
                 for f=1:nFeatures
                    sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                    {0,0, 0});   
                    sim2.setElementParameters({strcat('Feature_',num2str(f),'_Right')}, {'positionY','positionX','amplitude'},... 
                    {0,0, 0});     
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
    matTest.test_word = test_word{subject};
end