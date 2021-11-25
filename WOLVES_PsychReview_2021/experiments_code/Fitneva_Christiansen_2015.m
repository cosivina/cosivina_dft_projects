%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

nObjects = 10; %%total number of objects
maxFam = nObjects; % total number of familiarization trials
maxTR = 15; % 3 blocks of 5 trials each, total # of training trials
maxTRt = 15; % 3 blocks of 5 trials each total # of test trials
t_maxF = floor((3000+1000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
t_max = floor((4500+1000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
t_maxt = floor((1000+1000)/scale_factor); %specify simulation time  % scale with Experiment test trial Duration

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
    visuals_Off = floor(3000/scale_factor);
    words_On = floor(1000/scale_factor);
    words_Off = floor(2000/scale_factor);
    
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[]; 
    chosen_indices = randperm(size(Property{1},2),nObjects);%% Create objects
    for j=1:nFeatures     
        TempFe=[];
        for i=1:nObjects
            TempFe = [TempFe  Property{j}(chosen_indices(i))];
            %TempFe = [TempFe  Property{j}(i)];
        end
        Feature{j}=TempFe;
    end
    Words=[]; 
    chosen_indices = randperm(size(Names,2), nObjects);%% Create  pseudowords
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
    %savestate_historyL{subject}=zeros(maxFam,t_maxF);
    %savestate_historyR{subject}=zeros(maxFam,t_maxF);
    %training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%
    
    %% RUN FAMILIARIZATION TRIALS
    for tr=1:maxFam
        %sim2.init();
        sim2.t =sim2.tZero;
        targ = tr; % Finds position of a random tuple
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
        while sim2.t <= t_maxF

            t = sim2.t;

            % STEP 1 I have skipped the attention getter
            if t == visuals_On
                for f=1:nFeatures
                     sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                    {cell2mat(Feature{f}(targ)),spt_C, vstrength});                    
                end
            end
            if t == words_On            
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(targ)),wstrength});               
            end
            if t == words_Off            
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(targ)),0});               
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
%               savestate_historyLt{subject}(tr,:) = sim2.getComponent('history lookLt', 'output');
%               savestate_historyRt{subject}(tr,:) = sim2.getComponent('history lookRt', 'output');

                for f=1:nFeatures
                    sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                    {cell2mat(Feature{f}(targ)),spt_C, 0});    
                end
            end
            if (mode == 1)
                if ~gui.pauseSimulation
                    sim2.step();
                end
                if gui.quitSimulation
                    gui.close();
                    break;
                end
                if (mod(t,15)==0)
                    gui.step();
                end
                %gui.checkAndUpdateControls();
                %gui.updateVisualizations();                  
           else
                sim2.step();
           end
        end 
    end
%     OutName = [simName,num2str(subject),'_fam.mat'];
%     matFam=matfile(OutName,'writable',true);
% 
%     matFam.hcon_f = savestate_hcon_f{subject};
%     matFam.hwm_f = savestate_hwm_f{subject};
%     matFam.hwm_c = savestate_hwm_c{subject};
%     matFam.hwf = savestate_hwf{subject};
%     matFam.hcon_s = savestate_hcon_s{subject};
%     matFam.hwm_s = savestate_hwm_s{subject};
%     matFam.hword = savestate_hword{subject};
%   matFam.historyL = savestate_historyL{subject};
%   matFam.historyR = savestate_historyR{subject};
%   matFam.training_pair = training_pair{subject};
%     matFam.Feature1 = Feature{1};
%     matFam.Feature2 = Feature{2};
%     matFam.Words = Words;
    
    %% LEARNING PHASE
    visuals_On = 0;
    visuals_Off = floor(4500/scale_factor);
    word1_On  = floor(1000/scale_factor); 
    word1_Off = floor(2000/scale_factor);
    word2_On  = floor(2750/scale_factor);
    word2_Off = floor(3750/scale_factor);
    
    %% DESIGN TRAINING TRIALS
     Words_Fam = Words; 
    if mod(subject,2)== 0        
    %% Condition HIGH IA: from familiarized 10 objs.. reshuffle 4 pairings..
        nShuffle=4;%must be even %Take any four pairings to shuffle
        Sequence = randperm(nObjects, nShuffle); 
        seqF= 1:nShuffle;
        seqW= randperm(nShuffle);
        while (~prod(seqF-seqW))
            seqW=randperm(nShuffle); 
        end 
        for i=1:nShuffle
            Words(Sequence(i)) = Words_Fam(Sequence(seqW(i))); %shuffle words 
        end
    elseif mod(subject,2)== 1
    %% Condition LOW IA: from familiarized 10 objs.. reshuffle 6 pairings..
        nShuffle=6;
        Sequence = randperm(nObjects, nShuffle); %Take any six pairings to shuffle
        seqF= 1:nShuffle;
        seqW= randperm(nShuffle);
        while (~prod(seqF-seqW))
            seqW=randperm(nShuffle); 
        end
        for i=1:nShuffle
            Words(Sequence(i)) = Words_Fam(Sequence(seqW(i))); %shuffle words 
        end
    end    
   %create 15 train trials, 3 blocks of, 5 trials each... 
    seq0=randperm(nObjects); % plan a rand seq of 10 objects
    block1=randperm(nObjects/2); %1 ordr for first 5 objs, last 5 fixed
    block2=randperm(nObjects/2);%2nd order
    while (~prod(block1-block2))
        block2=randperm(nObjects/2);
    end
    block3=randperm(nObjects/2); %3rd ordr
    while (~prod(block1-block3) || ~prod(block2-block3))
        block3= randperm(nObjects/2); 
    end 
    blocks=[block1; block2; block3];
    object_pairs=[];
    for blok=1:3
        for i=1:nObjects/2
            temp=[];
            if (rand>0.5)
                temp = [seq0(i+nObjects/2) seq0(blocks(blok,i))];
            else
                temp = [seq0(blocks(blok,i)) seq0(i+nObjects/2)];
            end
            if (rand>0.5)
                temp = [temp seq0(i+nObjects/2) seq0(blocks(blok,i))];
            else
                temp = [temp seq0(blocks(blok,i)) seq0(i+nObjects/2)];
            end
            object_pairs = [object_pairs; temp];
        end
     end
          
    %% RUN TRAINING TRIALS
    savestate_historyL{subject}=zeros(maxTR,t_max);
    savestate_historyR{subject}=zeros(maxTR,t_max);
    for tr=1:maxTR
        %sim2.init();
        sim2.t =sim2.tZero;
        pos = tr; % Finds postion of a random tuple
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

        while sim2.t <= t_max

            t = sim2.t;

            % STEP 1 Skipped attention getters and delay follwing them
            if t == visuals_On
                for f=1:nFeatures               

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                        strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(object_pairs(pos,1))),cell2mat(Feature{f}(object_pairs(pos,2)))});                

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set spatial dimension value
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'positionX', 'positionX'}, {spt_L,spt_R});

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set left right stimuli up
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'amplitude', 'amplitude'}, vstrength * ones(1, 2));                     
                end
            end
            if t == word1_On 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(object_pairs(pos,3))),wstrength});
            end
            if t == word1_Off
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(object_pairs(pos,3))),0});
            end
            if t == word2_On
                sim2.setElementParameters({'Word_Label_2'}, {'position','amplitude'},... 
                    {cell2mat(Words(object_pairs(pos,4))),wstrength});
            end
            if t == word2_Off
                sim2.setElementParameters({'Word_Label_2'}, {'position','amplitude'},... 
                    {cell2mat(Words(object_pairs(pos,4))),0});
            end

            if t == visuals_Off
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

           if (mode == 1)
                if ~gui.pauseSimulation
                    sim2.step();
                end
                if gui.quitSimulation
                    gui.close();
                    break;
                end
                if (mod(t,15)==0)
                    gui.step();
                end
                %gui.checkAndUpdateControls();
                %gui.updateVisualizations();                  
           else
                sim2.step();
           end
        end % Time sim.t
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
    matTrain.historyL = savestate_historyL{subject};
    matTrain.historyR = savestate_historyR{subject};
%   matTrain.training_pair = training_pair{subject};
    matTrain.Feature1 = Feature{1};
    matTrain.Feature2 = Feature{2};
    matTrain.Words = Words;
    matTrain.Words_Fam = Words_Fam;
%     
    %% TEST PHASE
    %if (mode == 1), tts('Test Trials start now','Microsoft Hazel Desktop - English (Great Britain)'),end
    visuals_On = 0;
    visuals_Off = floor(1000/scale_factor);
    word1_On  = 1; %750 in experiment CHANGED!
    word1_Off = floor(1000/scale_factor);
    if (mode == 1); disp ('Training Phase Complete - Test phase begins now'); end
    %% DESIGN TEST TRIALS 
    seq0=randperm(nObjects); % plan a rand seq of 10 objects
    block1=randperm(nObjects/2); %1 ordr for foil 5 objs, last 5 fixed
    block2=randperm(nObjects/2);%2nd order
    while (~prod(block1-block2))
        block2=randperm(nObjects/2);
    end
    block3=randperm(nObjects/2); %3rd ordr
    while (~prod(block1-block3) || ~prod(block2-block3))
        block3= randperm(nObjects/2); 
    end
    blocks=[block1; block2; block3];
    test_pairs=[];
    for blok=1:3
        for i=1:nObjects/2
            temp=[];
            if (rand>0.5)
                temp = [seq0(i+nObjects/2) seq0(blocks(blok,i)) seq0(i+nObjects/2) double('L')];                
            else
                temp = [seq0(blocks(blok,i)) seq0(i+nObjects/2) seq0(i+nObjects/2)  double('R')];
            end
            
            test_pairs = [test_pairs; temp]; 
        end
    end
      
    %% RUN TEST TRIALS 
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);   
    test_vector{subject}=test_pairs;
    for trt=1:maxTRt
        %sim2.init();
        sim2.t =sim2.tZero;
        pos = trt; % Finds postion of a random tuple
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

            % STEP 1 Skipped attention getters and delay follwing them
            if t == visuals_On
                for f=1:nFeatures               

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                        strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(test_pairs(pos,1))),cell2mat(Feature{f}(test_pairs(pos,2)))});                

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set spatial dimension value
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'positionX', 'positionX'}, {spt_L,spt_R});

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set left right stimuli up
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'amplitude', 'amplitude'}, vstrength * ones(1, 2));                     
                end
            end
            if t == word1_On 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(test_pairs(pos,3))),wstrength});
            end
            if t == word1_Off
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(test_pairs(pos,3))),0});
            end

            if t == visuals_Off
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
            if isalmost(t, t_maxt,tolerance)
%                 for i = 1 : nFeatures
%                     n=num2str(i);
%                     savestate_hcon_f{subject}(i,:) = sim2.getComponent(['hcon_f' n], 'output');
%                     savestate_hwm_f{subject}(i,:) = sim2.getComponent(['hwm_f' n], 'output');
%                     savestate_hwf{subject}(i,:,:) = sim2.getComponent(['hwf' n], 'output');
%                     savestate_hwm_c{subject}(i,:,:) = sim2.getComponent(['hwm_c' n], 'output');
%                 end
%                 savestate_hcon_s{subject} = sim2.getComponent('hcon_s', 'output');
%                 savestate_hwm_s{subject} = sim2.getComponent('hwm_s', 'output');
%                 savestate_hword{subject} = sim2.getComponent('hword', 'output');
                savestate_historyLt{subject}(trt,:) = flipud(sim2.getComponent('history lookLt', 'output'));
                savestate_historyRt{subject}(trt,:) = flipud(sim2.getComponent('history lookRt', 'output'));
            end

           if (mode == 1)
                if ~gui.pauseSimulation
                    sim2.step();
                end
                if gui.quitSimulation
                    gui.close();
                    break;
                end
                if (mod(t,15)==0)
                    gui.step();
                end
                %gui.checkAndUpdateControls();
                %gui.updateVisualizations();                  
           else
                sim2.step();
           end
        end % Time sim.t
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
    matTest.test_vector = test_vector(subject);
%     if (mode == 1)
%          gui.close();
%     end
end