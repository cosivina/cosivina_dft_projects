%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

nObjects = 6; %%total number of objects
maxTR = 30; % total # of training trials
maxTRt = 2*nObjects; % total # of test trials
t_max = floor((4000+1000)/scale_factor); %specify simulation time  = scaled  Study(empirical) training trial Duration in millisecs + 1000 millisec break as gap between the trials   
t_maxt = floor((8000+1000)/scale_factor); %specify simulation time  = scaled  Study test trial Duration in millisecs + 1000 millisec break between the trials

History_Reset;
sim.init();
Gui_History_Reset;

for subject = 1:numSubjects % change to parfor for mode = 2 
%parfor subject = 1:numSubjects % change to parfor for mode = 2 
    
    sim2 = 5;
    if (mode == 1 )             
        sim2 = sim;
    elseif (mode == 0)
        sim2 = sim;
    elseif (mode == 2)
        sim2 = sim.copy();
    end
    
    visuals_On = 1;    visuals_Off = floor(4000/scale_factor); % visual display on off timing
    word1_On = floor([500, 2000]/scale_factor); word1_Off = floor([1500, 3000]/scale_factor); 
    
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[]; 
    Words=[]; 
    for j=1:nFeatures
        chosen_indices = randperm(size(Property{j},2),nObjects);
        TempFe=[];
        for i=1:nObjects
            TempFe = [TempFe  Property{j}(chosen_indices(i))];
            %TempFe = [TempFe  Property{j}(3*i)]; %diagonally aligned associations, use for an intuitive visualization  
        end
        Feature{j}=TempFe;
    end
    % Choose Words to associate on random
    chosen_indices = randperm(size(Names,2),nObjects);
    for i=1:nObjects
        Words = [Words  Names(chosen_indices(i))];
        %Words = [Words  Names(3*i)];    %diagonally aligned associations, use for an intuitive visualization 
    end
    %% DESIGN TRAINING TRIALS
    Object_pairs_named=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);
    Object_pairs_named_rev=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);
    county=1;
    for j=1:size(Feature{1},2)
        for k=1:size(Feature{1},2)
            if (cell2mat(Feature{1}(k))~= cell2mat(Feature{1}(j)))
                temp=[];
                for f=1:nFeatures
                  temp  =[temp cell2mat(Feature{f}(j)) cell2mat(Feature{f}(k))];
                end
                Object_pairs_named(county,:)= [temp cell2mat(Words(j)) cell2mat(Words(k)) 'P'];
                Object_pairs_named_rev(county,:)= [temp cell2mat(Words(k)) cell2mat(Words(j)) 'X'];
                county=county+1;
            end
        end
    end
    Paired_combs_all = [Object_pairs_named ; Object_pairs_named_rev]; %all possible Training tuples
    Sequence = randperm(size(Object_pairs_named,1));  % generate a random sequence of 30 numbers
   
    finSequence {subject}= zeros(maxTR,1);
    for tr=1:size(Object_pairs_named,1)  %shuffle the 30 tuples with the a random sequence
        finSequence {subject}(tr) = ((randi(2)-1)*size(Object_pairs_named,1))+ Sequence(tr); % Finds postion of a random tuple
    end
    
    % Handles to memory trace fields
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
    
    % Initialize variables to save the state of the hebbian fields 
    savestate_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
    savestate_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
    savestate_hcon_s{subject}=zeros(fieldSize_spt);
    savestate_hwm_s{subject}=zeros(fieldSize_spt);
    savestate_hword{subject}=zeros(fieldSize_wd);
    savestate_historyL{subject}=zeros(maxTR,t_max);
    savestate_historyR{subject}=zeros(maxTR,t_max);
    training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%
    
    %% RUN TRAINING TRIALS
    for tr=1:maxTR % trials loop
        sim2.t =sim2.tZero;
        pos = mod(finSequence{subject}(tr)-1,size(Object_pairs_named,1))+1;% sample trial tr with words in AB or BA order
        training_pair{subject}(tr,:)= Paired_combs_all(pos,:);
        
        if tr > 1 % Save the prior state of the hebbian fields
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
        % Run thorugh the simulation time of a trial
        while sim2.t <= t_max % time loop run

            t = sim2.t;

            if isalmost(t, visuals_On,tolerance) %present/display visual stimuli to the subject
                for f=1:nFeatures
                    cellNum=(2*f)-1;
                    
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                        strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY'}, ...
                        {(Paired_combs_all(pos,cellNum)),(Paired_combs_all(pos,cellNum+1))});

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
                    wrd1 = 2*nFeatures+www;
                    sim2.setElementParameters({['Word_Label_',num2str(www)]}, {'position','amplitude'},... 
                        {(Paired_combs_all(pos,wrd1)),wstrength});
                end
                if isalmost(t, word1_Off(www),tolerance)
                    wrd1 = 2*nFeatures+www;
                    sim2.setElementParameters({['Word_Label_',num2str(www)]}, {'position','amplitude'},... 
                        {0,0});
                end
            end

            if isalmost(t, visuals_Off,tolerance) %end visual stimuli 

                for f=1:nFeatures
                    cellNum=(2*f)-1;

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
            
            % Backup looking data and hebbian memory activity states
            if isalmost(t, t_max,tolerance) % save loking data and memory activity states
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

           if (mode == 1) && (mod(t,gui_speed)==0)% Gui and pause control
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
        % Training trials ended
    %Save the looking data, memory etc so far
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
        matTrain.training_pair = training_pair{subject};
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
        disp('Error saving a train file');
    end
    
    
%% TEST TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %if (mode == 1), tts('Test Trials start now','Microsoft Hazel Desktop - English (Great Britain)'),end
%% SY tests
    visuals_On = 1; visuals_Off = floor(8000/scale_factor);
    
    %SY2008 - CSWL word sequence
    word_On = floor([500 2000 4500 6000]/scale_factor); word_Off = floor([1500 3000 5500 7000]/scale_factor);%XSIT
    
    % SY2013 XTRAP word sequence
    %word_On =  floor([1 1800 3500 5200 6900]/scale_factor); word_Off = floor((745+[1 1800 3500 5200 6900])/scale_factor); %XTRAP   
    
%% 1 sec test   
    %visuals_On = 1;visuals_Off = floor(1000/scale_factor);  
    %word_On = floor([1]/scale_factor); word_Off =floor([1000]/scale_factor); 


%% Design Test Trials %%%%%%%%%%%%%%%%%%
    Test_trial_set=[]; %%blank ordered test set 12 trials for two test sets 6+6..
    % using rand(5) as index starts from 1, so 1+5 maxiumm=6
    tst1 = randi(nObjects-1); % a random shift to one object index to generate its distractor
    tst2 = randi(nObjects-1); % shift Left object to generate combins with Right
   while (~prod(tst1-tst2))
        tst2 = randi(nObjects-1); 
   end
    
    for i=1:nObjects
        if (i+tst1 > nObjects) % shift the object index to generate its distractor
            indnew = i+tst1-nObjects;
        else
            indnew = i+tst1;
        end
        TempFe=[];
        if rand >0.5
            for j=1:nFeatures  %  
                TempFe = [TempFe  cell2mat(Feature{j}(i)) cell2mat(Feature{j}(indnew))]; %add features of both
            end
            Test_trial_set(i,:) = [TempFe cell2mat(Words(i)) 'L'];
        else
            for j=1:nFeatures  % 
                TempFe = [TempFe  cell2mat(Feature{j}(indnew)) cell2mat(Feature{j}(i))]; %add features of both
            end
            Test_trial_set(i,:) = [TempFe cell2mat(Words(i)) 'R'];
        end

        if (i+tst2 > nObjects)  %  shift the object index to generate the corresp distractor
            indnew = i+tst2-nObjects;
        else
            indnew = i+tst2;
        end
        TempFe=[];
        if rand >0.5
            for j=1:nFeatures % creat tuple Right object with randomly chosen Left object 
                TempFe = [TempFe  cell2mat(Feature{j}(indnew)) cell2mat(Feature{j}(i))]; %add features of both
            end
            Test_trial_set(nObjects+i,:) = [TempFe cell2mat(Words(i)) 'R'];
        else
            for j=1:nFeatures % creat tuple Right object with randomly chosen Left object 
                TempFe = [TempFe  cell2mat(Feature{j}(i)) cell2mat(Feature{j}(indnew))]; %add features of both
            end
            Test_trial_set(nObjects+i,:) = [TempFe cell2mat(Words(i)) 'L'];            
        end
    end
    trial_order=[randperm(nObjects) nObjects+randperm(nObjects)];
    
    %% RUN TEST TRIALS
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    test_pair{subject}=zeros(2*nObjects,2*nFeatures+2);
    for trt=1:maxTRt
        pos=trial_order(trt);
        test_pair{subject}(trt,:)= Test_trial_set(pos,:); %test_pair{subject}(trt)= char(Test_trial_set(pos,2*nFeatures+2))
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
        
        % Run thorugh the simulation time of a test trial  
        while sim2.t <= t_maxt

            t = sim2.t;
            if isalmost(t, visuals_On,tolerance)

                for f=1:nFeatures
                    cellNum=(2*f)-1;

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                        strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY'}, ...
                        {(Test_trial_set(pos,cellNum)),(Test_trial_set(pos,cellNum+1))});                

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set spatial dimension value
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'positionX', 'positionX'}, {spt_L,spt_R});

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set left right stimuli up
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'amplitude', 'amplitude'}, vstrength * ones(1, 2)); 
                end
            end
            if isalmost(t, word_On,tolerance)
                wrd1 = 2*nFeatures+1;
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {(Test_trial_set(pos,wrd1)),wstrength});
            end
            if isalmost(t, word_Off,tolerance) 
                 wrd1 = 2*nFeatures+1;
                 sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {0,0});
            end

            if isalmost(t, visuals_Off,tolerance)
                 for f=1:nFeatures
                    cellNum=(2*f)-1;
                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                        strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY'}, {0,0});                                      

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set spatial dimension value
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'positionX', 'positionX'}, {0,0});

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set left right stimuli up
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'amplitude', 'amplitude'}, 0 * ones(1, 2));
                 end
            end
            
            % Backup looking data and hebbian memory activity states
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

           if (mode == 1) && (mod(t,1)==0)
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

        matTest.hcon_ft = savestate_hcon_f{subject};
        matTest.hwm_ft = savestate_hwm_f{subject};
        matTest.hwm_ct = savestate_hwm_c{subject};
        matTest.hwft = savestate_hwf{subject};
        matTest.hcon_st = savestate_hcon_s{subject};
        matTest.hwm_st = savestate_hwm_s{subject};
        matTest.hwordt = savestate_hword{subject};

        matTest.historyLt = savestate_historyLt{subject};
        matTest.historyRt = savestate_historyRt{subject};
        matTest.test_pair = test_pair{subject};

        matTest.visuals_On = visuals_On;
        matTest.visuals_Off = visuals_Off;
        matTest.word_On = word_On;
        matTest.word_Off = word_Off;
    catch
        disp('Error saving a test file');
    end
end