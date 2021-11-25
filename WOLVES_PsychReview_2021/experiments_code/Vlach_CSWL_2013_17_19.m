%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

nObjects = 12; %%total number of objects
maxTR = 36; % total # of training trials
maxTRt = nObjects; % total # of test trials
t_max = floor((4000+1000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
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
    
    visuals_On = 0;  visuals_Off = floor(4000/scale_factor);
    word1_On = floor([500, 2500]/scale_factor); word1_Off = floor([1500, 3500]/scale_factor); 
    
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
    Paired_combs_all = zeros(maxTR, 7);
    Sequence = 1:12;% randperm(nObjects); %Of this Sequence Use first 6 as Massed, Last 6 as Interleaved
    
    ranD= mod(randi(10),2); % for even odd location switching only
    for j=1:6
        rndx=randperm(12);
         for k=1:6
             indx= ((j-1)*6)+k;
             temp=[];
             if rndx(k) <= 2     % if rand smaller, we alter normal switiching locations to same-as-previous-loc            
                if mod(ranD,2) %normal switching of locations
                    for f=1:nFeatures
                        temp  =[temp cell2mat(Feature{f}(Sequence(j))) cell2mat(Feature{f}(Sequence(6+k)))];
                    end
                else
                     for f=1:nFeatures
                        temp  =[temp cell2mat(Feature{f}(Sequence(6+k))) cell2mat(Feature{f}(Sequence(j)))];
                     end 
                end
                ranD=ranD+2;
             else
                 if ~ mod(ranD,2) %normal switching of locations
                    for f=1:nFeatures
                        temp  =[temp cell2mat(Feature{f}(Sequence(j))) cell2mat(Feature{f}(Sequence(6+k)))];
                    end
                else
                     for f=1:nFeatures
                        temp  =[temp cell2mat(Feature{f}(Sequence(6+k))) cell2mat(Feature{f}(Sequence(j)))];
                     end
                 end
                 ranD=ranD+1;
             end
             if rndx(k) <= 6
                temp = [temp cell2mat(Words(Sequence(j))) cell2mat(Words(Sequence(6+k))) 'P'];
             else
                temp = [temp cell2mat(Words(Sequence(6+k))) cell2mat(Words(Sequence(j))) 'X'];  
             end
             Paired_combs_all(indx,:)= temp;             
         end
    end
    %%
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
    training_pair{subject}=zeros(36,(2*nFeatures)+3);%
    %% RUN TRAINING TRIALS
    for tr=1:maxTR

        sim2.t =sim2.tZero;
        pos=tr;
        training_pair{subject}(tr,:)= Paired_combs_all(pos,:);
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
             if isalmost(t, visuals_On,tolerance)
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

            if isalmost(t, visuals_Off,tolerance)

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
            
            if isalmost (t, t_max,tolerance)
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
                hL=sim2.getComponent('history lookL', 'output');
                hR=sim2.getComponent('history lookR', 'output');
                savestate_historyL{subject}(tr,:) = flipud(hL);
                savestate_historyR{subject}(tr,:) =flipud(hR);
            end

            if (mode == 1) && (mod(t,20)==0)
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
        matTrain.training_pair = training_pair{subject};
        matTrain.Feature1 = Feature{1};
        matTrain.Feature2 = Feature{2};
        matTrain.Words = Words;
        matTrain.Sequence = Sequence;

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
    %% Test Trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    if (mode == 1), tts('Test Trials start now','Microsoft Hazel Desktop - English (Great Britain)'),end

%% Vlach 8-sec Test
%     visuals_On = 1; %     visuals_Off = floor(8000/scale_factor);
%     word_On = floor([500 2500 4500 6500]/scale_factor);   
%     word_Off = floor([1500 3500 5500 7500]/scale_factor);

%% I-sec Test
   visuals_On = 1;    visuals_Off = floor(1000/scale_factor);
   word_On = 1;    word_Off = floor(1000/scale_factor);
   
    %% DESIGN TEST TRIALS %     
    Test_trial_set=[];
    county=1;
    for i=1:3
        TempFe = [];
        if rand > 0.5
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i))) cell2mat(Feature{j}(Sequence(i+3)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i))) 'L']; county=county+1;
        else
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i+3))) cell2mat(Feature{j}(Sequence(i)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i))) 'R']; county=county+1;
        end

        TempFe = [];
        if rand > 0.5
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i+3))) cell2mat(Feature{j}(Sequence(i+6)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i+3))) 'L']; county=county+1;
        else
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i+6))) cell2mat(Feature{j}(Sequence(i+3)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i+3))) 'R']; county=county+1;
        end

        TempFe = [];
        if rand > 0.5
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i+6))) cell2mat(Feature{j}(Sequence(i+9)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i+6))) 'L']; county=county+1;
        else
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i+9))) cell2mat(Feature{j}(Sequence(i+6)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i+6))) 'R']; county=county+1;
        end

        TempFe = [];
        if rand > 0.5
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i+9))) cell2mat(Feature{j}(Sequence(i)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i+9))) 'L']; county=county+1;
        else
            for j=1:nFeatures  
                TempFe = [TempFe  cell2mat(Feature{j}(Sequence(i))) cell2mat(Feature{j}(Sequence(i+9)))]; %add features of both
            end           
            Test_trial_set(county,:)=[TempFe cell2mat(Words(Sequence(i+9))) 'R']; county=county+1;
        end     
    end
    trial_order=randperm(12);
    %% RUN TEST TRIALS
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    test_pair{subject}=zeros(nObjects,2*nFeatures+2);
    for trt=1:maxTRt
        pos=trial_order(trt);
        test_pair{subject}(trt,:)= Test_trial_set(pos,:);
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

            if (mode == 1)
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
    if (mode == 1)
         gui.close();
    end