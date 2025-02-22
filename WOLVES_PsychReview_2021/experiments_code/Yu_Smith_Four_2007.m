%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

Prediction_On = false;% set to true for WOLVES prediction simulation
if Prediction_On == true
    scale_factor = scale_factor*3;
end
nObjects = 18; %%total number of objects
maxTR = 27; % total # of training trials
maxTRt = nObjects; % total # of test trials
t_max = floor((12000+1000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
t_maxt = floor((1000+1000)/scale_factor); %sAssumed % scale with Experiment test trial Duration

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
    visuals_Off = floor((12000)/(scale_factor));
    word1_On  = floor([1000 4000 7000 10000]/scale_factor);  
    word1_Off = floor([2000 5000 8000 11000]/scale_factor);
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
   
    Object_indA=1+mod(0:nObjects-1,nObjects);% A B C D...
    Object_indB=1+mod(3:nObjects+2,nObjects);  % B C D E...
    Object_indD=1+mod(5:nObjects+4,nObjects);% D E F G...
    Object_indH=1+mod(14:nObjects+13,nObjects);% H I J K...

    Indices_Set=[Object_indA' Object_indB' Object_indD' Object_indH'];
    x1=randperm(9);
    x2=randperm(9);while prod(abs(x1-x2))<1;x2=randperm(9);end 
    x3=randperm(9);
    x4=randperm(9);while prod(abs(x3-x4))<1;x4=randperm(9);end
    x3=x3+9; x4=x4+9; v=[x1' x2' x3' x4'];
    Indices_Set = [Indices_Set; v];
    Sequence = randperm(size(Indices_Set,1)); %generate a random permuted sequence of 30 tuples
    
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
    %training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%
    %% RUN TRAINING TRIALS
    for tr=1:maxTR
        sim2.init();
        pos = Sequence(tr); % Finds postion of a random tuple
        obj=randperm(4);
        wrd=randperm(4);
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
            if t == visuals_On                
                for f=1:nFeatures

                    n=num2str(f);                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(Indices_Set(pos,obj(1)))),cell2mat(Feature{f}(Indices_Set(pos,obj(2)))), ...
                        cell2mat(Feature{f}(Indices_Set(pos,obj(3)))),cell2mat(Feature{f}(Indices_Set(pos,obj(4))))});                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX'}, ...
                         {spt_L, spt_Lm, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
                          (vstrength) * ones(1, 4) );                     
                end
            end
            
             for www=1:length(word1_On)  
                
                if isalmost(t, word1_On(www),tolerance) 
                    sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                        {cell2mat(Words(Indices_Set(pos,wrd(www)))),wstrength});
                end
                if isalmost(t, word1_Off(www),tolerance)
                    sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                        {cell2mat(Words(Indices_Set(pos,wrd(www)))),0});
                end
             end
            
            if t == visuals_Off
                for f=1:nFeatures
                    cellNum=(2*f)-1;

                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {0, 0, 0, 0});                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX'}, ...
                         {spt_L, spt_Lm, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
                         { 0, 0, 0, 0});
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
                 savestate_historyR{subject}(tr,:) = flipud(sim2.getComponent('history lookR', 'output'));
                 savestate_historyLWB{subject}(tr,:) = flipud(sim2.getComponent('history lookLWB', 'output'));
                 savestate_historyRWB{subject}(tr,:) = flipud(sim2.getComponent('history lookRWB', 'output'));

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
    end; % Trial maxTR loop
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
    matTrain.historyLWB = savestate_historyLWB{subject};
    matTrain.historyRWB = savestate_historyRWB{subject};
    %matTrain.training_pair = training_pair{subject};
    matTrain.Feature1 = Feature{1};
    matTrain.Feature2 = Feature{2};
    matTrain.Words = Words;
    
    matTrain.Prediction_On = Prediction_On;
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
        disp('Unable to write to train MAT-file for subject number ');
        %continue;
    end;
   %% TEST TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    visuals_On = 0;
    visuals_Off = floor((1000)/scale_factor);
    word1_On  =1;% assumed
    word1_Off  = floor((1000)/scale_factor);
    Testseq=[];
    wrdList=randperm(nObjects);
    for tseq=1:maxTRt
     % for tseq_th word   %remove target to generate distractor list 
     disSet=wrdList(wrdList~=wrdList(tseq));
     dis3=disSet(randperm(length(disSet),3)); %pick 3 distractors
     indxr=randperm(4);
     for ds=1:3 Testseq(tseq,indxr(ds))=dis3(ds);end
     Testseq(tseq,indxr(4))=wrdList(tseq);%object
     Testseq(tseq,5)=indxr(4);%location of target   
    end
    Testset{subject}=Testseq;
    %% RUN TEST TRIALS
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyLWBt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRWBt{subject}=zeros(maxTRt,t_maxt);
    for trt=1:maxTRt
              
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
            
            if t== word1_On
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(wrdList(trt))),wstrength});
            end
            if t == visuals_On
               
                for f=1:nFeatures
                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(Testseq(trt,1))),cell2mat(Feature{f}(Testseq(trt,2))), ...
                        cell2mat(Feature{f}(Testseq(trt,3))),cell2mat(Feature{f}(Testseq(trt,4)))});                

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
           if t== word1_Off
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(wrdList(trt))),0});
            end

            if t == visuals_Off
               for f=1:nFeatures
                    n=num2str(f);
                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY'}, ...
                        0 * ones(1, 4));                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX'}, ...
                         {spt_L, spt_Lm, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
                         { 0, 0, 0, 0}); 
                    
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
    end;
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
    matTest.Testset = Testset{subject};
    
    matTest.visuals_On = visuals_On;
    matTest.visuals_Off = visuals_Off;
    matTest.word1_On = word1_On;
    matTest.word1_Off = word1_Off;
    catch
        disp('Unable to write to test MAT-file for subject number ');
        %continue;
    end;
end
          