%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

Prediction_On = false;% set to true for WOLVES prediction simulation
if Prediction_On == true
    scale_factor = scale_factor*3;
end
nObjects = 18; %%total number of objects
maxTR = 54; % total # of training trials
maxTRt = nObjects; % total # of test trials
t_max = floor((6000+1000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   %6000
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
    visuals_Off=floor(6000/(scale_factor));
    word1_On  = floor(1000/(scale_factor));% word length values assumed, i could not find in the article  
    word1_Off = floor(2000/(scale_factor));%
    word2_On  = floor(4000/(scale_factor));%
    word2_Off = floor(5000/(scale_factor));%
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
   %%% pair first 9 objects (in differnt orders) with the last 9 objects
    Object_ind_1=1+mod(0:(nObjects)/2-1,(nObjects)/2);% A B C D E F G H I
    Object_ind_2=1+mod(1:(nObjects)/2+0,(nObjects)/2);  % B C D E F G H I A
    Object_ind_3=1+mod(2:(nObjects)/2+1,(nObjects)/2);% C D E F...
    Object_ind_4=1+mod(3:(nObjects)/2+2,(nObjects)/2);% D E F G...
    Object_ind_5=1+mod(4:(nObjects)/2+3,(nObjects)/2);
    Object_ind_6=1+mod(5:(nObjects)/2+4,(nObjects)/2);
    Object_indX=(nObjects)/2+1:nObjects;%              % J K L M N O P Q R
    Indices_Set=[Object_ind_1' Object_indX'; Object_ind_2' Object_indX'; Object_ind_3' Object_indX'; ...
        Object_ind_4' Object_indX'; Object_ind_5' Object_indX'; Object_ind_6' Object_indX';];
    Sequence = randperm(size(Indices_Set,1)); %generate a random permuted sequence of 54 tuples
    
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
    training_pair{subject}=zeros(maxTR,(2*nFeatures)+3);
    %% RUN TRAINING TRIALS
    for tr=1:maxTR
        %sim2.init();
        sim2.t =sim2.tZero;
        pos = Sequence(tr); % Finds postion of a random tuple
        obj=randperm(2);
        wrd=randperm(2);
        if obj==wrd; px='P'; else; px='X'; end 
        training_pair{subject}(tr,:)=[cell2mat(Feature{1}(Indices_Set(pos,obj(1)))),cell2mat(Feature{1}(Indices_Set(pos,obj(2)))),...
            cell2mat(Feature{2}(Indices_Set(pos,obj(1)))),cell2mat(Feature{2}(Indices_Set(pos,obj(2)))),...
            cell2mat(Words(Indices_Set(pos,wrd(1)))), cell2mat(Words(Indices_Set(pos,wrd(2)))), px];
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
            if  isalmost(t, visuals_On,tolerance)              
                for f=1:nFeatures

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                        strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY'}, ...
                        {cell2mat(Feature{f}(Indices_Set(pos,obj(1)))),cell2mat(Feature{f}(Indices_Set(pos,obj(2))))});                

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set spatial dimension value
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'positionX', 'positionX'}, {spt_L,spt_R});

                    sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set left right stimuli up
                        strcat('Feature_',num2str(f),'_Right')}, ...
                         {'amplitude', 'amplitude'}, vstrength * ones(1, 2));                     
                end
            end
            if isalmost(t, word1_On,tolerance)
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(Indices_Set(pos,wrd(1)))),wstrength});
            end
            if isalmost(t, word1_Off,tolerance)
                wrd1 = 2*nFeatures+1;
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(Indices_Set(pos,wrd(1)))),0});
            end
            if isalmost(t, word2_On,tolerance)
                wrd2 = 2*nFeatures+2;
                sim2.setElementParameters({'Word_Label_2'}, {'position','amplitude'},... 
                    {cell2mat(Words(Indices_Set(pos,wrd(2)))),wstrength});
            end
            if isalmost(t, word2_Off,tolerance)
                wrd2 = 2*nFeatures+2;
                sim2.setElementParameters({'Word_Label_2'}, {'position','amplitude'},... 
                    {cell2mat(Words(Indices_Set(pos,wrd(2)))),0});
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
            end
            

           if (mode == 1) && (mod(t,gui_speed)==0)
                %if tr==50; waitforbuttonpress;end
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
    
    matTrain.Prediction_On = Prediction_On;
    matTrain.nObjects = nObjects;
    matTrain.maxTR = maxTR;
    matTrain.maxTRt = maxTRt;
    matTrain.t_max = t_max;
    matTrain.t_maxt = t_maxt;
    matTrain.visuals_On = visuals_On;
    matTrain.visuals_Off = visuals_Off
    matTrain.word1_On = word1_On;
    matTrain.word1_Off = word1_Off; 
    matTrain.word2_On = word2_On;
    matTrain.word2_Off = word2_Off;
    catch
        disp('Unable to write to train MAT-file for subject number ');
        %continue;
    end;
   %% TEST TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    visuals_On = 1;
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
            
            if isalmost(t, word1_On,tolerance) 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(wrdList(trt))),wstrength});
            end
            if isalmost(t, visuals_On,tolerance)
               
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
           if isalmost(t, word1_Off,tolerance)
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(wrdList(trt))),0});
            end

            if isalmost(t, visuals_Off,tolerance)                 
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
    matTest.historyLWBt = savestate_historyLWBt{subject};
    matTest.historyRWBt = savestate_historyRWBt{subject};
    matTest.Testset = Testset{subject};
    
    matTest.visuals_On = visuals_On;
    matTest.visuals_Off = visuals_Off
    matTest.word1_On = word1_On;
    matTest.word1_Off = word1_Off;
    catch
        disp('Unable to write to test MAT-file for subject number ');
        %continue;
     end;
end
          