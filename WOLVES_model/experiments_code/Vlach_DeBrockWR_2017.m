nObjects = 18; %%total number of objects % 12 trained + 6 novel at test
maxTR = 12; % 12 in exp actually total # of training trials
maxTRt = 12; % total # of test trials
t_max = floor(1200/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
t_maxt = floor(1200/scale_factor); %specify simulation time  % scale with Experiment test trial Duration

%sim.setElementParameters('history lookL', 'timeSlots', historyDuration ); % set spatial dimension value
% handle_historyL=sim.getElement('history lookL'); %reset lookingduration history variable                 
% handle_historyL.timeSlots = t_max; 
% handle_historyR=sim.getElement('history lookR');                       
% handle_historyR.timeSlots = t_max; 
% handle_historyLWB=sim.getElement('history lookLWB'); %reset lookingduration history variable                 
% handle_historyLWB.timeSlots = t_max; 
% handle_historyRWB=sim.getElement('history lookRWB');                       
% handle_historyRWB.timeSlots = t_max; 
% handle_historyC=sim.getElement('history lookC');                       
% handle_historyC.timeSlots = t_max; 
% 
% handle_historyLt=sim.getElement('history lookLt');                       
% handle_historyLt.timeSlots = t_maxt; 
% handle_historyRt=sim.getElement('history lookRt');                       
% handle_historyRt.timeSlots = t_maxt;
% handle_historyLWBt=sim.getElement('history lookLWBt'); %reset lookingduration history variable                 
% handle_historyLWBt.timeSlots = t_maxt; 
% handle_historyRWBt=sim.getElement('history lookRWBt');                       
% handle_historyRWBt.timeSlots = t_maxt; 
% handle_historyCt=sim.getElement('history lookCt');                       
% handle_historyCt.timeSlots = t_maxt; 

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

%% WR MEMORY TASK (Word Recognition)
%% TRAINING PHASE %skipped

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
    words_On = 0;
    words_Off = floor(1000/scale_factor);
    
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Words=[]; 
    chosen_indices = randperm(nObjects);%% Create all 18 words, 12 for learning, 6 as novel at tests
    for i=1:nObjects
        Words = [Words  Names(chosen_indices(i))];
        %Words = [Words  Names(2*i)];
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
%     savestate_historyL{subject}=zeros(maxTR,t_max);
%     savestate_historyR{subject}=zeros(maxTR,t_max);
    training_pair{subject}=zeros(nObjects*(nObjects-1),(2*nFeatures)+3);%
    %% RUN TRAINING TRIALS
     for tr=1:maxTR
        %sim2.init();
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

            % STEP 1 I have skipped the training part
            if t == words_On            
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos)),wstrength});               
            end
            if t == words_Off
                                
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos)),0});
            end
            if t == t_max
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
%                 savestate_historyLt{subject}(tr,:) = sim2.getComponent('history lookLt', 'output');
%                 savestate_historyRt{subject}(tr,:) = sim2.getComponent('history lookRt', 'output');
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
    OutName = [simName,num2str(subject),'_train.mat'];
    matTrain=matfile(OutName,'writable',true);

    matTrain.hcon_f = savestate_hcon_f{subject};
    matTrain.hwm_f = savestate_hwm_f{subject};
    matTrain.hwm_c = savestate_hwm_c{subject};
    matTrain.hwf = savestate_hwf{subject};
    matTrain.hcon_s = savestate_hcon_s{subject};
    matTrain.hwm_s = savestate_hwm_s{subject};
    matTrain.hword = savestate_hword{subject};

%     matTrain.historyL = savestate_historyL{subject};
%     matTrain.historyR = savestate_historyR{subject};
    matTrain.training_pair = training_pair{subject};
   % matTrain.Feature1 = Feature{1};
   % matTrain.Feature2 = Feature{2};
    matTrain.Words = Words;
%% TESTING PHASE
    words_On = 0;
    words_Off = floor(1000/scale_factor);
    %savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    %savestate_hword{subject}=zeros(fieldSize_wd);
    savestate_historyW{subject}=zeros(maxTRt,t_maxt,fieldSize_wd);
    test_word{subject}=zeros(maxTRt,2);
    seq= randperm(12);
    pos = randperm(12);
    postDis=randperm(6);
    pDis=1;
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
            %
            if t == words_On
                 if seq(trt)>6 % 50% times..half the number of size(seq) like if rand >0.5
                    sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos(trt))),wstrength-2.3}); 
                test_word{subject}(trt,:)= [ cell2mat(Words(pos(trt))) 'F'];
                 else
                     sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(12+postDis(pDis))),wstrength-2.3}); 
                 test_word{subject}(trt,:)= [ cell2mat(Words(12+postDis(pDis))) 'N']; pDis=pDis+1;
                 end
            end
            if t == words_Off
                              
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(pos(trt))),0}); 
            end
            if t == t_maxt
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
                %savestate_historyLt{subject}(tr,:) = sim2.getComponent('history lookLt', 'output');
            end
            savestate_historyW{subject}(trt,t+1,:) = flipud(sim2.getComponent('word', 'output'));
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
    OutName = [simName,num2str(subject),'_test.mat'];
    matTest=matfile(OutName,'writable',true);

    matTrain.hcon_ft = savestate_hcon_f{subject};
    matTest.hwm_ft = savestate_hwm_f{subject};
    matTest.hwm_ct = savestate_hwm_c{subject};
    matTest.hwft = savestate_hwf{subject};
    matTrain.hcon_st = savestate_hcon_s{subject};
    matTrain.hwm_st = savestate_hwm_s{subject};
    matTest.hwordt = savestate_hword{subject};
    matTest.historyW = savestate_historyW{subject};
    matTest.test_word = test_word{subject};
    if (mode == 1)
         gui.close();
    end
end