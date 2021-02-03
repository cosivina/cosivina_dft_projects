nObjects = 12; %%total number of objects
maxTR = 60; % total # of training trials
maxTRt = 0; % total # of test trials
t_max = floor((6000+1000)/scale_factor); %specify simulation time  % scale with Experiment training trial Duration   
t_maxt = floor((100+1000)/scale_factor); %sAssumed % scale with Experiment test trial Duration

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
    visuals_Off = floor((6000)/scale_factor);
     %word1_On  = 1;% floor(200/scale_factor); %values assumed, 
     %word1_Off = floor(1000/scale_factor);
        % word1_On  = floor([1 2150 4050 5950  7850]/scale_factor); %values assumed, authors provide time after word onset 
        % word1_Off = floor([1350 3150 5050 6950 9000]/scale_factor);
        
   word1_On  = floor([8 2500 4250 ]/scale_factor); %values assumed, authors provide time after word onset 
   word1_Off = floor([1500 3250 5000]/scale_factor);    
    
   %word1_On  = floor([1 1250 2500 3750]/scale_factor); 
    %word1_Off = floor([750 2000 3250 4500]/scale_factor);
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[]; 
    Words=[]; 
    for j=1:nFeatures
        chosen_indices = randperm(size(Property{j},2),nObjects);
        TempFe=[];
        for i=1:nObjects
            %TempFe = [TempFe  Property{j}(chosen_indices(i))];
            TempFe = [TempFe  Property{j}(i)];
        end
        Feature{j}=TempFe;
    end
    % Choose Words to associate on random
    chosen_indices = randperm(size(Names,2),nObjects);
    for i=1:nObjects
        %Words = [Words  Names(chosen_indices(i))];
        Words = [Words  Names(i)];
    end
    %% DESIGN TRAINING TRIALS
    Indices_Set{subject}=zeros(maxTR,7);
    repeat=1;% conduct the trial sequence finding process or reconduct if previously got stuck
    while(repeat)
        tic();
        Indices_Set{subject}=zeros(maxTR,7);
        pos=1;
        cross_correlations=zeros(12,12);% need to make sure no more than TWO cross-correlations in final set
       for word_rep=1:(maxTR/nObjects)

           for wordy=1:nObjects
               locs=randperm(5);%choose 5 random locations for objects to place
               Indices_Set{subject}(pos,locs(5)) = wordy;% place the target 
               Indices_Set{subject}(pos,6) = locs(5);% remember the target location
               Indices_Set{subject}(pos,7) = wordy;% add word for the target.. same as target here
               disty= randperm(nObjects,4);% choose four random objects as distractors
               bad_correl_flag=1;% assuming we have bad correlations; flags if we have more than 2 correlations
               while (ismember(wordy,disty) ||  bad_correl_flag) % keep finding new distractors until no bad correlations
                 if  toc() > 5; repeat=1; break; end %%looks like made bad choices and cannot find any new correlation-free pairs, so exit
                 if ismember(wordy,disty); disty= randperm(nObjects,4);end %checking if the random set does not have target
                 correlation_count=0;% to check if no more than two correlations will form
                   for kl=1:4
                       if (cross_correlations(wordy,disty(kl)) + 1 > 2); bad_correl_flag=1; else; correlation_count=correlation_count+1; end 
                   end
                   if correlation_count==4; bad_correl_flag=0; else; disty= randperm(nObjects,4); end
               end
                %found a good set of distractors
               for kl=1:4; Indices_Set{subject}(pos,locs(kl))=disty(kl); cross_correlations(wordy,disty(kl))=cross_correlations(wordy,disty(kl))+1; end
               pos=pos+1;
           end
       end      
       repeat=0;
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
    savestate_historyL{subject}=zeros(maxTR,t_max);
    savestate_historyR{subject}=zeros(maxTR,t_max);
    savestate_historyLWB{subject}=zeros(maxTR,t_max);
    savestate_historyRWB{subject}=zeros(maxTR,t_max);
    savestate_historyC{subject}=zeros(maxTR,t_max);
    training_pair{subject}=zeros(maxTR,7);% save trial information, objs on screen with word sequence
    %% RUN TRAINING TRIALS
    for tr=1:maxTR
         sim2.t =sim2.tZero;
        pos = tr; % Finds postion of a random tuple
        training_pair{subject}(tr,:) = Indices_Set{subject}(pos,:);
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
                        strcat('Feature_',n,'_Centre'), strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY','positionY'}, ...
                        {cell2mat(Feature{f}(Indices_Set{subject}(pos,1))),cell2mat(Feature{f}(Indices_Set{subject}(pos,2))), ...
                        cell2mat(Feature{f}(Indices_Set{subject}(pos,3))),cell2mat(Feature{f}(Indices_Set{subject}(pos,4))), cell2mat(Feature{f}(Indices_Set{subject}(pos,5)))});                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_Centre'),strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX','positionX'}, ...
                         {spt_L, spt_Lm, spt_C, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_Centre'),strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude','amplitude'}, ...
                          (vstrength) * ones(1, 5) );                     
                end
            end
            if isalmost(t, word1_On,tolerance) 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(Indices_Set{subject}(pos,7))),wstrength});
            end
            if isalmost(t, word1_Off,tolerance)
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(Indices_Set{subject}(pos,7))),0});
            end

            if t == visuals_Off
                for f=1:nFeatures

                    n=num2str(f);                    
                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_Centre'), strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionY', 'positionY','positionY', 'positionY','positionY'}, ...
                         (0) * zeros(1, 5) );                

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_Centre'),strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'positionX', 'positionX','positionX', 'positionX','positionX'}, ...
                         {spt_L, spt_Lm, spt_C, spt_Rm, spt_R});

                    sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
                        strcat('Feature_',n,'_Centre'),strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
                        {'amplitude', 'amplitude','amplitude', 'amplitude','amplitude'}, ...
                          (0) * zeros(1, 5) );  
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
                 savestate_historyL{subject}(tr,:) = flipud(sim2.getComponent('history lookL', 'output'));
                 savestate_historyR{subject}(tr,:) = flipud(sim2.getComponent('history lookR', 'output'));
                 savestate_historyLWB{subject}(tr,:) = flipud(sim2.getComponent('history lookLWB', 'output'));
                 savestate_historyRWB{subject}(tr,:) = flipud(sim2.getComponent('history lookRWB', 'output'));

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
    matTrain.historyC = savestate_historyC{subject};
    matTrain.historyLWB = savestate_historyLWB{subject};
    matTrain.historyRWB = savestate_historyRWB{subject};
    matTrain.training_pair = training_pair{subject};
    matTrain.Feature1 = Feature{1};
    matTrain.Feature2 = Feature{2};
    matTrain.Words = Words;
    catch
                disp('Error saving a train file');
    end
%    %% TEST TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     visuals_On = 0;
%     visuals_Off = floor((1000)/scale_factor);
%     word1_On  =1;% assumed
%     word1_Off  = floor((1000)/scale_factor);
%     Testseq=[];
%     wrdList=randperm(nObjects);
%     for tseq=1:maxTRt
%      % for tseq_th word   %remove target to generate distractor list 
%      disSet=wrdList(wrdList~=wrdList(tseq));
%      dis3=disSet(randperm(length(disSet),3)); %pick 3 distractors
%      indxr=randperm(4);
%      for ds=1:3 Testseq(tseq,indxr(ds))=dis3(ds);end
%      Testseq(tseq,indxr(4))=wrdList(tseq);%object
%      Testseq(tseq,5)=indxr(4);%location of target   
%     end
%     Testset{subject}=Testseq;
%     %% RUN TEST TRIALS
%     savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
%     savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
%     savestate_historyLWBt{subject}=zeros(maxTRt,t_maxt);
%     savestate_historyRWBt{subject}=zeros(maxTRt,t_maxt);
%     for trt=1:maxTRt
%               
%         sim2.init();  
%            
%         for i = 1 : nFeatures
%                 n=num2str(i);
%                 handle_hcon_f{subject}(i).output = squeeze(savestate_hcon_f{subject}(i,:));
%                 handle_hwm_f{subject}(i).output = squeeze(savestate_hwm_f{subject}(i,:));
%                 handle_hwf{subject}(i).output = squeeze(savestate_hwf{subject}(i,:,:));
%                 handle_hwm_c{subject}(i).output = squeeze(savestate_hwm_c{subject}(i,:,:));
%         end
%         handle_hcon_s{subject}.output = savestate_hcon_s{subject};
%         handle_hwm_s{subject}.output = savestate_hwm_s{subject};
%         handle_hword{subject}.output = savestate_hword{subject};
%      
%         while sim2.t <= t_maxt
% 
%             t = sim2.t;
%             
%             if t== word1_On
%                 sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
%                     {cell2mat(Words(wrdList(trt))),wstrength});
%             end
%             if t == visuals_On
%                
%                 for f=1:nFeatures
%                     n=num2str(f);
%                     
%                     sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
%                         strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
%                         {'positionY', 'positionY','positionY', 'positionY'}, ...
%                         {cell2mat(Feature{f}(Testseq(trt,1))),cell2mat(Feature{f}(Testseq(trt,2))), ...
%                         cell2mat(Feature{f}(Testseq(trt,3))),cell2mat(Feature{f}(Testseq(trt,4)))});                
% 
%                     sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
%                         strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
%                         {'positionX', 'positionX','positionX', 'positionX'}, ...
%                          {spt_L, spt_Lm, spt_Rm, spt_R});
% 
%                     sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
%                         strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
%                         {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
%                          { vstrength, vstrength, vstrength, vstrength }); 
%                     
%                 end              
%             end
%            if t== word1_Off
%                 sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
%                     {cell2mat(Words(wrdList(trt))),0});
%             end
% 
%             if t == visuals_Off
%                for f=1:nFeatures
%                     n=num2str(f);
%                     
%                     sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
%                         strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
%                         {'positionY', 'positionY','positionY', 'positionY'}, ...
%                         0 * ones(1, 4));                
% 
%                     sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
%                         strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
%                         {'positionX', 'positionX','positionX', 'positionX'}, ...
%                          {spt_L, spt_Lm, spt_Rm, spt_R});
% 
%                     sim2.setElementParameters({ strcat('Feature_',n,'_Left'), strcat('Feature_',n,'_LeftMid'), ... % set feature dimension value
%                         strcat('Feature_',n,'_RightMid'), strcat('Feature_',n,'_Right')}, ...
%                         {'amplitude', 'amplitude','amplitude', 'amplitude'}, ...
%                          { 0, 0, 0, 0}); 
%                     
%                 end
%             end
%             if isalmost(t, t_maxt,tolerance) 
%                for i = 1 : nFeatures
%                     n=num2str(i);
%                     savestate_hcon_f{subject}(i,:) = sim2.getComponent(['hcon_f' n], 'output');
%                     savestate_hwm_f{subject}(i,:) = sim2.getComponent(['hwm_f' n], 'output');
%                     savestate_hwf{subject}(i,:,:) = sim2.getComponent(['hwf' n], 'output');
%                     savestate_hwm_c{subject}(i,:,:) = sim2.getComponent(['hwm_c' n], 'output');
%                 end
%                 savestate_hcon_s{subject} = sim2.getComponent('hcon_s', 'output');
%                 savestate_hwm_s{subject} = sim2.getComponent('hwm_s', 'output');
%                 savestate_hword{subject} = sim2.getComponent('hword', 'output');
%                 savestate_historyLt{subject}(trt,:) = flipud(sim2.getComponent('history lookLt', 'output'));
%                 savestate_historyRt{subject}(trt,:) = flipud(sim2.getComponent('history lookRt', 'output'));
%                 savestate_historyLWBt{subject}(trt,:) = flipud(sim2.getComponent('history lookLWBt', 'output'));
%                 savestate_historyRWBt{subject}(trt,:) = flipud(sim2.getComponent('history lookRWBt', 'output'));
%                 
%             end
%             
% 
%            if (mode == 1)
%                 if ~gui.pauseSimulation
%                     sim2.step();
%                 end
%                 if gui.quitSimulation
%                     gui.close();
%                     break;
%                 end
%                 if (mod(t,15)==0)
%                     gui.step();
%                 end
%                 %gui.checkAndUpdateControls();
%                 %gui.updateVisualizations();                  
%            else
%                 sim2.step();
%            end
%         end
%     end
    try
    OutName = [simName,num2str(subject),'_test.mat'];
    matTest=matfile(OutName,'writable',true);

    %matTrain.hcon_ft = savestate_hcon_f{subject};
    %matTest.hwm_ft = savestate_hwm_f{subject};
    %matTest.hwm_ct = savestate_hwm_c{subject};
    matTest.hwft = savestate_hwf{subject};
    %matTrain.hcon_st = savestate_hcon_s{subject};
    %matTrain.hwm_st = savestate_hwm_s{subject};
    %matTest.hwordt = savestate_hword{subject};
    
%     matTest.historyLt = savestate_historyLt{subject};
%     matTest.historyRt = savestate_historyRt{subject};
%     matTest.historyLWBt = savestate_historyLWBt{subject};
%     matTest.historyRWBt = savestate_historyRWBt{subject};
%     matTest.Testset = Testset{subject};
    catch
                disp('Error saving a test file');
    end
end
          