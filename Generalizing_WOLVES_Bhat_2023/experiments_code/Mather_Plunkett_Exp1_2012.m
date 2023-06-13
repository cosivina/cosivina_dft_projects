%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn
%% Age group - 22 month olds.. 
%% pretrain... 2 known objects...  1 obje per trial at random of 3 locs
% possibly add 'wake up trials' for the known objects
%% familiarize 2 new objects ..... 1 obj per trial at rand of 3 locations
%% 2 test trials.. 3 objects per trial.. one known, one familiar, one novel.. at rand of 3 locs
% measure looking at test
nObjects = 6; %%total number of objects
nBlocks = 2;
maxTRwakeUp = 2; % 2 blocks of 2 trials
maxTR = 6; % total # of familiarization trials, 2 blocks of 6 trials each
maxTRt = 3; % total # of test trials, 3 for each novel label
t_maxWakeUp = floor((4000+1000)/scale_factor);
t_max = floor((4000+1000)/scale_factor); %specify simulation time  % scale with Experiment familiarization trial Duration   
t_maxt = floor((8000+1000)/scale_factor); %specify simulation time  % scale with Experiment test trial Duration

History_Reset;
sim.init();
Gui_History_Reset;

for subject = 1:numSubjects
                
    sim2 = 5;
    if (mode == 1 )             
        sim2 = sim;
    elseif (mode == 0)
        sim2 = sim;
    elseif (mode == 2)
        sim2 = sim.copy();
    end
    
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[]; 
    chosen_indices = randperm(size(Property{1},2),nObjects);%% Create objects
    for j=1:nFeatures     
        TempFe=[];
        for i=1:nObjects
            TempFe = [TempFe  Property{j}(chosen_indices(i))];
            %TempFe = [TempFe  Property{j}(2*i)];
        end
        Feature{j}=TempFe;
    end
    Words=[]; 
    chosen_indices = randperm(size(Names,2), nObjects);%% Create  pseudowords
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
    %savestate_historyL{subject}=zeros(maxFam,t_maxF);
    %savestate_historyR{subject}=zeros(maxFam,t_maxF);
    
    %% MAKE PRE_KNOWN WORDS INTO MEMORY
    association_strength=0.30;
    assoc_strength=association_strength;%%strength of known association
    inputMapping1=zeros(fieldSize_ftr,fieldSize_wd);
    inputMapping2=zeros(fieldSize_ftr,fieldSize_wd);
    
    inputSceneRep1=zeros(fieldSize_ftr,fieldSize_spt);
    inputSceneRep2=zeros(fieldSize_ftr,fieldSize_spt);
    
    vis_trace1=zeros(fieldSize_ftr,1);vis_trace2=zeros(fieldSize_ftr,1);
    for kk=1:2
        %wob's
        xx1=cell2mat(Feature{1}(kk)); xx2=cell2mat(Feature{2}(kk)); yy=cell2mat(Words(kk));
            inputMapping1(xx1,yy)= assoc_strength;
            inputMapping2(xx2,yy)= assoc_strength;
            vis_trace1(xx1)=1; vis_trace2(xx2)=1;
            for jj=1:8
                inputMapping1(xx1+jj-4,yy)= assoc_strength;
                inputMapping2(xx2+jj-4,yy)= assoc_strength;
                vis_trace1(xx1+jj-4)=1; vis_trace2(xx2+jj-4)=1;
            end
        %scene mems
        for kmk=1:fieldSize_spt
            inputSceneRep1(xx1,kmk)= assoc_strength;
            inputSceneRep2(xx2,kmk)= assoc_strength;
            for jj=1:8
                inputSceneRep1(xx1+jj-4,kmk)= assoc_strength;
                inputSceneRep2(xx2+jj-4,kmk)= assoc_strength;
            end
        end
       
    end
    sv_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
    sv_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
    sv_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
    sv_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
    sv_hcon_s{subject}=zeros(fieldSize_spt);
    sv_hwm_s{subject}=zeros(fieldSize_spt);
    sv_hword{subject}=zeros(fieldSize_wd);
    for i = 1 : nFeatures
        n=num2str(i);
        sv_hcon_f{subject}(i,:) = zeros(1,fieldSize_ftr); %squeeze(xsit_result.train(1).hcon_f(i,:))*0;
        sv_hwm_f{subject}(i,:) = zeros(1,fieldSize_ftr); %squeeze(xsit_result.train(1).hwm_f(i,:))*0;
        sv_hwf{subject}(i,:,:) = zeros(fieldSize_ftr,fieldSize_wd); %squeeze(xsit_result.train(1).hwf(i,:,:))*0;
        sv_hwm_c{subject}(i,:,:)= zeros(fieldSize_ftr,fieldSize_spt); %squeeze(xsit_result.train(1).hwm_c(i,:,:))*0;
    end
    sv_hcon_s{subject} = zeros(1,fieldSize_spt);% squeeze(xsit_result.train(1).hcon_s)*0;
    sv_hwm_s{subject}= zeros(1,fieldSize_spt);%squeeze(xsit_result.train(1).hwm_s)*0;
    sv_hword{subject}= zeros(1,fieldSize_wd);%squeeze(xsit_result.train(1).hword)*0;


    for i = 1 : nFeatures
        savestate_hcon_f{subject}(i,:) = sv_hcon_f{subject}(i,:);
        savestate_hwm_f{subject}(i,:) = sv_hwm_f{subject}(i,:);
        savestate_hwf{subject}(i,:,:) = sv_hwf{subject}(i,:,:);
        savestate_hwm_c{subject}(i,:,:)= sv_hwm_c{subject}(i,:,:);
    end
    savestate_hcon_s{subject} = sv_hcon_s{subject};
    savestate_hwm_s{subject} = sv_hwm_s{subject};
    savestate_hword{subject} = sv_hword{subject};
    %%% ADD HERE KNOWN WORDS 
    savestate_hwf{subject}(1,:,:) = inputMapping1(:,:);
    savestate_hwf{subject}(2,:,:)=  inputMapping2(:,:);
    %%% into scene memory
    %savestate_hwm_c{subject}(1,:,:) = inputSceneRep1(:,:);
    %savestate_hwm_c{subject}(2,:,:) = inputSceneRep2(:,:);
    %%% into contrast mem
    %savestate_hcon_f{subject}(1,:) = vis_trace1(:);
    %savestate_hcon_f{subject}(2,:) = vis_trace2(:);
    %% RUN WAKE-UP TRIALS 
    visuals_On = 0;
    visuals_Off = floor(4000/scale_factor);
    word_On =  floor([1000 2500]/scale_factor); word_Off = floor((600+[1000 2500])/scale_factor); %% used 600 instead of 550 & 650
    
%     visuals_Off = floor(32000/scale_factor);
%     word_On =  floor([0 16000]/scale_factor); word_Off = floor((16000+[0 16000])/scale_factor); %% used 600 instead of 550 & 650
  
    objSeq = randperm(2);
    locss = [spt_L spt_C spt_R]; locSeq= locss(randperm(length(locss)));
    for ib=1:nBlocks
        if ib>1; objSeq = flip(objSeq);end
        
        for tr=1:maxTRwakeUp
            sim2.t =sim2.tZero;            %sim2.init();
            targ = objSeq(tr); % Find the correct touple
             loc = locSeq(tr);
           % if (tr || ib) > 1
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
            %end
            while sim2.t <= t_maxWakeUp

                t = sim2.t;
                if isalmost(t, visuals_On,tolerance)
                    for f=1:nFeatures
                         sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(targ)),loc, vstrength});                    
                    end
                end
                if isalmost(t, word_On,tolerance)          
                    sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                        {cell2mat(Words(targ)),wstrength});               
                end
                if isalmost(t, word_Off,tolerance)           
                    sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                        {cell2mat(Words(targ)),0});               
                end
                if isalmost(t, visuals_Off,tolerance)
                    for f=1:nFeatures
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(targ)),loc, 0});    
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
    %               savestate_historyLt{subject}(tr,:) = sim2.getComponent('history lookLt', 'output');
    %               savestate_historyRt{subject}(tr,:) = sim2.getComponent('history lookRt', 'output');
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
        end% wakeUp ends

        % fill each object for 3 locs... ie. 2x3 mappings... then use order of these mappings
        objloc_pair =[];
        for i=1:length(objSeq);   for j=1:length(locSeq);  objloc_pair = [objloc_pair; objSeq(i)+2 locSeq(j)]; end ;  end
        objloc_pair = objloc_pair(randperm(size(objloc_pair,1)),:);
        %% RUN FAMILIARIZATION TRIALS
        for tr=1:maxTR
         savestate_hcon_f{subject}(1,:) = vis_trace1(:);
        savestate_hcon_f{subject}(2,:) = vis_trace2(:);
            sim2.t =sim2.tZero;            %sim2.init();
            targ = objloc_pair(tr, 1); % Find the correct touple
             loc = objloc_pair(tr, 2);
            
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
                if isalmost(t, visuals_On,tolerance)
                    for f=1:nFeatures
                         sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(targ)),loc, vstrength});                    
                    end
                end
                if isalmost(t, visuals_Off,tolerance)
                    for f=1:nFeatures
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(targ)),loc, 0});    
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
%                    savestate_historyL{subject}(tr,:) = sim2.getComponent('history lookL', 'output');
%                    savestate_historyR{subject}(tr,:) = sim2.getComponent('history lookR', 'output');
%                    savestate_historyC{subject}(tr,:) = sim2.getComponent('history lookC', 'output');
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
        end% familiarization ends       
    end %block ends
     try
        OutName = [simName,num2str(subject),'_train.mat'];
        matTrain=matfile([pwd,'\',OutName],'writable',true);
        matTrain.hcon_f = savestate_hcon_f{subject};
        matTrain.hwm_f = savestate_hwm_f{subject};
        matTrain.hwm_c = savestate_hwm_c{subject};
        matTrain.hwf = savestate_hwf{subject};
        matTrain.hcon_s = savestate_hcon_s{subject};
        matTrain.hwm_s = savestate_hwm_s{subject};
        matTrain.hword = savestate_hword{subject};
    %     matTrain.historyL = savestate_historyL{subject};
    %     matTrain.historyR = savestate_historyR{subject};
    %     matTrain.historyC = savestate_historyC{subject};
    %     %matTrain.training_pair = training_pair{subject};
        matTrain.Feature1 = Feature{1};
        matTrain.Feature2 = Feature{2};
        matTrain.Words = Words;
     catch
        disp('Unable to write to train MAT-file for subject number ');
        %continue;
     end;
%% novel trials
        savestate_historyLt{subject}=zeros(nBlocks*maxTRt,t_maxt);
        savestate_historyRt{subject}=zeros(nBlocks*maxTRt,t_maxt);
        savestate_historyCt{subject}=zeros(nBlocks*maxTRt,t_maxt);
        visuals_On = 0;
        visuals_Off = floor(8000/scale_factor);
        word_On =  floor([3633 5633]/scale_factor); word_Off = floor((600+[3633 5633])/scale_factor); %% used 600 instead of 550 & 650
        test_pair{subject}=zeros(nBlocks*maxTRt,6); 
            %% RUN NOVEL TRIALS
         for ib=1:nBlocks
            objloc_pair = [1 3 5]+ib-1;
            locT = randperm(3);
            locSeq(1,:) = locT;
            locSeq(2,:) = mod(locT,3)+1;
            locSeq(3,:) = mod(locT+1,3)+1;
            for tr=1:maxTRt
                    savestate_hcon_f{subject}(1,:) = vis_trace1(:);
                    savestate_hcon_f{subject}(2,:) = vis_trace2(:);
                test_pair{subject}(((ib-1)*maxTRt)+tr,:)=  [objloc_pair locSeq(tr,:)];
                sim2.t =sim2.tZero;            %sim2.init();            
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
                            sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set feature dimension value
                                strcat('Feature_',num2str(f),'_Centre'), strcat('Feature_',num2str(f),'_Right')},{'positionY', 'positionY', 'positionY'}, ...
                                {cell2mat(Feature{f}(objloc_pair(1))),...
                                cell2mat(Feature{f}(objloc_pair(2))),cell2mat(Feature{f}(objloc_pair(3)))});                

                            sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set spatial dimension value
                                strcat('Feature_',num2str(f),'_Centre'), strcat('Feature_',num2str(f),'_Right')}, ...
                                 {'positionX', 'positionX', 'positionX'}, ...
                                 {locss(locSeq(tr,1)),locss(locSeq(tr,2)),locss(locSeq(tr,3))});

                            sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... % set left right stimuli up
                                strcat('Feature_',num2str(f),'_Centre'), strcat('Feature_',num2str(f),'_Right')}, ...
                                 {'amplitude', 'amplitude', 'amplitude'}, vstrength * ones(1, 3));                       
                        end
                    end
                    if isalmost(t, word_On,tolerance)          
                        sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                            {cell2mat(Words(objloc_pair(3))),wstrength});               
                    end
                    if isalmost(t, word_Off,tolerance)           
                        sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                            {cell2mat(Words(objloc_pair(3))),0});               
                    end
                    if isalmost(t, visuals_Off,tolerance)
                        for f=1:nFeatures
                            cellNum=(2*f)-1;

                            sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... 
                             strcat('Feature_',num2str(f),'_Centre'),   strcat('Feature_',num2str(f),'_Right')}, ...
                                 {'amplitude', 'amplitude', 'amplitude'}, 0 * ones(1, 3));% stop left right stimuli

                            sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... 
                             strcat('Feature_',num2str(f),'_Centre'),   strcat('Feature_',num2str(f),'_Right')}, ...
                                {'positionY', 'positionY', 'positionY'},{0,0,0}); % Unset feature dimension value

                            sim2.setElementParameters({ strcat('Feature_',num2str(f),'_Left'), ... 
                             strcat('Feature_',num2str(f),'_Centre'),   strcat('Feature_',num2str(f),'_Right')}, ...
                                 {'positionX', 'positionX', 'positionX'}, {0,0,0});% Unset spatial dimension value
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

                     savestate_historyLt{subject}(((ib-1)*maxTRt)+tr,:) = flipud(sim2.getComponent('history lookLt', 'output'));
                     savestate_historyRt{subject}(((ib-1)*maxTRt)+tr,:) =flipud(sim2.getComponent('history lookRt', 'output'));
                     savestate_historyCt{subject}(((ib-1)*maxTRt)+tr,:) =flipud(sim2.getComponent('history lookCt', 'output'));

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
            end% a test word ends 
         end %  test ends
        try
            OutName = [simName,num2str(subject),'_test.mat'];
            matTest=matfile([pwd,'\',OutName],'writable',true);
            matTest.hcon_f = savestate_hcon_f{subject};
            matTest.hwm_f = savestate_hwm_f{subject};
            matTest.hwm_c = savestate_hwm_c{subject};
            matTest.hwf = savestate_hwf{subject};
            matTest.hcon_s = savestate_hcon_s{subject};
            matTest.hwm_s = savestate_hwm_s{subject};
            matTest.hword = savestate_hword{subject};
            matTest.historyLt = savestate_historyLt{subject};
            matTest.historyRt = savestate_historyRt{subject};
            matTest.historyCt = savestate_historyCt{subject};
            matTest.test_pair = test_pair{subject};
        catch
           disp('Unable to write to test MAT-file for subject number ');
           continue;
        end
end