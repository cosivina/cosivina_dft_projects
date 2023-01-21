%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn
% The study is run in  two conditions Silent (i.e. without any labels) and
% Labelling (i.e with a same label on all trials)... control_variable is Labelling_condition_ON
% We added a prediction case wherein the labelling condition, a new novel
% is presented on every trial. ... control_variable is prediction_mode_ON


Property{1} = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300};%%Shapes  %%more than or equal to nObjects
Property{2} = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300};%%Colors   %%more than or equal to nObjects
nObjects = 31; %%total number of objects
maxTRt = 30; % total # of test trials
t_max=floor((6000+1000)/scale_factor);
t_maxt = floor((6000+1000)/scale_factor); %specify simulation time  % scale with Experiment test trial Duration
%Labelling_condition_ON = 0;%%change to 0 for Silent condition
%prediction_mode_ON = 0;%% set to 1 for generating novel names on every new trial; change to 0 for same word as orginal experiment

History_Reset;
sim.init();
Gui_History_Reset;

%% 
for subject = 1:numSubjects
                
    sim2 = 5;
    if (mode == 1 )             
        sim2 = sim;
    elseif (mode == 0)
        sim2 = sim;
    elseif (mode == 2)
        sim2 = sim.copy();
    end
    
    visuals_On = 0;
    visuals_Off = floor(6000/scale_factor);
    words_On = floor([600 3200]/scale_factor);%[;
    words_Off = floor(([600 3200]+500)/scale_factor);
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[];
    for j=1:nFeatures
        chosen_indices = randperm(size(Property{1},2));%% Create all  objects
        TempFe=[];
        for i=1:nObjects
            TempFe = [TempFe  Property{j}(chosen_indices(i))];
            %TempFe = [TempFe  Property{j}(i)];
        end
        Feature{j}=TempFe;
    end
    Words=[];
    if prediction_mode_ON == 1
         chosen_index = randperm(size(Names,2));%% Create  all novel words for familiariztion
         chosen_index = [chosen_index chosen_index(1:12)];
        for i=1:maxTRt
            Words = [Words  Names(chosen_index(i))];
        end       
    else        
        chosen_index = randperm(size(Names,2), 1);%% Create  1 word for familiariztion
        for i=1:maxTRt
            Words = [Words  Names(chosen_index)];
        end
    end
    for i = 1 : nFeatures
        n=num2str(i);
        handle_hcon_f{subject}(i)=sim2.getElement(['hcon_f' n]);
        handle_hwm_f{subject}(i)=sim2.getElement(['hwm_f' n]);
        handle_hwf{subject}(i)=sim2.getElement(['hwf' n]);
        handle_hwm_c{subject}(i)=sim2.getElement(['hwm_c' n]); 
        
        handle_atn_f{subject}(i)=sim2.getElement(['atn_f' n]);
        handle_con_f{subject}(i)=sim2.getElement(['con_f' n]);
        handle_wm_f{subject}(i)=sim2.getElement(['wm_f' n]);
        handle_wf{subject}(i)=sim2.getElement(['wf' n]);
        handle_wm_c{subject}(i)=sim2.getElement(['wm_c' n]);
        handle_atn_c{subject}(i)=sim2.getElement(['atn_c' n]); 
    end
    handle_hcon_s{subject}=sim2.getElement('hcon_s');
    handle_hwm_s{subject}=sim2.getElement('hwm_s');
    handle_hword{subject}=sim2.getElement('hword');
    handle_atn_sa{subject}=sim2.getElement('atn_sa');
    handle_atn_sr{subject}=sim2.getElement('atn_sr');
    handle_wm_s{subject}=sim2.getElement('wm_s');
    handle_con_s{subject}=sim2.getElement('con_s');
    

%% TESTING PHASE
%     savestate2_atn_f{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr);
%     savestate2_con_f{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr);
%     savestate2_wm_f{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr);
%     savestate2_wf{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr,fieldSize_wd);
%     savestate2_atn_c{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr,fieldSize_spt);
%     savestate2_wm_c{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr,fieldSize_spt);
%     
%     savestate2_hcon_f{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr);
%     savestate2_hwm_f{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr);
%     savestate2_hwf{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr,fieldSize_wd);
%     savestate2_hwm_c{subject}=zeros(maxTRt,t_maxt,nFeatures,fieldSize_ftr,fieldSize_spt);
% 
%     savestate2_atn_sa{subject}=zeros(maxTRt,t_maxt,fieldSize_spt);
%     savestate2_atn_sr{subject}=zeros(maxTRt,t_maxt,fieldSize_spt);
%     savestate2_con_s{subject}=zeros(maxTRt,t_maxt,fieldSize_spt);
%     savestate2_wm_s{subject}=zeros(maxTRt,t_maxt,fieldSize_spt); 
%     savestate2_hcon_s{subject}=zeros(maxTRt,t_maxt,fieldSize_spt);
%     savestate2_hwm_s{subject}=zeros(maxTRt,t_maxt,fieldSize_spt);
%     savestate2_hword{subject}=zeros(maxTRt,t_maxt,fieldSize_wd);
    
    
    
%     savestate_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
%     savestate_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
%     savestate_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
%     savestate_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
%     savestate_hcon_s{subject}=zeros(fieldSize_spt);
%     savestate_hwm_s{subject}=zeros(fieldSize_spt);
%     savestate_hword{subject}=zeros(fieldSize_wd);
    
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    fam_word{subject}=zeros(maxTRt,1);
    
    side_balance=randperm(maxTRt);
    for trt=1:maxTRt
        pos=trt;
        %sim2.init(); 
        sim2.t =sim2.tZero;
%         if trt > 1
%             for i = 1 : nFeatures
%                     n=num2str(i);
%                     handle_hcon_f{subject}(i).output = squeeze(savestate_hcon_f{subject}(i,:));
%                     handle_hwm_f{subject}(i).output = squeeze(savestate_hwm_f{subject}(i,:));
%                     handle_hwf{subject}(i).output = squeeze(savestate_hwf{subject}(i,:,:));
%                     handle_hwm_c{subject}(i).output = squeeze(savestate_hwm_c{subject}(i,:,:));
%             end
%             handle_hcon_s{subject}.output = savestate_hcon_s{subject};
%             handle_hwm_s{subject}.output = savestate_hwm_s{subject};
%             handle_hword{subject}.output = savestate_hword{subject};
%         end
        while sim2.t <= t_maxt

            t = sim2.t;
            %
            if isalmost(t, visuals_On,tolerance) 
                if mod(side_balance(trt),2)
                    for f=1:nFeatures                   
                         sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(1)),spt_L, vstrength});   
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Right')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(pos+1)),spt_R, vstrength}); 
                    end
                    fam_word{subject}(trt,1)=  'L';
                else
                     for f=1:nFeatures
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(pos+1)),spt_L, vstrength});   
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Right')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(1)),spt_R, vstrength});
                     end
                     fam_word{subject}(trt,1)=  'R';
                 end
             end
            if Labelling_condition_ON == 1    
                if isalmost(t, words_On,tolerance)                 
                    sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                        {cell2mat(Words(trt)),wstrength});               
                end
                if isalmost(t, words_Off,tolerance)              
                    sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                        {cell2mat(Words(trt)),0}); 
                end
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
                savestate_historyLt{subject}(trt,:) = flipud(sim2.getComponent('history lookLt', 'output'));
                savestate_historyRt{subject}(trt,:) = flipud(sim2.getComponent('history lookRt', 'output'));
            end
%             if t > 1
%                 for i = 1 : nFeatures
%                     n=num2str(i);
%                     savestate2_atn_f{subject}(trt,t,i,:) = sim2.getComponent(['atn_f' n], 'output');
%                     savestate2_con_f{subject}(trt,t,i,:) = sim2.getComponent(['con_f' n], 'output');
%                     savestate2_wm_f{subject}(trt,t,i,:) = sim2.getComponent(['wm_f' n], 'output');
%                     savestate2_wf{subject}(trt,t,i,:,:) = sim2.getComponent(['wf' n], 'output');
%                     savestate2_atn_c{subject}(trt,t,i,:,:) = sim2.getComponent(['atn_c' n], 'output');
%                     savestate2_wm_c{subject}(trt,t,i,:,:) = sim2.getComponent(['wm_c' n], 'output');             
%                     
%                     savestate2_hcon_f{subject}(trt,t,i,:) = sim2.getComponent(['hcon_f' n], 'output');
%                     savestate2_hwm_f{subject}(trt,t,i,:) = sim2.getComponent(['hwm_f' n], 'output');
%                     savestate2_hwf{subject}(trt,t,i,:,:) = sim2.getComponent(['hwf' n], 'output');
%                     savestate2_hwm_c{subject}(trt,t,i,:,:) = sim2.getComponent(['hwm_c' n], 'output');
%                 end
%                 savestate2_atn_sa{subject}(trt,t,:) = sim2.getComponent('atn_sa', 'output');
%                 savestate2_atn_sr{subject}(trt,t,:) = sim2.getComponent('atn_sr', 'output');
%                 savestate2_con_s{subject}(trt,t,:) = sim2.getComponent('con_s', 'output');
%                 savestate2_wm_s{subject}(trt,t,:) = sim2.getComponent('wm_s', 'output');
%                 savestate2_hcon_s{subject}(trt,t,:) = sim2.getComponent('hcon_s', 'output');
%                 savestate2_hwm_s{subject}(trt,t,:) = sim2.getComponent('hwm_s', 'output');
%                 savestate2_hword{subject}(trt,t,:) = sim2.getComponent('hword', 'output');
%            end
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
    end
    try
    OutName = [simName,num2str(subject),'_test.mat'];
    matTest=matfile(OutName,'writable',true);
    
%     matTest.con_ft = savestate_con_f{subject};
%     matTest.hcon_ft = savestate_hcon_f{subject};
%     matTest.hwm_ft = savestate_hwm_f{subject};
%     matTest.hwm_ct = savestate_hwm_c{subject};
%     matTest.hwft = savestate_hwf{subject};
%     matTest.hcon_st = savestate_hcon_s{subject};
%     matTest.hwm_st = savestate_hwm_s{subject};
%     matTest.hwordt = savestate_hword{subject};
    matTest.historyLt = savestate_historyLt{subject};
    matTest.historyRt = savestate_historyRt{subject};
    matTest.fam_word = fam_word{subject};
    
    
%     matTest.atn_f = savestate2_atn_f{subject};
%     matTest.con_f = savestate2_con_f{subject};
%     matTest.wm_f = savestate2_wm_f{subject};
%     matTest.wf = savestate2_wf{subject};
%     matTest.wm_c = savestate2_wm_c{subject};
%     matTest.hcon_f = savestate2_hcon_f{subject};
%     matTest.hwm_f = savestate2_hwm_f{subject};
%     matTest.hwf = savestate2_hwf{subject};
%     matTest.hwm_c = savestate2_hwm_c{subject};
%     matTest.atn_sa = savestate2_atn_sa{subject};
%     matTest.atn_sr = savestate2_atn_sr{subject};
%     matTest.con_s = savestate2_con_s{subject};
%     matTest.wm_s = savestate2_wm_s{subject};
%     matTest.hcon_s = savestate2_hcon_s{subject};
%     matTest.hwm_s = savestate2_hwm_s{subject};
%     matTest.hword = savestate2_hword{subject};
    

    catch
    disp('Error on saving test file ');
    end
 end