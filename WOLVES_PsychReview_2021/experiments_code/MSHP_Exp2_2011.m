Property{1} = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300};%%Shapes  %%more than or equal to nObjects
Property{2} = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300};%%Colors   %%more than or equal to nObjects
nObjects = 28; %%total number of objects
maxTRt = 24; % total # of test trials
fam_objs = 4;
t_maxt = floor((5000+1000)/scale_factor); %specify simulation time  % scale with Experiment test trial Duration
t_max = t_maxt;

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
    visuals_Off = floor(5000/scale_factor);
    words_On = floor([1500 3000]/scale_factor);%[;
    words_Off = floor(([1500 3000]+700)/scale_factor);
    %% DESIGN OBJECTS, FEARTURES and WORDS through random permutations
    Feature=[];
    for j=1:nFeatures
        chosen_indices = randperm(size(Property{1},2),nObjects);%% Create all 28 objects
        %chosen_indices = [chosen_indices chosen_indices(5:14)]; %% add  10 more objects
        TempFe=[];
        for i=1:nObjects
            TempFe = [TempFe  Property{j}(chosen_indices(i))];
            %TempFe = [TempFe  Property{j}(i)];
        end
        Feature{j}=TempFe;
    end
    Words=[]; 
    chosen_indices = randperm(size(Names,2), fam_objs);%% Create  4 words for familiarizatiom
    %for i=1:fam_objs
        Words = Names(chosen_indices);
    %end
    fam_seq = [];for i=1:maxTRt/fam_objs ; fam_seq = [fam_seq randperm(fam_objs)]; end
    
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
  

%% TESTING PHASE
    savestate_hcon_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwm_f{subject}=zeros(nFeatures,fieldSize_ftr);
    savestate_hwf{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_wd);
    savestate_hwm_c{subject}=zeros(nFeatures,fieldSize_ftr,fieldSize_spt);
    savestate_hcon_s{subject}=zeros(fieldSize_spt);
    savestate_hwm_s{subject}=zeros(fieldSize_spt);
    savestate_hword{subject}=zeros(fieldSize_wd);
    savestate_historyLt{subject}=zeros(maxTRt,t_maxt);
    savestate_historyRt{subject}=zeros(maxTRt,t_maxt);
    fam_word{subject}=zeros(maxTRt,2);
    side_balance=randperm(maxTRt);
    for trt=1:maxTRt
        pos=trt;
        %sim2.init(); 
        sim2.t =sim2.tZero;
        if trt > 1
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
        while sim2.t <= t_maxt

            t = sim2.t;
            %
            if isalmost(t, visuals_On,tolerance) 
                if mod(side_balance(trt),2)
                    for f=1:nFeatures                   
                         sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(fam_seq(pos))),spt_L, vstrength});   
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Right')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(pos+fam_objs)),spt_R, vstrength}); 
                    end
                    fam_word{subject}(trt,1)=  'L';
                else
                     for f=1:nFeatures
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Left')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(pos+fam_objs)),spt_L, vstrength});   
                        sim2.setElementParameters({strcat('Feature_',num2str(f),'_Right')}, {'positionY','positionX','amplitude'},... 
                        {cell2mat(Feature{f}(fam_seq(pos))),spt_R, vstrength});
                     end
                     fam_word{subject}(trt,1)=  'R';
                end
                fam_word{subject}(trt,2)=  fam_seq(pos); 
             end
         if Labelling_condition_ON == 1  
            if isalmost(t, words_On,tolerance)                 
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {cell2mat(Words(fam_seq(pos))),wstrength});               
            end
            if isalmost(t, words_Off,tolerance)              
                sim2.setElementParameters({'Word_Label_1'}, {'position','amplitude'},... 
                    {0,0}); 
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
    matTest.fam_word = fam_word{subject};
    catch
    disp('Error on saving test subject number ');
    end
 end