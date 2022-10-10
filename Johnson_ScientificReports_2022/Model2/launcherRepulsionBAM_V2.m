%% setting up the simulator
close all; clear all;

HPC=0;

if HPC
    run('../COSIVINA/setpath.m') % add cosivina and jsoblab files to the matlab path
    addpath(genpath('../MemToolbox'));
    parpool('SlurmProfile1',96)
    %set different rand seed for each parallel process
    x = feature('getpid');
    s = RandStream('mt19937ar', 'seed', sum(x*clock));
    RandStream.setGlobalStream(s)
else
    randn('state',sum(100*clock));
    rand('state',sum(100*clock));
end

for exp=0:1 %0 = main exp; 1 = control exp
    
    headerflag = 1;
    SimName = 'M2-V40';
    if exp == 1
        SimName = [SimName '_control'];
    end
    OutName = [SimName, datestr(now, '_yyyy-mm-dd-THHMMSS') '.csv'];
    OutName2 = [SimName, '_errors', datestr(now, '_yyyy-mm-dd-THHMMSS') '.mat'];
    
    mode = 1; %0 for batch mode; 1 for auto mode with visualization; 2 for multicore (change for to parfor around line 47)
    
    if mode == 1
        n_reps = 1;
        nsubj = 1;
        finalplot = 0;
        storage = 0;
    else
        n_reps = 160; %160
        nsubj = 12;
        finalplot = 0;
        storage = 1;
    end
    gui_speed = 1; %display every X timesteps

    %% params that differ relative to Model1
    flatstr = 13; %12;
    wheelstr = 4; %4.5;
    visboost = 4;
    atnboost = 2.9; %2.8; %3
    wnoise = 3.5; %2.0;
  
    %% shared parameters
    fieldSize = 361;
    wheelSize = 2*fieldSize - 1;
    sigma_exc = 5;
    sigma_inh = 10;
    units = 1:fieldSize;    
    
    stimstr = 12;
    visnoise = 17.5;
    visnoisesig = 2;
    wheelnoise = 0.25;

    wnoise_sigma = 1.5;
    iorslowboost = -4;
    atnsrboost = 2.2;
    ulowboost = -12;
    wlowboost = -9; % -8;
    
    pdlow = -50;
    coslow = -50;
    
    ScaleTime=1; 
    
    if exp == 0
        presentation = fix(444*ScaleTime); %500; 800ms in exp so 1 timestep = 1.8ms
        delay = fix(556*ScaleTime); %626; 1000ms in exp   must be even
        attendwheel = fix(20*ScaleTime); %10
        stim_off = 0; %not used in main exp -- only used in control
        response = fix(334*ScaleTime); %376; %600ms in exp? until response; must be even
    else
        presentation = fix(278*ScaleTime); %500ms initial presentation; 1 timestep = 1.8ms
        delay = 0; %must be even
        attendwheel = fix(20*ScaleTime); %10
        stim_off = fix(166*ScaleTime); %300ms
        response = fix(334*ScaleTime); %376; %600ms in exp? until response; must be even
    end
    
    flat_pattern = zeros(1, fieldSize);
    
    %% run the simulator
    createRepulsionSim;
    sim.loadSettings('repulsion2020-2b.json'); 
    sim.init();
    
    createRepulsionGUI;
    if mode == 1
        gui.init();
    end
    
    ncond = 5; %% cond 5 runs the 'extra' unique trials
    target = [90, 120, 270, 270, 270];
    target_space = [180, 90, 270, 270, 270];
    
    errors_subj = zeros(nsubj, ncond, n_reps);
    
    for subj=1:nsubj
        
        errors = zeros(ncond, n_reps);
        wfinal = zeros(ncond, n_reps, fieldSize);
        finalDist = zeros(ncond, n_reps, fieldSize);
        
        tic
        
        for j = 1:ncond %run 160 in each condition (1-4) then run the extra trials to unique (5)...
            
            fprintf('TrialType %i\n',j);
            
            for i = 1:n_reps
            %parfor i = 1:n_reps
                
                sim2 = 5;
                if (mode == 0) || (mode == 1)
                    sim2 = sim;
                elseif mode == 2
                    sim2 = sim.copy();
                end
                                
                sim2.init();
                                
                if j ~= 4
                    sim2.setElementParameters('vis stim 1', 'amplitude', stimstr);
                    sim2.setElementParameters('vis stim 2', 'amplitude', stimstr);
                    sim2.setElementParameters('vis stim 1', 'positionX', target_space(1));
                    sim2.setElementParameters('vis stim 2', 'positionX', target_space(2));
                    sim2.setElementParameters('vis stim 1', 'positionY', target(1));
                    sim2.setElementParameters('vis stim 2', 'positionY', target(2));
                else
                    sim2.setElementParameters('vis stim 1', 'amplitude', 0);
                    sim2.setElementParameters('vis stim 2', 'amplitude', 0);
                end
                sim2.setElementParameters('vis stim 3', 'amplitude', stimstr);
                sim2.setElementParameters('vis stim 3', 'positionX', target_space(3));
                sim2.setElementParameters('vis stim 3', 'positionY', target(3));
                sim2.setElementParameters('noise kernel w', 'amplitude', wnoise);
                sim2.setElementParameters('noise kernel w', 'sigma', wnoise_sigma);
                sim2.setElementParameters('noise kernel vis', 'amplitude', visnoise);
                sim2.setElementParameters('noise kernel vis', 'sigmaX', visnoisesig);
                sim2.setElementParameters('noise kernel vis', 'sigmaY', visnoisesig);
                
                while sim2.t <= presentation + delay + attendwheel + response
                    
                    t = sim2.t;
                    
                    if (t == presentation) && (exp == 0)
                        sim2.setElementParameters('vis stim 1', 'amplitude', 0);
                        sim2.setElementParameters('vis stim 2', 'amplitude', 0);
                        sim2.setElementParameters('vis stim 3', 'amplitude', 0);
                    elseif t == presentation + delay
                        sim2.setElementParameters('field vis', 'h', visboost-5);
                        sim2.setElementParameters('field ior_s', 'h', iorslowboost-5);
                        sim2.setElementParameters('field atn_sr', 'h', atnsrboost-5);
                        sim2.setElementParameters('field atn', 'h', atnboost-7);
                        sim2.setElementParameters('field u', 'h', ulowboost-7);
                        sim2.setElementParameters('field w', 'h', wlowboost-4);
                        sim2.setElementParameters('cos', 'h', coslow-5);
                        sim2.setElementParameters('pd_atnc', 'h', pdlow-5);
                        switch j
                            case 1
                                sim2.setElementParameters('scale flat boost', 'amplitude', flatstr);
                                if exp == 1
                                    sim2.setElementParameters('vis stim 2', 'amplitude', stimstr*0.3);
                                    sim2.setElementParameters('vis stim 3', 'amplitude', stimstr*0.3);
                                end
                            case 2
                                sim2.setElementParameters('scale flat boost 2', 'amplitude', flatstr);
                                if exp == 1
                                    sim2.setElementParameters('vis stim 1', 'amplitude', stimstr*0.3);
                                    sim2.setElementParameters('vis stim 3', 'amplitude', stimstr*0.3);
                                end
                            case 3
                                sim2.setElementParameters('scale flat boost 3', 'amplitude', flatstr);
                                if exp == 1
                                    sim2.setElementParameters('vis stim 1', 'amplitude', stimstr*0.3);
                                    sim2.setElementParameters('vis stim 2', 'amplitude', stimstr*0.3);
                                end
                            case 4
                                sim2.setElementParameters('scale flat boost 3', 'amplitude', flatstr);
                                if exp == 1
                                    sim2.setElementParameters('vis stim 1', 'amplitude', stimstr*0.3);
                                    sim2.setElementParameters('vis stim 2', 'amplitude', stimstr*0.3);
                                end
                            case 5
                                sim2.setElementParameters('scale flat boost 3', 'amplitude', flatstr);
                                if exp == 1
                                    sim2.setElementParameters('vis stim 1', 'amplitude', stimstr*0.3);
                                    sim2.setElementParameters('vis stim 2', 'amplitude', stimstr*0.3);
                                end
                        end
                    elseif t == presentation + delay + attendwheel
                        sim2.setElementParameters('wheel pattern', 'amplitude', wheelstr);
                        sim2.setElementParameters('wheel noise', 'amplitude', wheelnoise); %was 1/wheelstr
                    elseif (t == presentation + delay + stim_off) && (exp == 1)
                        sim2.setElementParameters('vis stim 1', 'amplitude', 0);
                        sim2.setElementParameters('vis stim 2', 'amplitude', 0);
                        sim2.setElementParameters('vis stim 3', 'amplitude', 0);
                    elseif t == presentation + delay + attendwheel + response
                        sim2.setElementParameters('field vis', 'h', -5);
                        sim2.setElementParameters('field ior_s', 'h', -5);
                        sim2.setElementParameters('field atn_sr', 'h', -5);
                        sim2.setElementParameters('field atn', 'h', -7);
                        sim2.setElementParameters('field u', 'h', -7);
                        sim2.setElementParameters('field w', 'h', -4);
                        sim2.setElementParameters('cos', 'h', -5);
                        sim2.setElementParameters('pd_atnc', 'h', -5);
                        sim2.setElementParameters('scale flat boost', 'amplitude', 0);
                        sim2.setElementParameters('scale flat boost 2', 'amplitude', 0);
                        sim2.setElementParameters('scale flat boost 3', 'amplitude', 0);
                        sim2.setElementParameters('wheel pattern', 'amplitude', 0);
                        sim2.setElementParameters('wheel noise', 'amplitude', 0);
                    end
                    
                    if mode == 1 & (mod(t,gui_speed) == 0)
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
                
                wfinal(j,i,:) = sim2.getComponent('field w', 'activation');
                
                X = sim2.getComponent('sum field vis', 'verticalSum');
                peakloc = fieldSize - X*units'/sum(X) + 1; %need to add 1 since matlab matrices start at 1
                
                Y = sim2.getComponent('sum field vis', 'horizontalSum');
                finalDist(j,i,:) = Y;
                
                errors(j, i) = target(j)-peakloc;
                
            end
            
            if finalplot == 1
                for i = 1:n_reps
                    figure
                    plot(squeeze(wfinal(j,i,:))');
                    hold;
                    plot(squeeze(finalDist(j,i,:))','r');
                    hold;
                end
            end
        end
        
        toc
        
        if storage
            
            if ncond == 5
                errors_unique = zeros(1,n_reps*2);
                errors_unique(1,1:n_reps) = errors(3,:);
                errors_unique(1,n_reps+1:n_reps*2) = errors(5,:);
            else
                errors_unique = zeros(1,n_reps);
                errors_unique = errors(3,:);
            end
            
            model = WithBias(StandardMixtureModel);
            SS1 = MLE(errors(4,:), model);
            unique = MLE(errors_unique(1,:), model);
            CW = MLE(errors(1,:), model);
            CCW = MLE(errors(2,:), model);
            
            SS1(2) = 1-SS1(2);
            unique(2) = 1-unique(2);
            CW(2) = 1-CW(2);
            CCW(2) = 1-CCW(2);
            
            OutFile = fopen(OutName,'a');
            if(headerflag)
                headerflag = 0;
                fprintf(OutFile,'Participant,Condition,mean,P_m,StDev\n');
            end
            fprintf(OutFile,'%i', subj);
            fprintf(OutFile,',SS1');
            fprintf(OutFile,',%6.4f', SS1);
            fprintf(OutFile,'\n');
            fprintf(OutFile,'%i', subj);
            fprintf(OutFile,',Unique');
            fprintf(OutFile,',%6.4f', unique);
            fprintf(OutFile,'\n');
            fprintf(OutFile,'%i', subj);
            fprintf(OutFile,',CW');
            fprintf(OutFile,',%6.4f', CW);
            fprintf(OutFile,'\n');
            fprintf(OutFile,'%i', subj);
            fprintf(OutFile,',CCW');
            fprintf(OutFile,',%6.4f', CCW);
            fprintf(OutFile,'\n');
            fclose(OutFile);
            
        end
        
        errors_subj(subj,:,:) = errors;
        
    end %subj loop
    
    if storage
        save(OutName2,'errors_subj');
    end
    
end %exp loop