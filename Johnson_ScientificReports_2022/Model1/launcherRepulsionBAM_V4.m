%% setting up the simulator
close all; clear all;

HPC=1;

if HPC
    run('../COSIVINA/setpath.m') % add cosivina and jsoblab files to the matlab path
    addpath(genpath('../MemToolbox'));
    parpool('SlurmProfile1',96)
end

headerflag = 1;
SimName = 'M1-V9';
OutName = [SimName, datestr(now, '_yyyy-mm-dd-THHMMSS') '.csv'];
OutName2 = [SimName, '_errors', datestr(now, '_yyyy-mm-dd-THHMMSS') '.mat'];

mode = 2; %0 for batch mode; 1 for auto mode with visualization; 2 for multicore (change for to parfor around line 47)

if mode == 1
    n_reps = 1;
    nsubj = 1;
    finalplot = 0;
    storage = 0;
else
    n_reps = 160;
    nsubj = 12;
    finalplot = 0;
    storage = 1;
end
gui_speed = 50; %display every X timesteps

%% params that differ relative to Model2
flatstr = 4.5;
wheelstr = 4.5;
visboost = 1.5;
atnboost = 3;
wnoise = 4;

%% shared parameters
fieldSize = 361; %401;
wheelSize = 2*fieldSize - 1;
sigma_exc = 5;
sigma_inh = 10;
units = 1:fieldSize;

stimstr = 12;
visnoise = 17.5;
visnoisesig = 2;
wheelnoise = 0.25;

%% timing -- scale using tau ratio of 40/12 = 3.333
presentation = fix(444*3.333); %800ms in exp -- 1 timestep = .54ms
delay = fix(556*3.333); %1000ms
attendwheel = fix(20*3.333);
response = fix(334*3.333); %600ms

flat_pattern = zeros(1, fieldSize);

%% run the simulator
createRepulsionSim2;
sim.loadSettings('repulsion3a_V3.json');
sim.init();

createRepulsionGUI_V2;
if (mode == 1)
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
    
    for j = 1:ncond
        
        fprintf('TrialType %i\n',j);
        
        %for i = 1:n_reps
        parfor i = 1:n_reps
            
            sim2 = 5;
            if (mode == 0) || (mode == 1)
                sim2 = sim;
            elseif mode == 2
                sim2 = sim.copy();
            end
            
            sim2.init();
            
            sim2.setElementParameters('atn -> vis', 'sigma', 30);
            sim2.setElementParameters('atn -> vis', 'amplitude', 5);
            sim2.setElementParameters('vis -> vis (global)', 'amplitude', -0.0175); %-0.005 in param file
            
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
            sim2.setElementParameters('noise kernel vis', 'amplitude', visnoise);
            sim2.setElementParameters('noise kernel vis', 'sigmaX', visnoisesig);
            sim2.setElementParameters('noise kernel vis', 'sigmaY', visnoisesig);
            
            while sim2.t <= presentation + delay + attendwheel + response
                
                t = sim2.t;
                
                if t == presentation
                    sim2.setElementParameters('vis stim 1', 'amplitude', 0);
                    sim2.setElementParameters('vis stim 2', 'amplitude', 0);
                    sim2.setElementParameters('vis stim 3', 'amplitude', 0);
                elseif t == presentation + delay
                    sim2.setElementParameters('field vis', 'h', visboost-5);
                    sim2.setElementParameters('field atn', 'h', atnboost-7);
                    switch j
                        case 1
                            sim2.setElementParameters('scale flat boost', 'amplitude', flatstr);
                        case 2
                            sim2.setElementParameters('scale flat boost 2', 'amplitude', flatstr);
                        case 3
                            sim2.setElementParameters('scale flat boost 3', 'amplitude', flatstr);
                        case 4
                            sim2.setElementParameters('scale flat boost 3', 'amplitude', flatstr);
                    end
                elseif t == presentation + delay + attendwheel
                    sim2.setElementParameters('wheel pattern', 'amplitude', wheelstr);
                    sim2.setElementParameters('wheel noise', 'amplitude', wheelnoise);
                elseif t == presentation + delay + attendwheel + response
                    sim2.setElementParameters('field vis', 'h', -5);
                    sim2.setElementParameters('field atn', 'h', -7);
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
        unique = MLE(errors(3,:), model);
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

end %subj

if storage
    save(OutName2,'errors_subj');
end
