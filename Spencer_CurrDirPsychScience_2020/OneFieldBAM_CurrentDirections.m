%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COSIVINA BAM File for model used in:
%
% Spencer, J. P. (2020). The development of working memory.
% Current Directions in Psychological Science, doi/10.1177/0963721420959835.
%
% For supporting documentation and videos, see www.dynamicfieldtheory.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

% shared parameters
fieldSize = 100;
sigma_exc = 5;
sigma_inh = 12.5;

StimOn=100; %turn stimulus 'on' 100 time steps into trial
StimOff=300; %turn stimulus 'off' 300 time steps into trial
t_max=800; %the duration of each trial
StimPos=[20, 50, 80]; %stimulus positions
StimSep=40; %present targets at random separation from StimPos above
StimStr=5; %stimulus strength

% mode control
mode = 1; % 0 for batch mode, 1 for auto mode with visualization, 2 for multicore (also need to switch to 'parfor' below)
n_reps = 1; % repetitions for batch mode, overwritten to 1 for auto
n_trials = 30; %1 for Fig1 simulations; 30 for Fig2 simulations 

% simulator options
plotResults = 0; %set to 1 if final mesh shadow plot desired
saveResults = 0; %set to 1 if output/results file is desired
resultBaseFilename = 'OneFieldResult';
trial_check = 10; %set the stimulus positions to their 'default' values every X trials to probe WM
samplingHistory = 5; %sample the field activity every X timesteps
historyDuration = ceil((t_max*n_trials)/samplingHistory);

% create model structure and initialize
OneFieldSim_CurrentDirections;
sim.init();

% if auto mode, create GUI for visualization
% If batch mode, GUI excluded for faster computation
if mode == 1
    n_reps = 1;
    OneFieldGUI_CurrentDirections;
    gui.init(); %just initialize this once
end

% specify the seeding of the random number generator as desired
% this can be useful when comparing simulations with memory traces...
if 0
    seed = 1;
    s = RandStream('mt19937ar', 'Seed', seed);
    RandStream.setGlobalStream(s);
else
    rand('state', sum(100*clock));
    randn('state', sum(100*clock));
end

% name the output file
if saveResults
    resultFilename = [resultBaseFilename, datestr(now, 'yyyy-mm-dd-THHMMSS') '.mat'];
end

%% setting up the simulator
sim.loadSettings('presetOneLayerField_stabilized.json');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to reproduce Fig1, run a single trial, set mem u -> u to 0, and set
% u -> u to desired developmental value
%
% to reproduce Fig2, run 30 trials, set mem u -> u to 1, and set
% u -> u to 24
sim.setElementParameters('u -> u', 'amplitudeExc', 24); %24 young, 25 middle, 26 old
sim.setElementParameters('mem u -> u', 'amplitude', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%other parameters I modified relative to the default parameter file above
sim.setElementParameters('noise kernel', 'amplitude', 2);
sim.setElementParameters('u -> u', 'amplitudeGlobal', -0.05);
sim.setElementParameters('u -> u', 'amplitudeInh', 17);

%initialize a data structure to store the results
results=zeros(n_reps, historyDuration, fieldSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toggle the next lines if multi-core mode (mode 2) desired

for rep = 1 : n_reps
%parfor rep = 1 : n_reps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sim2 = 5;
    
    if (mode == 0) || (mode == 1)
        sim2 = sim;
    elseif mode == 2
        sim2 = sim.copy();
    end
    
    % reset model at the start of each rep
    sim2.init();
       
    %embed a sequence of trials here...allows build-up of memory trace.
    for tct = 1: n_trials
        
        sim2.t=1; 
        StimPosNew = [StimPos(1)+((rand-0.5)*StimSep), StimPos(2)+((rand-0.5)*StimSep), StimPos(3)+((rand-0.5)*StimSep)];
        
        %use default stimulus positions every X trials
        if ((mod(tct,trial_check) == 0) | (tct == 1))
            StimPosNew = StimPos;
        end
        
        while sim2.t <= t_max
            
            % reset the field (wipe out all peaks) for trials 2 onward
            if (tct > 1) & (sim2.t < 50)
                sim2.setElementParameters('stimulus s_1', 'amplitude', -50);
            else
                sim2.setElementParameters('stimulus s_1', 'amplitude', 0);
            end
            
            %turn the stimuli on and off
            if (sim2.t >= StimOn) & (sim2.t < StimOff)
                sim2.setElementParameters({'stimulus 1', 'stimulus 2', 'stimulus 3'}, ...
                    {'amplitude','amplitude','amplitude'} , {StimStr, StimStr, StimStr});
                sim2.setElementParameters({'stimulus 1', 'stimulus 2', 'stimulus 3'}, ...
                    {'position','position','position'} , ...
                    {StimPosNew(1), StimPosNew(2), StimPosNew(3)});
            elseif (sim2.t >= StimOff)
                sim2.setElementParameters({'stimulus 1', 'stimulus 2', 'stimulus 3'}, ...
                    {'amplitude','amplitude','amplitude'} , {0, 0, 0});
            end
            
            % if in auto mode, monitor GUI buttons
            if mode == 1
                if ~gui.pauseSimulation
                    sim2.step();
                end
                if gui.quitSimulation
                    gui.close();
                    break;
                end
                gui.step();%so that parameter buttons work                
            else
                sim2.step();
            end
            
        end % time loop        
        
    end % trial loop

    %store results
    results(rep,:,:) = sim2.getComponent('history', 'output');
    
end % reps loop

if saveResults
    save(resultFilename,'results');
end

% plot results
if plotResults
    for rep = 1 : n_reps
        figure('color','w')
        mesh(squeeze(results(rep,:,:)));
        shadowplot x
        shadowplot y
        set(gca,'zlim',[-10, 20],'view',[-122,11],'fontsize',18,'color','w')
        xlabel('feature', 'Rotation', -18)
        ylabel('time', 'Rotation', 10)
        zlabel('activation u(x,t)')
        colormap default
    end
end

% if in auto mode, close the GUI window after the trial is over
% or comment this out if you want to keep it open...
if mode == 1
%    gui.close();
end
