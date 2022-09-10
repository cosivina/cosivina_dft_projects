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
trial_check = 10; %set the stimulus positions to their 'default' values every X trials to probe WM
samplingHistory = 5; %sample the field activity every X timesteps
historyDuration = ceil((t_max*n_trials)/samplingHistory);

% create model structure and initialize
OneFieldSim_Simple_CurrentDirections;
sim.init();

% if auto mode, create GUI for visualization
% If batch mode, GUI excluded for faster computation
if mode == 1
    n_reps = 1;
    OneFieldGUI_CurrentDirections;
    gui.init(); %just initialize this once
end

%% setting up the simulator
sim.loadSettings('presetOneLayerField_stabilized_simple.json');

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

for rep = 1 : n_reps
    sim2 = sim;
    
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
            
            % monitor GUI buttons
            if ~gui.pauseSimulation
                sim2.step();
            end
            if gui.quitSimulation
                gui.close();
                break;
            end
            gui.step();%so that parameter buttons work
            
        end % time loop        
        
    end % trial loop

    
end % reps loop

