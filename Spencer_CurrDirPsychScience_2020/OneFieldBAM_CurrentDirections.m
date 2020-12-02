%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COSIVINA BAM File for OneField Figs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;


% shared parameters
fieldSize = 100;
sigma_exc = 5;
sigma_inh = 12.5;

StimOn=100;
StimOff=300;
StimPos=[20, 50, 80];
StimStr=5;
StimStr2=0; %1.5;
mode = 1; % 0 for batch mode, 1 for auto mode with visualization, 2 for multicore (also need to switch to 'parfor' below)
n_reps = 1; % repetitions for batch mode, overwritten to 1 for auto
n_trials = 30;
trial_check = 10; 

t_max=800;
samplingHistory = 5;
historyDuration = ceil((t_max*n_trials)/samplingHistory);


% create model structure and initialize
OneFieldSim_CurrentDirections;
sim.init();

% if auto mode, create GUI for visualization
% If batch mode, GUI excluded for faster computation
% if mode == 1
    OneFieldGUI_CurrentDirections;
    %    gui.addVisualization(TimeDisplay(), [4, 2], [1, 1], 'control');
    %    n_reps = 1;
% end

if 0
    seed = 1;
    s = RandStream('mt19937ar', 'Seed', seed);
    RandStream.setGlobalStream(s);
else
    rand('state', sum(100*clock));
    randn('state', sum(100*clock));
end


resultBaseFilename = 'OneFieldResult';
if 1 %saveResults
    resultFilename = [resultBaseFilename, datestr(now, 'yyyy-mm-dd-THHMMSS') '.mat'];
end


%% setting up the simulator
sim.loadSettings('presetOneLayerField_stabilized.json');

%%toggle the line below to turn lateral interactions off
sim.setElementParameters('u -> u', 'amplitudeGlobal', -0.05);
sim.setElementParameters('u -> u', 'amplitudeExc', 24); %24 young, 25 middle, 26 old
sim.setElementParameters('u -> u', 'amplitudeInh', 17);

% % % "amplitudeExc": 17.5 v 30
% % % "sigmaInh": 10
% % % "amplitudeInh": 15 v 27.5
% % % "amplitudeGlobal": 0

%%adjust parameters
sim.setElementParameters('mem u -> u', 'amplitude', 1);
sim.setElementParameters('noise kernel', 'amplitude', 2);

% create a set of stimulus positions
firstcase = 1;

for params = 1:1
    
    for rep = 1 : n_reps
        %parfor rep = 1 : n_reps
        
        sim2 = 5;
        
        if (mode == 0) || (mode == 1)
            sim2 = sim;
        elseif mode == 2
            sim2 = sim.copy();
        end
        
        % reset model at the start of each rep and initialize variables
        sim2.init();
        if (mode == 1) && (firstcase == 1) && (rep == 1)
            gui.init(); %just initialize this once
            firstcase = 0;
        end
        sim2.t=1; %initialize time to 1--otherwise, get matrix index errors below
        
        for tct = 1: n_trials
            
            %embed a sequence of trials here...allows build-up of memory trace.
            sim2.t=1; %initialize time to 1--otherwise, get matrix index errors below
            StimPosNew = [StimPos(1)+((rand-0.5)*40), StimPos(2)+((rand-0.5)*40), StimPos(3)+((rand-0.5)*40)];

            if (mode == 0) && ((mod(tct,trial_check) == 0) | (tct == 1))
                StimPosNew = StimPos;
            end
            
            while sim2.t <= t_max
                                
                % code below is only used to track activation at one site and
                % display in the GUI as stimulus s_1...
% %                 unow = sim2.getComponent('field u', 'activation');
% %                 sim2.setElementParameters('stimulus s_1', 'amplitude', unow(1,StimPos(1)));
                if (tct > 1) & (sim2.t < 50)
                    sim2.setElementParameters('stimulus s_1', 'amplitude', -50);
                else 
                    sim2.setElementParameters('stimulus s_1', 'amplitude', 0);
                end
                if (sim2.t >= StimOn) & (sim2.t < StimOff)
                    sim2.setElementParameters({'stimulus 1', 'stimulus 2', 'stimulus 3'}, ...
                        {'amplitude','amplitude','amplitude'} , {StimStr, StimStr, StimStr});
                    sim2.setElementParameters({'stimulus 1', 'stimulus 2', 'stimulus 3'}, ...
                        {'position','position','position'} , ...
                        {StimPosNew(1), StimPosNew(2), StimPosNew(3)});
                elseif (sim2.t >= StimOff)
                    sim2.setElementParameters({'stimulus 1', 'stimulus 2', 'stimulus 3'}, ...
                        {'amplitude','amplitude','amplitude'} , {StimStr2, StimStr2, StimStr2});
                end
                
                % if in auto mode, monitor pause and quit buttons
                if mode == 1
                    if ~gui.pauseSimulation
                        sim2.step();
                    end
                    if gui.quitSimulation
                        gui.close();
                        break;
                    end
                    
                    %gui.checkAndUpdateControls();
                    %gui.updateVisualizations();
                    % if in batch mode, there are no buttons to monitor
                    gui.step();%so that parameter buttons work
                    
                else
                    sim2.step();
                end
                
            end % time loop
            
            if (mode == 0) && ((mod(tct,trial_check) == 0) | (tct == 1))
                gui.init()
                gui.updateVisualizations();
            end

            
        end % trial loop
        
    end % reps loop
    
end %params loop

% analysis of results
% % data2mesh = sim2.getComponent('history', 'output');
% % figure('color','w')
% % mesh(data2mesh);
% % shadowplot x
% % shadowplot y
% % set(gca,'zlim',[-10, 10],'view',[-122,11],'fontsize',18,'color','w')
% % xlabel('feature', 'Rotation', -18)
% % ylabel('time', 'Rotation', 10)
% % zlabel('activation u(x,t)')
% % colormap default

% if in auto mode, close the GUI window after the trial is over
if mode == 1
    %gui.close();
else
    %display gui at end
    %gui.init();
end
