close all; clear all; 

%% To run a batch of simulations on a HPC/ multicore pc, set variable mode = 2 (below) and use 'parfor' subjects loop in the corresponding experiment file
%To run a gui simulation only, set variable mode = 1 (below) and use 'for' subjects loop in the experiment file
mode = 1; % set to 1 for auto/gui mode, 0 for singlebatch non-gui, and 2 for multicore batch mode (also need to switch to 'parfor' in the experiment file for mode=2)

%% sample code needed to run simulator on a high performance cluster
HPC = 0;
if HPC
    run('../../../COSIVINA/setpath.m') % add cosivina and jsoblab files to the matlab path
    addpath(genpath('../../WOLVES_PsychReview_2021'));
    cd('../') % move from the wolves_core folder where you'll launch the job to the main folder
    parpool('SlurmProfile1',96) %this will be HPC-specific
end

gui_speed=15; %update gui after every n iterations: ideal values from 1 to 20.
notes = ['any notes you want to take specific to this simulation'];% notes about any variable changes  
simNamePar = ['WPR_']; % give a name to your simulation.
%sim_seed = rng('shuffle'); %get seed for reproducibility, use rng(sim_seed) for reproduction
createComboSim; %%%create the model, for Kachergis et al task (taskvar==7) change to, createComboSimKachergis; 
if sim.loadSettings('wolvesPaperPR.json','changeable') == 0; disp('json file load ERROR!!'); end; % loads the parameters file
createComboGUI;% create and initialize GUI
createComboControls;% create and initialize GUI controls
if (mode == 2), numSubjects = 300; tolerance = 0; % specify the number of simulations/subjects to run. 300 default 
else,   numSubjects = 1; tolerance = 3; simNamePar = ['guiTest']; end % tolerance is used for gui's only to ensure they don't skip equality(==)conditionals

%% Update Memory Build and Decay Parameters
parBuild = 1200; parDecay = 8000;  
setMemoryTraceTimescales(sim, parBuild, parDecay);

%% TURN ON below Noise Parameters for MSHF/Mather novelty tasks (taskvars 14,15,16 below: CD target journal)
%noise_ior_s = 1.3; noise_wm_f = 2;
%setNoveltySpecificNoiseParams(sim, noise_ior_s, noise_wm_f);

%% Choose and Run the Experiment
taskvar = 1; %update taskvar value to simulate correspoding task/experiment from below: default Smith & Yu (2008,11)

if (taskvar==1)
    %% Task - Smith & Yu, Dev Sci, 2008 - Yu & Smith, Dev Sci, 2011 - Infants - standard cross-sit
    taskName = 'Smith_Yu_2008_2011';
    simName = [simNamePar,taskName,'_'];
    Smith_Yu_2008_2011;   
elseif (taskvar==3)
    %% Trueswell Medina Hafri & Gleitman Cognitive Psychology(2013) - Adults - Eyetracking
    taskName = 'Trueswell_Medina_Hafri_Gleitman_2013';
    simName = [simNamePar,taskName,'_'];
    Trueswell_Medina_Hafri_Gleitman_2013;
elseif (taskvar==402)
    %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 2x2 condition
     taskName = 'Yu_Smith_Two_2007';
     simName = [simNamePar,taskName,'_'];
     Yu_Smith_Two_2007;
elseif (taskvar==403)
    %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 3x3 condition    
    taskName = 'Yu_Smith_Three_2007';
    simName = [simNamePar,taskName,'_'];
    Yu_Smith_Three_2007;
 elseif (taskvar==404)
    %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 4x4 condition         
    taskName = 'Yu_Smith_Four_2007';
    simName = [simNamePar,taskName,'_'];
    Yu_Smith_Four_2007;
elseif (taskvar==5)
    %% Yu Zhong & Fricker, Frontiers in Psychology(2012) - Adults - Eyetracking
    taskName = 'Yu_Zhong_Fricker_2012';
    simName = [simNamePar,taskName,'_'];
    Yu_Zhong_Fricker_2012;  
elseif (taskvar==6)
    %% Task - Yurovsky_Yu_Smith, Cognitive Science (2013) - Adults - Competitive Processes
    taskName = 'Yurovsky_Yu_Smith_2013';
    simName = [simNamePar,taskName,'_'];
    Yurovsky_Yu_Smith_2013; % use 300x3 conditions numSubjects 
elseif (taskvar==7)
    %% Task - Kachergis_Yu_Shiffrin, Psychon Bull Rev (2012) - Adults- Prior Learning, ME
    taskName = 'Kachergis_Yu_Shiffrin_2012_11AFC';
    simName = [simNamePar,taskName,'_'];
    Kachergis_Yu_Shiffrin_2012_11AFC; % use 300x12 conditions numSubjects 
elseif (taskvar==8)
    %% Task - Smith & Yu, Lang Learn Dev, (2013) - Infants - novelty trap
    taskName = 'Smith_Yu_2013';
    simName = [simNamePar,taskName,'_'];
    Smith_Yu_2013;
elseif (taskvar==9)
    %% Tasks - Vlach & Johnson, Cognition (2013); Vlach & DeBrock, JML, (2017)
    %% Vlach & DeBrock, JEP:LMC (2019) - Infants - massed vs interleaved
    taskName = 'Vlach_CSWL_2013_17_19';
    simName = [simNamePar,taskName,'_'];
    Vlach_CSWL_2013_17_19;
elseif (taskvar==11)
    %% Task - Vlach & DeBrock, JML, 2017 - Infants - CSWL vs Memory tests
    taskName = 'Vlach_DeBrock_WOB_2017';
    simName = [simNamePar,taskName,'_'];
    Vlach_DeBrock_WOB_2017;% run this & Taskvar 9 on multiple tau_Decay values
elseif (taskvar==12)
    %% Task - Suanda_Mugwanya_Namy, Jrl Exp Child Psy (2014) - 6 year olds - Referential Ambiguity
    taskName = 'Suanda_Mugwanya_Namy_2014';
    simName = [simNamePar,taskName,'_'];
    Suanda_Mugwanya_Namy_2014; % use 300x3 conditions numSubjects 
elseif (taskvar==13)
    %% Task - Fitneva & Christiansen, Cognitive Science (2015) - Partial Knowledge, Initial Accuracy
    taskName = 'Fitneva_Christiansen_2015';
    simName = [simNamePar,taskName,'_'];
    Fitneva_Christiansen_2015; % use 300x2 conditions numSubjects
elseif (taskvar==14)
    %% Mather, Schafer & Houston_Price,British Journal of Dev Psy (2011) - Infants - Eyetracking - Fam vs Novel objects/label vs not
    taskName = 'MSHP_Exp1_2011';
    Labelling_condition_ON = 0;%%set to 0 for Silent condition, 1 for Labelling 
    prediction_mode_ON = 0;%% set to 0 for same word repetition as in the orginal experiment; set to 1 for generating novel names on every new trial
    %sim naming 
    if Labelling_condition_ON, simName = [simNamePar,'Labelling_', taskName,'_'];
    else, simName = [simNamePar,'Silent_', taskName,'_'];end
    MSHP_Exp1_2011;
elseif (taskvar==15)
    %% Mather & Plunkett,Cognitive Science (2012) - Infants - Eyetracking - Novel word Novel Object 
    taskName = 'Mather_Plunkett_Exp1_2012';
    simName = [simNamePar,taskName,'_'];
    Mather_Plunkett_Exp1_2012;
 elseif (taskvar==16)
    %% Mather, Schafer & Houston_Price Exp2, British Journal of Dev Psy (2011) - Infants - Eyetracking - Fam vs Novel objects/label vs not
    taskName = 'MSHP_Exp2_2011';  
    Labelling_condition_ON = 0;%%set to 0 for Silent condition, 1 for Labelling 
    %sim naming 
    if Labelling_condition_ON, simName = [simNamePar,'Labelling_', taskName,'_'];
    else, simName = [simNamePar,'Silent_', taskName,'_'];end
    MSHP_Exp2_2011;
end

%% Save Simulation Results 
Results_Saving;