close all; clear all; 
%%%%%%run('..\COSIVINA\setpath.m') % add cosivina and jsoblab files to the matlab path
%pc = parcluster('local'); poolobj = parpool(pc, 8);
mode =2; %1 = auto/gui, 0 = single batch; 2 for multicore batch mode (also need to switch to 'parfor' in the experiment file)
gui_speed=50; %update gui after every n iterations: ideal values from 1 to 20.
notes = ['wolvesPaperPR.json, no other changes'];% notes for experimenting simulations
simNamePar = ['Modified11AFC_batch1']; % give a name to your simulation.

createComboSimKachergis; % create the model
if sim.loadSettings('wolvesPaperPR.json','changeable') == 0; disp('json file load ERROR!!'); end; % load the parameters file
sim.init(); %initialize the model
createComboGUI;% create and initialize GUI
createComboControls;% create and initialize GUI controls

%% To run a batch of simulations on a HPC/ multicore pc, change variable mode =2 (above). 
if (mode == 2), numSubjects = 480; tolerance = 0; % specify the number of simulations/subjects to run. 300 default 
else,   numSubjects = 1; tolerance = 3; end % tolerance is used for gui's only to ensure they don't skip equality(==)conditionals
 
%% Code to Modify Memory Build and Decay Parameters %% Uncomment to use
           parBuild=1000; parDecay=15000;
%           sim.setElementParameters({'hwm_s'}, {'tauBuild'}, parBuild); 
%           sim.setElementParameters({'hcon_s'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hword'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hcon_f1'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hwm_f1'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hwm_c1'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hwf1'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hcon_f2'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hwm_f2'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hwm_c2'}, {'tauBuild'}, parBuild);
%           sim.setElementParameters({'hwf2'}, {'tauBuild'}, parBuild);
% 
%           sim.setElementParameters({'hwm_s'}, {'tauDecay'}, parDecay); 
%           sim.setElementParameters({'hcon_s'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hword'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hcon_f1'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hwm_f1'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hwm_c1'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hwf1'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hcon_f2'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hwm_f2'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hwm_c2'}, {'tauDecay'}, parDecay);
%           sim.setElementParameters({'hwf2'}, {'tauDecay'}, parDecay);
%%     

%% Choose and Run the Experiment
taskvar = 6; %choose the task/experiment (taskvar value) to simulate: default Smith & Yu (2008,11)

if (taskvar==1)
    %% Task - Smith & Yu, Dev Sci, 2008 - Yu & Smith, Dev Sci, 2011 - Infants - standard cross-sit
    simName = [simNamePar '_Smith_Yu_2008_'];
    Trials90_Smith_Yu_2008;          
elseif (taskvar==2)
    %% Task - Smith & Yu, Lang Learn Dev, 2013 - Infants - novelty trap
    simName = [simNamePar '_Smith_Yu_2013_'];
    Smith_Yu_2013;
elseif (taskvar==3)
    %% Task - Vlach & Johnson, Cognition, 2013 - Infants - massed vs interleaved
    simName = [simNamePar '_Vlach_Johnson_2013_'];
    Vlach_Johnson_2013;
elseif (taskvar==401)
    %% Task - Vlach & DeBrock, JML, 2017 - Infants - CSWL vs Memory tests
    simName = [simNamePar '_Vlach_DeBrockOR_2017_'];
    Vlach_DeBrockOR_2017;
elseif (taskvar==402)
    %% Task - Vlach & DeBrock, JML, 2017 - Infants - CSWL vs Memory tests
    simName = [simNamePar '_Vlach_DeBrockWR_2017_'];
    Vlach_DeBrockWR_2017;
elseif (taskvar==403)
    %% Task - Vlach & DeBrock, JML, 2017 - Infants - CSWL vs Memory tests
    simName = [simNamePar '_Vlach_DeBrockWOB_2017_'];
    Vlach_DeBrockWOB_2017;
elseif (taskvar==5)
    %% Task - Fitneva & Christiansen, Cognitive Science (2015) - Partial Knowledge, Initial Accuracy
    simName = [simNamePar '_Fitneva_Christiansen_2015_'];
    Fitneva_Christiansen_2015; % use 300x2 conditions numSubjects 
elseif (taskvar==6)
    %% Task - Kachergis_Yu_Shiffrin, Psychon Bull Rev (2012) - Adults- Prior Learning, ME
    simName = [simNamePar '_Kachergis_Yu_Shiffrin_2012_'];
    Kachergis_Yu_Shiffrin_2012_Modified; % use 300x12 conditions numSubjects 
elseif (taskvar==7)
    %% Task - Yurovsky_Yu_Smith, Cognitive Science (2013) - Adults - Competitive Processes
    simName = [simNamePar '_Yurovsky_Yu_Smith_2013_'];
    Yurovsky_Yu_Smith_2013; % use 300x3 conditions numSubjects 
elseif (taskvar==8)
    %% Task - Suanda_Mugwanya_Namy, Jrl Exp Child Psy (2014) - 6 year olds - Referential Ambiguity
    simName = [simNamePar '_Suanda_Mugwanya_Namy_2014_'];
    Suanda_Mugwanya_Namy_2014; % use 300x3 conditions numSubjects 
elseif (taskvar==902)
    %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 2x2
     simName = [simNamePar '_Yu_Smith_Two_2007_'];
     Yu_Smith_Two_2007;
elseif (taskvar==903)
    %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 3x3     
    simName = [simNamePar '_Yu_Smith_Three_2007_'];
    Yu_Smith_Three_2007;
 elseif (taskvar==904)
    %% Task - Yu & Smith, Psychological Science (2007) - Adults - Referrent frequency 4x4           
    simName = [simNamePar '_Yu_Smith_Four_2007_'];
    Yu_Smith_Four_2007;
elseif (taskvar==11)
    %% Yu Zhong & Fricker, Frontiers in Psychology(2012) - Adults - Eyetracking
    simName = [simNamePar '_Yu_Zhong_Fricker_2012_'];
    Yu_Zhong_Fricker_2012;    
elseif (taskvar==12)
    %% Trueswell Medina Hafri & Gleitman Cognitive Psychology(2013) - Adults - Eyetracking
    simName = [simNamePar '_Trueswell_Medina_Hafri_Gleitman_2013_'];
    Trueswell_Medina_Hafri_Gleitman_2013; 
end
     
%% Save Simulated Results 
train=[];test=[];
for subject=1:numSubjects
   try
       OutName1 = [simName,num2str(subject),'_train.mat'];              
       OutName2 = [simName,num2str(subject),'_test.mat'];
       tempTrn=load(OutName1);
       tempTst=load(OutName2);
       [tempTrn(:).subject] = subject;
       train = [train; tempTrn];
       [tempTst(:).subject] = subject;
       test = [test; tempTst];
       delete(OutName1);
       delete(OutName2);
    catch
       disp('Error on concatenating subject number ');
       disp(subject);
       continue;
   end
end
OutName = [simName,'results.mat'];
save(OutName,'train','test','sim','notes');
 % delete parallel pool
%delete(poolobj)