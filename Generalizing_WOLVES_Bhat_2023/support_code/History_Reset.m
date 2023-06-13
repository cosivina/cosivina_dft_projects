%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

%% History_Variable  / Look_Duration / Task_Specific_Simulation_Time Resetting 
%Note: should probably add 'history motor' and 'history inpn' for completeness
handle_historyL=sim.getElement('history lookL'); %reset lookingduration history variable                 
handle_historyL.timeSlots = t_max; 
handle_historyR=sim.getElement('history lookR');                       
handle_historyR.timeSlots = t_max; 
handle_historyLWB=sim.getElement('history lookLWB'); %reset lookingduration history variable                 
handle_historyLWB.timeSlots = t_max; 
handle_historyRWB=sim.getElement('history lookRWB');                       
handle_historyRWB.timeSlots = t_max; 
handle_historyC=sim.getElement('history lookC');                       
handle_historyC.timeSlots = t_max; 

handle_historyLt=sim.getElement('history lookLt');                       
handle_historyLt.timeSlots = t_maxt; 
handle_historyRt=sim.getElement('history lookRt');                       
handle_historyRt.timeSlots = t_maxt;
handle_historyLWBt=sim.getElement('history lookLWBt'); %reset lookingduration history variable                 
handle_historyLWBt.timeSlots = t_maxt; 
handle_historyRWBt=sim.getElement('history lookRWBt');                       
handle_historyRWBt.timeSlots = t_maxt; 
handle_historyCt=sim.getElement('history lookCt');                       
handle_historyCt.timeSlots = t_maxt;