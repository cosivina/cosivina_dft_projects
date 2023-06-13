spt_L = 15; %spatial locations for stimuli
spt_Lm = 35;
spt_C = 50;
spt_Rm = 65;
spt_R = 85;
% sim.addElement(WeightMatrix('atn_sr -> lookL', gauss(1:fieldSize_spt, spt_L, 10)'), 'atn_sr', [], 'lookL');
% sim.addElement(WeightMatrix('atn_sr -> lookLWB', gauss(1:fieldSize_spt, spt_Lm, 10)'), 'atn_sr', [], 'lookLWB');
% sim.addElement(WeightMatrix('atn_sr -> lookC', gauss(1:fieldSize_spt, spt_C, 4)'), 'atn_sr', [], 'lookC');
% sim.addElement(WeightMatrix('atn_sr -> lookR', gauss(1:fieldSize_spt, spt_R, 10)'), 'atn_sr', [], 'lookR');
% sim.addElement(WeightMatrix('atn_sr -> lookRWB', gauss(1:fieldSize_spt, spt_Rm, 10)'), 'atn_sr', [], 'lookRWB');

handle_WM_L=sim.getElement('atn_sr -> lookL');
handle_WM_LWB=sim.getElement('atn_sr -> lookLWB');
handle_WM_C=sim.getElement('atn_sr -> lookC');
handle_WM_R=sim.getElement('atn_sr -> lookR');
handle_WM_RWB=sim.getElement('atn_sr -> lookRWB');

handle_WM_L.weights    = gauss(1:fieldSize_spt, spt_L, 10)'; 
handle_WM_LWB.weights  = gauss(1:fieldSize_spt, spt_Lm, 10)'; 
handle_WM_C.weights    = gauss(1:fieldSize_spt, spt_C, 4)'; 
handle_WM_R.weights    = gauss(1:fieldSize_spt, spt_R, 10)';
handle_WM_RWB.weights  = gauss(1:fieldSize_spt, spt_Rm, 10)'; 
