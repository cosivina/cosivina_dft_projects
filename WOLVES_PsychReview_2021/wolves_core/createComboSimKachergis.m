%% setting up the architecture (fields, interactions, and inputs)

scale_factor = 8; % Each simulation timestep equals 8 real-time milliseconds 
historyDuration = floor((8000+1000)/scale_factor);%Simulation time training; + a gap of 1 sec between every two trials%
historyDuration2 = floor((1000+1000)/scale_factor); %Simulation time test trials + a gap of 1 sec between every two trials

% create simulator object
rng('shuffle'); %a new seed for the random generator
sim = Simulator();

% parameters shared by multiple fields
fieldSize_spt = 300;
fieldSize_ftr = 306;
fieldSize_wd = 20;

sigma_exc = 5;
sigma_inh = 10;
vstrength =  5.25; % visual stimulus strength
wstrength = 8;% word stimulus strength

spt_L = 10; %spatial locations for stimuli
spt_Lm = 30;
spt_C = 50;
spt_Rm = 70;
spt_R = 90;

% spt_L = 15; %spatial locations for stimuli
% spt_Lm = 35;
% spt_C = 50;
% spt_Rm = 65;
% spt_R = 85;

nFeatures = 2;

kcf = 3; % kernel cutoff factor
tau = 5;
tauMotor = 5;
tauBuild = 1000;
tauDecay = 15000;

% add spatial fields
sim.addElement(NeuralField('atn_sr', fieldSize_spt, tau));
sim.addElement(NeuralField('ior_s', fieldSize_spt, tau));
sim.addElement(NeuralField('cos', 1, tau));
sim.addElement(NeuralField('cos_m', 1, tau));
sim.addElement(NeuralField('mrn', 1, tauMotor));
sim.addElement(NeuralField('inpn', 1, tau));
sim.addElement(NeuralField('lookL', 1, tau, -5));
sim.addElement(NeuralField('lookLWB', 1, tau, -5));
sim.addElement(NeuralField('lookC', 1, tau, -5));
sim.addElement(NeuralField('lookR', 1, tau, -5));
sim.addElement(NeuralField('lookRWB', 1, tau, -5));
sim.addElement(NeuralField('atn_sa', fieldSize_spt, tau));
sim.addElement(NeuralField('con_s', fieldSize_spt, tau));
sim.addElement(MemoryTrace('hcon_s', fieldSize_spt, tauBuild, tauDecay, 0), 'con_s','output');%new
sim.addElement(NeuralField('wm_s', fieldSize_spt, tau));
sim.addElement(MemoryTrace('hwm_s', fieldSize_spt, tauBuild, tauDecay, 0), 'wm_s','output');%new
sim.addElement(NeuralField('word', fieldSize_wd, tau));
sim.addElement(MemoryTrace('hword', fieldSize_wd, tauBuild, tauDecay, 0), 'word','output');
%sim.addElement(NeuralField('mwf', fieldSize_wd, tau));
%sim.addElement(NeuralField('wwf', fieldSize_wd, tau));
sim.addElement(NeuralField('motor', fieldSize_spt, tau));
%sim.addElement(NeuralField('motorL', fieldSize_spt, tau));
%sim.addElement(NeuralField('motorR', fieldSize_spt, tau));
%sim.addElement(NeuralField('ww', [fieldSize_wd, fieldSize_wd], tau));
%sim.addElement(MemoryTrace('hww', [fieldSize_wd, fieldSize_wd], tauBuild, tauDecay, 0),'ww','output');
for ikch=1:11
    n = num2str(ikch);
    sim.addElement(NeuralField(['look_K' n], 1, tau, -5));
    sim.addElement(ScaleInput(['look_K' n ' -> look_K' n], 1, 3), ['look_K' n], [], ['look_K' n]);
    sim.addElement(WeightMatrix(['atn_sr -> look_K' n], gauss(1:fieldSize_spt, ikch*25, 10)'), 'atn_sr', [], ['look_K' n]);
    sim.addElement(RunningHistory(['history look_K' n], [1, 1], historyDuration2, 1), ['look_K' n], 'output');
end

% add space-feature and feature fields
for i = 1 : nFeatures
    n = num2str(i);
    sim.addElement(NeuralField(['vis_f' n], [fieldSize_ftr, fieldSize_spt], tau));
    sim.addElement(NeuralField(['atn_f' n], fieldSize_ftr, tau));
    sim.addElement(NeuralField(['con_f' n], fieldSize_ftr, tau));
    sim.addElement(MemoryTrace(['hcon_f' n], fieldSize_ftr, tauBuild, tauDecay, 0), ['con_f' n],'output');%%new
    sim.addElement(NeuralField(['wm_f' n], fieldSize_ftr, tau));
    sim.addElement(MemoryTrace(['hwm_f' n], fieldSize_ftr, tauBuild, tauDecay, 0), ['wm_f' n],'output');%%newNo
    sim.addElement(NeuralField(['atn_c' n], [fieldSize_ftr, fieldSize_spt], tau));
    sim.addElement(NeuralField(['wm_c' n], [fieldSize_ftr, fieldSize_spt], tau));
    sim.addElement(MemoryTrace(['hwm_c' n], [fieldSize_ftr, fieldSize_spt], tauBuild, tauDecay, 0), ['wm_c' n],'output');
    sim.addElement(NeuralField(['pd_c' n], 1, tau));
    sim.addElement(NeuralField(['wf' n], [fieldSize_ftr, fieldSize_wd], tau));
    sim.addElement(MemoryTrace(['hwf' n], [fieldSize_ftr, fieldSize_wd], tauBuild, tauDecay, 0), ['wf' n],'output');
end

% add lateral connections for purely spatial fields
sim.addElement(LateralInteractions1D('atn_sr -> atn_sr', fieldSize_spt, sigma_exc, 0, sigma_inh, 0, 0, true, true, kcf), ...
    'atn_sr', [], 'atn_sr');
sim.addElement(LateralInteractions1D('ior_s -> ior_s', fieldSize_spt, sigma_exc, 0, sigma_inh, 0, 0, true, true, kcf), ...
    'ior_s', [], 'ior_s');
sim.addElement(ScaleInput('cos -> cos', 1, 0), 'cos', [], 'cos');
sim.addElement(ScaleInput('cos_m -> cos_m', 1, 0), 'cos_m', [], 'cos_m');
sim.addElement(ScaleInput('mrn -> mrn', 1, 0), 'mrn', [], 'mrn');
sim.addElement(LateralInteractions1D('atn_sa -> atn_sa', fieldSize_spt, sigma_exc, 0, sigma_inh, 0, 0, true, true, kcf), ...
    'atn_sa', [], 'atn_sa');
sim.addElement(LateralInteractions1D('con_s -> con_s', fieldSize_spt, sigma_exc, 0, sigma_inh, 0, 0, true, true, kcf), ...
    'con_s', [], 'con_s');
sim.addElement(LateralInteractions1D('wm_s -> wm_s', fieldSize_spt, sigma_exc, 0, sigma_inh, 0, 0, true, true, kcf), ...
    'wm_s', [], 'wm_s');

sim.addElement(ScaleInput('lookL -> lookL', 1, 3), 'lookL', [], 'lookL');
sim.addElement(ScaleInput('lookLWB -> lookLWB', 1, 3), 'lookLWB', [], 'lookLWB');
sim.addElement(ScaleInput('lookC -> lookC', 1, 3), 'lookC', [], 'lookC');
sim.addElement(ScaleInput('lookR -> lookR', 1, 3), 'lookR', [], 'lookR');
sim.addElement(ScaleInput('lookRWB -> lookRWB', 1, 3), 'lookRWB', [], 'lookRWB');
%sim.addElement(WeightMatrix('atn_sr -> motorL', gauss(1:fieldSize_spt, 25, 30)), 'atn_sr', [], 'motorL');
%sim.addElement(WeightMatrix('atn_sr -> motorR', gauss(1:fieldSize_spt, 75, 30)), 'atn_sr', [], 'motorR');
sim.addElement(WeightMatrix('atn_sr -> lookL', gauss(1:fieldSize_spt, spt_L, 10)'), 'atn_sr', [], 'lookL');
sim.addElement(WeightMatrix('atn_sr -> lookLWB', gauss(1:fieldSize_spt, spt_Lm, 10)'), 'atn_sr', [], 'lookLWB');
sim.addElement(WeightMatrix('atn_sr -> lookC', gauss(1:fieldSize_spt, spt_C, 4)'), 'atn_sr', [], 'lookC');
sim.addElement(WeightMatrix('atn_sr -> lookR', gauss(1:fieldSize_spt, spt_R, 10)'), 'atn_sr', [], 'lookR');
sim.addElement(WeightMatrix('atn_sr -> lookRWB', gauss(1:fieldSize_spt, spt_Rm, 10)'), 'atn_sr', [], 'lookRWB');
%sim.addElement(SumDimension('motorR -> lookR', 2, [1, 1], 0.75), 'motorR', [], 'lookR');
%sim.addElement(SumDimension('motorL -> lookL', 2, [1, 1], 0.75), 'motorL', [], 'lookL');

% lateral connections for DCCS fields
sim.addElement(LateralInteractions1D('word -> word', fieldSize_wd, 0, 0, sigma_inh, 0, 0, true, true, kcf), ...
    'word', [], 'word');
%sim.addElement(LateralInteractions1D('mwf -> mwf', fieldSize_wd, 0, 0, sigma_inh, 0, 0, true, true, kcf), ...
%    'mwf', [], 'mwf');
%sim.addElement(LateralInteractions1D('wwf -> wwf', fieldSize_wd, 0, 0, sigma_inh, 0, 0, true, true, kcf), ...
%    'wwf', [], 'wwf');
sim.addElement(LateralInteractions1D('motor -> motor', fieldSize_spt, sigma_exc, 0, sigma_inh, 0, 0, true, true, kcf), ...
    'motor', [], 'motor');
%sim.addElement(LateralInteractions2D('ww -> ww', [fieldSize_wd, fieldSize_wd], 0, 0, 1.5, 8, 20, 1, -0.03), ...
%    'ww', [], 'ww');
%sim.addElement(SumAllDimensions(['sum ww'], [fieldSize_wd, fieldSize_wd]), ['ww']);

% lateral connections for feature fields
for i = 1 : nFeatures
    n = num2str(i);
    
    % lateral connections and output sums for visual field
    sim.addElement(GaussKernel2D(['vis_f' n ' -> vis_f' n ' (exc.)'], [fieldSize_ftr, fieldSize_spt], ...
        sigma_exc, sigma_exc, 0, true, true, true, kcf), ['vis_f' n], [], ['vis_f' n]);
    sim.addElement(GaussKernel2D(['vis_f' n ' -> vis_f' n ' (inh.)'], [fieldSize_ftr, fieldSize_spt], ...
        sigma_inh, sigma_inh, 0, true, true, true, kcf), ['vis_f' n], [], ['vis_f' n]);
    sim.addElement(SumAllDimensions(['sum vis_f' n], [fieldSize_ftr, fieldSize_spt]), ['vis_f', n]);
    sim.addElement(ScaleInput(['vis_f' n ' -> vis_f' n ' (global)'], 1, 0), ['sum vis_f' n], 'fullSum', ['vis_f' n]);
    
    % lateral connections and output sums for 2d selection field
    sim.addElement(GaussKernel2D(['atn_c' n ' -> atn_c' n ' (exc.)'], [fieldSize_ftr, fieldSize_spt], ...
        sigma_exc, sigma_exc, 0, true, true, true, kcf), ['atn_c' n], [], ['atn_c' n]);
    %   sim.addElement(GaussKernel2D(['atn_c' n ' -> atn_c' n ' (inh.)'], [fieldSize_ftr, fieldSize_spt], ...
    %     sigma_inh, sigma_inh, 0, true, true, true, kcf), ['atn_c' n], [], ['atn_c' n]);
    sim.addElement(SumAllDimensions(['sum atn_c' n], [fieldSize_ftr, fieldSize_spt]), ['atn_c', n]);
    sim.addElement(ScaleInput(['atn_c' n ' -> atn_c' n ' (global)'], 1, 0), ['sum atn_c' n], 'fullSum', ['atn_c' n]);
    
    % lateral connections and output sums for association memory fields
    sim.addElement(GaussKernel2D(['wm_c' n ' -> wm_c' n ' (exc.)'], [fieldSize_ftr, fieldSize_spt], ...
        sigma_exc, sigma_exc, 0, true, true, true, kcf), ['wm_c' n], [], ['wm_c' n]);
    sim.addElement(GaussKernel2D(['wm_c' n ' -> wm_c' n ' (inh.)'], [fieldSize_ftr, fieldSize_spt], ...
        sigma_inh, sigma_inh, 0, true, true, true, kcf), ['wm_c' n], [], ['wm_c' n]);
    sim.addElement(SumAllDimensions(['sum wm_c' n], [fieldSize_ftr, fieldSize_spt]), ['wm_c', n]);
    sim.addElement(ScaleInput(['wm_c' n ' -> wm_c' n ' (global/feature)'], [1, fieldSize_spt], 0), ...
        ['sum wm_c' n], 'verticalSum');
    sim.addElement(ExpandDimension2D(['expand wm_c' n ' -> wm_c' n ' (global/feature)'], ...
        1, [fieldSize_ftr, fieldSize_spt]), ['wm_c' n ' -> wm_c' n ' (global/feature)'], [], ['wm_c' n]);
    sim.addElement(ScaleInput(['wm_c' n ' -> wm_c' n ' (global)'], 1, 0), ['sum wm_c' n], 'fullSum', ['wm_c' n]);
    
    % lateral connections for 1D fields
    sim.addElement(LateralInteractions1D(['atn_f' n ' -> atn_f' n], fieldSize_ftr, sigma_exc, 0, sigma_inh, 0, 0, ...
        true, true, kcf), ['atn_f' n], [], ['atn_f' n]);
    sim.addElement(LateralInteractions1D(['con_f' n ' -> con_f' n], fieldSize_ftr, sigma_exc, 0, sigma_inh, 0, 0, ...
        true, true, kcf), ['con_f' n], [], ['con_f' n]);
    sim.addElement(LateralInteractions1D(['wm_f' n ' -> wm_f' n], fieldSize_ftr, sigma_exc, 0, sigma_inh, 0, 0, ...
        true, true, kcf), ['wm_f' n], [], ['wm_f' n]);
    
    % lateral connections and such for 2D word fields
    sim.addElement(LateralInteractions2D(['wf' n ' -> wf' n], [fieldSize_ftr, fieldSize_wd], 5, 0, 1.5, 8, 20, 1, -0.03), ...
        ['wf' n], [], ['wf' n]);
    sim.addElement(SumAllDimensions(['sum wf' n], [fieldSize_ftr, fieldSize_wd]), ['wf', n]);
    
    % self-excitation for peak detector nodes
    sim.addElement(ScaleInput(['pd_c' n ' -> pd_c' n], 1, 0), ['pd_c' n], [], ['pd_c' n]);
end


% add connections between purely spatial fields
% from spatial attention (retinal) field
sim.addElement(GaussKernel1D('atn_sr -> atn_sa', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sr');
sim.addElement(ScaleInput('scale atn_sr -> atn_sa', fieldSize_spt, 1), 'atn_sr -> atn_sa', [], 'atn_sa');
sim.addElement(GaussKernel1D('atn_sr -> vis_f', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sr');
sim.addElement(ExpandDimension2D('expand atn_sr -> vis_f', 1, [fieldSize_ftr, fieldSize_spt]), ...
    'atn_sr -> vis_f', [], cellstr([repmat('vis_f', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(GaussKernel1D('atn_sr -> ior_s', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sr', [], 'ior_s');

% from inhibition of return field
sim.addElement(GaussKernel1D('ior_s -> atn_sr', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'ior_s', [], 'atn_sr');
sim.addElement(GaussKernel1D('ior_s -> atn_sa', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'ior_s');
sim.addElement(ScaleInput('scale ior_s -> atn_sa', fieldSize_spt, 1), 'ior_s -> atn_sa', [], 'atn_sa');
sim.addElement(SumDimension('motor -> mrn', 2, [1,1], 0), 'motor', [], 'mrn');
sim.addElement(SumDimension('motor -> atn_f', 2, [1,1], 0), 'motor', [], ...
    cellstr([repmat('atn_f', [nFeatures, 1]), num2str((1:nFeatures)')]));

% from condition of satisfaction node
sim.addElement(ScaleInput('cos -> ior_s', 1, 0), 'cos', [], 'ior_s');
sim.addElement(ScaleInput('cos -> atn_c', 1, 0), 'cos', [], ...
    cellstr([repmat('atn_c', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(ScaleInput('cos -> wf', 1, 0), 'cos', [], ...
    cellstr([repmat('wf', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(ScaleInput('cos -> atn_sr', 1, 0), 'cos', [], 'atn_sa');
sim.addElement(ScaleInput('cos -> atn_sa', 1, 0), 'cos', [], 'atn_sa');
sim.addElement(ScaleInput('cos -> atn_f', 1, 0), 'cos', [], ...
    cellstr([repmat('atn_f', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(ScaleInput('cos -> cos_m', 1, 0), 'cos', [], 'cos_m');

% from motor condition of satisfaction node
sim.addElement(ScaleInput('mrn -> cos_m', 1, 0), 'mrn', [], 'cos_m');
sim.addElement(ScaleInput('cos_m -> motor', 1, 0), 'cos_m', [], 'motor');
sim.addElement(ScaleInput('cos_m -> mrn', 1, 0), 'cos_m', [], 'mrn');

% from motor response node
sim.addElement(ScaleInput('mrn -> motor', 1, 0), 'mrn', [], 'motor');

% from spatial attention (allocentric) field
sim.addElement(GaussKernel1D('atn_sa -> atn_sr', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sa');
sim.addElement(ScaleInput('scale atn_sa -> atn_sr', fieldSize_spt, 1), 'atn_sa -> atn_sr', [], 'atn_sr');
sim.addElement(SumDimension('sum atn_sa', 2, [1,1], 1), 'atn_sa');
sim.addElement(ScaleInput('atn_sa -> atn_sr (inh.)', 1, 0), 'sum atn_sa', [], 'atn_sr');
sim.addElement(GaussKernel1D('atn_sa -> con_s', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sa', [], 'con_s');
sim.addElement(GaussKernel1D('atn_sa -> wm_s', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sa', [], 'wm_s');
sim.addElement(GaussKernel1D('atn_sa -> motor', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sa', [], 'motor');
sim.addElement(GaussKernel1D('motor -> atn_sa', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'motor', [], 'atn_sa');

sim.addElement(GaussKernel1D('atn_sa -> wm_c', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sa');
sim.addElement(ExpandDimension2D('expand atn_sa -> wm_c', 1, [fieldSize_ftr, fieldSize_spt]), ...
    'atn_sa -> wm_c', [], cellstr([repmat('wm_c', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(GaussKernel1D('atn_sa -> atn_c', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'atn_sa');
sim.addElement(ExpandDimension2D('expand atn_sa -> atn_c', 1, [fieldSize_ftr, fieldSize_spt]), ...
    'atn_sa -> atn_c', [], cellstr([repmat('atn_c', [nFeatures, 1]), num2str((1:nFeatures)')]));

% from contrast (spatial) field
sim.addElement(GaussKernel1D('con_s -> atn_sa', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'con_s', [], 'atn_sa');
sim.addElement(GaussKernel1D('con_s -> wm_s', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'con_s', [], 'wm_s');
sim.addElement(GaussKernel1D('con_s -> atn_c', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'con_s');
sim.addElement(ExpandDimension2D('expand con_s -> atn_c', 1, [fieldSize_ftr, fieldSize_spt]), ...
    'con_s -> atn_c', [], cellstr([repmat('atn_c', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(GaussKernel1D('hcon_s -> con_s', fieldSize_spt, sigma_exc, 0.2, true, true, kcf), 'hcon_s', [], 'con_s');%new

% from working memory (spatial)
sim.addElement(GaussKernel1D('wm_s -> con_s', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'wm_s', [], 'con_s');
sim.addElement(GaussKernel1D('wm_s -> atn_sa', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'wm_s', [], 'atn_sa');
sim.addElement(GaussKernel1D('wm_s -> wm_c', fieldSize_spt, sigma_exc, 0, true, true, kcf), 'wm_s');
sim.addElement(ExpandDimension2D('expand wm_s -> wm_c', 1, [fieldSize_ftr, fieldSize_spt]), ...
    'wm_s -> wm_c', [], cellstr([repmat('wm_c', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(GaussKernel1D('hwm_s -> wm_s', fieldSize_spt, sigma_exc, 0.2, true, true, kcf), 'hwm_s', [], 'wm_s');%new

% from 1D word field
sim.addElement(GaussKernel1D('word -> wf', fieldSize_wd, 0, 0, true, true, kcf), 'word');
sim.addElement(ExpandDimension2D('expand word -> wf', 1, [fieldSize_ftr, fieldSize_wd]), ...
    'word -> wf', [], cellstr([repmat('wf', [nFeatures, 1]), num2str((1:nFeatures)')]));

sim.addElement(GaussKernel1D('hword -> word', fieldSize_wd, 0, 0.2, true, true, kcf), 'hword', [], 'word');

%sim.addElement(GaussKernel1D('word -> ww', fieldSize_wd, 0, 0, true, true, kcf), 'word');
%sim.addElement(ExpandDimension2D('expand word -> ww', 1, [fieldSize_wd, fieldSize_wd]), ...
%    'word -> ww', [], 'ww');
%sim.addElement(GaussKernel1D('word -> mwf', fieldSize_wd, 0, 0, true, true, kcf), 'word');
%sim.addElement(ScaleInput('scale word -> mwf', fieldSize_wd, 1), 'word -> mwf', [], 'mwf');
%sim.addElement(GaussKernel1D('word -> wwf', fieldSize_wd, 0, 0, true, true, kcf), 'word');
%sim.addElement(ScaleInput('scale word -> wwf', fieldSize_wd, 1), 'word -> wwf', [], 'wwf');
sim.addElement(SumDimension('word -> mrn', 2, [1,1], 1), 'word', [], 'mrn');

% from mwf field
%sim.addElement(GaussKernel1D('mwf -> word', fieldSize_wd, 0, 0, true, true, kcf), 'mwf');
%sim.addElement(ScaleInput('scale mwf -> word', fieldSize_wd, 1), 'mwf -> word', [], 'word');
%sim.addElement(GaussKernel1D('mwf -> wwf', fieldSize_wd, 0, 0, true, true, kcf), 'mwf');
%sim.addElement(ScaleInput('scale mwf -> wwf', fieldSize_wd, 1), 'mwf -> wwf', [], 'wwf');
%sim.addElement(GaussKernel1D('wwf -> ww', fieldSize_wd, 0, 0, true, true, kcf), 'wwf');
%sim.addElement(ScaleInput('scale wwf -> ww', fieldSize_wd, 1), 'wwf -> ww');
%sim.addElement(ExpandDimension2D('expand wwf -> ww', 2, [fieldSize_wd, fieldSize_wd]), ...
%    'scale wwf -> ww', [], 'ww');

%sim.addElement(GaussKernel2D('ww -> hww', [fieldSize_wd, fieldSize_wd], 0, 0, 1, ...
%    true, true, true, kcf), 'ww', 'output', 'hww');
%sim.addElement(GaussKernel2D('hww -> ww', [fieldSize_wd, fieldSize_wd], 1, 0, 0.8), 'hww', [], 'ww');
%sim.addElement(ScaleInput('scale hww', [fieldSize_wd, fieldSize_wd], 100), 'hww');



% from motor field, silly way to sum for activity on a 1D field
sim.addElement(SumDimension('motor -> cos_m', 2, [1,1], 0), 'motor', [], 'cos_m');
sim.addElement(RunningHistory('history motor', [1, 1], historyDuration, 1), 'motor -> cos_m', 'output');
sim.addElement(RunningHistory('history inpn', [1, 1], historyDuration, 1), 'inpn', 'h');

sim.addElement(RunningHistory('history lookL', [1, 1], historyDuration, 1), 'lookL', 'output');
sim.addElement(RunningHistory('history lookLWB', [1, 1], historyDuration, 1), 'lookLWB', 'output');
sim.addElement(RunningHistory('history lookC', [1, 1], historyDuration, 1), 'lookC', 'output');
sim.addElement(RunningHistory('history lookR', [1, 1], historyDuration, 1), 'lookR', 'output');
sim.addElement(RunningHistory('history lookRWB', [1, 1], historyDuration, 1), 'lookRWB', 'output');
sim.addElement(RunningHistory('history lookLt', [1, 1], historyDuration2, 1), 'lookL', 'output');
sim.addElement(RunningHistory('history lookLWBt', [1, 1], historyDuration2, 1), 'lookLWB', 'output');
sim.addElement(RunningHistory('history lookCt', [1, 1], historyDuration2, 1), 'lookC', 'output');
sim.addElement(RunningHistory('history lookRt', [1, 1], historyDuration2, 1), 'lookR', 'output');
sim.addElement(RunningHistory('history lookRWBt', [1, 1], historyDuration2, 1), 'lookRWB', 'output');

% known words - trace projected across corresponding word fields
%sig_exc=0.5*sigma_exc;
%sim.addElement(GaussStimulus1D('colorTrace', fieldSize_wd, sig_exc, 0, 25, true, false));
%sim.addElement(ExpandDimension2D('expand colorTrace -> wf1', 1, [fieldSize_ftr, fieldSize_wd]), ...
%    'colorTrace', [], 'wf1');
%sim.addElement(GaussStimulus1D('shapeTrace', fieldSize_wd, sig_exc, 0, 75, true, false));
%sim.addElement(ExpandDimension2D('expand shapeTrace -> wf2', 1, [fieldSize_ftr, fieldSize_wd]), ...
%    'shapeTrace', [], 'wf2');

for i = 1 : nFeatures
    n = num2str(i);
    
    % from visual field
    sim.addElement(GaussKernel1D(['vis_f' n ' -> atn_sr'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum vis_f' n], 'verticalSum', 'atn_sr');
    sim.addElement(GaussKernel1D(['vis_f' n ' -> ior_s'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum vis_f' n], 'verticalSum', 'ior_s');
    sim.addElement(GaussKernel1D(['vis_f' n ' -> con_s'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum vis_f' n], 'verticalSum', 'con_s');
    sim.addElement(GaussKernel1D(['vis_f' n ' -> wm_s'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum vis_f' n], 'verticalSum', 'wm_s');
    sim.addElement(GaussKernel1D(['vis_f' n ' -> atn_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum vis_f' n], 'horizontalSum', ['atn_f' n]);
    sim.addElement(GaussKernel1D(['vis_f' n ' -> con_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum vis_f' n], 'horizontalSum', ['con_f' n]);
    sim.addElement(GaussKernel1D(['vis_f' n ' -> wm_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum vis_f' n], 'horizontalSum', ['wm_f' n]);
    
    % from feature attention field
    sim.addElement(GaussKernel1D(['atn_f' n ' -> vis_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ['atn_f' n]);
    sim.addElement(ExpandDimension2D(['expand atn_f' n ' -> vis_f' n], 2, [fieldSize_ftr, fieldSize_spt]), ...
        ['atn_f' n ' -> vis_f' n], [], ['vis_f' n]);
    
    sim.addElement(GaussKernel1D(['atn_f' n ' -> con_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['atn_f' n], [], ['con_f' n]);
    sim.addElement(GaussKernel1D(['atn_f' n ' -> wm_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['atn_f' n], [], ['wm_f' n]);
    
    sim.addElement(GaussKernel1D(['atn_f' n ' -> atn_c' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ['atn_f' n]);
    sim.addElement(ExpandDimension2D(['expand atn_f' n ' -> atn_c' n], 2, [fieldSize_ftr, fieldSize_spt]), ...
        ['atn_f' n ' -> atn_c' n], [], ['atn_c' n]);
    
    sim.addElement(GaussKernel1D(['atn_f' n ' -> wf' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ['atn_f' n]);
    sim.addElement(ExpandDimension2D(['expand atn_f' n ' -> wf' n], 2, [fieldSize_ftr, fieldSize_wd]), ...
        ['atn_f' n ' -> wf' n], [], ['wf' n]);
    
    sim.addElement(GaussKernel1D(['atn_f' n ' -> wm_c' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ['atn_f' n]);
    sim.addElement(ExpandDimension2D(['expand atn_f' n ' -> wm_c' n], 2, [fieldSize_ftr, fieldSize_spt]), ...
        ['atn_f' n ' -> wm_c' n], [], ['wm_c' n]);
    
    % from perceptual intention (feature) field
    sim.addElement(GaussKernel1D(['con_f' n ' -> atn_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['con_f' n], [], ['atn_f' n]);
    sim.addElement(GaussKernel1D(['con_f' n ' -> wm_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['con_f' n], [], ['wm_f' n]);
    
    sim.addElement(GaussKernel1D(['con_f' n ' -> atn_c' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ['con_f' n]);
    sim.addElement(ExpandDimension2D(['expand con_f' n ' -> atn_c' n], 2, [fieldSize_ftr, fieldSize_spt]), ...
        ['con_f' n ' -> atn_c' n], [], ['atn_c' n]);
    
    sim.addElement(GaussKernel1D(['hcon_f' n ' -> con_f' n], fieldSize_ftr, sigma_exc, 0.2, true, true, kcf), ['hcon_f' n], [], ['con_f' n]);

    % from working memory memory (feature) field
    sim.addElement(GaussKernel1D(['wm_f' n ' -> con_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['wm_f' n], [], ['con_f' n]);
    sim.addElement(GaussKernel1D(['wm_f' n ' -> atn_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['wm_f' n], [], ['atn_f' n]);
    
    sim.addElement(GaussKernel1D(['wm_f' n ' -> wm_c' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ['wm_f' n]);
    sim.addElement(ExpandDimension2D(['expand wm_f' n ' -> wm_c' n], 2, [fieldSize_ftr, fieldSize_spt]), ...
        ['wm_f' n ' -> wm_c' n], [], ['wm_c' n]);

    sim.addElement(GaussKernel1D(['hwm_f' n ' -> wm_f' n], fieldSize_ftr, sigma_exc, 0.2, true, true, kcf), ['hwm_f' n], [], ['wm_f' n]);
    
    % from feature selection field
    sim.addElement(GaussKernel1D(['atn_c' n ' -> atn_sa'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum atn_c' n], 'verticalSum', 'atn_sa');
    sim.addElement(GaussKernel1D(['atn_c' n ' -> con_s'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum atn_c' n], 'verticalSum', 'con_s');
    sim.addElement(GaussKernel1D(['atn_c' n ' -> wm_s'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum atn_c' n], 'verticalSum', 'wm_s');
    
    sim.addElement(GaussKernel1D(['atn_c' n ' -> con_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum atn_c' n], 'horizontalSum', ['con_f' n]);
    sim.addElement(GaussKernel1D(['atn_c' n ' -> wm_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum atn_c' n], 'horizontalSum', ['wm_f' n]);
    
    sim.addElement(GaussKernel2D(['atn_c' n ' -> wm_c' n], [fieldSize_ftr, fieldSize_spt], sigma_exc, sigma_exc, 0, ...
        true, true, true, kcf), ['atn_c' n], [], ['wm_c' n]);
    sim.addElement(GaussKernel2D(['atn_c' n ' -> wm_c' n ' (inh.)'], [fieldSize_ftr, fieldSize_spt], sigma_exc, sigma_exc, 0, ...
        true, true, true, kcf), ['atn_c' n], [], ['wm_c' n]);
    sim.addElement(GaussKernel1D(['atn_c' n ' -> wf' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum atn_c' n], 'horizontalSum');
    sim.addElement(ExpandDimension2D(['expand atn_c' n ' -> wf' n], 2, [fieldSize_ftr, fieldSize_wd]), ...
        ['atn_c' n ' -> wf' n], [], ['wf' n]);
    %sim.addElement(GaussKernel2D(['atn_c' n ' -> hatn_c' n], [fieldSize_ftr, fieldSize_spt], sigma_exc, sigma_exc, 1, ...
    %    true, true, true, kcf), ['atn_c' n], 'output', ['hatn_c' n]);
    % sim.addElement(GaussKernel2D(['wm_c' n ' -> hatn_c' n], [fieldSize_ftr, fieldSize_spt], sigma_exc, sigma_exc, 1, ...
    %   true, true, true, kcf), ['wm_c' n], 'output', ['hatn_c' n]);
    sim.addElement(GaussKernel2D(['hwm_c' n ' -> atn_c' n], [fieldSize_ftr, fieldSize_spt], sigma_exc, sigma_exc, 20, ...
        true, true, true, kcf), ['hwm_c' n], [], ['atn_c' n]);
    sim.addElement(GaussKernel2D(['hwm_c' n ' -> wm_c' n], [fieldSize_ftr, fieldSize_spt], sigma_exc, sigma_exc, 3, ...
        true, true, true, kcf), ['hwm_c' n], [], ['wm_c' n]);
    %sim.addElement(ScaleInput(['scale hatn_c', n], [fieldSize_ftr, fieldSize_spt], 100), ['hatn_c' n]);
    %sim.addElement(ScaleInput(['hatn_c' n '-> atn_c' n], [fieldSize_ftr, fieldSize_spt], 20), ['hatn_c' n], [], ['atn_c' n]);
    %sim.addElement(GaussKernel2D(['wf' n ' -> hwf' n], [fieldSize_ftr, fieldSize_wd], sigma_exc, sigma_exc, 1, ...
    %  true, true, true, kcf), ['wf' n], 'output', ['hwf' n]);
    sim.addElement(GaussKernel2D(['hwf' n ' -> wf' n], [fieldSize_ftr, fieldSize_wd], 1, 0, 0.8), ['hwf' n], [], ['wf' n]);
    %sim.addElement(ScaleInput(['scale hwf', n], [fieldSize_ftr, fieldSize_wd], 100), ['hwf' n]);
    
    sim.addElement(ScaleInput(['atn_c' n ' -> pd_c' n], 1, 0), ['sum atn_c' n], 'fullSum', ['pd_c' n]);
   
    % from peak detector node
    sim.addElement(ScaleInput(['pd_c' n ' -> cos'], 1, 0), ['pd_c' n], [], 'cos');
    
    % from integrated working memory field
    sim.addElement(GaussKernel1D(['wm_c' n ' -> wm_s'], fieldSize_spt, sigma_exc, 0, true, true, kcf), ...
        ['sum wm_c' n], 'verticalSum', 'wm_s');
    sim.addElement(GaussKernel1D(['wm_c' n ' -> wm_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum wm_c' n], 'horizontalSum', ['wm_f' n]);
    sim.addElement(GaussKernel2D(['wm_c' n ' -> atn_c' n], [fieldSize_ftr, fieldSize_spt], sigma_exc, sigma_exc, 0, ...
        true, true, true, kcf), ['wm_c' n], [], ['atn_c' n]);
    
    % from 2D word fields
    % sum across word dim, project across into atn field
    sim.addElement(GaussKernel1D(['wf' n ' -> atn_c' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum wf' n], 'horizontalSum');
    sim.addElement(ExpandDimension2D(['expand wf' n ' -> atn_c' n], 2, [fieldSize_ftr, fieldSize_spt]), ...
        ['wf' n ' -> atn_c' n], [], ['atn_c' n]);
    % 2D word field project back into word 1D
    sim.addElement(GaussKernel1D(['wf' n ' -> word'], fieldSize_wd, 0, 0, true, true, kcf), ...
        ['sum wf' n], 'verticalSum', 'word');
    sim.addElement(GaussKernel1D(['wf' n ' -> atn_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum wf' n], 'horizontalSum', ['atn_f' n]);
    
    sim.addElement(GaussKernel1D(['wf' n ' -> con_f' n], fieldSize_ftr, sigma_exc, 0, true, true, kcf), ...
        ['sum wf' n], 'horizontalSum', ['con_f' n]);
    
    % 2D word field exciting motor response node
    % (depends on previous horizontal sum)
    %sim.addElement(ScaleInput(['wf' n ' -> motor' n], 1, 0), ['wf' n ' -> atn_c' n], [], 'motor');
    
 
end

%% CREATE STIMULI

% nObjects = 6; %%total number of objects
% %nFeatures = 2; %% Total number of features 
% %%CAUTION..if nFeatures< size of Property{}, re-number up the desired Property{i}
% Property{1} = {15, 30, 45, 60, 75, 90};%%Shapes  %%more than or equal to nObjects
% Property{2} = {11, 22, 33, 44, 55, 66};%%Colors   %%more than or equal to nObjects
% Names = {1, 2, 3, 4, 5, 6};  %%more than or equal to nObjects

%nFeatures = 2; %% Total number of features 
%%CAUTION..if nFeatures< size of Property{}, re-number up the desired Property{i}
Property{1} = {10, 27, 44, 61, 78, 95, 112, 129, 146, 163, 180, 197, 214, 231, 248, 265, 282, 299};%%Shapes  %%more than or equal to nObjects
Property{2} = {10, 27, 44, 61, 78, 95, 112, 129, 146, 163, 180, 197, 214, 231, 248, 265, 282, 299};%%Colors   %%more than or equal to nObjects
%Property{1} = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120};%%Shapes  %%more than or equal to nObjects
%Property{2} = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120};%%Colors   %%more than or equal to nObjects
%Property{1} = {10, 22, 34, 46, 58, 70, 82, 94, 106, 118, 130, 142, 154, 166, 178, 190, 202, 214, 226, 238, 250, 262, 274, 286, 298, 310, 322, 334};%%Shapes  %%more than or equal to nObjects
%Property{2} = {10, 22, 34, 46, 58, 70, 82, 94, 106, 118, 130, 142, 154, 166, 178, 190, 202, 214, 226, 238, 250, 262, 274, 286, 298, 310, 322, 334};

Names = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};  %%more than or equal to nObjects

for i= 1:nFeatures  
    n =num2str(i);
    sim.addElement(GaussStimulus2D(strcat('Feature_',n,'_Left'), [fieldSize_ftr, fieldSize_spt],... 
        sigma_exc, sigma_exc, 0, 0, 0, true, true, false), ...
        [], [], ['vis_f' n]);
    sim.addElement(GaussStimulus2D( strcat('Feature_',n,'_Right'), [fieldSize_ftr, fieldSize_spt],... 
        sigma_exc, sigma_exc, 0, 0, 0, true, true, false), ...
        [], [], ['vis_f' n]);
    sim.addElement(GaussStimulus2D(strcat('Feature_',n,'_LeftMid'), [fieldSize_ftr, fieldSize_spt],... 
        sigma_exc, sigma_exc, 0, 0, 0, true, true, false), ...
        [], [], ['vis_f' n]);
    sim.addElement(GaussStimulus2D( strcat('Feature_',n,'_RightMid'), [fieldSize_ftr, fieldSize_spt],... 
        sigma_exc, sigma_exc, 0, 0, 0, true, true, false), ...
        [], [], ['vis_f' n]);
    sim.addElement(GaussStimulus2D( strcat('Feature_',n,'_Centre'), [fieldSize_ftr, fieldSize_spt],... 
        sigma_exc, sigma_exc, 0, 0, 0, true, true, false), ...
        [], [], ['vis_f' n]);
    
end
sim.addElement(GaussStimulus1D('Word_Label_1', fieldSize_wd, 0, 0, 0, true), [], [], 'word');
sim.addElement(GaussStimulus1D('Word_Label_2', fieldSize_wd, 0, 0, 0, true), [], [], 'word');
sim.addElement(GaussStimulus1D('Word_Label_3', fieldSize_wd, 0, 0, 0, true), [], [], 'word');
sim.addElement(GaussStimulus1D('Word_Label_4', fieldSize_wd, 0, 0, 0, true), [], [], 'word');
for i= 1:nFeatures  
    n =num2str(i);
    for ikch=1:12
        mx =num2str(ikch);
    sim.addElement(GaussStimulus2D(strcat('Feature_',n,'_',mx), [fieldSize_ftr, fieldSize_spt],... 
        sigma_exc, sigma_exc, 0, 0, 0, true, true, false), ...
        [], [], ['vis_f' n]);
    end
end

% for i= 1:size(Feature,2)
%     for j=1:size(Feature{i},2)
%         sim.addElement(GaussStimulus2D(char(Feature{i}(j)), [fieldSize_ftr, fieldSize_spt],... 
%             sigma_exc, sigma_exc, 0, 15, 20, true, true, false), ...
%             [], [], 'vis_f1');
%     end
% end
% for i= 1:size(Names,2)
%     sim.addElement(GaussStimulus1D(char(Names(i)), fieldSize_wd, 0, 0, 90, true), [], [], 'word');
% end
%% create stimuli
% sim.addElement(GaussStimulus2D('input 1Al', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 15, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Am', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 15, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Ar', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 15, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Bl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 30, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Bm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 30, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Br', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 30, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Cl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 45, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Cm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 45, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Cr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 45, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Dl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 60, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Dm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 60, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Dr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 60, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1El', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Em', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Er', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Fl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 90, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Fm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 90, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Fr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 90, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Gl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 37, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Gm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 37, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Gr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 37, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Hl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 82, 20, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Hm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 82, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 1Hr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 82, 80, true, true, false), ...
%     [], [], 'vis_f1');
% 
% sim.addElement(GaussStimulus2D('input 2Al', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 45, 20, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Am', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 45, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Ar', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 45, 80, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Bl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 20, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Bm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Br', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 80, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Cl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 15, 20, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Cm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 15, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Cr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 15, 80, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Dl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 60, 20, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Dm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 60, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Dr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 60, 80, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2El', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 30, 20, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Em', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 30, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Er', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 30, 80, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Fl', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 90, 20, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Fm', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 90, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 2Fr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 90, 80, true, true, false), ...
%     [], [], 'vis_f2');

% sim.addElement(GaussStimulus2D('input 1Err', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 90, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Crr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 75, 90, true, true, false), ...
%     [], [], 'vis_f2');
% 
% %Known objects
% sim.addElement(GaussStimulus2D('input 1Cup', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 7, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Cup', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 7, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 1Ball', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 22, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Ball', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 22, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 1Cat', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 37, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Cat', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 37, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 1Dog', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 52, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Dog', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 52, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 1Dogr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 52, 80, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Dogr', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 52, 80, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 1Bottle', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 67, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Bottle', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 67, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 1Book', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 82, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Book', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 82, 50, true, true, false), ...
%     [], [], 'vis_f2');
% sim.addElement(GaussStimulus2D('input 1Car', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 95, 50, true, true, false), ...
%     [], [], 'vis_f1');
% sim.addElement(GaussStimulus2D('input 2Car', [fieldSize_ftr, fieldSize_spt], ...
%     sigma_exc, sigma_exc, 0, 95, 50, true, true, false), ...
%     [], [], 'vis_f2');
% 
% 
% % stimuli for word field
% sim.addElement(GaussStimulus1D('inputColor', fieldSize_wd, 0, 0, 90, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputShape', fieldSize_wd, 0, 0, 80, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputCup', fieldSize_wd, 0, 0, 5, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputBall', fieldSize_wd, 0, 0, 15, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputCat', fieldSize_wd, 0, 0, 25, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputDog', fieldSize_wd, 0, 0, 35, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputBottle', fieldSize_wd, 0, 0, 45, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputBook', fieldSize_wd, 0, 0, 55, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputCar', fieldSize_wd, 0, 0, 65, true), [], [], 'word');
% 
% sim.addElement(GaussStimulus1D('inputBosa', fieldSize_wd, 0, 0, 10, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputGasser', fieldSize_wd, 0, 0, 20, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputManu', fieldSize_wd, 0, 0, 30, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputColat', fieldSize_wd, 0, 0, 40, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputKaki', fieldSize_wd, 0, 0, 50, true), [], [], 'word');
% sim.addElement(GaussStimulus1D('inputRegli', fieldSize_wd, 0, 0, 60, true), [], [], 'word');

sim.addElement(GaussStimulus1D('inputSpacel', fieldSize_spt, sigma_exc, 0, 20, true), [], [], 'atn_sr');
sim.addElement(GaussStimulus1D('inputSpacem', fieldSize_spt, sigma_exc, 0, 50, true), [], [], 'atn_sr');
sim.addElement(GaussStimulus1D('inputSpacer', fieldSize_spt, sigma_exc, 0, 80, true), [], [], 'atn_sr');

sim.addElement(GaussStimulus1D('inputCONF1B', fieldSize_ftr, sigma_exc, 0, 30, true), [], [], 'con_f1');
sim.addElement(GaussStimulus1D('inputCONF1E', fieldSize_ftr, sigma_exc, 0, 75, true), [], [], 'con_f1');
sim.addElement(GaussStimulus1D('inputCONF2B', fieldSize_ftr, sigma_exc, 0, 75, true), [], [], 'con_f2');
sim.addElement(GaussStimulus1D('inputCONF2E', fieldSize_ftr, sigma_exc, 0, 30, true), [], [], 'con_f2');
sim.addElement(GaussStimulus1D('inputWMF1B', fieldSize_ftr, sigma_exc, 0, 30, true), [], [], 'wm_f1');
sim.addElement(GaussStimulus1D('inputWMF1E', fieldSize_ftr, sigma_exc, 0, 75, true), [], [], 'wm_f1');
sim.addElement(GaussStimulus1D('inputWMF2B', fieldSize_ftr, sigma_exc, 0, 75, true), [], [], 'wm_f2');
sim.addElement(GaussStimulus1D('inputWMF2E', fieldSize_ftr, sigma_exc, 0, 30, true), [], [], 'wm_f2');

% boost stimuli for all fields
sim.addElement(BoostStimulus('boost ior_s', 0), [], [], 'ior_s');
sim.addElement(BoostStimulus('boost atn_sr', 0), [], [], 'atn_sr');
sim.addElement(BoostStimulus('boost atn_sa', 0), [], [], 'atn_sa');
sim.addElement(BoostStimulus('boost con_s', 0), [], [], 'con_s');
sim.addElement(BoostStimulus('boost wm_s', 0), [], [], 'wm_s');

sim.addElement(BoostStimulus('boost vis_f', 0), [], [], cellstr([repmat('vis_f', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(BoostStimulus('boost atn_f', 0), [], [], cellstr([repmat('atn_f', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(BoostStimulus('boost con_f', 0), [], [], cellstr([repmat('con_f', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(BoostStimulus('boost wm_f', 0), [], [], cellstr([repmat('wm_f', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(BoostStimulus('boost atn_c', 0), [], [], cellstr([repmat('atn_c', [nFeatures, 1]), num2str((1:nFeatures)')]));
sim.addElement(BoostStimulus('boost atn_c2', 0), [], [], 'atn_c2');
sim.addElement(BoostStimulus('boost wm_c', 0), [], [], cellstr([repmat('wm_c', [nFeatures, 1]), num2str((1:nFeatures)')]));

% add correlated noise
sim.addElement(NormalNoise('noise atn_sr', fieldSize_spt, 1));
sim.addElement(GaussKernel1D('noise kernel atn_sr', fieldSize_spt, 1, 1, true, true, kcf), ...
    'noise atn_sr', [], 'atn_sr');
sim.addElement(NormalNoise('noise ior_s', fieldSize_spt, 1));
sim.addElement(GaussKernel1D('noise kernel ior_s', fieldSize_spt, 1, 1, true, true, kcf), ...
    'noise ior_s', [], 'ior_s');
sim.addElement(NormalNoise('noise atn_sa', fieldSize_spt, 1));
sim.addElement(GaussKernel1D('noise kernel atn_sa', fieldSize_spt, 1, 1, true, true, kcf), ...
    'noise atn_sa', [], 'atn_sa');

%sim.addElement(NormalNoise('noise ww', [fieldSize_wd, fieldSize_wd], 1), [], [], 'ww');
sim.addElement(NormalNoise('noise word', fieldSize_wd, 1), [], [], 'word');
%sim.addElement(NormalNoise('noise mwf', fieldSize_wd, 1), [], [], 'mwf');
%sim.addElement(NormalNoise('noise wwf', fieldSize_wd, 1), [], [], 'wwf');

sim.addElement(NormalNoise('noise cos', 1, 1), [], [], 'cos');


sim.addElement(NormalNoise('noise con_s', fieldSize_spt, 1));
sim.addElement(GaussKernel1D('noise kernel con_s', fieldSize_spt, 0, 0, true, true, kcf), ...
    'noise con_s', [], 'con_s');
sim.addElement(NormalNoise('noise wm_s', fieldSize_spt, 1));
sim.addElement(GaussKernel1D('noise kernel wm_s', fieldSize_spt, 0, 0, true, true, kcf), ...
    'noise wm_s', [], 'wm_s');

for i = 1 : nFeatures
    n = num2str(i);
    
    sim.addElement(NormalNoise(['noise pd' n], 1, 1), [], [], ['pd_c' n]);
    
    sim.addElement(NormalNoise(['noise wf' n], [fieldSize_ftr, fieldSize_wd], 1), [], [], ['wf' n]);
    
    sim.addElement(NormalNoise(['noise vis_f' n], [fieldSize_ftr, fieldSize_spt], 1));
    sim.addElement(GaussKernel2D(['noise kernel vis_f' n], [fieldSize_ftr, fieldSize_spt], 1, 1, 1, ...
        true, true, true, kcf), ['noise vis_f' n], [], ['vis_f' n]);
    
    sim.addElement(NormalNoise(['noise atn_f' n], fieldSize_ftr, 1));
    sim.addElement(GaussKernel1D(['noise kernel atn_f' n], fieldSize_ftr, 1, 1, true, true, kcf), ...
        ['noise atn_f' n], [], ['atn_f' n]);
    sim.addElement(NormalNoise(['noise con_f' n], fieldSize_ftr, 1));
    sim.addElement(GaussKernel1D(['noise kernel con_f' n], fieldSize_ftr, 1, 1, true, true, kcf), ...
        ['noise con_f' n], [], ['con_f' n]);
    sim.addElement(NormalNoise(['noise wm_f' n], fieldSize_ftr, 1));
    sim.addElement(GaussKernel1D(['noise kernel wm_f' n], fieldSize_ftr, 1, 1, true, true, kcf), ...
        ['noise wm_f' n], [], ['wm_f' n]);
    
    sim.addElement(NormalNoise(['noise atn_c' n], [fieldSize_ftr, fieldSize_spt], 1));
    sim.addElement(GaussKernel2D(['noise kernel atn_c' n], [fieldSize_ftr, fieldSize_spt], 1, 1, 1, ...
        true, true, true, kcf), ['noise atn_c' n], [], ['atn_c' n]);
    sim.addElement(NormalNoise(['noise wm_c' n], [fieldSize_ftr, fieldSize_spt], 1));
    sim.addElement(GaussKernel2D(['noise kernel wm_c' n], [fieldSize_ftr, fieldSize_spt], 1, 1, 1, ...
        true, true, true, kcf), ['noise wm_c' n], [], ['wm_c' n]);
end


%% settings for element list in advanced parameter panel

elementGroups = sim.elementLabels;
elementsInGroup = sim.elementLabels;

% exclude elements that should not appear in the param panel
excludeLabels = {'sum', 'expand', 'boost', 'scale'};
excludeIndices = zeros(size(elementGroups));
for i = 1 : length(excludeLabels)
    excludeIndices = excludeIndices | strncmp(excludeLabels{i}, elementGroups, length(excludeLabels{i}));
end
elementGroups(excludeIndices) = [];
elementsInGroup(excludeIndices) = [];

% group together elements that only differ in index (and remove the index
for i = 1 : length(elementGroups)
    if elementGroups{i}(1) ~= 'i' % include the inputs
        nonumber = true(1, length(elementGroups{i}));
        nonumber([strfind(elementGroups{i}, '1'), strfind(elementGroups{i}, '2')]) = false;
        elementGroups{i} = elementGroups{i}(nonumber);
    end
end

i = 1;
while i <= length(elementGroups)
    %if ~isempty(elementGroups{i})
    same = strcmp(elementGroups{i}, elementGroups);
    elementsInGroup{i} = elementsInGroup(same);
    same(find(same, 1)) = 0;
    elementsInGroup(same) = [];
    elementGroups(same) = [];
    i = i + 1;
end


