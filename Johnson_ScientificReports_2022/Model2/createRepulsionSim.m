% create simulator object
sim = Simulator();

% create neural field
sim.addElement(NeuralField('field atn_sr', fieldSize, 20));
sim.addElement(NeuralField('field ior_s', fieldSize, 20));
sim.addElement(NeuralField('field vis', [fieldSize, fieldSize], 20, -5, 5));
sim.addElement(NeuralField('field atn', fieldSize, 20, -5, 5));
sim.addElement(NeuralField('field u', fieldSize, 20, -5, 5));
sim.addElement(NeuralField('field v', fieldSize, 5, -5, 5));
sim.addElement(NeuralField('field w', fieldSize, 20, -5, 5));
sim.addElement(NeuralField('field atn_c', fieldSize, 20, -5, 5));
sim.addElement(NeuralField('cos', 1, 20));
sim.addElement(NeuralField('pd_atnc', 1, 20));

% create within field (or within layer) interactions
sim.addElement(LateralInteractions1D('atn_sr -> atn_sr', fieldSize, sigma_exc, 0, sigma_inh, 0, 0, true, true), ...
  'field atn_sr', [], 'field atn_sr');
sim.addElement(LateralInteractions1D('ior_s -> ior_s', fieldSize, sigma_exc, 0, sigma_inh, 0, 0, true, true), ...
  'field ior_s', [], 'field ior_s');

sim.addElement(GaussKernel2D('vis -> vis (exh)', [fieldSize, fieldSize], ...
    sigma_exc, sigma_exc, 0, true, true, true), 'field vis', [], 'field vis');
sim.addElement(GaussKernel2D('vis -> vis (inh)', [fieldSize, fieldSize], ...
    sigma_exc, sigma_exc, 0, true, true, true), 'field vis', [], 'field vis');
sim.addElement(SumAllDimensions('sum vis', [fieldSize, fieldSize]), 'field vis');
sim.addElement(ScaleInput('vis -> vis (global)', 1, 0), 'sum vis', 'fullSum', 'field vis');

sim.addElement(LateralInteractions1D('atn -> atn', fieldSize, sigma_exc, 0, sigma_inh, 0, 0, ...
    true, true), 'field atn', [], 'field atn');

sim.addElement(GaussKernel1D('u -> u', fieldSize, sigma_exc, 0, true, true), 'field u', 'output', 'field u');
sim.addElement(GaussKernel1D('u -> v', fieldSize, sigma_exc, 0, true, true), 'field u', 'output', 'field v');
%sim.addElement(GaussKernel1D('u -> w', fieldSize, sigma_exc, 0, true, true), 'field u', 'output', 'field w');
sim.addElement(SumDimension('sum v', 2, 1, 1), 'field v', 'output');
sim.addElement(GaussKernel1D('v -> u (local)', fieldSize, sigma_inh, 0, true, true), 'field v', 'output', 'field u');
sim.addElement(ScaleInput('v -> u (global)', fieldSize, 0), 'sum v', 'output', 'field u');
sim.addElement(GaussKernel1D('v -> w (local)', fieldSize, sigma_inh, 0, true, true), 'field v', 'output', 'field w');
sim.addElement(ScaleInput('v -> w (global)', fieldSize, 0), 'sum v', 'output', 'field w');
%sim.addElement(GaussKernel1D('w -> u', fieldSize, sigma_exc, 0, true, true), 'field w', 'output', 'field u');
sim.addElement(GaussKernel1D('w -> v', fieldSize, sigma_exc, 0, true, true), 'field w', 'output', 'field v');
sim.addElement(GaussKernel1D('w -> w', fieldSize, sigma_exc, 0, true, true), 'field w', 'output', 'field w');

sim.addElement(LateralInteractions1D('atn_c -> atn_c', fieldSize, sigma_exc, 0, sigma_inh, 0, 0, ...
    true, true), 'field atn_c', [], 'field atn_c');
sim.addElement(ScaleInput('cos -> cos', 1, 0), 'cos', [], 'cos');
sim.addElement(ScaleInput('pd_atnc -> pd_atnc', 1, 0), 'pd_atnc', [], 'pd_atnc');

%add connections between purely spatial fields -- spatial attention (retinal) field
sim.addElement(GaussKernel1D('vis -> atn_sr', fieldSize, sigma_exc, 0, true, true), ...
    'sum vis', 'verticalSum', 'field atn_sr');
sim.addElement(GaussKernel1D('vis -> ior_s', fieldSize, sigma_exc, 0, true, true), ...
    'sum vis', 'verticalSum', 'field ior_s');
sim.addElement(GaussKernel1D('atn_sr -> vis', fieldSize, sigma_exc, 0, true, true), 'field atn_sr');
sim.addElement(ExpandDimension2D('expand atn_sr -> vis', 1, [fieldSize, fieldSize]), ...
    'atn_sr -> vis', [], 'field vis');
sim.addElement(GaussKernel1D('atn_sr -> ior_s', fieldSize, sigma_exc, 0, true, true), 'field atn_sr', [], 'field ior_s');
sim.addElement(GaussKernel1D('ior_s -> atn_sr', fieldSize, sigma_exc, 0, true, true), 'field ior_s', [], 'field atn_sr');

%add projections between feature fields (beyond 3-layer)
sim.addElement(GaussKernel1D('vis -> atn', fieldSize, sigma_exc, 0, true, true), ...
    'sum vis', 'horizontalSum', 'field atn');
sim.addElement(GaussKernel1D('vis -> u', fieldSize, sigma_exc, 0, true, true), ...
    'sum vis', 'horizontalSum', 'field u');
sim.addElement(GaussKernel1D('vis -> w', fieldSize, sigma_exc, 0, true, true), ...
    'sum vis', 'horizontalSum', 'field w');

sim.addElement(GaussKernel1D('atn -> vis', fieldSize, sigma_exc, 0, true, true), 'field atn');
sim.addElement(ExpandDimension2D('expand atn -> vis', 2, [fieldSize, fieldSize]), ...
    'atn -> vis', [], 'field vis');
sim.addElement(GaussKernel1D('atn -> u', fieldSize, sigma_exc, 0, true, true), 'field atn', 'output', 'field u');
sim.addElement(GaussKernel1D('atn -> w', fieldSize, sigma_exc, 0, true, true), 'field atn', 'output', 'field w');
sim.addElement(GaussKernel1D('atn -> atn_c', fieldSize, sigma_exc, 0, true, true), 'field atn', 'output', 'field atn_c');

sim.addElement(GaussKernel1D('u -> atn', fieldSize, sigma_exc, 0, true, true), 'field u', 'output', 'field atn');
%sim.addElement(GaussKernel1D('u -> vis', fieldSize, sigma_exc, 0, true, true), 'field u');
%sim.addElement(ExpandDimension2D('expand u -> vis', 2, [fieldSize, fieldSize]), ...
%    'u -> vis', [], 'field vis');

sim.addElement(GaussKernel1D('w -> atn', fieldSize, sigma_exc, 0, true, true), 'field w', 'output', 'field atn');
sim.addElement(GaussKernel1D('w -> atn_c', fieldSize, sigma_exc, 0, true, true), 'field w', 'output', 'field atn_c');
sim.addElement(GaussKernel1D('atn_c -> w', fieldSize, sigma_exc, 0, true, true), 'field atn_c', 'output', 'field w');
%sim.addElement(GaussKernel1D('w -> vis', fieldSize, sigma_exc, 0, true, true), 'field w');
%sim.addElement(ExpandDimension2D('expand w -> vis', 2, [fieldSize, fieldSize]), ...
%    'w -> vis', [], 'field vis');

%add cos and pd
sim.addElement(ScaleInput('cos -> ior_s', 1, 0), 'cos', [], 'field ior_s');
sim.addElement(ScaleInput('cos -> atn_sr', 1, 0), 'cos', [], 'field atn_sr');
sim.addElement(ScaleInput('cos -> atn', 1, 0), 'cos', [], 'field atn');
sim.addElement(ScaleInput('cos -> atn_c', 1, 0), 'cos', [], 'field atn_c');

sim.addElement(SumAllDimensions('sum atn_c', fieldSize), 'field atn_c', 'output');
sim.addElement(ScaleInput('atn_c -> pd_atnc', 1, 0), 'sum atn_c', [], 'pd_atnc');
sim.addElement(ScaleInput('pd_atnc -> cos', 1, 0), 'pd_atnc', [], 'cos');

% create noise stimulus and noise kernel
sim.addElement(NormalNoise('noise atn_sr', fieldSize, 1));
sim.addElement(GaussKernel1D('noise kernel atn_sr', fieldSize, 0, 0, true, true), 'noise atn_sr', 'output', 'field atn_sr');
sim.addElement(NormalNoise('noise ior_s', fieldSize, 1));
sim.addElement(GaussKernel1D('noise kernel ior_s', fieldSize, 0, 0, true, true), 'noise ior_s', 'output', 'field ior_s');
sim.addElement(NormalNoise('noise vis', [fieldSize, fieldSize], 1));
sim.addElement(GaussKernel2D('noise kernel vis', [fieldSize, fieldSize], 0, 0, 0, ...
    true, true, true), 'noise vis', [], 'field vis');
sim.addElement(NormalNoise('noise atn', fieldSize, 1));
sim.addElement(GaussKernel1D('noise kernel atn', fieldSize, 0, 0, true, true), 'noise atn', 'output', 'field atn');
sim.addElement(NormalNoise('noise u', fieldSize, 1));
sim.addElement(GaussKernel1D('noise kernel u', fieldSize, 0, 0, true, true), 'noise u', 'output', 'field u');
sim.addElement(NormalNoise('noise v', fieldSize, 1));
sim.addElement(GaussKernel1D('noise kernel v', fieldSize, 0, 0, true, true), 'noise v', 'output', 'field v');
sim.addElement(NormalNoise('noise w', fieldSize, 1));
sim.addElement(GaussKernel1D('noise kernel w', fieldSize, 0, 0, true, true), 'noise w', 'output', 'field w');
sim.addElement(NormalNoise('noise atn_c', fieldSize, 1));
sim.addElement(GaussKernel1D('noise kernel atn_c', fieldSize, 0, 0, true, true), 'noise atn_c', 'output', 'field atn_c');
sim.addElement(NormalNoise('noise cos', 1, 1), [], [], 'cos');
sim.addElement(NormalNoise('noise pd_atnc', 1, 1), [], [], 'pd_atnc');

% flat boost for recall
sim.addElement(CustomStimulus('flat boost', flat_pattern));
sim.addElement(ScaleInput('scale flat boost', fieldSize, 0), 'flat boost', [], 'field atn_c');

sim.addElement(CustomStimulus('flat boost 2', flat_pattern));
sim.addElement(ScaleInput('scale flat boost 2', fieldSize, 0), 'flat boost 2', [], 'field atn_c');

sim.addElement(CustomStimulus('flat boost 3', flat_pattern));
sim.addElement(ScaleInput('scale flat boost 3', fieldSize, 0), 'flat boost 3', [], 'field atn_c');

% stimuli
sim.addElement(GaussStimulus2D('vis stim 1', [fieldSize, fieldSize], ...
    sigma_exc, sigma_exc, 0, 90, 20, true, true, false), ...
    [], [], 'field vis');
sim.addElement(GaussStimulus2D('vis stim 2', [fieldSize, fieldSize], ...
    sigma_exc, sigma_exc, 0, 120, 120, true, true, false), ...
    [], [], 'field vis');
sim.addElement(GaussStimulus2D('vis stim 3', [fieldSize, fieldSize], ...
    sigma_exc, sigma_exc, 0, 270, 40, true, true, false), ...
    [], [], 'field vis');

%wheel
sim.addElement(NormalNoise('wheel noise', 1, .25));
sim.addElement(NeuralField('wheel noise shift', 1, 100, 1, 5));
sim.addElement(SumInputs('shifted wheel noise', 1), {'wheel noise', 'wheel noise shift'}, {'output', 'output'});
sim.addElement(GaussStimulus1D('wheel pattern', wheelSize, sigma_exc, 0, fieldSize, true, false));
sim.addElement(PointwiseProduct('noisy wheel pattern', wheelSize), ...
    {'shifted wheel noise', 'wheel pattern'}, {'output', 'output'});
sim.addElement(DiagonalExpansion('wheel stimulus', [1,wheelSize], 1), 'noisy wheel pattern', [], 'field vis');

% diagnostics
sim.addElement(RunningHistory('history w', [1, fieldSize], 2450, 1), 'field w', 'activation');

sim.addElement(SumAllDimensions('sum field vis', [fieldSize, fieldSize]), 'field vis', 'output');