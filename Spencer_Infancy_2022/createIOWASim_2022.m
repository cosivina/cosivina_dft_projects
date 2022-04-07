%% create simulator object
sim = Simulator();

%% create neural field
sim.addElement(NeuralField('field a', pa.spatialFieldSize, pa.tau_a, pa.h_a, pa.beta_a));
sim.addElement(NeuralField('field s', pa.spatialFieldSize, pa.tau_s, pa.h_s, pa.beta_s));
sim.addElement(NeuralField('r', 1, pa.tau_s, pa.h_r, pa.beta_r));
sim.addElement(NeuralField('x', 1, pa.tau_a, pa.h_x, pa.beta_x));
sim.addElement(NeuralField('g', 1, pa.tau_a, pa.h_g, pa.beta_g));

%% create lateral interactions
sim.addElement(LateralInteractions1D('input_aa', pa.spatialFieldSize, pa.sigma_aa_exc, pa.c_aa_exc, pa.sigma_aa_inh, pa.c_aa_inh, -pa.c_aa_gi, false, true, pa.kernelWidthMultiplier), ...
    'field a', 'output', 'field a');
sim.addElement(LateralInteractions1D('input_ss', pa.spatialFieldSize, pa.sigma_ss_exc, pa.c_ss_exc, pa.sigma_ss_inh, pa.c_ss_inh, -pa.c_ss_gi, false, true, pa.kernelWidthMultiplier), ...
    'field s', 'output', 'field s');
sim.addElement(ScaleInput('input_rr', 1, pa.c_rr), 'r', 'output', 'r');
sim.addElement(ScaleInput('input_xx', 1, pa.c_xx), 'x', 'output', 'x');
sim.addElement(ScaleInput('input_gg', 1, pa.c_gg), 'g', 'output', 'g');

%% create inputs and projections

%into field a
sim.addElement(CustomStimulus('input_v', zeros(1,pa.spatialFieldSize)), [], [], 'field a'); 
sim.addElement(GaussKernel1D('input_as', pa.spatialFieldSize, pa.sigma_as, pa.c_as, false, true, pa.kernelWidthMultiplier), 'field s', 'output', 'field a');
sim.addElement(GaussStimulus1D('w_ax', pa.spatialFieldSize, pa.sigma_ax_exc, pa.c_ax, spatialHalfSize, false, false));
sim.addElement(ExpandDimension('expand x', 1, [1, pa.spatialFieldSize], 1), 'x', 'output', []);
sim.addElement(PointwiseProduct('input_ax', pa.spatialFieldSize), {'w_ax', 'expand x'}, {'output', 'output'}, 'field a');
sim.addElement(GaussStimulus1D('w_ag_init', pa.spatialFieldSize, pa.sigma_ag_exc, -pa.c_ag_inh, spatialHalfSize, false, false));
sim.addElement(CustomStimulus('w_ag_offset', pa.c_ag_exc*ones(1,pa.spatialFieldSize))); 
sim.addElement(SumInputs('w_ag', pa.spatialFieldSize), {'w_ag_init', 'w_ag_offset'}, [], []);
sim.addElement(ExpandDimension('expand g', 1, [1, pa.spatialFieldSize], 1), 'g', 'output', []);
sim.addElement(PointwiseProduct('input_ag', pa.spatialFieldSize), {'w_ag', 'expand g'}, {'output', 'output'}, 'field a');
sim.addElement(ScaleInput('input_ar', 1, -pa.c_ar_gi), 'r', 'output', 'field a');

%into field s
suppressionPattern = 1 - exp(-0.5 * ((-spatialHalfSize:spatialHalfSize)-0).^2 / pa.sigma_sa^2);
sim.addElement(CustomStimulus('fovPattern', suppressionPattern));
sim.addElement(PointwiseProduct('fovSup', pa.spatialFieldSize), {'fovPattern', 'field a'}, {'output', 'output'});
sim.addElement(GaussKernel1D('input_sa', pa.spatialFieldSize, pa.sigma_sa, pa.c_sa, false, true, pa.kernelWidthMultiplier), 'fovSup', 'output', 'field s');
sim.addElement(ScaleInput('input_sr', 1, -pa.c_sr_gi), 'r', 'output', 'field s');

%into r
sim.addElement(SumDimension('field s -> node', 2, [1, 1], 1), 'field s', 'output');
sim.addElement(ScaleInput('input_rs', 1, pa.c_rs), 'field s -> node', 'output', 'r');

%into x
sim.addElement(CustomStimulus('input_x', 1), [], [], 'x'); 
sim.addElement(ScaleInput('input_xg', 1, -pa.c_gx_inh), 'g', 'output', 'x');
sim.addElement(ScaleInput('input_xr', 1, -pa.c_xr_inh), 'r', 'output', 'x');

%into g
sim.addElement(CustomStimulus('input_g', 1), [], [], 'g'); 
sim.addElement(ScaleInput('input_gx', 1, -pa.c_gx_inh), 'x', 'output', 'g');
sim.addElement(ScaleInput('input_gr', 1, -pa.c_gr_inh), 'r', 'output', 'g');

%% create noise stimulus and noise kernel
sim.addElement(NormalNoise('noise a', pa.spatialFieldSize, 1)); 
sim.addElement(GaussKernel1D('noise kernel a', pa.spatialFieldSize, pa.sigma_q, pa.q_a, false, true, pa.kernelWidthMultiplier), 'noise a', 'output', 'field a');
sim.addElement(NormalNoise('noise s', pa.spatialFieldSize, 1)); 
sim.addElement(GaussKernel1D('noise kernel s', pa.spatialFieldSize, pa.sigma_q, pa.q_s, false, true, pa.kernelWidthMultiplier), 'noise s', 'output', 'field s');
sim.addElement(NormalNoise('noise r', 1, pa.q_r), [], [], 'r');
sim.addElement(NormalNoise('noise x', 1, pa.q_x), [], [], 'x');
sim.addElement(NormalNoise('noise g', 1, pa.q_x), [], [], 'g'); 



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
    same = strcmp(elementGroups{i}, elementGroups);
    elementsInGroup{i} = elementsInGroup(same);
    same(find(same, 1)) = 0;
    elementsInGroup(same) = [];
    elementGroups(same) = [];
    i = i + 1;
end