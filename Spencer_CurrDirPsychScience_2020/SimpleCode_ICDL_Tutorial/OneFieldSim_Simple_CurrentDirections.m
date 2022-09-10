%% setting up the simulator
sim = Simulator();

% create inputs (and sum for visualization)
sim.addElement(GaussStimulus1D('stimulus 1', fieldSize, sigma_exc, 0, round(1/4*fieldSize), true, false));
sim.addElement(GaussStimulus1D('stimulus 2', fieldSize, sigma_exc, 0, round(1/2*fieldSize), true, false));
sim.addElement(GaussStimulus1D('stimulus 3', fieldSize, sigma_exc, 0, round(3/4*fieldSize), true, false));
sim.addElement(SumInputs('stimulus sum', fieldSize), {'stimulus 1', 'stimulus 2', 'stimulus 3'});

% create neural field
sim.addElement(NeuralField('field u', fieldSize, 20, -5, 4), 'stimulus sum');
sim.addElement(MemoryTrace('mem u', fieldSize, 1000, 100000, 0.5), 'field u', 'output');

% create interactions
sim.addElement(LateralInteractions1D('u -> u', fieldSize, sigma_exc, 0, sigma_inh, 0, 0, true, true), ...
  'field u', 'output', 'field u');
sim.addElement(GaussKernel1D('mem u -> u', fieldSize, sigma_exc, 0, true, true), 'mem u', 'output', 'field u');

% create noise stimulus and noise kernel
sim.addElement(NormalNoise('noise', fieldSize, 1.0));
sim.addElement(GaussKernel1D('noise kernel', fieldSize, 0, 1.0, true, true), 'noise', 'output', 'field u');

% extra input needed to reset after each trial
sim.addElement(BoostStimulus('stimulus s_1', 0), 'stimulus s_1', 'output', 'field u');

% store activation history
sim.addElement(RunningHistory('history trial', fieldSize, t_max, 1), 'field u', 'activation');









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