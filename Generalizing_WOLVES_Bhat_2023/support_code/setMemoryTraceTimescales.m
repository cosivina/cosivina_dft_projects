function  setMemoryTraceTimescales(sim, parBuild, parDecay )
%Update the memory trace field time sclaes
    sim.setElementParameters({'hwm_s'}, {'tauBuild'}, parBuild); 
    sim.setElementParameters({'hcon_s'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hword'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hcon_f1'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hwm_f1'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hwm_c1'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hwf1'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hcon_f2'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hwm_f2'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hwm_c2'}, {'tauBuild'}, parBuild);
    sim.setElementParameters({'hwf2'}, {'tauBuild'}, parBuild);

    sim.setElementParameters({'hwm_s'}, {'tauDecay'}, parDecay); 
    sim.setElementParameters({'hcon_s'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hword'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hcon_f1'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hwm_f1'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hwm_c1'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hwf1'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hcon_f2'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hwm_f2'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hwm_c2'}, {'tauDecay'}, parDecay);
    sim.setElementParameters({'hwf2'}, {'tauDecay'}, parDecay);
    
    %% Update other model parameters
    %sim.setElementParameters({'atn_sa -> atn_c'}, {'amplitude'}, 4.9);

end

