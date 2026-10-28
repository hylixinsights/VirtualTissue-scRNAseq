"""Units and evidence for constants outside the editable CD8 model."""
from .tissue import RADII
MECHANICS_UNITS={
 'dt_hours':'hours', 'packing_fraction':'area fraction', 'immune_speed_um_h':'micrometers/hour',
 'contact_gap_um':'micrometers', 'perception_um':'micrometers', 'cycle_hours':'hours',
 'attack_damage_per_hour':'normalized damage/hour', 'fatigue_per_hour':'normalized fatigue/hour',
 'recovery_per_hour':'normalized fatigue/hour','corpse_hours':'hours','cue_diffusion_um2_h':'micrometers^2/hour',
 'cue_decay_per_hour':'1/hour','cue_source_per_hour':'arbitrary cue units/hour per tumor cell'}
FIXED_CONSTANTS={
 'radii_um':RADII, 'cue_grid_side':32, 'fatigue_attack_migration_gate':.85,
 'tumor_cycle_baseline_fraction':.2,'tumor_cycle_RNA_fraction':.8,
 'division_baseline_fraction':.2,'division_antigen_fraction':.8,
 'effector_transition_baseline_fraction':.1,'effector_transition_antigen_fraction':.9,
 'dysfunction_transition_baseline':.15,'dysfunction_tumor_fraction_coefficient':1.,
 'dysfunction_antigen_coefficient':1.,'movement_dysfunction_penalty':.4,
 'exploratory_direction_normal_sd':.65,'migration_fatigue_recovery_fraction':.5,
 'minimal_domain_side_um':100.}

def parameter_ledger(config):
    return {'model':config['model']['parameter_provenance'],
     'mechanics':{k:{'value':config['mechanics'][k],'unit':u,'evidence_class':'assumed','source':'A05 software scenario, not measured or fitted'} for k,u in MECHANICS_UNITS.items()},
     'fixed_constants':{'values':FIXED_CONSTANTS,'evidence_class':'assumed',
        'units':'radii/minimal side in micrometers; grid in bins; other coefficients dimensionless',
        'source':'Versioned implementation constants; no numerical biological estimates'},
     'output_precision':{'position_um':1e-5,'state_time_hours':1e-8,'kind':'recording precision, not measurement accuracy'}}
