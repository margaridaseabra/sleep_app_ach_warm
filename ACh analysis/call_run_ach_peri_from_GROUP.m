% recompute peri matrices using the new transition rules
GROUP_TRANS = run_ach_peri_from_GROUP('ACh_GROUP_all.mat');

% group-level plots
GROUP_STATE_NREM = ach_group_state_plot('ACh_periTransitions_all.mat','NREM');
GROUP_STATE_Wake = ach_group_state_plot('ACh_periTransitions_all.mat','Wake');
GROUP_STATE_REM  = ach_group_state_plot('ACh_periTransitions_all.mat','REM');

%%
% NREM onset heatmaps
ach_group_state_heatmap('ACh_periTransitions_all.mat','NREM','WT','baseline');
ach_group_state_heatmap('ACh_periTransitions_all.mat','NREM','APP','baseline');
ach_group_state_heatmap('ACh_periTransitions_all.mat','NREM','WT','ambtemp');
ach_group_state_heatmap('ACh_periTransitions_all.mat','NREM','APP','ambtemp');
ach_group_state_heatmap('ACh_periTransitions_all.mat','NREM','WT','drugs');
ach_group_state_heatmap('ACh_periTransitions_all.mat','NREM','APP','drugs');

%%
% Wake onset
ach_group_state_heatmap('ACh_periTransitions_all.mat','Wake','WT','drugs');
ach_group_state_heatmap('ACh_periTransitions_all.mat','Wake','APP','drugs');
%%
% REM onset
ach_group_state_heatmap('ACh_periTransitions_all.mat','REM','WT','baseline');
ach_group_state_heatmap('ACh_periTransitions_all.mat','REM','APP','baseline');
ach_group_state_heatmap('ACh_periTransitions_all.mat','REM','WT','ambtemp');
ach_group_state_heatmap('ACh_periTransitions_all.mat','REM','APP','ambtemp');
ach_group_state_heatmap('ACh_periTransitions_all.mat','REM','WT','drugs');
ach_group_state_heatmap('ACh_periTransitions_all.mat','REM','APP','drugs');
