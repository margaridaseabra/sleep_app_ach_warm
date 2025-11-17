GROUP = run_ach_batch_auto( ...
    '/Users/margaridaseabra/signalnotscored', ...
    '/Users/margaridaseabra/15.11scores');
%% 2) Run group comparison plots (APP vs WT across conditions)
run_ach_group_plots(GROUP);
%%
group_plot_all_transitions_paper(GROUP);
