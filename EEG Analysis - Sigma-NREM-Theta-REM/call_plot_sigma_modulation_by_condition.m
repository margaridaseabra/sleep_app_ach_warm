eegDir   = '/Users/margaridaseabra/24.11notchedsignal';
scoreDir = '/Users/margaridaseabra/24.11scores';

run_sigma_theta_batch(eegDir, scoreDir, 'SigmaTheta_ModAnalysis');
%%
plot_sigma_modulation_by_condition('SigmaTheta_ModAnalysis');
%%
plot_theta_modulation_by_condition('SigmaTheta_ModAnalysis');
%% 
plot_rem_theta_power_by_condition('SigmaTheta_ModAnalysis');
