eegDir   = '/Users/margaridaseabra/signalnotscored';
scoreDir = '/Users/margaridaseabra/15.11scores';

run_sigma_theta_batch(eegDir, scoreDir, 'SigmaTheta_ModAnalysis');
%%
plot_sigma_modulation_by_condition('SigmaTheta_ModAnalysis');

plot_theta_modulation_by_condition('SigmaTheta_ModAnalysis');
