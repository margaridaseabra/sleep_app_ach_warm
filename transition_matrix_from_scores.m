function T = transition_matrix_from_scores(scores)
% scores: integer labels (0=Wake,1=NREM,2=REM)
S = scores(:); S = S(isfinite(S));
states = 0:2; n = numel(states);
C = zeros(n,n);
for i=1:numel(S)-1
    a = S(i)+1; b = S(i+1)+1;
    if a>=1 && a<=n && b>=1 && b<=n
        C(a,b) = C(a,b)+1;
    end
end
P = C ./ max(sum(C,2),1); % row-normalized transition probabilities
T = array2table(P,'VariableNames',{'to_Wake','to_NREM','to_REM'}, ...
                  'RowNames',{'Wake','NREM','REM'});
end
