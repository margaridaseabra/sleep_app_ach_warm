function plot_bandpower_change_by_genotype(csvFile, state, bandsUse, outPng)

    T = readtable(csvFile);
    T.MouseID   = string(T.MouseID);
    T.Genotype  = string(T.Genotype);
    T.Condition = string(T.Condition);
    T.State     = string(T.State);
    T.Band      = string(T.Band);

    % Only baseline + ambtemp for this state
    T = T(ismember(lower(T.Condition),["baseline","ambtemp"]) & T.State==state, :);

    mice = unique(T.MouseID);
    figure(...)

    for each band in bandsUse
        for each mouse
            take baseline Power_dB, ambtemp Power_dB
            delta = amb - base
            store (MouseID, Genotype, Delta)
        end

        % two groups: WT, APP
        valsWT  = Delta(Genotype=="WT");
        valsAPP = Delta(Genotype=="APP");

        % box/scatter plot WT vs APP
        % ttest(valsWT,0) -> star above WT
        % ttest(valsAPP,0) -> star above APP
        % ttest2(valsWT,valsAPP) -> star between groups
    end
end
