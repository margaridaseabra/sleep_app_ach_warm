function bin_out = runlength_prune(bin, min_dur, min_gap)
% Keep "1" runs >= min_dur; merge gaps < min_gap; output logical row
bin = bin(:)'; n = numel(bin);
% Merge small gaps
i = 1;
while i <= n
    if bin(i)==1
        j = i;
        while j<=n && bin(j)==1, j=j+1; end
        % j is first 0 after run
        % check following gap
        g0 = j;
        while g0<=n && bin(g0)==0, g0=g0+1; end
        gap = g0-j;
        if gap>0 && gap < min_gap
            bin(j:g0-1) = 1;  % merge gap
            j = g0;           % continue run
            while j<=n && bin(j)==1, j=j+1; end
        end
        i = j;
    else
        i = i+1;
    end
end
% Remove short runs
bin_out = false(1,n);
i = 1;
while i <= n
    if bin(i)==1
        j = i;
        while j<=n && bin(j)==1, j=j+1; end
        if (j-i) >= min_dur
            bin_out(i:j-1) = true;
        end
        i = j;
    else
        i = i+1;
    end
end
end


