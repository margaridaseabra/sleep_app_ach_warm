function [P,F] = quick_ps(x, fs, method, FR)
x = double(x(:));
switch lower(method)
    case 'spectrogram'
        wlen = max(64, round(0.5*fs));
        nover = round(0.9*wlen);
        nfft = 2^nextpow2(max(wlen, 4*wlen));
        [Sxx,F,~] = spectrogram(x, wlen, nover, nfft, fs, 'yaxis');
        P = abs(Sxx).^2;
    otherwise % 'cwt'
        if ~(exist('cwt','file')==2) % fallback
            wlen = max(64, round(0.5*fs));
            nover = round(0.9*wlen);
            nfft = 2^nextpow2(max(wlen, 4*wlen));
            [Sxx,F,~] = spectrogram(x, wlen, nover, nfft, fs, 'yaxis');
            P = abs(Sxx).^2;
        else
            [cfs,F,~] = cwt(x, fs, 'FrequencyLimits', FR);
            P = abs(cfs).^2;
        end
end
mask = F>=FR(1) & F<=FR(2); F = F(mask); P = P(mask,:);
end

