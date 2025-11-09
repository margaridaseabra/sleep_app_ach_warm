function y = bandpass_zero(x, fs, fband)
% Zero-phase FIR bandpass with safe order; falls back to filtfilt(IIR) if needed
x = double(x(:));
nyq = fs/2;
f1 = max(0.1, fband(1)); f2 = min(fband(2), nyq-0.1);
if f2 <= f1, y = x; return; end
try
    % ~2 s window FIR
    ord = max(100, 2*round(fs));
    b = fir1(ord, [f1 f2]/nyq, 'bandpass', blackman(ord+1));
    y = filtfilt(b,1,x);
catch
    d = designfilt('bandpassiir','FilterOrder',8, ...
        'HalfPowerFrequency1',f1,'HalfPowerFrequency2',f2, ...
        'SampleRate',fs);
    y = filtfilt(d, x);
end
end

