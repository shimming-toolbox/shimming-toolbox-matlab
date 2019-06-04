function y = bpfilt(signal, f1, f2, fs, isplot)
%% Bandpass filtering
%
% Syntax:
%   y = bpfilt(signal, f1, f2, [options])
%
% Description:
%   This function performs bandpass filtering of a time series 
%   with rectangle window.
%
% Input Arguments:
%   signal 	- a column vector of time series.
%   f1 		- the lower bound of frequencies (in Hz).
%   f2 		- the upper bound of frequencies (in Hz).
%
% Options:
%   fs      - the sampling frequency in Hz. Default is 1 Hz.
%   isplot  - whether to produce plots.
%
% Output Arguments:
%   y 		- the filtered time series.
%
% Examples:
%   fs = 100;
%   t  = 1:1/fs:10;
%   x  = sin(t);
%   y  = bpfilt(x,20,30);
%
% See also 
%
% References:
%
% History:
%   07/13/2016 - Initial script.
%   04/02/2019 by ryan.topfer@polymtl.ca
%
%__________________________________________________________________________
% Wonsang You(wsgyou@gmail.com)
% 07/13/2016
% Copyright (c) 2016 Wonsang You.

%% getting options
if nargin < 4 || isempty(fs)
	fs = 1;
end

if nargin < 5 || isempty(isplot)
	isplot = 1;
end

%% define variables
if isrow(signal)
    signal = signal';
end
N  = length(signal);
dF = fs/N;
f  = (-fs/2:dF:fs/2-dF)';

%% Band-Pass Filter:
if isempty(f1) || f1==-Inf
    BPF = (abs(f) < f2);
elseif isempty(f2) || f2==Inf
    BPF = (f1 < abs(f));
else
    BPF = ((f1 < abs(f)) & (abs(f) < f2));
end

% Power spectrum of the original signal
% signal 	 = signal-mean(signal);
spectrum = fftshift(fft(signal))/N;

%% Power spectrum of the band-pass filtered signal
filteredSpectrum = BPF.*spectrum;

%% The band-pass filtered time series
y = N*ifft(ifftshift(filteredSpectrum)); %inverse ifft
y = real(y);

%% Power spectrum of the band-pass filter
if isplot
    figure;
    subplot(311);
    plot(f,abs(spectrum));
    title('Power spectrum of the original signal');
    
    subplot(312);
    plot(f,BPF);
    title(sprintf('Power spectrum of the band-pass filter in (%.3f, %.3f) Hz',f1,f2));
    
    subplot(313);
    plot(f,abs(filteredSpectrum));
    title(sprintf('Power spectrum of the band-pass filtered signal in (%.3f, %.3f) Hz',f1,f2));
    
    time = 1/fs*(0:N-1)';
	
    figure;
    subplot(2,1,1);
    plot(time,signal);
    axis tight
    yl = ylim();
    title('The original time series');
    subplot(2,1,2);
    plot(time,y);
    ylim(yl);
    title(sprintf('The band-pass filtered time series in (%.3f, %.3f) Hz',f1,f2));
end

end
