%% General function for generating filters for use in extracting LFP and MUA activity from raw neural data
function [MUAfilt, LFPfilt, NOTCHfilt] = generate_Filters(FS,MUA_len,MUA_BW,LFP_len,LFP_BW)
%% Instructions
% Takes two parameters and returns a filter. _len represents the sequence
% length in seconds, and _BW represents the effective bandwidth.
% FS is the sampling frequency of the raw datafile.
% Requires the presence of dpsschk and mtfilt in MATLAB's path.

% tapers = [sequence length (seconds), effective bandwidth (Hz)]
%% Initialisation of parameters
if (nargin <= 5)
    NOTCH_len = 0.5;
    NOTCH_BW = 2;
end
if (nargin <= 3)
    LFP_len = 0.0035;%0.0025
    LFP_BW = 300;%500
end
if (nargin <= 1)
    MUA_len = 0.002;%0.002
    MUA_BW = 3000;%2900 3700
end
if (nargin == 0)
    FS = 30000;
end
%% MUAFILTER
tapers = [MUA_len,MUA_BW];
n = tapers(1);
w = tapers(2);
p = n*w;
k = floor(2*p-1);
tapers = [n,p,k];
tapers(1) = tapers(1).*FS;
% These are the Slepian tapers
tapers = dpsschk(tapers);
MUAfilt = 2*real(mtfilt(tapers, FS, 3300)); %3400 4000

%% LFPFILTER
tapers = [LFP_len,LFP_BW];
n = tapers(1);
w = tapers(2);
p = n*w;
k = floor(2*p-1);
tapers = [n,p,k];
tapers(1) = tapers(1).*FS;
tapers =  dpsschk(tapers);
filt = mtfilt(tapers, FS, 0);
LFPfilt = single(filt./sum(filt));

%% NOTCHFILTER
tapers = [NOTCH_len,NOTCH_BW];
n = tapers(1);
w = tapers(2);
p = n*w;
k = floor(2*p-1);
tapers = [n,p,k];
tapers(1) = tapers(1).*FS;
% These are the Slepian tapers
tapers = dpsschk(tapers);
NOTCHfilt = 2*real(mtfilt(tapers, FS, 50));
end