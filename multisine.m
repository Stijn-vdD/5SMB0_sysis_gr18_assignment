function y = multisine( frequencyLimits, fs, Ns, varargin )
%%% Multisine: a function to generate a multi-sine signal with various
%%% properties, most notably phase distribution to minimise crest-factor.

% Author: Ben Holmes, adapted from a pseudonoise signal by Maarten van
% Walstijn.
% Date: 2019/01/09
% License: All rights reserved. (until the code is cleaned up)

% Inputs
% Required:
%   - frequencyLimits: the boundaries between which all sine components
%   will fall. Each sine will fall at multiples of f0=fs/Ns, which the
%   frequency limits will be rounded to.

%   - fs: sampling frequency.

%   - Ns: signal length in samples. It is recommended for multiple periods
%   of the signal to use repmat instead of a high value of Ns as the
%   Schroeder phases are slow to calculate, and increasing Ns will increase
%   the density of sine wave components.

% Other parameters:
%   - MagnitudeResponse: the amplitude of the sine wave components. Zero
%   values should be used for all magnitudes outside of the frequency
%   limits. Default is a flat response.

%   - PhaseResponse: Either a string to select "Schroeder",
%   "NormallyDistributed", or "ZeroValued" for the different phase options,
%   or a vector of the desired phase values. Default is "Schroeder".

%   - StartAtZero: a boolean flag to indicate whether to wrap the signal
%   such that it starts at the lowest gradient zero crossing. Default true.

%   - Normalise: a boolean flag that indicates whether to normalise the
%   signal to a peak value of 1. Default true.

%   - TimeDomain: a boolean flag that indicates whether to generate the
%   signal in the time domain or frequency domain. Default false. Used for
%   debugging the IFFT.

%   - InitialPhase: a scalar element that sets the first value of the
%   Schroeder phases, ignored for other settings. Default 0.

% Output
%   y: the multi-sine output signal.

%% Input parsing

p = inputParser;

is2ElementPositiveVector =@(x) isnumeric(x) && any(size(x) == 1)...
                                            && any(size(x) == 2)...
                                            && ~any(x < 0);
                                        
addRequired(p, 'frequencyLimits', is2ElementPositiveVector);

isPositiveScalarInteger =@(x) isnumeric(x) && isscalar(x) && (round(x) == x) && x > 0;

addRequired(p, 'fs', isPositiveScalarInteger);
addRequired(p, 'Ns', isPositiveScalarInteger);

addParameter(p, 'MagnitudeResponse', false, @(x) isnumeric(x) && any(size(x) == 1) && ~any(x < 0) && sum(x)>0)

addParameter(p, 'PhaseResponse', 'Schroeder', @(x) ischar(x) || isnumeric(x))

addParameter(p, 'StartAtZero', true, @(x) islogical(x) && isscalar(x));

addParameter(p, 'Normalise', true, @(x) islogical(x) && isscalar(x));

addParameter(p, 'TimeDomain', false, @(x) islogical(x) && isscalar(x));

addParameter(p, 'InitialPhase', 0, @(x) isnumeric(x) && isscalar(x));
         
parse(p, frequencyLimits, fs, Ns, varargin{:});

%% Find frequency indices

f0 = fs/Ns;

% DC is in bin 1 so everything must start at 2
fInds = 1 + round(frequencyLimits./f0);

if any(fInds > Ns/2)
    error('Frequency limits must be beneath Nyquist.');
end

indexVector = fInds(1):fInds(2);

NN = length(indexVector);

%% Find amplitude response

if ~any(p.Results.MagnitudeResponse)
    mag = zeros(1, Ns);
    mag(indexVector) = 1./length(indexVector);
else
    if length(p.Results.MagnitudeResponse) ~= Ns
        error('Magnitude response must be the same length as the desired signal.');
    end
    
    mag = p.Results.MagnitudeResponse.^2;
    
    % Find indices at which no components should be present.
    fullIndices = (1:Ns);
    zeroValueIndices = fullIndices;
    zeroValueIndices(indexVector) = [];
    
    if any(mag(zeroValueIndices))
        warning('Non-zero magnitude values present outside of frequency limits.');
        mag(zeroValueIndices) = zeros(size(zeroValueIndices));
    end
    
    
    if sum(mag(indexVector)) ~= 1
        mag(indexVector) = mag(indexVector)./sum(mag(indexVector));
    end
end

%% Find phase response

if ischar(p.Results.PhaseResponse)
    switch p.Results.PhaseResponse
        case 'Schroeder'
            phase = schroederPhases(NN, Ns, indexVector, mag, p.Results.InitialPhase);
        case 'ZeroPhase'
            phase = zeros(1, Ns);
        case 'NormalDistribution'
            phase = randn(1, Ns);
        otherwise
            error('Phase Response string must be Schroeder, ZeroPhase, or NormalDistribution.');
    end
else
    phase = p.Results.PhaseResponse;
    if ~any(size(phase) == 1) || ~any(size(phase) == Ns)
        error('Phase response must be Ns x 1.');
    end
end


%% Generate multisine signal

% Switch between time and frequency domain generation
if p.Results.TimeDomain
    y = zeros(1, Ns);
    t = (0:Ns-1)./fs;
    for nn=1:NN
        mm = indexVector(nn);
        y = y + sqrt(mag(mm)/2)*cos(2*pi*f0*(mm-1)*t + phase(mm));
    end
else
    % Frequency domain representation
    Y = sqrt(mag/2).*exp(sqrt(-1).*phase);

    % IFFT to time domain
    y = ifft(forceFFTSymmetry(Y))*(Ns/2);
end

%% Normalise peak abs value to 1

if p.Results.Normalise
    y = y./max(abs(y));
end

%% Find zero crossing closest to zero 

% Heuristic method of finding minimum gradient zero crossing
if p.Results.StartAtZero
    ySign = y>0;

    zeroInds = find((ySign(2:end) ~= ySign(1:end-1)));

    % Find the index with the smallest gradient around the zero crossing.
    zeroGrad = zeros(1,length(zeroInds));
    for nn=1:length(zeroInds)
        zeroGrad(nn) = abs(y(zeroInds(nn)) - y(zeroInds(nn)+1));
    end
    [~, minInd] = min(zeroGrad);

    yWrapped = [y(zeroInds(minInd):end) y(1:zeroInds(minInd)-1)];

    y = yWrapped;
end

end

function phase = schroederPhases(nComponents, Ns, indexVector, magnitude, p1)
% Generate phases as defined in "Synthesis of Low-Peak-Factor Signals and
% Binary Sequences With Low Autocorrelation" by M. R. Schroeder

% Preallocate vector
phase = zeros(1, Ns);

% Bin 1 is DC, so bin 2 is the phase of the first component
phase(2) = p1;

% Iterate over phase values using Schroeder's algorithm
for nn=2:nComponents
    ll=1:(nn-1);
    phase(indexVector(nn)) = phase(2) -2*pi*sum((nn-ll).*magnitude(indexVector(ll)));
end

end

function Y = forceFFTSymmetry(X)
% forceFFTSymmetry  A function to force conjugate symmetry on an FFT such that when an
% IFFT is performed the result is a real signal.

% The function has been written to replace MATLAB's ifft(X,'symmetric'), as this function
% is not compatible with MATLAB Coder.
  
Y = X;
XStartFlipped = fliplr(X(2:floor(end/2)));
Y(ceil(end/2)+2:end) = real(XStartFlipped) - sqrt(complex(-1))*imag(XStartFlipped);

end

