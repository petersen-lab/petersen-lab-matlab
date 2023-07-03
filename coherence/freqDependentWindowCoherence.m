% performs coherence/power analysis on spike trains
% Inputs: spk1, spk2, dt - spike counts in windows of length dt seconds (as column vectors)
%         spk3 - optional, if provided computes partical coherence (of spk1 and spk2 given spk3)
%         opt: structure with the following fields
%              .typespk1 - 'pb' (default) or 'c'
%              .typespk2 - 'pb' (default) or 'c'
%              .typespk3 - 'pb' (default) or 'c'
%              .winfactor - 5 by default (each window is at least this many times than 1/(highest frequency) we estimate using it)
%              .freqfactor - 2 by default (has to be > 1, each window is at least opt.winfactor/opt.freqfacor times than 1/(lowest frequency) we estimate using it)
%              .tapers - 3 by default (one number!)
%              .debug - if true will print information on sizes of windows used
%              .maxFreq - by default it is 0.5/dt
%              .decimate - true by default, if false will not decimate the signals (for low frequencies) to reduce runtime
%              .monotoneFreq - false by default, if true the results won't have the same frequency twice (possible at transition across different window,
%                              and serves as a control for the effect of window size)
%              .jack - to use jackknife error estimates
%              .pad - a scalar representing
%                     FFT padding (can take values -1,0,1,2...). -1 corresponds to no
%                     padding, 0 corresponds to padding to the next highest power of 2 etc.
%                     e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad
%                     the FFT to 512 points, if pad=1, we pad to 1024 points etc. Setting
%                     pad to larger values, gives denser frequency grids. Defaults to 0.
%              .minFreq - by default it is 0.
% Outputs: freq - frequencies
%          coh  - coherence (or PSD if input is a single signal)
%          phase- phase (of spk2 wrt spk1, e.g. a small positive (linear) phase means spk2 precedes spk1
%          conf - 95% confidence interval for coherence (or PSD)
%          phaseCu, phaseCl - upper and lower confidence intervals for phase
% For spk2 = [] will compute the PSD (in coh output, phi is meaningless).
% For PSD estimation in the continuous case, provided the entire bandwith is used, it should hold that trapz(freq, coh) equals to var(spk1)/2
%
% The function was originally written by Michael Okun and later modified by
% Martynas Dervinis (martynas.dervinis@gmail.com).
function [freq, coh, phase, conf, phaseCu, phaseCl] = freqDependentWindowCoherence(spk1, spk2, dt, spk3, opt)

if nargin < 4
  spk3 = [];
end
if nargin < 5
  opt = [];
end
if ~isfield(opt, 'typespk1') || isempty(opt.typespk1)
  opt.typespk1 = 'pb';
end
if ~isfield(opt, 'typespk2') || isempty(opt.typespk2)
  opt.typespk2 = 'pb';
end
if ~isfield(opt, 'typespk3') || isempty(opt.typespk3)
  opt.typespk3 = 'pb';
end
if ~isfield(opt, 'winfactor') || isempty(opt.winfactor)
  opt.winfactor = 5;
end
if ~isfield(opt, 'freqfactor') || isempty(opt.freqfactor)
  opt.freqfactor = 2;
end
assert(opt.freqfactor > 1)
if ~isfield(opt, 'tapers') || isempty(opt.tapers)
  opt.tapers = 3;
elseif numel(opt.tapers) > 1
  error('this format of opt.tapers not supported')
end
if ~isfield(opt, 'debug') || isempty(opt.debug)
  opt.debug = false;
end
if ~isfield(opt, 'monotoneFreq') || isempty(opt.monotoneFreq)
  opt.monotoneFreq = false;
end
if isfield(opt, 'ag')
  error('no longer supported')
end
if ~isfield(opt, 'maxFreq') || isempty(opt.maxFreq)
  opt.maxFreq = 0.5/dt;
elseif opt.maxFreq > 0.5/dt
  error('opt.maxFreq > Nyquist?!')
end
if ~isfield(opt, 'decimate') || isempty(opt.decimate)
  opt.decimate = true;
end
if ~isfield(opt, 'pad') || isempty(opt.pad)
  opt.pad = 0;
end
if ~isfield(opt, 'minFreq') || isempty(opt.minFreq)
  opt.fpass = 0;
end

assert(opt.winfactor >= 2, 'winfactor of less than 2 not allowed')

% There are two issues:
% * To have an estimate for some particular frequency f, we need at least
% 30 intervals of duration 1/f. So the length of available data provides an
% automatic constraint on how low we can go.
%
% * On the other hand for (high) frequency f, no point using intervals
% longer than 20/f.

spk1 = torow(spk1)';
spk2 = torow(spk2)';
spk3 = torow(spk3)';
if ~isempty(spk2)
  assert(numel(spk1) == numel(spk2), 'inputs of different length?!')
end

params.Fs = 1/dt;

if ~isfield(opt, 'jack') || isempty(opt.jack) || ~opt.jack
  opt.jack = false;
end
params.err = [1 0.05]; % will be switched to [2 0.05] where appropriate in the code below...

%startFreq = 10.1*(opt.winfactor / (numel(spk1)*dt));
startFreq = opt.maxFreq;

freq = [];
coh = [];
phase = [];
conf = [];
phaseCu = []; phaseCl = [];
while true
  L = opt.winfactor/startFreq;
  if opt.winfactor*L > numel(spk1)*dt
    if opt.debug
      fprintf('Cannot use windows of %2.3f sec. (total recording length %2.1f sec) ==> done\n', ...
        L, numel(spk1)*dt);
    end
    break % We cannot have even 10 repeats with windows of this size
  end
  params.fpass = [max([startFreq/opt.freqfactor opt.minFreq]) startFreq];
  if params.fpass(1) >= params.fpass(2)
    break
  end
  params.tapers = [opt.tapers/L L opt.tapers];
  if params.Fs > 500*params.fpass(1) && opt.decimate
    % we're now estimating frequencies that are way lower than signals'
    % sampling rate, so we will decimate the signals, to make the
    % computations more efficient
    spk1 = spk1(1:floor(numel(spk1)/10)*10);
    spk2 = spk2(1:floor(numel(spk2)/10)*10);
    spk3 = spk3(1:floor(numel(spk3)/10)*10);
    if strcmpi(opt.typespk1, 'pb')
      spk1 = sum(reshape(spk1, 10, []))';
    else
      spk1 = decimate(spk1, 10);
      spk1 = spk1 - mean(spk1);
    end
    if ~isempty(spk2) && strcmpi(opt.typespk2, 'pb')
      spk2 = sum(reshape(spk2, 10, []))';
    elseif ~isempty(spk2)
      spk2 = decimate(spk2, 10);
      spk2 = spk2 - mean(spk2);
    end
    if ~isempty(spk3) && strcmpi(opt.typespk3, 'pb')
      spk3 = sum(reshape(spk3, 10, []))';
    elseif ~isempty(spk3)
      spk3 = decimate(spk3, 10);
    end
    if opt.debug
      fprintf('Signal decimation (from %2.2f Hz to %2.2f Hz)\n', params.Fs, params.Fs/10)
    end
    params.Fs = params.Fs/10;
    dt = 1/params.Fs;
  end % signal decimation

  if isempty(spk2) % PSD computation
    switch opt.typespk1
      case 'pb'
        [C, f, ~, ~, ~, ~, confC] = mtspectrumsegpb(spk1, L, params, true, true);
      case 'c'
        [C, f, ~, ~, confC] = mtspectrumsegc(spk1, L, params, true);
      otherwise
        error('unknown type of spk1')
    end
    phi = []; phiCu = []; phiCl = [];
    confC = abs(confC(1,:)'-C);
  elseif ~isempty(spk3) % partial coherence computation
    assert(max(spk1) - min(spk1) > 0 && max(spk2) - min(spk2) > 0, 'spk1 or/and spk2 is constant?!')
    if strcmpi(opt.typespk1, 'pb') && strcmpi(opt.typespk2, 'pb') && strcmpi(opt.typespk3, 'pb')
      [f, C] = partialCoherenceChronux(spk1, spk2, spk3, 'segpb', params, L);
    elseif strcmpi(opt.typespk1, 'c') && strcmpi(opt.typespk2, 'c') && strcmpi(opt.typespk3, 'c')
      [f, C] = partialCoherenceChronux(spk1, spk2, spk3, 'segc', params, L);
    else
      error('such partial coherence not supported')
    end
    phi = []; phiCu = []; phiCl = [];
    confC = [];
  else % coherence computation
    if ~sum(spk1) || ~sum(spk2)
      f = mean(params.fpass); C = NaN; phi = NaN; confC = NaN; phiCu = NaN; phiCl = NaN;
    else
      if opt.jack && opt.tapers * ceil(numel(spk1)*dt / L) < 1000 % indication for degree of freedoms relevant for jackknife
        % when the above product is too high, jackkinife is prohibitively long...
        params.err(1) = 2;
      end
      assert(max(spk1) - min(spk1) > 0 && max(spk2) - min(spk2) > 0, 'spk1 or/and spk2 is constant?!')
      if strcmpi(opt.typespk1, 'pb') && strcmpi(opt.typespk2, 'pb')
        [~,phi,~,~,~,f,zerosp] = coherencysegpb(spk1, spk2, L, params, false);
        if params.err(1) == 2
          [C, ~, ~, ~, ~, ~, ~, ~, ~, confC] = coherencysegpb(spk1, spk2, L, params, true); % recompute with averaging, to get confidence for coherence
        else
          [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegpb(spk1, spk2, L, params, true); % recompute with averaging, to get confidence for coherence
        end
        %[~,~,S12sh,S1sh,S2sh] = coherencysegpb(spk1, circshift(spk2, round(numel(spk2)/2)), L, params, false);
      elseif strcmpi(opt.typespk1, 'c') && strcmpi(opt.typespk2, 'c')
        assert(abs(mean(spk1)) < 1e-3 && abs(mean(spk2)) < 1e-3, 'continuous data must be mean subtracted')
        [~,phi,~,~,~,f] = coherencysegc(spk1, spk2, L, params, false);
        if params.err(1) == 2
          [C, ~, ~, ~, ~, ~, ~, ~, ~, confC] = coherencysegc(spk1, spk2, L, params, true);
        else
          [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegc(spk1, spk2, L, params, true);
        end
      elseif strcmpi(opt.typespk1, 'c') && strcmpi(opt.typespk2, 'pb')
        assert(abs(mean(spk1)) < 1e-3, 'continuous data must be mean subtracted')
        [~,phi,~,~,~,f,zerosp] = coherencysegcpb(spk1, spk2, L, params, false);
        if params.err(1) == 2
          [C, ~, ~, ~, ~, ~, ~, ~, ~, confC] = coherencysegcpb(spk1, spk2, L, params, true);
        else
          [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegcpb(spk1, spk2, L, params, true);
        end
      elseif strcmpi(opt.typespk1, 'pb') && strcmpi(opt.typespk2, 'c')
        assert(abs(mean(spk2)) < 1e-3, ['continuous data must be mean subtracted' num2str(abs(mean(spk2)))])
        [~,phi,~,~,~,f,zerosp] = coherencysegcpb(spk2, spk1, L, params, false);
        phi = -phi;
        if params.err(1) == 2
          [C, ~, ~, ~, ~, ~, ~, ~, ~, confC] = coherencysegcpb(spk2, spk1, L, params, true);
        else
          [C, ~, ~, ~, ~, ~, ~, confC] = coherencysegcpb(spk2, spk1, L, params, true);
        end
      else
        error('unknown type of spk1 & spk2')
      end
      % C = abs(mean(S12, 2) ./sqrt(mean(S1,2).*mean(S2,2))); % this produces the very same result as the above...
      if numel(confC) == 1 % confC is just a number which we turn into a vector corresponding to frequencies...
        confC = confC*ones(size(f));
      elseif params.err(1) == 2 % for Jackknife confidence intervals...
        assert(size(confC, 2) == numel(f) && size(confC, 1) == 2)
        confC = confC(2,:)-torow(C);
      end
      % for coherency between two continuous signals confC is actually a vector (for each frequency)

      % In some intervals there may be 0 spikes, and so C & phi are NaN there
      if exist('zerosp', 'var')
        [phi, phiCu, phiCl] = circ_mean(phi(:, ~zerosp), [], 2);
      else % when both signals are continuous
        [phi, phiCu, phiCl] = circ_mean(phi, [], 2);
      end
      % In addition, I've encountered rare cases (where the mean is on just 2 values?) where confidence intervals were complex
      phiCu(abs(phiCu - real(phiCu)) > 0) = NaN;
      phiCl(abs(phiCl - real(phiCl)) > 0) = NaN;
    end

    if opt.debug
      fprintf('Using windows of %2.3f sec. for frequencies in the %2.3f - %2.3f Hz range (%d tapers gives resolution of %2.3f Hz), conf: %d \n', ...
        L, params.fpass, opt.tapers, params.tapers(1), params.err(1));
    end
  end

  freq = [freq torow(f)];
  coh = [coh torow(C)];
  phase = [phase torow(phi)];
  conf = [conf torow(confC)];
  phaseCu = [phaseCu torow(phiCu)]; phaseCl = [phaseCl torow(phiCl)];
  startFreq = startFreq/opt.freqfactor;
end
[~,i] = sort(freq);
if opt.monotoneFreq
  I = [true diff(freq(i)) > 0];
  i = i(I);
end
freq = freq(i);
coh = coh(i);
if ~isempty(phase)
  phase = phase(i);
  phaseCu = phaseCu(i); phaseCl = phaseCl(i);
end
if ~isempty(conf)
  conf = conf(i);

  % I have seen that it can happen that the function finds a phase (with a quite small confidence interval!)
  % where there is none by construction (and the coherence was correctly identified as 0).
  % Therefore I do the following:
  phaseCu(coh - conf <= 0) = NaN;
  phaseCl(coh - conf <= 0) = NaN;
end

%------ For debugging only:
% figure;
% semilogx(freq, coh);
% hold on, semilogx(freq, coh+conf, '--k');
% hold on, semilogx(freq, coh-conf, '--k');

% % if x is a vector returns a matrix where each column is k sequential values of x
% % if x is a matrix, returns a new matrix where each row is k original rows
% % Note that elements from the tail of x can be dropped (if number of elements/rows is not divisable by k)
% function x = agregate(x, k)
% if isvector(x)
%   x = x(1:floor(numel(x)/k)*k);
%   x = reshape(x, k, []);
%   return
% end
% % we're here ==> x is a matrix
% x = x(1:floor(size(x,1)/k)*k,:);
% x = reshape(x', [], size(x,1)/k)';

