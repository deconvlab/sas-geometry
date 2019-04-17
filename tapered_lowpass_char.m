function [ mu ] = tapered_lowpass_char(p, m, taper, bwres, trials, pltchar)
if nargin < 3 || isempty(taper);    taper = .05;      end
if nargin < 4 || isempty(bwres);    bwres = 100;      end
if nargin < 5 || isempty(trials);   trials = 20;      end
if nargin < 6 || isempty(pltchar);  pltchar = true;   end

if numel(bwres) == 1
  bw = linspace(0.01, 1, bwres);
elseif numel(bwres) == 3
  bw = linspace(bwres(1), bwres(2), bwres(3));
else
  bw = bwres;
end

mu = zeros(numel(bw),trials);
for i = 1:numel(bw)
  for j = 1:trials
    a = tapered_lowpass(p, bw(i), taper);
    mu(i,j) = shift_coherence(a, [], m);
  end
end

if pltchar
  plot(bw, mean(mu, 2));
  xlabel('Bandwidth');  xlim([min(bw) max(bw)]);
  ylabel('Shift Coherence');
end
end