clc; clear;
addpath('sbd_tools')

%% Settings
p = 100;
m = 10*p;
bw = 0.2;

theta = 0.2;
lda = 0.01/sqrt(p*theta);

nsamp = 200;
s1 = 0;  s2 = 5;

%% Create objective landscape
[a0, win] = tapered_lowpass(p, bw, 0.05);
a1 = circshift(a0, s1);  a2 = circshift(a0, s2);
x0 = (rand(m,1) <= theta) .* randn(m,1);
solver = lasso_fista();  solver.y = cconvfft(a0,x0);

alph = linspace(0., 2, nsamp);
phiDQ = NaN(nsamp);
for i = 1:nsamp
  for j = 1:nsamp
    if abs(alph(i)) + abs(alph(j)) >= 0.2
      a = alph(i)*a1 + alph(j)*a2;
      phiDQ(i,j) = solver.evaluate(a, lda);
    end
  end
end

%%
imagesc(phiDQ);
cax = caxis;  caxis([cax(1) 0.1]);
colormap jet;
colorbar;
axis equal;