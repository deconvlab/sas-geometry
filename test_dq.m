clc;  run('sbd_tools/initpkg.m');
set(0,'DefaultFigureWindowStyle','docked')

fidx = 1:3;  
fidx = 1*numel(fidx)+fidx;

p2s = @(a) a/norm(a);
p2t = @(x,a) x - dot(x,a)*a;

% TODO:
%  - Find rgrad for norm(DQ rgrad)^2
%  - Test convergence of phi_step using data init

%% Problem setup
p0 = 5e3;
m = p0*1e2;
p = 3*p0-2;

bw = .7;  taper = .05;
a0 = tapered_lowpass(p0, bw, taper);
%a0 = randn(p0,1);  a0 = a0/norm(a0);
fprintf('Shift coherence of a0 = %.4f.\n', ...
  shift_coherence(a0, [], m));

s0 = [a0 ; zeros(p-p0,1)];
ell = randi(p-p0);
s1 = circshift(s0, ell);

theta = p0^(-3/4);
x0 = randn(m,1) .* (rand(m,1) <= theta);

y = cconvfft({a0 x0});
dq = dropped_quadratic(y);
lambda = .8/sqrt(p0*theta);

%% Plot objective and gradient between two shifts
alph = linspace(-.5, 1.5, 100);
obj = NaN(size(alph));
rgnorm = NaN(size(alph));
for i = 1:numel(alph)
  a = (1-alph(i))*s0 + alph(i)*s1;
  a = a/norm(a);
  obj(i) = dq.phi_val(a, lambda);
  rg = dq.phi_grad(a, lambda, 'r');
  rgnorm(i) = norm(dot(rg, p2t(s1-s0, a)));
end

figure(fidx(1)); clf;
yyaxis left; plot(alph, obj);
yyaxis right; plot(alph, rgnorm);
drawnow;
% magnitude of objective plot increases w.r.t. mu

%% Test gradient and Hessian in random projections
t = .51;  
a_0 = (1-t)*s0 + t*s1;  a_0 = a_0/norm(a_0);
ts = linspace(-1,1,50) * 1e-2;
del = s1-s0; %randn(p,1);
del = del/norm(del);

args = { 
  @(t) p2s(a_0 + t*del), ts, ...
  @(x0, xi) acos(dot(x0,xi)) * p2s(p2t(xi,x0)), ...
  @(a) dq.phi_val(a, lambda), ...
  dq.phi_grad(a_0, lambda, 'r'), ...
  dq.phi_hess(a_0, lambda, 'r') 
};
%figure(fidx(2));
%[phi, phi_hat] = test_grad_hess(args{:});
%plot(ts, [phi; phi_hat]'); drawnow;

%% Check Hessian at shifts and balanced points
t = 0.5;
a = (1-t)*s0 + t*s1;  a = a/norm(a);
rH = dq.phi_hess(a, lambda, 'r');

pwrmthd = power_method(rH, -1e3);
pwrmthd.mmtm = 0.1;
ncomps = 5;  niter = 1e2;
sig = NaN(ncomps, niter);  j = 1;
disp(' ');
while (j<=ncomps) && pwrmthd.init(randn(p,1))
  for i = 1:niter;  [u, sig(j,i)] = pwrmthd.step();  end
  fprintf('Recovered component %d.\n', j);  j=j+1;
end
disp(' ');  figure(fidx(2));  clf; 
plot(5:niter, sig(:,5:end)'); drawnow;

%% Try SaS-BD using phi_step
a = 0.6;
a = (1-a)*s0 + a*s1;  a = a/norm(a);
a_prev = a;

niter = 20;  costs = NaN(1,niter);
btinfo = zeros(2, niter);
for i = 1:niter
  [a_new, costs(i), x_opt, bt] = dq.phi_step(a, lambda, a_prev);
  a_prev = a;  a = a_new;
  btinfo(1,i) = bt.success;  btinfo(2,i) = bt.stepsize;
end

figure(fidx(3)); clf;
stem(a0, '.'); hold on;
stem(a, '.'); hold off;

ax1 = gca;
color2 = [0 .8 0];
ax2 = axes('Position', ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none', 'XColor',color2, 'YColor', color2);
line(1:niter, costs,'Parent',ax2,'Color',color2,'LineWidth',1.5)

[mu,ell_hat,xr] = shift_coherence(a0, a, m);
fprintf('Shift coherence between a0, a = %.2f.\n', mu); drawnow;

%% Finish.
disp(' ');
