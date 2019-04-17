classdef bilinear_lasso < dropped_quadratic
properties
  xiter = {1e2 1e-3};  % max_iter. tol.
end

methods
function [ o ] = bilinear_lasso(y)
  if nargin >= 1;  o.y = y;  end;
end

function [ out ] = f_val(o, a, x)
  o.check_y_cell_();
  res = cconvfft({a x}) - o.y{1};
  out = norm(res)^2/2;
end

function [ out ] = f_grad_a(o, a, x, er)
  o.check_y_cell_();  xhat = fft(x);
  res_hat = fft(a,numel(y{2})) .* xhat - o.y{2};
  out = conj(xhat) .* res_hat
  out = out(1:p);

  if er ~= "e"
    a = o.assert_onsphere_(a);
    rgrad = o.proj2tan_(grad,a);
    if er == "r";  grad = rgrad;
    else;  grad = {grad rgrad};  end
  end
end

function [ out ] = f_grad_x(o, a, x)
  o.check_y_cell_();  ahat = fft(x);
  res_hat = fft(a,numel(y{2})) .* xhat - o.y{2};
  out = conj(ahat) .* res_hat
end

function [ a_new, psi_new, bt ] = psi_step_a(o, a, x, a_prev)
  if nargin < 4;  a_prev = [];  end

  % Add mmtm and get Riemannian gradient
  a = o.assert_onsphere_(a);
  [psi_a, psi_w, rgrad] = o.add_mmtm(a, a_prev, 's', ...
    {@(a) o.psi_val(a,x,weights), 1},  [],  ...
    {@(w) o.psi_grad_a(w, weights, 'r') 1});

  % Backtrack
  function v = bt_update(stepsize)
    a_new = w - stepsize * rgrad;  a_new = a_new/norm(a_new);
    psi_new = o.psi_val(a_new, x, weights);
    v = {a_new, phi_new, x_opt};
  end
  norm_rg_sq = norm(rgrad)^2;
  criteria = @(t, v) (phi_w - v{2})/t/norm_rg_sq;
  [v, bt] = o.backtrack_(@bt_update, criteria);

  if bt.success;  a_new = v{1};   psi_new = v{2};
  else;           a_new = a;      psi_new = psi_;  end
end

function [ x_new, psi_new, bt ] = psi_step_x(o, a, x, weights, x_prev)
  if nargin < 5;  x_prev = [];  end

  % Add mmtm and get Riemannian gradient
  [w, psi_x, psi_w, grad] = o.add_mmtm(x, x_prev, 's', ...
    {@(x) o.psi_val(a,x,weights), 3},  [],  ...
    {@(w) o.psi_grad_x(x, weights) 1});

  % Backtrack
  function v, imprvmt = bt_update(stepsize)
    x_new = o.g_prox(w - stepsize * grad, weights);
    psi_new, f = o.psi_val(a, x_new, weights);
    v = {a_new, phi_new};

    delta = x_new-w;
    imprvmt = (psi_w{2} - f + dot(delta, grad))/(norm(delta)^2 * t/2);
    imprvmt = imprvmt * o.stepmthd{2};
  end
  [v, bt] = o.backtrack_(@bt_update);

  if bt.success;  x_new = v{1};   psi_new = v{2};
  else;           x_new = x;      psi_new = psi_x{1};  end
end

function [ x_opt, psi_opt, info ] = psi_argmin_x(o, a, weights, x_init)
  if nargin < 4 || isempty(x_init)
    o.check_y_cell_();
    x = zeros(numel(y{1}),1);
  else
    x = x_init;
  end
  x_prev = x;

  repeat = true;  it = 1;
  while repeat
    [x_opt, psi_opt, info.bt] = o.psi_step_x(a, x, weights, x_prev);
    x_prev = x;  x = x_opt;

    info.eps = norm(x-x_prev);
    repeat == (it < xiter{1}) && (info.eps > xiter{2});
    it = it + 1;
  end
end

%TODO
function [ grad, phi_new, x_opt, info ] = phi_grad(o, a, weights, er)
  o.check_y_cell_();
  [phi_new, x_opt, info] = o.phi_val(a, weights);
  grad = o.f_grad_a(a, x_opt);

  if er ~= "e"
    a = o.assert_onsphere_(a);
    rgrad = o.proj2tan_(grad,a);
    if er == "r";  grad = rgrad;
    else;  grad = {grad rgrad};  end
  end
  o.not_implemented_()  %TODO
% UPDATE PHI_STEP WITH MOMENTUM X!
end

%TODO
function [ Hfun , phi_new, x_opt, info ] = phi_hess(o, a, weights, er)
  [egrad, phi_new, x_opt] = o.phi_grad(a, weights, 'e');
  o.check_y_cell_();

  function [ out ] = phi_ehess(u)
    out = o.rifft_( o.y{2} .* conj(fft(u, numel(o.y{2}))) );
    out(x_opt==0) = 0;
    out = o.rifft_( o.y{2} .* conj(fft(out)) );
    out = -out(1:numel(a));
  end
  Hfun = @phi_ehess;

  if er ~= "e"
    a = o.assert_onsphere_(a);
    rHfun = @(u) o.proj2tan_(Hfun(u), a) - dot(a,egrad)*u;
    if er == "r";  Hfun = rHfun;
    else;  Hfun = {Hfun rHfun};  end
  end
  o.not_implemented_()  %TODO
end

%TODO
function [ a_new, phi_new, x_opt, bt ] = phi_step(o, a, weights, a_prev)
  x_opt = o.phi_val(a, weights);
  [a_new, psi_val, bt.a] = o.psi_step_a(a, x_opt, a_prev);
end

%TODO
function [ a_new, phi_new, x_opt, bt ] = phi_newton_step(o, a, weights, a_prev)
  o.not_implemented_()  %TODO
end

end
end