classdef dropped_quadratic < handle
properties
  y;
  stepmthd = {'backtrack' 1e-1 1e-1 1e-10};  % mthd. crit. dec. tol.
  mmtm = 0;
end

methods
function [ o ] = dropped_quadratic(y)
  if nargin >= 1;  o.y = y;  end;
end

function [ out ] = f_val(o, a, x)
  o.check_y_cell_();
  out = norm(x)^2/2 - dot(cconvfft({a x}), o.y{1}) + norm(o.y{1})^2/2;
end

% We never use this anywhere
function [ out ] = f_grad_a(o, a, x)
  o.not_implemented_();  out = [];
end

% We never use this anywhere
function [ out ] = f_grad_x(o, a, x)
  o.not_implemented_();  out = [];
end

function [ out ] = g_val(o, x, weights)
  out = sum(x(:).*weights(:));
end

function [ out ] = g_prox(o, x, weights)
  out = sign(x) .* max(abs(x)-weights,0);
end

function [ out, f, g ] = psi_val(o, a, x, weights)
  f = o.f_val(a,x);  g = o.g_val(x, weights);
  out = f + g;
end

%TODO
function [ a_new, psi_new ] = psi_step_a(o, a, x, a_prev)
  o.not_implemented_()
end

function [ x_new, psi_new ] = psi_step_x(o, a, x, weights)
  o.check_y_cell_();
  x_new = o.rifft_(conj(fft(a, numel(o.y{1}))) .* o.y{2});
  x_new = o.g_prox(x_new, weights);
  psi_new = -norm(x_new)^2/2;
end

function [ x_opt, psi_opt ] = psi_argmin_x(o, a, weights)
  [x_opt, psi_opt] = o.psi_step_x(a, [], weights);
end

function [ out, x_opt ] = phi_val(o, a, weights)
  [x_opt, out] = o.psi_argmin_x(a, weights);
end

function [ grad, phi_val, x_opt ] = phi_grad(o, a, weights, er)
  [phi_val, x_opt] = o.phi_val(a, weights);

  o.check_y_cell_();
  grad = o.rifft_( conj(fft(x_opt)).*o.y{2} );
  grad = -grad(1:numel(a));

  if er ~= "e"
    a = o.assert_onsphere_(a);
    rgrad = o.proj2tan_(grad,a);
    if er == "r";  grad = rgrad;
    else;  grad = {grad rgrad};  end
  end
end

function [ Hfun , phi_val, x_opt ] = phi_hess(o, a, weights, er)
  [egrad, phi_val, x_opt] = o.phi_grad(a, weights, 'e');
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
end

function [ a_new, phi_val, x_opt, bt ] = phi_step(o, a, weights, a_prev)
  o.check_y_cell_();  a = o.assert_onsphere_(a);

  % Add mmtm and get Riemannian gradient
  [w, phi_a, phi_w, rgrad] = o.add_mmtm(a, a_prev, 's', ...
    {@(a) o.phi_grad(a,x,weights,'r'), 3}, [], []);
  rgrad = phi_w{1};  phi_w = phi_w{2};

  % Backtrack
  norm_rg_sq = norm(rgrad)^2;
  function v, imprvmt = bt_update(stepsize)
    a_new = w - stepsize * rgrad;  a_new = a_new/norm(a_new);
    [ phi_new, x_opt ] = o.phi_val(a_new, weights);
    v = {a_new, phi_new, x_opt};
    imprvmt = (phi_w - phi_new)/t/norm_rg_sq;
  end
  [v, bt] = o.backtrack_(@bt_update);

  if bt.success
    a_new = v{1};  phi_new = v{2}; x_opt = v{3};
  else
    a_new = a;  phi_new= phi_a{2}; x_opt = phi_a{3};
  end
end

%TODO
function [ a_new, phi_val, x_opt ] = phi_newton_step(o, a, weights, a_prev)
  o.not_implemented_()
end

% Helpers
function [ w,valx,valw,gradw ] = add_mmtm_(o, x, x_, es, valfx, valfw, gradf)
  if ~isempty(valfx);   valfx = @(x) [];      end
  if ~isempty(gradf);   gradf = @(x) [];      end
  if ~iscell(valfx);    valfx = {valfx 1};    end
  if ~iscell(gradf);    gradf = {gradf 1};    end

  % if valfw is empty we simply use valfx
  if ~isempty(valfw) && ~iscell(valfw)
    valfw = {valfw 1};
  end

  valx = cell(1,valfx{2});  valw = cell(1,valfw{2});
  gradw = cell(1,gradf{2});

  if es == 's';  x = o.assert_onsphere_(x);  end
  [valx{:}] = valfx{1}(x);

  if ~isempty(x_) && (o.mmtm > 0)
    if es == 's'
      x_ = o.assert_onsphere_(x_);
      w = (1 + o.mmtm/dot(x, x_))*x - x_;
      w = w/norm(w);
    else;
      w = x + o.mmtm*(x - x_);
    end
    [valw{:}] = valfw{1}(w);
  elseif isempty(valfw)
    w = x;  valw = valx;
  else
    w = x;  [valw{:}] = valfw{1}(w);
  end
  [gradw{:}] = gradf{1}(w);

  if valfx{2}==1;  valx = valx{1};    end
  if valfw{2}==1;  valw = valw{1};    end
  if gradf{2}==1;  gradw = gradw{1};  end
end

function [ out ] = assert_onsphere_(o, in)
  if(abs(1-norm(in))>1e-5)
    warning("Forcing input onto the sphere.");
    out = in/norm(in);
  else;  out = in;  end
end

function [ v, bt ] = backtrack_(o, update)
  bt = true;  stepsize = 1;
  while bt
    % update, then check criteria
    v, imprvmt = update(stepsize);
    success = imprvmt >= o.stepmthd{2};
    if success;  bt = false;  else
      stepsize = o.stepmthd{3}*stepsize;
      bt = stepsize > o.stepmthd{4};
    end
  end
  bt.success = success;  bt.stepsize = stepsize;  bt.imprvmt = imprvmt
end

function check_y_cell_(o)
  assert(~isempty(o.y), "Observation  y  needs to be initialized.");
  if ~iscell(o.y);  o.y = o.y(:);  o.y = {o.y fft(o.y)};  else
    assert(iscell(o.y), "y must either by a vector or a cell.")
    assert(all(size(o.y{2}) == size(o.y{1})), ...
     "y{2} should be the fft of y{1}.");
  end
end

function not_implemented_(o)
  warning('This function is not implemented!')
end

function [ out ] = proj2tan_(o, in, a)
  a = o.assert_onsphere_(a);
  out = in - a*(a'*in);
end

function [ out ] = rifft_(o, in)
  out = real(ifft(in));
end

end
end

%#ok<*ALIGN>
%#ok<*INUSD>
%#ok<*INUSL>
%#ok<*MANU>
%#ok<*STOUT>
