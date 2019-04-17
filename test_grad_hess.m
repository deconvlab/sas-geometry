function [fx, fx_hat] = test_grad_hess(xt, t, delta, f, f0grad, f0hess)
  x0 = xt(0);  f0 = f(x0);

  if nargin < 6 || isempty(f0hess)
    f0hess = @(x) zeros(size(x0));
  else
    if  isnumeric(f0hess);  f0hess = @(x) f0hess*x;  end
  end

  fx = NaN(size(t));  fx_hat = NaN(size(t));
  for i = 1:numel(t)
    xi = xt(t(i));  del = delta(x0, xi);
    fx(i) = f(xi);
    fx_hat(i) = f0 + dot(del, f0grad) + f0hess(del)'*del;
  end
end