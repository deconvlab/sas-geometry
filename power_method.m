classdef power_method < handle
properties
  H;  U;  Sig;  bias = 0;  tsh = 1e-1;
  mmtm = 0;
end
properties(Access=private)
  u_;  U_ = 0;  Sig_ = 0;
end

methods
function o = power_method(H, bias, tsh, mmtm)
  o.H = H;
  if nargin >= 2 && ~isempty(bias);  o.bias = bias;   end
  if nargin >= 3 && ~isempty(tsh);   o.tsh = tsh;     end
  if nargin >= 3 && ~isempty(mmtm);  o.mmtm = mmtm;   end
end

function success = init(o, u)
  u = u(:);  u = u/norm(u);   % place on the sphere

  % project to perp(U)
  if numel(o.U) > 0;  u = u-o.U*(o.U'*u);  end

  if norm(u) > o.tsh       % add what is left to U
    if numel(o.U) > 0
      o.U_ = o.U;
      o.Sig_ = o.Sig(:);
    end

    o.u_ = u;
    o.U = [o.U u/norm(u)];
    o.Sig = [o.Sig 0];
    success = true;
  else
    warning('Initialization lies in nullspace, ignoring.');
    success = false;
  end
end

function [u, sig] = step(o)
  assert(numel(o.U)>0, "Run the init() method first.");
  o.check_H_fun_();

  u = o.U(:,end);
  if o.mmtm > 0
    u = u + o.mmtm*(u/dot(u,o.u_) - o.u_);
    u = u/norm(u);
  end
  u_proj = (o.U_'*u);
  Hu_bias = o.bias * (u - o.U_*u_proj);
  Hu = o.H(u) - o.U_*(o.Sig_.*u_proj);

  sig = dot(u, Hu);
  u = Hu + Hu_bias;  u = u/norm(u);
  o.u_ = o.U(:,end);  o.U(:,end) = u;  o.Sig(end) = sig;
end

% Helpers
function check_H_fun_(o)
  if isnumeric(o.H)
    H = o.H;  o.H = @(u) H*u;
  end
end

end
end
%#ok<*ALIGN>
%#ok<*PROP>
