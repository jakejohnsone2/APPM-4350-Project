function jp = circ_roots(nmax, smax)
%CIRC_ROOTS  Positive zeros of J_n'(x) for free-edge circular-plate modes.
%   jp = circ_roots(nmax, smax) returns jp(n+1, s) = s-th positive zero of
%   J_n'(x), n = 0..nmax, s = 1..smax. For n = 0 these are the zeros of J_1
%   (the rigid-body zero at x = 0 is not included).
dJ = @(nu, x) 0.5*(besselj(nu-1, x) - besselj(nu+1, x));   % J_n' for any n >= 0
jp = zeros(nmax+1, smax);
for n = 0:nmax
    xmax = (smax + n/2 + 1)*pi + n;              % generous upper bound
    xs   = linspace(0.5, xmax, ceil(200*xmax/pi));
    v    = dJ(n, xs);
    idx  = find(v(1:end-1).*v(2:end) < 0);
    if numel(idx) < smax
        error('circ_roots:notEnoughRoots', 'Only %d zeros found for n = %d.', numel(idx), n);
    end
    for s = 1:smax
        jp(n+1, s) = fzero(@(t) dJ(n, t), [xs(idx(s)), xs(idx(s)+1)]);
    end
end
end
