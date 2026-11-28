function T = closest_mode_triangle(f_meas, L, cf, x0, y0, nShow)
%CLOSEST_MODE_TRIANGLE  Nearest triangle-plate resonances that a force at (x0,y0) can drive.
%   T = closest_mode_triangle(f_meas, L, cf, x0, y0)          nearest mode for each f
%   T = closest_mode_triangle(f_meas, L, cf, x0, y0, nShow)   nShow candidates per f
%
%   f_meas : measured frequency or vector of frequencies [Hz]
%   L      : side length [m]
%   cf     : sqrt(D/(rho*h)), so f = cf*k^2/(2*pi)
%   x0, y0 : drive point [m] (centroid: L/2, sqrt(3)L/6)
%
% Model: f_mn = cf*k1^2*(m^2 + m*n + n^2)/(2*pi),  k1 = 4*pi/(3L).
% A mode is driven only if phi_mn(x0,y0) ~= 0 (the projection coefficient is
% proportional to it). At the centroid this keeps modes with m - n a multiple of 3.
% Listed with m >= n; the rigid mode (0,0) is excluded. Needs tri_mode.m.

if nargin < 6 || isempty(nShow), nShow = 1; end

k1 = 4*pi/(3*L);
rows = [];
for f = f_meas(:).'
    Nmax = ceil(2*sqrt(2*pi*f/cf)/k1) + 2;      % search well past the drive k
    cand = [];
    for m = 0:Nmax
        for n = 0:m
            if m + n == 0, continue; end
            if abs(tri_mode(m, n, k1, x0, y0)) < 1e-6, continue; end   % not driven
            cand = [cand; m, n, cf*k1^2*(m^2 + m*n + n^2)/(2*pi)];   %#ok<AGROW>
        end
    end
    err = 100*(cand(:, 3) - f)/f;
    [~, idx] = sort(abs(err));
    idx = idx(1:min(nShow, numel(idx)));

    rows = [rows; repmat(f, numel(idx), 1), cand(idx, :), err(idx)]; %#ok<AGROW>
end

T = array2table(rows, 'VariableNames', {'f_meas_Hz', 'm', 'n', 'f_pred_Hz', 'error_pct'});
end
