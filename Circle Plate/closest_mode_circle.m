function T = closest_mode_circle(f, cf, R, nShow, doPlot, rc)
%CLOSEST_MODE_CIRCLE  Nearest circular-plate model resonances (n,s) to a measured frequency.
%   T = closest_mode_circle(f, cf, R)                f [Hz], cf = sqrt(D/(rho*h)), R [m]
%   T = closest_mode_circle(f, cf, R, nShow, doPlot, rc)
%
% Model: f_ns = cf*(j'_ns/R)^2/(2*pi), where j'_ns is the s-th positive zero of
% J_n'(x) (n nodal diameters, s-th root). Candidates are sorted by percent error.
% center_node = true for n >= 1 (W = 0 at the center, as a screw-pinned plate needs).
% Plots show J_n(k r) cos(n*theta); for n >= 1 the pattern can appear rotated
% (cos and sin versions are degenerate).
%
% Needs circ_roots.m.  Example:  closest_mode_circle(3201, cf, 0.0667, 6)

if nargin < 4 || isempty(nShow),  nShow  = 5;      end
if nargin < 5 || isempty(doPlot), doPlot = true;   end
if nargin < 6 || isempty(rc),     rc     = 0.0026; end

xneed = R*sqrt(2*pi*f/cf);                 % x = j'_ns at the target frequency
xcap  = 1.5*xneed + 3;
nmax  = ceil(xcap);                        % j'_{n,1} > n
smax  = ceil(xcap/pi) + 2;
jp = circ_roots(nmax, smax);

[nn, ss] = ndgrid(0:nmax, 1:smax);
keep  = jp <= xcap;
n_all = nn(keep);  s_all = ss(keep);  x_all = jp(keep);
fpred = cf*(x_all/R).^2/(2*pi);
errq  = 100*(fpred - f)/f;
[~, ord] = sort(abs(errq));
ord = ord(1:min(nShow, numel(ord)));

T = table(n_all(ord), s_all(ord), fpred(ord), errq(ord), n_all(ord) >= 1, ...
          'VariableNames', {'n', 's', 'f_pred_Hz', 'error_pct', 'center_node'});
fprintf('Measured f = %.0f Hz   (cf = %.4f)\n', f, cf);
disp(T);

if doPlot
    nP = min(3, height(T));
    xo = x_all(ord);
    [X, Y] = meshgrid(linspace(-R, R, 400));
    Rr = hypot(X, Y);  Th = atan2(Y, X);
    figure('Position', [100 100 330*nP 360]);
    for i = 1:nP
        ni = T.n(i);
        Wm = besselj(ni, xo(i)*Rr/R).*cos(ni*Th);
        Wm = abs(Wm)/max(abs(Wm(Rr <= R)));
        Wm(Rr > R | Rr < rc) = NaN;
        subplot(1, nP, i);
        imagesc([-R R], [-R R], Wm, 'AlphaData', ~isnan(Wm)); axis image xy;
        colormap(flipud(gray)); caxis([0 0.15]);
        title(sprintf('(n=%d, s=%d), f_{pred} = %.0f Hz (%+.1f%%)', ...
              ni, T.s(i), T.f_pred_Hz(i), T.error_pct(i)));
    end
end
end
