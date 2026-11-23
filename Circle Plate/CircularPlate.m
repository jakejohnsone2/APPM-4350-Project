clear
close all
clc
%% chladni_circle.m -- Chladni figure on a free-edge circular plate, center drive
% Solves  lap(W) + (k + i*gamma)^2 W = -delta(r)  with  dW/dr = 0 at r = R
% using Bessel modes J_n(k r) cos(n*theta), free edge J_n'(kR) = 0 (Report Sec. 2.5).
% A center drive excites only n = 0 (rings):  W(r) = sum_s a_s J_0(k_s r).
% Light regions in the plot = |W| ~ 0 = where the sand collects.
% Needs circ_roots.m and closest_mode_circle.m in the same folder.
%% 1. Plate properties and frequency -> wave number (Sec. 2.2, Eq. 4)
E = 69e9;  nu = 0.33;  rho = 2700;  h = .64e-3;   % aluminum -- EDIT to your plate
R = 0.0745;                                        % plate radius [m] -- EDIT to your plate
D  = E*h^3/(12*(1 - nu^2));                        % bending stiffness
cf = sqrt(D/(rho*h));                              % omega = cf*k^2  (from k^4 = rho*h*omega^2/D)
% jp(n+1, s) = s-th positive zero of J_n'(x)  (n = 0..12, s = 1..40)
% Mode (n,s): n nodal diameters, s-th root.  For n = 0 these are the zeros of J_1.
jp  = circ_roots(12, 40);
jpf = @(n,s) jp(sub2ind(size(jp), n+1, s));
our_n = 7;                                         % azimuthal index (center drive: n = 0 only)
our_m = 1;                                         % radial index
% Predicted resonances of the model: k_ns = j'_ns / R
fnm = @(n,s) cf*(jpf(n,s)/R).^2/(2*pi);            % [Hz]
f     = fnm(our_n,our_m);              % <-- drive frequency [Hz]; replace with a lab value
omega = 2*pi*f;
k     = (rho*h*omega^2/D)^(1/4);
%% 2. Numerical settings (Sec. 2.7)
N    = 40;                     % radial modes kept (zeros of J_1), N <= 40
Ng   = 500;                    % grid points per side
gam  = 0.01*pi/R;              % damping gamma [1/m]; keeps resonance finite
rc   = 0.0026;                 % screw/clamp radius [m] -- measure yours
if jp(1,N)/R < 2*k
    warning('N = %d may be too small for k = %.1f 1/m. Increase N.', N, k);
end
%% 3. Series solution (n = 0 modes, finite clamp at the center)
kj = [0, jp(1,1:N)]/R;                           % radial wavenumbers incl. rigid mode k = 0
xd = kj*rc;                                      % finite clamp: disk-average of the source
Fd = ones(size(xd));  nz = xd > 0;               % mean-value factor 2*J1(x)/x for Helmholtz modes
Fd(nz) = 2*besselj(1, xd(nz))./xd(nz);
Nj = pi*R^2*besselj(0, kj*R).^2;                 % norm: int J0(kj r)^2 dA = pi R^2 J0(kj R)^2
a  = Fd ./ (Nj .* (kj.^2 - (k + 1i*gam)^2));     % coefficients
r  = linspace(0, R, Ng);
Wr = besselj(0, r(:)*kj) * a(:);                 % radial profile W(r)
%% 4. Plot |W| (light = no motion = sand)
[X, Y] = meshgrid(linspace(-R, R, Ng));
Rr = hypot(X, Y);
W  = interp1(r, Wr, Rr);                         % NaN outside the disk
Wabs = abs(W)/max(abs(Wr));
Wabs(Rr < rc) = NaN;                             % hide the screw region
figure;
imagesc([-R R], [-R R], Wabs, 'AlphaData', ~isnan(Wabs)); axis image xy;
colormap(flipud(gray)); caxis([0 0.15]);         % saturate so nodal lines stand out
title(sprintf('Circular plate (center drive), f = %.1f Hz, Mode: n = %.0f s = %.0f', f, our_n, our_m));
xlabel('x [m]'); ylabel('y [m]');
saveas(gcf, ["circle_pattern_n"+ num2str(our_n) + "_m" + num2str(our_m) + ".png"]);

% Mode shape of the selected (n,s):  J_n(k r) cos(n*theta).  Use this for n >= 1
% (a center point force cannot excite those, but a screw-pinned plate can resonate in them).
% Th = atan2(Y, X);
% Wms = besselj(our_n, jpf(our_n,our_s)*Rr/R).*cos(our_n*Th);
% Wms = abs(Wms)/max(abs(Wms(Rr <= R)));
% Wms(Rr > R | Rr < rc) = NaN;
% figure;
% imagesc([-R R], [-R R], Wms, 'AlphaData', ~isnan(Wms)); axis image xy;
% colormap(flipud(gray)); caxis([0 0.15]);
% title(sprintf('Mode shape (n = %.0f, s = %.0f), f_{pred} = %.1f Hz', our_n, our_s, f));
% xlabel('x [m]'); ylabel('y [m]');
% saveas(gcf, 'circle_mode_shape.png');

%% 5. Finding Closest Mode
f_meas = 2250;                     % <-- measured frequency [Hz] from the lab
T = closest_mode_circle(f_meas, cf, R, 5, 0, rc);