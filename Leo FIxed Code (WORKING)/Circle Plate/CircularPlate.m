%% CircularPlate.m -- Chladni figures on a free-edge circular plate, center drive
% Model (Report Sec. 2.3, 2.5):
%   lap(W) + (k + i*gamma)^2 W = -delta(r),   dW/dr = 0 at r = R
% with Bessel modes J_n(k r) cos(n*theta), free edge J_n'(kR) = 0.
%
% The center screw is modeled as a point force (averaged over the screw head).
% J_n(0) = 0 for every n >= 1, so the center force only drives the n = 0 (ring)
% modes:  W(r) = sum_s a_s J_0(k_s r). Every pattern below is computed at a MEASURED
% lab frequency and labeled with the nearest ring mode.
%
% Light regions in the plot = |W| ~ 0 = where the sand collects.
% Needs circ_roots.m and closest_mode_circle.m in the same folder.
clear; clc; close all

%% 1. Plate properties and frequency -> wave number (Sec. 2.2)
E   = 69e9;  nu = 0.33;  rho = 2700;  h = 0.64e-3;   % aluminum
R   = 0.0745;                                        % plate radius [m]
D   = E*h^3/(12*(1 - nu^2));                         % bending stiffness
cf  = sqrt(D/(rho*h));                               % omega = cf*k^2

%% 2. Numerical settings (Sec. 2.7)
P.N   = 40;                 % ring modes kept (zeros of J_0' = zeros of J_1)
P.Ng  = 500;                % grid points per side
P.gam = 0.01*pi/R;          % damping gamma [1/m]; keeps resonance finite
P.rc  = 0.0026;             % screw head radius [m]

%% 3. Measured frequencies from the lab [Hz]
f_lab = [160 780 795 1785 2250 3580];

%% 4. Pattern at each measured frequency
for f = f_lab
    [Wr, r] = circle_response(f, R, cf, P);
    M = closest_mode_circle(f, cf, R, 1);            % nearest ring mode

    [X, Y] = meshgrid(linspace(-R, R, P.Ng));
    Rr   = hypot(X, Y);
    W    = interp1(r, Wr, Rr);                       % NaN outside the disk
    Wabs = abs(W)/max(abs(Wr));
    Wabs(Rr < P.rc) = NaN;                           % hide the screw region

    figure;
    imagesc([-R R], [-R R], Wabs, 'AlphaData', ~isnan(Wabs)); axis image xy;
    colormap(flipud(gray)); caxis([0 0.15]);         % saturate so nodal lines stand out
    title(sprintf('Circle, f = %.0f Hz (nearest driven mode (0,%d), %.0f Hz)', ...
          f, M.s, M.f_pred_Hz));
    xlabel('x [m]'); ylabel('y [m]');
    saveas(gcf, sprintf('circle_pattern_f%dHz.png', round(f)));
end

%% 5. Table for the report (Table 1)
T = closest_mode_circle(f_lab, cf, R, 1);
disp(T);

%% ------------------------------------------------------------------------
function [Wr, r] = circle_response(f, R, cf, P)
%CIRCLE_RESPONSE  Radial profile W(r) for a center force at frequency f [Hz].
k  = sqrt(2*pi*f/cf);                            % from omega = cf*k^2
jp = circ_roots(0, P.N);                         % zeros of J_0' (= zeros of J_1)
if jp(1, P.N)/R < 2*k
    warning('N = %d may be too small for k = %.1f 1/m. Increase N.', P.N, k);
end
kj = [0, jp(1, :)]/R;                            % ring wavenumbers, incl. rigid mode k = 0

% Finite screw: average the point force over a disk of radius rc.
xd = kj*P.rc;
Fd = ones(size(xd));  nz = xd > 0;
Fd(nz) = 2*besselj(1, xd(nz))./xd(nz);

% <delta, phi_s> = J_0(0)*Fd = Fd;  <phi_s, phi_s> = pi R^2 J_0(j'_s)^2
nrm = pi*R^2*besselj(0, kj*R).^2;
a   = Fd ./ (nrm .* (kj.^2 - (k + 1i*P.gam)^2)); % projection coefficients

r  = linspace(0, R, P.Ng);
Wr = besselj(0, r(:)*kj) * a(:);                 % W(r)
end
