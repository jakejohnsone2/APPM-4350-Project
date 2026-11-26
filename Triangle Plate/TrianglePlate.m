clear
close all
clc
%% chladni_triangle.m

%% 1. Plate properties and frequency to wave number
E = 69e9;  nu = 0.33;  rho = 2700;  h = .52e-3;    % aluminum
H = 0.1268;                                        % triangle height [m]
L = 2*H/sqrt(3);                                   % side length [m]
D  = E*h^3/(12*(1 - nu^2));                        % bending stiffness
cf = sqrt(D/(rho*h));                              % omega = cf*k^2  (from k^4 = rho*h*omega^2/D)
our_m = 3;
our_n = 0;
% Predicted resonances of the model: k_mn^2 = (4*pi/(3L))^2 (m^2 + m*n + n^2)
k1  = 4*pi/(3*L);                                  % smallest wave number of the lattice
fmn = @(m,n) cf*k1^2*(m.^2 + m.*n + n.^2)/(2*pi);  % [Hz]
f   = fmn(our_m,our_n);                % drive frequency [Hz]

%% Series solution and plots
chladni_triangle(f, L, cf);

%% Save
for f_lab = [230 600 1950 2250 3390 4835], chladni_triangle(f_lab, L, cf); end


%% ------------------------------------------------------------------------
function chladni_triangle(f, L, cf)

%% Numerical settings
H    = sqrt(3)/2*L;            % height [m]
x0   = L/2;  y0 = H/3;         % center drive (centroid of the triangle)
Ng   = 1000;                   % grid points per side
gam  = 0.01*pi/L;              % damping gamma [1/m]; keeps resonance finite
rc   = 0.0026;                 % screw/clamp radius [m] -- measure yours
k    = sqrt(2*pi*f/cf);        % wave number of the drive
k1   = 4*pi/(3*L);             % smallest wave number of the lattice
N    = ceil(2.5*k/k1) + 5;     % modes kept (m,n = 0..N); enough to pass 2k

%% 3. Series solution: sum over Lame modes (m >= n >= 0)
x = linspace(x0 - 1.25*L, x0 + 1.25*L, Ng);      % window holds the plate + its mirror images
y = linspace(y0 - 1.25*L, y0 + 1.25*L, Ng);
[X, Y] = meshgrid(x, y);
inside = (Y >= 0) & (Y <= sqrt(3)*X) & (Y <= sqrt(3)*(L - X));   % points on the plate
dA = (x(2) - x(1))*(y(2) - y(1));                % area of one grid cell

W = zeros(size(X));
fbest = inf;
for m = 0:N
    for n = 0:m
        kv  = k1*[m + n/2, n*sqrt(3)/2];         % lattice wave vector
        kmn = norm(kv);
        % Mode = sum of 6 plane waves: kv turned by 0, 120, 240 degrees,
        % plus the same three for kv mirrored across the x axis
        phi  = zeros(size(X));
        phi0 = 0;                                % value of the mode at the screw
        for a = [0, 2*pi/3, 4*pi/3]
            for s = [1, -1]
                kx = cos(a)*kv(1) - sin(a)*s*kv(2);
                ky = sin(a)*kv(1) + cos(a)*s*kv(2);
                phi  = phi  + cos(kx*X  + ky*Y);
                phi0 = phi0 + cos(kx*x0 + ky*y0);
            end
        end
        Fd = 1;                                  % finite clamp: disk-average of the source
        if kmn > 0
            Fd = 2*besselj(1, kmn*rc)/(kmn*rc);
        end
        nrm = sum(phi(inside).^2)*dA;            % int phi^2 dA over the plate
        A = Fd*phi0/(nrm*(kmn^2 - (k + 1i*gam)^2));
        W = W + A*phi;
        % remember the closest mode the center screw can drive (for the title)
        fm = cf*kmn^2/(2*pi);
        if abs(phi0) > 1e-6 && abs(fm - f) < abs(fbest - f)
            fbest = fm;  mbest = m;  nbest = n;
        end
    end
end

%% 4. Plot |W| (light = no motion = sand)
Wabs = abs(W)/max(abs(W(inside)));
Wabs(hypot(X - x0, Y - y0) < rc) = NaN;          % hide the screw region
xt = [0, L, L/2, 0];  yt = [0, 0, H, 0];         % triangle outline

% Figure 1: full mirrored pattern
figure;
imagesc(x, y, Wabs, 'AlphaData', ~isnan(Wabs)); axis image xy;
colormap(flipud(gray)); caxis([0 0.15]);         % saturate so nodal lines stand out
hold on; plot(xt, yt, 'r-', 'LineWidth', 1.5); hold off;   % the real plate
title(sprintf('Triangle plate (mirrored), f = %.1f Hz, Mode: m = %.0f n = %.0f', f, mbest, nbest));
xlabel('x [m]'); ylabel('y [m]');
saveas(gcf, "triangle_full_pattern_f" + num2str(round(f)) + "Hz.png");

% Figure 2: just the single triangle
Wtri = Wabs;
Wtri(~inside) = NaN;                             % hide everything off the plate
figure;
imagesc(x, y, Wtri, 'AlphaData', ~isnan(Wtri)); axis image xy;
colormap(flipud(gray)); caxis([0 0.15]);
xlim([0 L]); ylim([0 H]);
title(sprintf('Triangle plate, f = %.1f Hz, Mode: m = %.0f n = %.0f', f, mbest, nbest));
xlabel('x [m]'); ylabel('y [m]');
saveas(gcf, "triangle_pattern_f" + num2str(round(f)) + "Hz.png");
end
