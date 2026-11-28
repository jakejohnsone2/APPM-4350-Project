function [phi, kmn] = tri_mode(m, n, k1, X, Y)
%TRI_MODE  Lame free-edge mode (m,n) of the equilateral triangle with vertices
%   (0,0), (L,0), (L/2, sqrt(3)L/2), evaluated at points (X,Y).
%
%   k1 = 4*pi/(3L) is the smallest lattice wave number. The mode is a sum of six
%   plane waves: the lattice vector kv turned by 0, 120 and 240 degrees, plus the
%   same three mirrored across the x axis. Each has |k| = kmn, so phi solves the
%   Helmholtz equation, and the symmetric sum has zero normal slope on all three
%   edges.  kmn^2 = k1^2 (m^2 + m*n + n^2).

kv  = k1*[m + n/2, n*sqrt(3)/2];                 % lattice wave vector
kmn = norm(kv);
phi = zeros(size(X));
for a = [0, 2*pi/3, 4*pi/3]
    for s = [1, -1]
        kx  = cos(a)*kv(1) - sin(a)*s*kv(2);
        ky  = sin(a)*kv(1) + cos(a)*s*kv(2);
        phi = phi + cos(kx*X + ky*Y);
    end
end
end
