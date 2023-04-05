function [a, e, incl, w, W, nu] = cart2oe(x, y, z, vx, vy, vz, mu)
% Convert Cartesian state to Keplerian orbital elements

% Compute specific angular momentum
r = [x; y; z];
v = [vx; vy; vz];
h = cross(r, v);

% Compute magnitude of position and velocity
R = norm(r);
V = norm(v);

% Compute energy
E = V^2/2 - mu/R;

% Compute semi-major axis
a = -mu/(2*E);

% Compute eccentricity
e_vec = cross(v, h)/mu - r/R;
e = norm(e_vec);

% Compute inclination
incl = acosd(h(3)/norm(h));

% Compute right ascension of ascending node
n = cross([0;0;1], h);
W = mod(atan2d(n(2), n(1)), 360);

% Compute argument of perigee
w = mod(atan2d(dot(n, cross(e_vec, h)), dot(e_vec, n)), 360);

% Compute true anomaly
nu = mod(atan2d(dot(r, cross(v, h))/(R*V), dot(r, v)/(R*V)), 360);
end
