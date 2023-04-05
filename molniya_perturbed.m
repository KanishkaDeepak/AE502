function dXdt = molniya_perturbed(~, X, mu, J2, ~)
% MOLNIYA_PERTURBED Calculates the derivatives of the state vector for a
% perturbed Molniya orbit.
%
% Inputs:
%   - t     : current time [s]
%   - X     : current state vector [km] and [km/s]
%   - mu    : gravitational parameter of the central body [km^3/s^2]
%   - J2    : J2 perturbation coefficient of the central body
%   - R     : radius of the central body [km]
%
% Outputs:
%   - dXdt  : derivative of the state vector [km/s] and [km/s^2]
%
% Author: Siddharth Kherada (https://github.com/siddharthkherada)

% Extract position and velocity from state vector
r = X(1:3);
v = X(4:6);

% Calculate acceleration due to J2 perturbation
J2_factor = (3/2)*J2*(mu/norm(r))^2;
a_J2 = J2_factor * [(r(1)/norm(r))*(5*(r(3)/norm(r))^2 - 1);
                    (r(2)/norm(r))*(5*(r(3)/norm(r))^2 - 1);
                    (r(3)/norm(r))*(5*(r(3)/norm(r))^2 - 3)];

% Calculate acceleration due to central body gravity
a_gravity = -mu/norm(r)^3*r;

% Calculate acceleration due to third body gravity
a_third_body = [0; 0; 0]; % Assumes no third body perturbations

% Calculate total acceleration
a_total = a_gravity + a_J2 + a_third_body;

% Derivative of state vector
dXdt = [v; a_total];

end
