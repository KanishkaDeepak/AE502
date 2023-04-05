clc
clear all
% Orbital elements
a = 26600;              % semimajor axis [km]
i = 1.10654;            % inclination [rad]
e = 0.74;               % eccentricity
omega = 5;              % argument of perigee [deg]
Omega = 90;             % longitude of ascending node [deg]
M0 = 10;                % mean anomaly [deg]

% Convert angles to radians
omega = deg2rad(omega);
Omega = deg2rad(Omega);
M0 = deg2rad(M0);

% Earth parameters
mu = 3.986e5;           % gravitational parameter [km^3/s^2]
J2 = 1.0826e-3;         % Earth's J2 perturbation
R = 6378.137;           % Earth's radius [km]

% Time span
tspan = linspace(0, 100*24*3600, 1000); % 100 days

% Initial state vector
r0 = a*(1 - e^2)/(1 + e*cos(M0));
v0 = sqrt(mu/a)*(1 + e*cos(M0))/(1 - e^2)*[-sin(M0); sqrt(1 - e^2)*cos(M0); 0];
X0 = [r0*cos(M0); r0*sin(M0); 0; v0(1); v0(2); v0(3)];

% Set tolerances and maximum step size
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', 1500);

% Solve equations of motion
[t, X] = ode45(@(t,X) molniya_perturbed(t, X, mu, J2, R), tspan, X0, options);

% Extract orbital elements from state vector
[~, ~, ~, w, W, nu] = cart2oe(X(:,1), X(:,2), X(:,3), X(:,4), X(:,5), X(:,6), mu);

% Convert angles to degrees
w = rad2deg(w);
W = rad2deg(W);
nu = rad2deg(nu);

% Plot results
figure;
subplot(3,2,1);
plot(t/(24*3600), X(:,1));
ylabel('r [km]');
title('Molniya Orbit');
grid on;
subplot(3,2,2);
plot(t/(24*3600), X(:,2));
ylabel('v_r [km/s]');
title('Molniya Orbit');
grid on;
subplot(3,2,3);
plot(t/(24*3600), X(:,3));
ylabel('v_\theta [km/s]');
title('Molniya Orbit');
grid on;
subplot(3,2,4);
plot(t/(24*3600), w);
ylabel('\omega [deg]');
title('Molniya Orbit');
grid on;
subplot(3,2,5);
plot(t/(24*3600), W);
ylabel('\Omega [deg]');
title('Molniya Orbit');
grid on;
subplot(3,2,6);
plot(t/(24*3600), nu);
ylabel('\nu [deg]');
xlabel('Time [days]');
title('Molniya Orbit');
grid on;
