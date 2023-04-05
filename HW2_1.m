% Constants
GM = 3.986e5; % km^3/s^2
J2 = 0.00108;
R = 6370; % km
T = 12 * 60 * 60; % 12 hours in seconds
rp = 6970; % km

% Calculate semimajor axis
a = (GM * (T / (2 * pi))^2)^(1/3);
disp(['Semimajor axis: ', num2str(round(a, 2)), ' km']);

% Calculate eccentricity
e = 1 - rp / a;
disp(['Eccentricity: ', num2str(round(e, 3))]);

% Calculate inclination
i = pi; % 180 degrees in radians
disp(['Inclination: ', num2str(round(rad2deg(i), 2)), ' degrees']);

% Calculate nodal precession rate
n = 2 * pi / T;
dOdt = -1.5 * n * J2 * (R / a)^2 / (1 - e^2)^2;
disp(['Nodal precession rate: ', num2str(round(rad2deg(dOdt), 3)), ' degrees/day']);
