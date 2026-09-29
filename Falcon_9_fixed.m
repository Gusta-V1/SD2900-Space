clear; clc; close all;
%% 1. ENVIRONMENT AND LAUNCH PARAMETERS
p.g0 = 9.80665;                          % Gravitational acceleration [m/s^2]
p.my = 3.986004418 * 10^14;              % Earth gravitational parameter [m^3/s^2]
p.wE = 7.292115 * 10^-5;                 % Earth rotation rate [rad/s]
p.R_E = 6378137;                         % Earth radius [m]
p.P_SL = 101325;                         % Sea level pressure [Pa]

Lat = deg2rad(28.5623);                  % Latitude (Radians)
Long = deg2rad(-80.5774);                % Longitude (Negative for West, Radians)

i = deg2rad(74);                         % Target inclination (Radians)
p.beta = asin(cos(i)/cos(Lat));          % Launch azimuth for the launch-window calc (Radians)

% Inertial thrust azimuth used by the stage-1 pitch program [rad].
% Earth's rotation gives ~408 m/s of eastward velocity, so the thrust azimuth that
% ends up at 74 deg inclination is smaller than p.beta. 7.5 deg was tuned to give
% i ~ 73.9 deg with the gains below. Retune if you change payload or pitch program.
p.beta_thrust = deg2rad(7.5);

%% --- LAUNCH TIME CALCULATION (Two Daily Windows) ---
% Define targeting parameters
target_RAAN_deg = 45.0;            % Target Right Ascension of Ascending Node [deg]
launch_date_str = '2026-10-15';    % Planned launch date (YYYY-MM-DD)

% 1. Calculate Right Ascension Offset (Delta Alpha)
sin_delta_alpha = tan(Lat) / tan(i);
delta_alpha = asin(sin_delta_alpha); % [radians]

% 2. Calculate Required Local Sidereal Time (LST) for BOTH windows
% Ascending pass (Northeast launch)
theta_LST_asc = deg2rad(target_RAAN_deg) + delta_alpha;
% Descending pass (Southeast launch)
theta_LST_desc = deg2rad(target_RAAN_deg) + pi - delta_alpha;

% 3. Calculate Required Greenwich Sidereal Time (GST)
theta_GST_asc_deg = mod(rad2deg(theta_LST_asc - Long), 360);
theta_GST_desc_deg = mod(rad2deg(theta_LST_desc - Long), 360);

% 4. Calculate GMST at 00:00 UTC on the launch date
d_0h = datetime(launch_date_str, 'TimeZone', 'UTC');
jd_0h = juliandate(d_0h);
T_cen = (jd_0h - 2451545.0) / 36525; % Julian centuries since J2000.0

% IAU-83 formula for GMST at 00:00 UTC [degrees]
gmst_0h_deg = mod(100.4606184 + 36000.7700536 * T_cen + ...
                  0.000387933 * T_cen^2 - 2.583e-8 * T_cen^3, 360);

% 5. Calculate Time of Flight (UTC hours past midnight)
hours_utc_asc = mod(theta_GST_asc_deg - gmst_0h_deg, 360) / 15.0410686386;
hours_utc_desc = mod(theta_GST_desc_deg - gmst_0h_deg, 360) / 15.0410686386;

% 6. Final Launch Times
launch_time_asc = d_0h + hours(hours_utc_asc);
launch_time_desc = d_0h + hours(hours_utc_desc);

% 7. Calculate Launch Azimuths (Beta)
beta_asc_deg = rad2deg(p.beta);
beta_desc_deg = 180 - beta_asc_deg;

% Display the launch windows in the command window
fprintf('=== LAUNCH WINDOW TARGETING ===\n');
fprintf('Target RAAN: %.2f deg | Target Inclination: %.2f deg\n\n', target_RAAN_deg, rad2deg(i));

fprintf('--- WINDOW 1 (Ascending Pass) ---\n');
fprintf('Launch Azimuth (Beta): %.2f deg (Northeast)\n', beta_asc_deg);
fprintf('Launch Time (UTC): %s\n\n', datestr(launch_time_asc, 'yyyy-mm-dd HH:MM:SS'));

fprintf('--- WINDOW 2 (Descending Pass) ---\n');
fprintf('Launch Azimuth (Beta): %.2f deg (Southeast)\n', beta_desc_deg);
fprintf('Launch Time (UTC): %s\n\n', datestr(launch_time_desc, 'yyyy-mm-dd HH:MM:SS'));
%% 2. ROCKET PARAMETERS (Falcon 9 Block 5)
m_pay = 15000;                           % Payload [kg]
fairing = 1900;                          % Fairing [kg]
diam_falcon = 3.7;                       % Diameter [m]
p.A = pi * (diam_falcon / 2)^2;          % Cross-sectional area [m^2]
p.C_D = 0.3;                             % Drag coefficient (simplified)

m_struc_1 = 28100 + 107500 + 25000 + m_pay + fairing; 
p.m_prop_1 = 386000;                     % Propellant mass for Stage 1 [kg]
p.m0 = m_struc_1 + p.m_prop_1;           % Gross liftoff mass [kg]

p.Isp_SL = 283;                          % Isp sea level [s]
p.Isp_vac = 312;                         % Isp vacuum [s]
p.T_vac = 8227000;                       % Vacuum thrust for 9x Merlin [N]

% Stage 1 burn time = propellant / mass flow (~143.6 s). Stage 1 is cut off by an
% event when the propellant is gone, NOT at a hard-coded 164 s.
p.t_burn1 = p.m_prop_1 / (p.T_vac / (p.g0 * p.Isp_vac));

% Stage 2
p.Isp_vac_2 = 348;                       % Vacuum Isp Stage 2 (Merlin Vac)
p.T_vac_2 = 981000;                      % Vacuum thrust Stage 2 [N]
m_struc2 = 4000;                         % Stage 2 dry mass [kg]
m_prop2 = 107500;                        % Stage 2 propellant [kg]
p.m_dry2 = m_struc2 + m_pay;             % Stage 2 dry + payload (engine cutoff mass) [kg]
p.m02 = m_prop2 + m_struc2 + fairing + m_pay;

% Target orbit (circular). 200 km fits in one burn with a 15 t payload.
% 719 km does NOT fit: stage 2 runs out of propellant (needs less payload or a 2nd burn).
p.r_target = p.R_E + 200e3;              % Target orbit radius [m]

% Stage 2 linear-tangent steering: elevation = atan(a - b*(t - t_s2)) above local horizon.
% Tuned for 200 km / 15 t payload. Retune if you change payload or target altitude.
p.lt_a = 0.755956;
p.lt_b = 0.00211479;

%% 3. INITIAL CONDITIONS (t = 0)
x0 = p.R_E * cos(Lat) * cos(Long);
y0 = p.R_E * cos(Lat) * sin(Long);
z0 = p.R_E * sin(Lat);
vx0 = -p.wE * y0;
vy0 = p.wE * x0;
vz0 = 0;

Y0 = [x0; y0; z0; vx0; vy0; vz0; p.m0];

%% 4. RUN SIMULATION
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-6);

% --- Stage 1: runs until propellant is depleted (event) ---
options_s1 = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'Events', @(t, Y) stage1_event(t, Y, p));
[T1, Y1] = ode45(@(t, Y) rocket_dynamics(t, Y, p), [0 300], Y0, options_s1);

m_struc_1_total = p.m0 - p.m_prop_1;
final_mass = Y1(end, 7);
remaining_propellant = final_mass - m_struc_1_total;
remaining_percentage = (remaining_propellant / p.m_prop_1) * 100;

fprintf('--- STEG 1 MECO SAMMANSTÄLLNING ---\n');
fprintf('MECO vid t = %.1f s\n', T1(end));
fprintf('Höjd vid MECO: %.1f km | Hastighet (ECI): %.0f m/s\n', ...
        (norm(Y1(end,1:3)) - p.R_E)/1000, norm(Y1(end,4:6)));
fprintf('Kvarvarande bränsle: %.1f kg\n', remaining_propellant);
fprintf('Kvarvarande proportion: %.1f%% av ursprungligt bränsle\n', remaining_percentage);

%% ----------------------------- Second Stage -----------------------------
p.t_s2 = T1(end);                        % Stage 2 ignition time (used by steering)

Y0_2 = Y1(end, :)';
Y0_2(7) = p.m02;                         % Stage 1 separated

% Phase 2: with fairing (24 s)
tspan2 = [T1(end), T1(end) + 24];
[T2, Y2] = ode45(@(t, Y) stage2_dynamics(t, Y, p), tspan2, Y0_2, options);

%% --------------------- Second Stage W/O fairing -------------------------
Y0_3 = Y2(end, :)';
Y0_3(7) = Y0_3(7) - fairing;

tspan3 = [T2(end), T2(end) + 800];       % Event function stops it early

% SECO when perigee reaches the target radius (or propellant is depleted)
options_meco = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'Events', @(t, Y) meco_event(t, Y, p));
[T3, Y3] = ode45(@(t, Y) stage2_dynamics(t, Y, p), tspan3, Y0_3, options_meco);

% Orbit at SECO
r_s = Y3(end,1:3)'; v_s = Y3(end,4:6)';
rn_s = norm(r_s);
a_s = -p.my / (2*(dot(v_s,v_s)/2 - p.my/rn_s));
e_s = norm(cross(v_s, cross(r_s,v_s))/p.my - r_s/rn_s);
h_s = cross(r_s, v_s);
inc_s = acosd(h_s(3)/norm(h_s));

fprintf('\n--- STEG 2 SECO ---\n');
fprintf('SECO inträffade vid t = %.1f s\n', T3(end));
fprintf('Perigeum: %.1f km | Apogeum: %.1f km | e = %.4f\n', ...
        (a_s*(1-e_s) - p.R_E)/1000, (a_s*(1+e_s) - p.R_E)/1000, e_s);
fprintf('Inklination: %.2f deg (mål %.1f)\n', inc_s, rad2deg(i));
fprintf('Kvarvarande bränsle steg 2: %.0f kg\n', Y3(end,7) - p.m_dry2);

%% --------------------- Coasting i Rymden (Fas 4) -------------------------
Y0_4 = Y3(end, :)';
tspan4 = [T3(end), T3(end) + 6000];      % ~ one orbit
[T4, Y4] = ode45(@(t, Y) coast_dynamics(t, Y, p), tspan4, Y0_4, options);

%% 5. CALCULATE RESULTS FOR PLOTTING
T_out = [T1; T2; T3; T4];
Y_out = [Y1; Y2; Y3; Y4];

r_out = Y_out(:, 1:3);
v_out = Y_out(:, 4:6);
m_out = Y_out(:, 7);

altitudes = vecnorm(r_out, 2, 2) - p.R_E; % Altitude in meters
velocities = vecnorm(v_out, 2, 2);        % Absolute velocity in m/s

fprintf('\nMin altitude during coast: %.2f km | Max: %.2f km\n', ...
        min(altitudes(numel([T1;T2;T3])+1:end))/1000, max(altitudes(numel([T1;T2;T3])+1:end))/1000);

%% 6. PLOTS
figure('Name', 'Falcon 9 - First Stage Trajectory', 'Position', [100, 100, 900, 600]);

% Plot 1: Altitude over time
subplot(2, 2, 1);
plot(T_out, altitudes / 1000, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Altitude (km)');
title('Altitude Profile'); grid on;

% Plot 2: Velocity over time
subplot(2, 2, 2);
plot(T_out, velocities / 1000, 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Velocity (km/s)');
title('Velocity Profile (ECI)'); grid on;

% Plot 3: Mass over time
subplot(2, 2, 3);
plot(T_out, m_out / 1000, 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Total Mass (tons)');
title('Mass Depletion'); grid on;
subplot(2, 2, 4);

start_vec = r_out(1, :)'; % Transponera till stående kolumnvektor
downrange = zeros(length(T_out), 1);

for k = 1:length(T_out)
    t_current = T_out(k);
    theta = p.wE * t_current; % Vinkeln jorden har hunnit rotera
    
    % Rotationsmatris från ECI till markfixerat (ECEF) system
    Rz = [cos(theta),  sin(theta), 0;
         -sin(theta),  cos(theta), 0;
          0,           0,          1];
          
    % Raketens position i förhållande till startplattan (snurrar med jorden)
    r_ecef = Rz * r_out(k, :)';
    
    % Beräkna vinkeln längs med jordytan
    cos_angle = dot(start_vec, r_ecef) / (norm(start_vec) * norm(r_ecef));
    
    % Säkerhetsspärr för numeriska fel i MATLAB (så att cos_angle inte blir > 1)
    cos_angle = min(1, max(-1, cos_angle)); 
    
    angle = acos(cos_angle);
    downrange(k) = p.R_E * angle; 
end

plot(downrange / 1000, altitudes / 1000, 'g', 'LineWidth', 2);
xlabel('Ground distance (km)'); ylabel('Altitude (km)');
title('Flight path rellative to ground'); grid on;
%{
% Plot 4: 2D Trajectory Profile
subplot(2, 2, 4);
% Calculate downrange distance (Simplified via angle difference)
start_vec = r_out(1, :);
downrange = zeros(length(T_out), 1);
for k = 1:length(T_out)
    angle = acos(dot(start_vec, r_out(k, :)) / (norm(start_vec) * norm(r_out(k, :))));
    downrange(k) = p.R_E * angle; 
end
plot(downrange / 1000, altitudes / 1000, 'g', 'LineWidth', 2);
xlabel('Downrange Distance (km)'); ylabel('Altitude (km)');
title('Flight Trajectory'); grid on;
%}

%% 8. DYNAMISKT TRYCK OCH MAX-Q
q_array = zeros(length(T_out), 1);

for k = 1:length(T_out)
    % Hämta data för aktuellt tidssteg
    r_k = r_out(k, :)';
    v_k = v_out(k, :)';
    h_k = altitudes(k);
    
    % Hämta luftdensitet (rho)
    if h_k < 84000
        [~, ~, ~, rho_k] = atmosisa(h_k);
    else
        rho_k = 0;
    end

    % Beräkna relativ hastighet mot vinden
    omega_vec = [0; 0; p.wE];
    v_wind = cross(omega_vec, r_k);
    v_rel = v_k - v_wind;
    v_rel_mag = norm(v_rel);
    
    % Beräkna dynamiskt tryck (q) i Pascal
    q_array(k) = 0.5 * rho_k * v_rel_mag^2;
end

% Hitta värdet och tidpunkten för Max-Q
[max_q, max_idx] = max(q_array);
t_max_q = T_out(max_idx);
h_max_q = altitudes(max_idx);

% Skapa en plot för Dynamiskt tryck
figure('Name', 'Dynamic Pressure (Max-Q)', 'Position', [200, 200, 800, 400]);
plot(T_out, q_array / 1000, 'm', 'LineWidth', 2);
hold on;

% Markera Max-Q med en röd prick
plot(t_max_q, max_q / 1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Begränsa x-axeln till de första 90 sekunderna
xlim([0 90]);

% Formatera grafen
xlabel('Tid (s)');
ylabel('Dynamiskt tryck (kPa)');
title(sprintf('Max-Q: %.1f kPa vid t = %.1f s (Höjd: %.1f km)', max_q/1000, t_max_q, h_max_q/1000));
grid on;

%% 9. 3D OMLOPPSBANA (ECI)
figure('Name', '3D Trajectory (ECI)', 'Position', [250, 150, 800, 800]);

% Skapa en 3D-sfär som representerar jorden
R_E_km = p.R_E / 1000;
[X_sphere, Y_sphere, Z_sphere] = sphere(50);

% Plotta jorden (med lite transparens så man ser banan på andra sidan)
surf(X_sphere * R_E_km, Y_sphere * R_E_km, Z_sphere * R_E_km, ...
    'EdgeColor', 'none', 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.6);
hold on;

% Plotta raketens hela färdväg (röd linje)
x_km = r_out(:, 1) / 1000;
y_km = r_out(:, 2) / 1000;
z_km = r_out(:, 3) / 1000;
plot3(x_km, y_km, z_km, 'r-', 'LineWidth', 2.5);

% Markera startpunkten (grön) och slutpunkten (svart)
plot3(x_km(1), y_km(1), z_km(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot3(x_km(end), y_km(end), z_km(end), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% Inställningar för grafens utseende
axis equal; % Extremt viktigt, annars ser jorden ut som ett ägg
grid on;
xlabel('X-axel ECI (km)');
ylabel('Y-axel ECI (km)');
zlabel('Z-axel ECI (km)');
title('Falcon 9 - Omloppsbana i rymden (ECI)');
legend('Jorden', 'Raketens bana', 'Start (Florida)', 'Slutposition', 'Location', 'best');
view(3); % Sätter vyn i standard 3D
hold off;


%% -----------------------------------------------------------------------
% LOCAL FUNCTIONS (must be at the end of a script file)
%% -----------------------------------------------------------------------

% ---------------------------- STAGE 1 ------------------------------------
function dY = rocket_dynamics(t, Y, p)
    r = Y(1:3);
    v = Y(4:6);
    m = Y(7);
    h = norm(r) - p.R_E;

    % Atmosphere (atmosisa is only valid below ~86 km)
    if h < 84000
        [~, ~, Patm, rho] = atmosisa(h);
    else
        Patm = 0;
        rho = 0;
    end

    % Relative velocity (air rotates with Earth)
    omega_vec = [0; 0; p.wE];
    v_wind = cross(omega_vec, r);
    v_rel = v - v_wind;
    v_rel_mag = norm(v_rel);

    % Dynamic pressure and Max-Q throttling
    % NOTE: q_limit = 355400 Pa is far above real Max-Q (~30-35 kPa), so this
    % never triggers. If you lower it, retune p.t_burn1 / pitch / p.lt_a, p.lt_b.
    q = 0.5 * rho * v_rel_mag^2;
    q_limit = 355400;
    q_warning = q_limit * 0.8;
    throttle_q = 1.0;
    if q > q_warning && h < 84000
        throttle_q = 1.0 - 0.2 * ((q - q_warning) / (q_limit - q_warning));
        throttle_q = max(throttle_q, 0.60);
    end

    % Engine performance
    m_flow_rate = p.T_vac / (p.g0 * p.Isp_vac) * throttle_q;
    Isp = p.Isp_vac - (p.Isp_vac - p.Isp_SL) * (Patm / p.P_SL);
    T_current = m_flow_rate * p.g0 * Isp;

    % Local unit vectors
    Z_loc = r / norm(r);
    k_pole = [0; 0; 1];
    E_loc = cross(k_pole, Z_loc);
    E_loc = E_loc / norm(E_loc);
    N_loc = cross(Z_loc, E_loc);
    H_dir = cos(p.beta_thrust) * N_loc + sin(p.beta_thrust) * E_loc;

    % Steering: vertical for 10 s, then pitch over as f^1.5 to 85 deg from vertical
    if t < 10
        u = Z_loc;
    else
        f = min((t - 10) / (p.t_burn1 - 10), 1);
        pitch_angle = deg2rad(85 * f^1.5);
        u = cos(pitch_angle) * Z_loc + sin(pitch_angle) * H_dir;
    end

    % Forces
    Fg = -(p.my * m / norm(r)^3) * r;
    if v_rel_mag > 0
        FD = -0.5 * rho * v_rel_mag^2 * p.A * p.C_D * (v_rel / v_rel_mag);
    else
        FD = [0; 0; 0];
    end

    % Engine cannot burn structure
    if m <= (p.m0 - p.m_prop_1)
        T_current = 0;
        m_flow_rate = 0;
    end

    FT = T_current * u;

    dY = [v; (Fg + FD + FT) / m; -m_flow_rate];
end

% Stops stage 1 when the propellant is used up
function [value, isterminal, direction] = stage1_event(t, Y, p)
    value = Y(7) - (p.m0 - p.m_prop_1) - 1;   % 1 kg margin
    isterminal = 1;
    direction = -1;
end

% ---------------------------- STAGE 2 ------------------------------------
function dY = stage2_dynamics(t, Y, p)
    r = Y(1:3);
    v = Y(4:6);
    m = Y(7);

    m_flow_rate = p.T_vac_2 / (p.g0 * p.Isp_vac_2);
    T_current = p.T_vac_2;

    % Engine off when propellant is gone
    if m <= p.m_dry2
        T_current = 0;
        m_flow_rate = 0;
    end

    % Steering: linear tangent law, elevation above the local horizontal in the
    % current orbital plane
    r_hat = r / norm(r);
    h_vec = cross(r, v);
    if norm(h_vec) > 0
        horiz_dir = cross(h_vec, r);
        horiz_dir = horiz_dir / norm(horiz_dir);
        theta = atan(p.lt_a - p.lt_b * (t - p.t_s2));
        u = cos(theta) * horiz_dir + sin(theta) * r_hat;
    else
        u = v / norm(v);
    end

    Fg = -(p.my * m / norm(r)^3) * r;
    FT = T_current * u;

    dY = [v; (Fg + FT) / m; -m_flow_rate];
end

% SECO event: (1) perigee reaches target radius, (2) stage 2 propellant depleted
function [value, isterminal, direction] = meco_event(t, Y, p)
    r = Y(1:3);
    v = Y(4:6);
    rn = norm(r);

    energy = dot(v, v) / 2 - p.my / rn;
    a = -p.my / (2 * energy);
    e = norm(cross(v, cross(r, v)) / p.my - r / rn);
    r_perigee = a * (1 - e);

    value = [r_perigee - p.r_target;  Y(7) - p.m_dry2 - 1];
    isterminal = [1; 1];
    direction = [1; -1];
end

%% -----------------------------------------------------------------------
% DIFFERENTIAL EQUATION FUNCTION - PHASE 4 (COASTING)
%% -----------------------------------------------------------------------
function dY = coast_dynamics(t, Y, p)
    % Ren gravitation, ingen motor, inget luftmotstånd
    r = Y(1:3);
    v = Y(4:6);
    m = Y(7);
    
    Fg = -(p.my * m / norm(r)^3) * r;
    
    r_dot = v;
    v_dot = Fg / m;
    m_dot = 0; % Massan förändras inte
        
    dY = [r_dot; v_dot; m_dot];
end