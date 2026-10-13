clear; clc; close all;

%% 1. ENVIRONMENT AND LAUNCH PARAMETERS
p.g0 = 9.80665;                          % Gravitational acceleration [m/s^2]
p.my = 3.986004418 * 10^14;              % Earth gravitational parameter [m^3/s^2]
p.wE = 7.292115 * 10^-5;                 % Earth rotation rate [rad/s]
p.R_E = 6378137;                         % Earth radius [m]
p.P_SL = 101325;                         % Sea level pressure [Pa]

Lat = deg2rad(28.5623);                  % Latitude (Cape Canaveral)
Long = deg2rad(-80.5774);                % Longitude
i = deg2rad(82.95426);                   % Target inclination (Radians)
p.beta = asin(cos(i)/cos(Lat));          % INERTIAL azimuth of the target orbit plane (Radians)
p.beta_steer = p.beta;                   % CHANGED: thrust azimuth actually flown (solved in section 4)

%% 2. ROCKET PARAMETERS, Falcon 9
m_pay = 15000;                           % Payload [kg]
fairing = 1900;                          % Fairing [kg]
diam_falcon = 3.7;                       % Diameter [m]
p.A = pi * (diam_falcon / 2)^2;          % Cross-sectional area [m^2]
p.C_D = 0.3;                             % Drag coefficient, stage 1
p.A_2 = p.A;                             % CHANGED: Stage 2 reference area [m^2] (same 3.7 m diameter)
p.C_D_2 = 0.3;                           % CHANGED: Stage 2 drag coefficient (fairing off, same value assumed)

% Mass configuration
p.m_prop_1 = 395700;                     % Stage 1 Propellant [kg]
m_struc_1 = 28100;                       % Dry mass [kg]
m_prop2 = 107500;                        % Propellant mass [kg]
m_prop_return = 25000;                   % Propellant mass for return [kg]

m_struc_1 = m_struc_1 + m_prop2 + m_prop_return + m_pay + fairing; 
p.m0 = m_struc_1 + p.m_prop_1;           % Gross liftoff mass [kg]

m_struc2 = 4000;                         % Stage 2 dry mass [kg]
p.m_dry2 = m_struc2 + m_pay;             % Stage 2 mass at burnout (dry + payload) [kg]

% Engine performance
p.Isp_SL = 285;                          % Stage 1 Isp sea level [s]
p.Isp_vac = 313;                         % Stage 1 Isp vacuum [s]
p.T_vac = 8227000;                       % Stage 1 Vacuum thrust [N]
p.Isp_vac_2 = 348;                       % Stage 2 Vacuum Isp [s]
p.T_vac_2 = 981000;                      % Stage 2 Thrust [N]

% Target Orbit
p.target_radius = 7351.8*10^3;               % Target orbital radius (a0) [m]
target_RAAN_deg = 18.80579;                  % Target Right Ascension of Ascending Node [deg]
launch_date_str = '2026-10-15';              % Planned launch date (YYYY-MM-DD)

% ---- CHANGED: Stage 1 pitch program (angle from local vertical) ----
% pitch(t) = pitch_final * ((t - t_kick)/(t_meco_est - t_kick))^pitch_exp
p.t_kick      = 10;                      % Vertical rise time [s]
p.t_meco_est  = 147;                     % Estimated MECO time [s] (pitch program time base)
p.pitch_final = 80;                      % Pitch from vertical at MECO [deg]   (was 60)
p.pitch_exp   = 1.0;                     % Pitch curve exponent [-]            (was 2, i.e. quadratic)

% ---- CHANGED: Stage 2 transfer-burn guidance (altitude/vertical-speed feedback) ----
p.h_ref  = 150e3;                        % Altitude the transfer orbit perigee is steered towards [m]
p.Kh     = 0.01;                         % Commanded vertical speed per metre of altitude error [1/s]
p.Kv     = 0.02;                         % Thrust elevation per m/s vertical-speed error [rad/(m/s)]
p.th_max = deg2rad(70);                  % Max thrust elevation above local horizontal [rad]

%% 2.5 Satellite and debris parameters 

% Satellite parameters
% ISP
% mass
% Thrust

% Debris parameters
% Target RAAN
% Target inclination
% Target Alt
% mass

%% 3. INITIAL CONDITIONS (t = 0)
x0 = p.R_E * cos(Lat) * cos(Long);
y0 = p.R_E * cos(Lat) * sin(Long);
z0 = p.R_E * sin(Lat);
vx0 = -p.wE * y0;
vy0 = p.wE * x0;
vz0 = 0;

Y0 = [x0; y0; z0; vx0; vy0; vz0; p.m0];
options_standard = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

%% 4. AZIMUTH TARGETING (CHANGED)
inc_target_deg = rad2deg(i);
f_inc = @(b) meco_plane(b, p, Y0) - inc_target_deg;
p.beta_steer = fzero(f_inc, [deg2rad(-10), p.beta], optimset('TolX', 1e-9));

%% 5. SIMULATION PHASES

% ================= PHASE 1: STAGE 1 (ASCENT) =================
fprintf('--- PHASE 1: STAGE 1 ASCENT ---\n');
options_meco = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @(t, Y) event_meco(t, Y, p));

tspan1 = [0 200]; 
[T1, Y1] = ode45(@(t, Y) stage1_dynamics(t, Y, p), tspan1, Y0, options_meco);

Y_meco = Y1(end, :)';
final_mass = Y_meco(end);

remaining_propellant = final_mass - m_struc_1;
remaining_percentage = (remaining_propellant / p.m_prop_1) * 100;

fprintf('--- STEG 1 MECO SAMMANSTÄLLNING ---\n');
fprintf('MECO at t = %.f s\n', T1(end));
fprintf('MECO altitude: %.1f km | inertial speed: %.0f m/s\n', (norm(Y_meco(1:3)) - p.R_E)/1000, norm(Y_meco(4:6)));
fprintf('Remaining fuel: %.1f kg\n', remaining_propellant);

% ================= PHASE 2: STAGE 2 (TRANSFER BURN) =================
fprintf('--- PHASE 2: STAGE 2 TRANSFER BURN ---\n');
Y0_2 = Y_meco;
Y0_2(7) = m_prop2 + m_struc2 + m_pay; 

options_transfer = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @(t, Y) event_target_apogee(t, Y, p));
tspan2 = [T1(end), T1(end) + 500]; 
[T2, Y2] = ode45(@(t, Y) stage2_transfer(t, Y, p), tspan2, Y0_2, options_transfer);
fprintf('SECO-1 (Transfer Orbit Achieved) at t = %.1f s\n', T2(end));
fprintf('Remaining Stage 2 Mass: %.1f kg (Fuel left: %.1f kg)\n', Y2(end, 7), Y2(end, 7) - p.m_dry2);

% ================= PHASE 3: COAST TO APOGEE =================
fprintf('--- PHASE 3: COASTING TO APOGEE ---\n');
Y0_3 = Y2(end, :)';
options_coast = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @(t, Y) event_reach_apogee(t, Y));
tspan3 = [T2(end), T2(end) + 6000];
[T3, Y3] = ode45(@(t, Y) coast_dynamics(t, Y, p), tspan3, Y0_3, options_coast);
fprintf('Reached Apogee at t = %.1f s\n', T3(end));

% ================= PHASE 4: STAGE 2 (CIRCULARIZATION BURN) =================
fprintf('--- PHASE 4: CIRCULARIZATION BURN ---\n');
Y0_4 = Y3(end, :)';
options_circ = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @(t, Y) event_circularize(t, Y, p));
tspan4 = [T3(end), T3(end) + 300]; 
[T4, Y4] = ode45(@(t, Y) stage2_circularize(t, Y, p), tspan4, Y0_4, options_circ);
fprintf('SECO-2 (Orbit Circularized) at t = %.1f s\n', T4(end));
fprintf('Final Remaining Fuel: %.1f kg\n', Y4(end, 7) - p.m_dry2);

% ================= PHASE 5: FINAL ORBIT COAST =================
fprintf('--- PHASE 5: COASTING IN ORBIT ---\n');
Y0_5 = Y4(end, :)';
tspan5 = [T4(end), T4(end) + 6000]; 
[T5, Y5] = ode45(@(t, Y) coast_dynamics(t, Y, p), tspan5, Y0_5, options_standard);

% ================= PHASE 6: RENDEZOUS ===========================
% fprintf('--- PHASE 7: Rendezvous with debris---\n');
% Y0_6 = Y5(end, :)';
% options_rendezvous = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @(t, Y) event_rendezvous(t, Y, p));
% tspan6 = [T5(end), T5(end) + 6000]; % change 60000 so that it has engough
% time to catch up
% [T6, Y6] = ode45(@(t, Y) rendezvous_dynamics(t, Y, p), tspan6, Y0_6, options_rendezvous);

% ================= PHASE 5: DEORBIT =============================
% fprintf('--- PHASE 5: De-Orbit ---\n');
% Y0_7 = Y6(end, :)';
% tspan7 = [T6(end), T6(end) + 6000]; % Change 6000 to when it is de.orbited
% [T7, Y7] = ode45(@(t, Y) deorbit_dynamics(t, Y, p), tspan7, Y0_7, options_standard);

%% 6. CALCULATE RESULTS
T_out_1 = [T1; T2; T3; T4; T5]; % Use for falcon 9 plots (Figure 1)
% T_out_2 = [T6; T7]            % Use for satellite plots (Figure 2)
% T_out_tot = [T1; T2; T3; T4; T5; T6; T7]; % Use for the 3D plot of the orbits

Y_out_1 = [Y1; Y2; Y3; Y4; Y5]; % Use for falcon 9 plots (Figure 1)
% Y_out_2 = [Y6; Y7]            % Use for satellite plots (Figure 2)
% Y_out_tot = [Y1; Y2; Y3; Y4; Y5; Y6; Y7]; % Use for the 3D plot of the orbits


r_out_1 = Y_out_1(:, 1:3);
% r_out_2 = Y_out_2(:, 1:3);
% r_out_tot = Y_out_tot(:, 1:3);

v_out_1 = Y_out_1(:, 4:6);
% v_out_2 = Y_out_2(:, 4:6);
% v_out_tot = Y_out_tot(:, 4:6);

m_out_1 = Y_out_1(:, 7);
% m_out_2 = Y_out_1(:, 7);
% m_out_tot = Y_out_1(:, 7);


radii = vecnorm(r_out_1, 2, 2); 
altitudes = radii - p.R_E; 
velocities_1 = vecnorm(v_out_1, 2, 2);  
% velocities_2 = vecnorm(v_out_2, 2, 2);  



[rp_f, ra_f, inc_f, raan_f, ecc_f] = orbit_summary(Y5(end, 1:3)', Y5(end, 4:6)', p.my);

%% 7. LAUNCH WINDOW CALCULATION (CHANGED)
% GMST at 00:00 UTC on the launch date
d_0h = datetime(launch_date_str, 'TimeZone', 'UTC');
jd_0h = juliandate(d_0h);
T_cen = (jd_0h - 2451545.0) / 36525;     % Julian centuries since J2000.0
gmst_0h_deg = mod(100.4606184 + 36000.7700536 * T_cen + ...
                  0.000387933 * T_cen^2 - 2.583e-8 * T_cen^3, 360);

% Window 1 (ascending pass, northeast launch): required GST at liftoff
theta_GST_asc_deg = mod(target_RAAN_deg - raan_f, 360);

% Window 2 (descending pass, southeast launch): analytic estimate only (not simulated)
sin_delta_alpha = tan(Lat) / tan(i);
delta_alpha = asin(sin_delta_alpha);
theta_LST_desc = deg2rad(target_RAAN_deg) + pi - delta_alpha;
theta_GST_desc_deg = mod(rad2deg(theta_LST_desc - Long), 360);

% GST -> UTC
hours_utc_asc = mod(theta_GST_asc_deg - gmst_0h_deg, 360) / 15.0410686386;
hours_utc_desc = mod(theta_GST_desc_deg - gmst_0h_deg, 360) / 15.0410686386;
launch_time_asc = d_0h + hours(hours_utc_asc);
launch_time_desc = d_0h + hours(hours_utc_desc);

beta_asc_deg = rad2deg(p.beta);
beta_desc_deg = 180 - beta_asc_deg;

fprintf('\n=== LAUNCH WINDOW TARGETING ===\n');
fprintf('Target RAAN: %.2f deg | Target Inclination: %.2f deg\n\n', target_RAAN_deg, rad2deg(i));

fprintf('--- WINDOW 1 (Ascending Pass, calculated from the simulated final orbit) ---\n');
fprintf('Inertial plane azimuth: %.2f deg | Thrust steering azimuth: %.2f deg\n', beta_asc_deg, rad2deg(p.beta_steer));
fprintf('Required GST at liftoff: %.4f deg\n', theta_GST_asc_deg);
fprintf('Launch Time (UTC): %s\n\n', datestr(launch_time_asc, 'yyyy-mm-dd HH:MM:SS'));

fprintf('--- WINDOW 2 (Descending Pass, analytic estimate) ---\n');
fprintf('Launch Azimuth (Beta): %.2f deg (Southeast)\n', beta_desc_deg);
fprintf('Launch Time (UTC): %s\n\n', datestr(launch_time_desc, 'yyyy-mm-dd HH:MM:SS'));

fprintf('\n--- FINAL ORBIT STATS ---\n');
fprintf('Target Radius: %.2f km\n', p.target_radius / 1000);
fprintf('Actual Radius: %.2f km\n', radii(end) / 1000);
fprintf('Perigee radius: %.2f km | Apogee radius: %.2f km | Eccentricity: %.5f\n', rp_f/1000, ra_f/1000, ecc_f);
fprintf('Inclination: %.3f deg (target %.3f) | RAAN (inertial, at launch epoch): %.3f deg\n', inc_f, inc_target_deg, mod(raan_f + theta_GST_asc_deg, 360));
fprintf('Final Velocity: %.2f m/s (circular speed at target: %.2f m/s)\n', velocities_1(end), sqrt(p.my/p.target_radius));
fprintf('Final Altitude: %.2f km\n', altitudes(end) / 1000);

%% 8. PLOTTING

% --- FIGURE 1: 2x2 Telemetry Plots ---
figure('Name', 'Falcon 9 - Telemetry', 'Position', [100, 100, 1000, 600]);

subplot(2, 2, 1);
plot(T_out_1, radii / 1000, 'b', 'LineWidth', 2); hold on;
yline(p.target_radius / 1000, 'r--', 'Target Orbit');
xlabel('Time (s)'); ylabel('Radius (km)'); title('Orbital Radius'); grid on;

subplot(2, 2, 2);
plot(T_out_1, velocities_1 / 1000, 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Velocity (km/s)'); title('Absolute Velocity'); grid on;

subplot(2, 2, 3);
plot(T_out_1, m_out_1 / 1000, 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Mass (tons)'); title('Mass Depletion'); grid on;

subplot(2, 2, 4);
start_vec = r_out_1(1, :)';
downrange = zeros(length(T_out_1), 1);
for k = 1:length(T_out_1)
    t_current = T_out_1(k);
    theta = p.wE * t_current;
    Rz = [cos(theta),  sin(theta), 0;
         -sin(theta),  cos(theta), 0;
          0,           0,          1];
    r_ecef = Rz * r_out_1(k, :)';
    cos_angle = dot(start_vec, r_ecef) / (norm(start_vec) * norm(r_ecef));
    cos_angle = min(1, max(-1, cos_angle));
    downrange(k) = p.R_E * acos(cos_angle);
end
plot(downrange / 1000, altitudes / 1000, 'g', 'LineWidth', 2);
xlabel('Ground Distance (km)'); ylabel('Altitude (km)');
title('Flight Path Relative to Ground'); grid on;

% ----FIGURE 2: Same as figure 1 but with the satelite operations data-----



% --- FIGURE 3: Separate 3D Orbital Plot with Phase Markers ---
figure('Name', '3D Orbital Trajectory with Phase Markers', 'Position', [150, 150, 800, 800]);

[X_sph, Y_sph, Z_sph] = sphere(40);
surf(X_sph * (p.R_E/1000), Y_sph * (p.R_E/1000), Z_sph * (p.R_E/1000), ...
    'EdgeColor', 'none', 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.6); hold on;

plot3(r_out_1(:,1)/1000, r_out_1(:,2)/1000, r_out_1(:,3)/1000, 'r-', 'LineWidth', 2);

r_meco   = Y1(end, 1:3) / 1000;
r_seco1  = Y2(end, 1:3) / 1000;
r_seco2  = Y4(end, 1:3) / 1000;

plot3(r_meco(1), r_meco(2), r_meco(3), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
text(r_meco(1), r_meco(2), r_meco(3), '  MECO', 'Color', 'm', 'FontWeight', 'bold');

plot3(r_seco1(1), r_seco1(2), r_seco1(3), 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'c');
text(r_seco1(1), r_seco1(2), r_seco1(3), '  SECO-1', 'Color', 'c', 'FontWeight', 'bold');

plot3(r_seco2(1), r_seco2(2), r_seco2(3), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(r_seco2(1), r_seco2(2), r_seco2(3), '  SECO-2', 'Color', 'k', 'FontWeight', 'bold');

axis equal; grid on; view(3);
title('Falcon 9 - 3D Orbital Trajectory & Phase Markers'); 
xlabel('X ECI (km)'); ylabel('Y ECI (km)'); zlabel('Z ECI (km)');
legend('Earth', 'Orbit', 'MECO', 'SECO-1', 'SECO-2', 'Location', 'best');
hold off;

% -----------------Figure 3 

%% ==================== DYNAMICS FUNCTIONS ====================

function dY = stage1_dynamics(t, Y, p)
    r = Y(1:3); v = Y(4:6); m = Y(7);
    [FD, Patm] = aero_drag(r, v, p, p.A, p.C_D);

    m_flow_rate = p.T_vac / (p.g0 * p.Isp_vac); 
    Isp = p.Isp_vac - (p.Isp_vac - p.Isp_SL) * (Patm / p.P_SL);
    T_current = m_flow_rate * p.g0 * Isp;
    
    Z_loc = r / norm(r);                        
    E_loc = cross([0;0;1], Z_loc); E_loc = E_loc / norm(E_loc);          
    N_loc = cross(Z_loc, E_loc);            
    H_dir = cos(p.beta_steer) * N_loc + sin(p.beta_steer) * E_loc;      
    
    if t < p.t_kick
        u = Z_loc;
    else
        time_frac = min((t - p.t_kick) / (p.t_meco_est - p.t_kick), 1);
        pitch_angle = deg2rad(p.pitch_final * (time_frac^p.pitch_exp)); 
        u = cos(pitch_angle) * Z_loc + sin(pitch_angle) * H_dir;        
    end
    
    Fg = -(p.my * m / norm(r)^3) * r;
    FT = T_current * u;
    
    dY = [v; (Fg + FD + FT) / m; -m_flow_rate];
end

function dY = stage2_transfer(~, Y, p)
    r = Y(1:3); v = Y(4:6); m = Y(7);
    m_flow_rate = p.T_vac_2 / (p.g0 * p.Isp_vac_2); 
    T_current = p.T_vac_2;
    if m <= p.m_dry2 + 1                
        T_current = 0; m_flow_rate = 0;
    end

    Z_loc = r / norm(r);
    v_r = dot(v, Z_loc);                        
    v_h = v - v_r * Z_loc;                      
    H_dir = v_h / norm(v_h);                    
    h = norm(r) - p.R_E;
    v_r_des = p.Kh * (p.h_ref - h);
    theta = min(max(p.Kv * (v_r_des - v_r), -p.th_max), p.th_max);
    u = cos(theta) * H_dir + sin(theta) * Z_loc;

    Fg = -(p.my * m / norm(r)^3) * r;
    FD = aero_drag(r, v, p, p.A_2, p.C_D_2);    
    FT = T_current * u;
    dY = [v; (Fg + FD + FT) / m; -m_flow_rate];
end

function dY = stage2_circularize(~, Y, p)
    r = Y(1:3); v = Y(4:6); m = Y(7);
    m_flow_rate = p.T_vac_2 / (p.g0 * p.Isp_vac_2); 
    T_current = p.T_vac_2;
    
    Z_loc = r / norm(r);
    v_r = dot(v, Z_loc);
    v_h = v - v_r * Z_loc;
    H_dir = v_h / norm(v_h);
    theta = min(max(-0.02 * v_r, -0.3), 0.3);
    u = cos(theta) * H_dir + sin(theta) * Z_loc;

    if m <= p.m_dry2 + 1                
        T_current = 0; m_flow_rate = 0;
    end
    Fg = -(p.my * m / norm(r)^3) * r;
    FD = aero_drag(r, v, p, p.A_2, p.C_D_2);    
    FT = T_current * u;
    dY = [v; (Fg + FD + FT) / m; -m_flow_rate];
end

function dY = coast_dynamics(~, Y, p)
    r = Y(1:3); v = Y(4:6); m = Y(7);
    Fg = -(p.my * m / norm(r)^3) * r;
    FD = aero_drag(r, v, p, p.A_2, p.C_D_2);    
    dY = [v; (Fg + FD) / m; 0];
end

% ----------- INSERT RENDEZOUS FUNCTION HERE ----------------
% function dY = rendezvous_dynamics(~, Y, p)

% Pick one of the debris and target that orbit
% Might be better to add all the variables outside the function. Look at
% the other functions to see how I did it wit .p

% ----------- INSERT DECENT FUNCTION HERE -------------------
% % function dY = deorbit_dynamics(~, Y, p)

% Add the mass of the debris and de-orbit it


%% ==================== HELPER FUNCTIONS ====================

function [FD, Patm] = aero_drag(r, v, p, A, C_D)
    h = norm(r) - p.R_E;
    if h < 20000            % atmosisa is only valid to h<20000 m
        [~, ~, Patm, rho] = atmosisa(h); 
    elseif h >= 20000 && h < 32000      % At 32000 m, rho & patm ≈ 1% of sea level
        [~, ~, Patm0, rho0] = atmosisa(20000); 
        Patm = (32000 - h) / 12000 * Patm0;
        rho = (32000 - h) / 12000 * rho0;
    else
        Patm = 0; rho = 0;
    end
    v_wind = cross([0; 0; p.wE], r);
    v_rel = v - v_wind; v_rel_mag = norm(v_rel);
    FD = (v_rel_mag > 0) * (-0.5 * rho * v_rel_mag^2 * A * C_D * (v_rel / (v_rel_mag + 1e-6)));
end

function [inc_deg, raan_deg] = meco_plane(beta_s, p, Y0)
    p.beta_steer = beta_s;
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @(t, Y) event_meco(t, Y, p));
    [~, Yq] = ode45(@(t, Y) stage1_dynamics(t, Y, p), [0 200], Y0, opts);
    [~, ~, inc_deg, raan_deg] = orbit_summary(Yq(end, 1:3)', Yq(end, 4:6)', p.my);
end

function [rp, ra, inc_deg, raan_deg, ecc] = orbit_summary(r, v, mu)
    h_vec = cross(r, v); h_mag = norm(h_vec);
    epsilon = norm(v)^2 / 2 - mu / norm(r);
    a = -mu / (2 * epsilon);
    ecc = sqrt(max(0, 1 + 2 * epsilon * h_mag^2 / mu^2));
    rp = a * (1 - ecc);
    ra = a * (1 + ecc);
    inc_deg = acosd(h_vec(3) / h_mag);
    raan_deg = mod(atan2d(h_vec(1), -h_vec(2)), 360);
end

%% ==================== EVENT FUNCTIONS ====================

function [value, isterminal, direction] = event_meco(~, Y, p)
    m = Y(7);
    m_dry_stage1 = p.m0 - p.m_prop_1;
    
    value = m - m_dry_stage1; 
    isterminal = 1;           
    direction = -1;           
end

function [value, isterminal, direction] = event_target_apogee(~, Y, p)
    r = Y(1:3); v = Y(4:6);
    epsilon = (norm(v)^2)/2 - p.my/norm(r);
    h_mag = norm(cross(r, v));
    e = sqrt(1 + (2 * epsilon * h_mag^2) / (p.my^2));
    a = -p.my / (2 * epsilon);
    r_apogee = a * (1 + e);
    
    value = r_apogee - p.target_radius; 
    isterminal = 1;
    direction = 1; 
end

function [value, isterminal, direction] = event_reach_apogee(~, Y)
    r = Y(1:3); v = Y(4:6);
    value = dot(r, v);
    isterminal = 1;
    direction = -1; 
end

function [value, isterminal, direction] = event_circularize(~, Y, p)
    r = Y(1:3); v = Y(4:6);
    v_mag = norm(v);
    v_circ = sqrt(p.my / norm(r));
    
    value = v_mag - v_circ; 
    isterminal = 1;
    direction = 1;
end

% ---------INSTERT RENDEZOUS EVENT FUNCTION HERE -------------------
% function [value, isterminal, direction] = event_rendezvous(~, Y, p)

% Function that stops the ODE solver at rendezvous with the debris