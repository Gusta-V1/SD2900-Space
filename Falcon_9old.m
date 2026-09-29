
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
p.beta = asin(cos(i)/cos(Lat));          % Launch azimuth (Radians)

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
p.m_prop_1 = 386000;                       % Propellant mass for Stage 1 [kg]
p.m0 = m_struc_1 + p.m_prop_1;               % Gross liftoff mass [kg]

p.Isp_SL = 283;                          % Isp sea level [s]
p.Isp_vac = 312;                         % Isp vacuum [s]
p.T_vac = 8227000;                       % Vacuum thrust for 9x Merlin [N]

%% 3. INITIAL CONDITIONS (t = 0)
x0 = p.R_E * cos(Lat) * cos(Long);
y0 = p.R_E * cos(Lat) * sin(Long);
z0 = p.R_E * sin(Lat);
vx0 = -p.wE * y0;
vy0 = p.wE * x0;
vz0 = 0;

Y0 = [x0; y0; z0; vx0; vy0; vz0; p.m0];
tspan = [0 164];                         % Stage 1 burn time MECO target=162



%% 4. RUN SIMULATION
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
% Pass 'p' using an anonymous function
[T1, Y1] = ode45(@(t, Y) rocket_dynamics(t, Y, p), tspan, Y0, options);


m_struc_1_total = p.m0 - p.m_prop_1;

% Raketens totala massa vid exakt slutdatum (MECO)
Y0_2 = Y1(end, :)';
final_mass = Y0_2(end);

% Kvarvarande bränsle är den totala slutmassan minus strukturmassan
remaining_propellant = final_mass - m_struc_1_total;
remaining_percentage = (remaining_propellant / p.m_prop_1) * 100;

% Skriv ut resultatet i kommandofönstret
fprintf('--- STEG 1 MECO SAMMANSTÄLLNING ---\n');
fprintf('Kvarvarande bränsle: %.1f kg\n', remaining_propellant);
fprintf('Kvarvarande proportion: %.1f%% av ursprungligt bränsle\n', remaining_percentage);


%% ----------------------------- Second Stage -----------------------------

p.Isp_vac_2 = 348;       % Vakuum Isp för Steg 2 (Merlin Vac)
p.T_vac_2 = 981000;      % Thrust i vakuum för Steg 2

% 1. Hämta SISTA tillståndet från Steg 1
Y0_2 = Y1(end, :)'; 

% 2. Byt ut massan till Steg 2:s startmassa
m_struc2 = 4000;
m_prop2 = 107500;
p.m02 = m_prop2 + m_struc2 + fairing + m_pay;
Y0_2(7) = p.m02; 

% 3. Kör simuleringen (notera att vi anropar stage2_dynamics)
tspan2 = [164, 188];
[T2, Y2] = ode45(@(t, Y) stage2_dynamics(t, Y, p), tspan2, Y0_2, options);


%% --------------------- Second Stage W/O fairing -------------------------
Y0_3 = Y2(end, :)';
Y0_3(7) = Y0_3(7) - fairing;

tspan3 = [188, 1000]; % Ge den gott om tid, event-funktionen stoppar den i förtid!

% Skapa options som säger åt ode45 att lyssna på vår MECO-event
options_meco = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @(t, Y) meco_event(t, Y, p));

% Kör fas 3 med den nya options-inställningen
[T3, Y3] = ode45(@(t, Y) stage2_dynamics(t, Y, p), tspan3, Y0_3, options_meco);

fprintf('SECO inträffade vid t = %.1f s\n', T3(end));

%% --------------------- Coasting i Rymden (Fas 4) -------------------------
% Hämta tillståndet där motorn stängdes av
Y0_4 = Y3(end, :)';

% Simulera i 10 000 sekunder i fritt fall
tspan4 = [T3(end), T3(end) + 10000]; 
[T4, Y4] = ode45(@(t, Y) coast_dynamics(t, Y, p), tspan4, Y0_4, options);

%% 5. CALCULATE RESULTS FOR PLOTTING
T_out = [T1; T2; T3; T4];
Y_out = [Y1; Y2; Y3; Y4];

r_out = Y_out(:, 1:3);
v_out = Y_out(:, 4:6);
m_out = Y_out(:, 7);

altitudes = vecnorm(r_out, 2, 2) - p.R_E; % Altitude in meters
velocities = vecnorm(v_out, 2, 2);        % Absolute velocity in m/s
target_vel = sqrt(p.my/(altitudes(end)+p.R_E));


fprintf('Final altitude: %.2f km\n', altitudes(end) / 1000);
fprintf('Final velocity: %.2f m/s\n', velocities(end));
fprintf('Needed velocity at the alt: %.2f m/s\n', target_vel)


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
%% -----------------------------------------------------------------------
% DIFFERENTIAL EQUATION FUNCTION
%% -----------------------------------------------------------------------
function dY = rocket_dynamics(t, Y, p)
    

    % 1. Extract from Y
    r = Y(1:3);
    v = Y(4:6);
    m = Y(7);
    
    h = norm(r) - p.R_E;
    

    % 2. Handle Atmosphere (atmosisa crashes above 86km)
    if h < 84000
        [~, ~, Patm, rho] = atmosisa(h);
    else
        Patm = 0;
        rho = 0;
    end
    
    
    % 3. Relative velocity (Wind)
    omega_vec = [0; 0; p.wE];
    v_wind = cross(omega_vec, r);
    v_rel = v - v_wind;
    v_rel_mag = norm(v_rel);
    
    % 3.5 Dynamic pressure
    q = 0.5 * rho * v_rel_mag^2;

    % --- NY KOD: Throttling för Max-Q ---
    q_limit = 355400; % Smärtgräns för strukturen: 30 kPa
    q_warning = q_limit * 0.8; % Börja strypa vid 24 kPa
    
    throttle_q = 1.0; % Som standard gasar vi 100%
    
    if q > q_warning && h < 84000
        % Procentuell nedstrypning: om q = 30000, blir kvoten 1, och vi drar av 40% gas
        throttle_q = 1.0 - 0.2 * ((q - q_warning) / (q_limit - q_warning));
        
        % Motorerna kan inte gå lägre än 60% i verkligheten
        if throttle_q < 0.60
            throttle_q = 0.60;
        end
    end

    % 4. Engine performance 
    m_flow_rate = p.T_vac / (p.g0 * p.Isp_vac) * throttle_q; 
    Isp = p.Isp_vac - (p.Isp_vac - p.Isp_SL) * (Patm / p.P_SL);
    T_current = m_flow_rate * p.g0 * Isp; % * throttle_q;
    

    % 5. Create unit vectors for direction
    Z_loc = r / norm(r);                  
    k_pole = [0; 0; 1];                   
    E_loc = cross(k_pole, Z_loc);
    E_loc = E_loc / norm(E_loc);          
    N_loc = cross(Z_loc, E_loc);          
    
    H_dir = cos(p.beta) * N_loc + sin(p.beta) * E_loc;
    
    % 6. Steering logic
    if t < 10
        u = Z_loc;
    elseif t >= 10 && t <= 15
        pitch_angle = deg2rad(0.9); 
        u = cos(pitch_angle) * Z_loc + sin(pitch_angle) * H_dir;
    else
        % Avoid division by zero right at the launch pad
        if v_rel_mag > 0.1
            u = v_rel / v_rel_mag;
        else
            u = Z_loc;
        end
    end
    
    % 7. Forces
    Fg = -(p.my * m / norm(r)^3) * r;
    
    if v_rel_mag > 0
        FD = -0.5 * rho * v_rel_mag^2 * p.A * p.C_D * (v_rel / v_rel_mag);
    else
        FD = [0; 0; 0];
    end


    % Saftey so that the engin cant burn m_struc
    if m <= (p.m0 - p.m_prop_1)
        T_current = 0;
        m_flow_rate = 0;
    end

    FT = T_current * u;
    
    % 8. Derivatives
    r_dot = v;
    v_dot = (Fg + FD + FT) / m;
    m_dot = -m_flow_rate;
        
    dY = [r_dot; v_dot; m_dot];
end


%% -----------------------------------------------------------------------
% DIFFERENTIAL EQUATION FUNCTION - STAGE 2
%% -----------------------------------------------------------------------
function dY = stage2_dynamics(t, Y, p)
    % 1. Extrahera variabler
    r = Y(1:3);
    v = Y(4:6);
    m = Y(7);
    
    % 2. Motorprestanda i Vakuum
    m_flow_rate = p.T_vac_2 / (p.g0 * p.Isp_vac_2); 
    T_current = p.T_vac_2*3;

    if m <= (4000 + 15000) % Steg 2 torrvikt + payload
        T_current = 0;
        m_flow_rate = 0;
    end

    % 3. Styrning (Prograde: Peka i hastighetens riktning)
    v_mag = norm(v);
    if v_mag > 0
        u = v / v_mag;
    else
        u = r / norm(r); % Fallback
    end
    
    % 3. Styrning (Lokal horisontell styrning)
    % Räkna ut vad som är rakt upp från jorden (radial)
    r_hat = r / norm(r); 
    
    % Beräkna rörelsemängdsmomentet (en vektor som står vinkelrätt mot planet raketen flyger i)
    h_vec = cross(r, v); 
    
    if v_mag >= sqrt(p.my/norm(r))
        T_current=0;
        m_flow_rate=0;
    end
    
    
    % Beräkna den lokala horisontella riktningen (rakt fram parallellt med marken)
    if norm(h_vec) > 0
        horiz_dir = cross(h_vec, r); 
        horiz_dir = horiz_dir / norm(horiz_dir);
        
        % Tvinga motorn att peka perfekt horisontellt
        u = horiz_dir; 
    else
        % Fallback om något går fel
        u = v / norm(v); 
    end
    

    % 4. Krafter
    Fg = -(p.my * m / norm(r)^3) * r;
    FD = [0; 0; 0]; % Inget luftmotstånd i rymden
    FT = T_current * u;
    
    % 5. Derivator
    r_dot = v;
    v_dot = (Fg + FD + FT) / m;
    m_dot = -m_flow_rate;
        
    dY = [r_dot; v_dot; m_dot];
end


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
% EVENT FUNCTION - STOPPA NÄR OMLOPPSHASTIGHET NÅS
%% -----------------------------------------------------------------------
function [value, isterminal, direction] = meco_event(t, Y, p)
    r = Y(1:3);
    v = Y(4:6);
    v_mag = norm(v);
    v_target = sqrt(p.my / norm(r));
    
    value = v_mag - v_target; % Blir 0 exakt när vi når omloppshastighet
    isterminal = 1;           % 1 = Stoppa simuleringen (MECO)
    direction = 1;            % Triggern sker när vi går från för långsam till snabb
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