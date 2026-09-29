clear;
clc;
close all;

%% ================================================================
% FALCON 9 ORBITAL SIMULATION
%
% The simulation includes:
% 1. Launch-window calculation
% 2. First-stage ascent
% 3. Stage separation
% 4. Second-stage ascent
% 5. Fairing separation
% 6. Propellant reservation for circularization
% 7. Target-apogee cutoff
% 8. Coast to apogee
% 9. Circularization burn
% 10. Orbital-element calculation
% 11. Ground-impact protection
% 12. Trajectory and dynamic-pressure plots
%
% All distances are expressed in meters internally.
% All velocities are expressed in meters per second internally.
%% ================================================================


%% 1. ENVIRONMENT PARAMETERS

p.g0 = 9.80665;                         % Standard gravity [m/s^2]
p.mu = 3.986004418e14;                  % Earth gravitational parameter [m^3/s^2]
p.wE = 7.292115e-5;                     % Earth rotation rate [rad/s]
p.R_E = 6378137;                        % Earth radius [m]
p.P_SL = 101325;                        % Sea-level pressure [Pa]

launch_latitude = deg2rad(28.5623);     % Launch latitude [rad]
launch_longitude = deg2rad(-80.5774);   % Launch longitude [rad]

target_inclination = deg2rad(74);       % Target inclination [rad]


%% 2. TARGET ORBIT

% The target orbit altitude is approximately 719.3 km.
p.target_altitude = 719318;             % Target altitude [m]
p.target_radius = p.R_E + p.target_altitude;

p.target_circular_velocity = ...
    sqrt(p.mu / p.target_radius);

fprintf('============================================================\n');
fprintf('TARGET ORBIT\n');
fprintf('============================================================\n');
fprintf('Target altitude:             %.3f km\n', ...
    p.target_altitude / 1000);
fprintf('Target radius:               %.3f km\n', ...
    p.target_radius / 1000);
fprintf('Target circular velocity:    %.3f m/s\n', ...
    p.target_circular_velocity);
fprintf('Target inclination:          %.3f deg\n\n', ...
    rad2deg(target_inclination));


%% 3. LAUNCH AZIMUTH

% Calculate the launch azimuth for a direct launch into the target
% inclination.
sin_launch_azimuth = ...
    cos(target_inclination) / cos(launch_latitude);

sin_launch_azimuth = ...
    max(-1, min(1, sin_launch_azimuth));

p.launch_azimuth = asin(sin_launch_azimuth);


%% 4. LAUNCH WINDOW CALCULATION

target_raan_deg = 45.0;
launch_date_string = '2026-10-15';

sin_delta_alpha = ...
    tan(launch_latitude) / tan(target_inclination);

sin_delta_alpha = ...
    max(-1, min(1, sin_delta_alpha));

delta_alpha = asin(sin_delta_alpha);

local_sidereal_time_ascending = ...
    deg2rad(target_raan_deg) + delta_alpha;

local_sidereal_time_descending = ...
    deg2rad(target_raan_deg) + pi - delta_alpha;

greenwich_sidereal_time_ascending_deg = mod( ...
    rad2deg(local_sidereal_time_ascending - launch_longitude), 360);

greenwich_sidereal_time_descending_deg = mod( ...
    rad2deg(local_sidereal_time_descending - launch_longitude), 360);

launch_date_utc = datetime( ...
    launch_date_string, ...
    'TimeZone', 'UTC');

julian_date = juliandate(launch_date_utc);
julian_centuries = ...
    (julian_date - 2451545.0) / 36525;

gmst_at_midnight_deg = mod( ...
    100.4606184 + ...
    36000.7700536 * julian_centuries + ...
    0.000387933 * julian_centuries^2 - ...
    2.583e-8 * julian_centuries^3, 360);

ascending_launch_hours = mod( ...
    greenwich_sidereal_time_ascending_deg - ...
    gmst_at_midnight_deg, 360) / 15.0410686386;

descending_launch_hours = mod( ...
    greenwich_sidereal_time_descending_deg - ...
    gmst_at_midnight_deg, 360) / 15.0410686386;

ascending_launch_time = ...
    launch_date_utc + hours(ascending_launch_hours);

descending_launch_time = ...
    launch_date_utc + hours(descending_launch_hours);

ascending_launch_azimuth_deg = ...
    rad2deg(p.launch_azimuth);

descending_launch_azimuth_deg = ...
    180 - ascending_launch_azimuth_deg;

fprintf('============================================================\n');
fprintf('LAUNCH WINDOW TARGETING\n');
fprintf('============================================================\n');
fprintf('Target RAAN:                 %.2f deg\n', target_raan_deg);
fprintf('Target inclination:          %.2f deg\n\n', ...
    rad2deg(target_inclination));

fprintf('Ascending launch window:\n');
fprintf('  Launch azimuth:            %.2f deg\n', ...
    ascending_launch_azimuth_deg);
fprintf('  Launch time UTC:           %s\n\n', ...
    datestr(ascending_launch_time, 'yyyy-mm-dd HH:MM:SS'));

fprintf('Descending launch window:\n');
fprintf('  Launch azimuth:            %.2f deg\n', ...
    descending_launch_azimuth_deg);
fprintf('  Launch time UTC:           %s\n\n', ...
    datestr(descending_launch_time, 'yyyy-mm-dd HH:MM:SS'));


%% 5. VEHICLE PARAMETERS

payload_mass = 15000;                    % Payload mass [kg]
fairing_mass = 1900;                     % Fairing mass [kg]

vehicle_diameter = 3.7;                  % Vehicle diameter [m]

p.reference_area = ...
    pi * (vehicle_diameter / 2)^2;

p.drag_coefficient = 0.3;

% First-stage propellant.
p.first_stage_propellant_mass = 386000;

% First-stage dry mass includes:
% First-stage structure
% Second-stage structure
% Second-stage propellant
% Payload
% Fairing
p.first_stage_dry_mass = ...
    28100 + 107500 + 25000 + payload_mass + fairing_mass;

p.first_stage_initial_mass = ...
    p.first_stage_dry_mass + ...
    p.first_stage_propellant_mass;

p.first_stage_isp_sea_level = 283;       % Sea-level Isp [s]
p.first_stage_isp_vacuum = 312;          % Vacuum Isp [s]
p.first_stage_vacuum_thrust = 8227000;   % Vacuum thrust [N]

% Second-stage masses.
p.second_stage_structural_mass = 4000;
p.second_stage_propellant_mass = 107500;

p.second_stage_dry_mass = ...
    p.second_stage_structural_mass + payload_mass;

p.second_stage_initial_mass = ...
    p.second_stage_structural_mass + ...
    p.second_stage_propellant_mass + ...
    fairing_mass + payload_mass;

% Reserve propellant for the circularization burn.
%
% This value can be tuned. Increasing it gives more propellant for
% circularization but reduces the propellant available during ascent.
p.circularization_reserve_mass = 18000;

p.second_stage_ascent_cutoff_mass = ...
    p.second_stage_dry_mass + ...
    p.circularization_reserve_mass;

p.second_stage_isp_vacuum = 348;         % Vacuum Isp [s]
p.second_stage_vacuum_thrust = 981000;   % Vacuum thrust [N]


%% 6. INITIAL CONDITIONS

initial_position = [
    p.R_E * cos(launch_latitude) * cos(launch_longitude);
    p.R_E * cos(launch_latitude) * sin(launch_longitude);
    p.R_E * sin(launch_latitude)
];

initial_velocity = [
    -p.wE * initial_position(2);
     p.wE * initial_position(1);
     0
];

initial_state = [
    initial_position;
    initial_velocity;
    p.first_stage_initial_mass
];


%% 7. NUMERICAL OPTIONS

base_solver_options = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', 1e-8);


%% 8. FIRST-STAGE ASCENT

first_stage_time_span = [0 164];

[first_stage_time, first_stage_state] = ode45( ...
    @(time, state) first_stage_dynamics(time, state, p), ...
    first_stage_time_span, ...
    initial_state, ...
    base_solver_options);

first_stage_final_state = ...
    first_stage_state(end, :)';

first_stage_final_mass = ...
    first_stage_final_state(7);

remaining_first_stage_propellant = ...
    first_stage_final_mass - p.first_stage_dry_mass;

fprintf('============================================================\n');
fprintf('FIRST-STAGE SUMMARY\n');
fprintf('============================================================\n');
fprintf('MECO time:                  %.3f s\n', ...
    first_stage_time(end));
fprintf('Mass at MECO:               %.3f kg\n', ...
    first_stage_final_mass);
fprintf('Remaining first-stage propellant: %.3f kg\n\n', ...
    remaining_first_stage_propellant);


%% 9. STAGE SEPARATION

% Discard the first stage and initialize the second-stage mass.
second_stage_initial_state = ...
    first_stage_final_state;

second_stage_initial_state(7) = ...
    p.second_stage_initial_mass;


%% 10. SECOND-STAGE INITIAL BURN

% The fairing remains attached during this segment.
second_stage_initial_time_span = [164 188];

[second_stage_initial_time, ...
 second_stage_initial_state_history] = ode45( ...
    @(time, state) second_stage_dynamics(time, state, p), ...
    second_stage_initial_time_span, ...
    second_stage_initial_state, ...
    base_solver_options);

fairing_separation_state = ...
    second_stage_initial_state_history(end, :)';


%% 11. FAIRING SEPARATION

fairing_separation_state(7) = ...
    fairing_separation_state(7) - fairing_mass;


%% 12. SECOND-STAGE ASCENT TO TARGET APOGEE

% The second stage stops when:
% 1. The predicted apogee reaches the target altitude, or
% 2. The ascent propellant limit is reached.
second_stage_options = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', 1e-8, ...
    'Events', @(time, state) ...
        second_stage_cutoff_events(time, state, p));

second_stage_time_span = [188 2500];

[second_stage_time, second_stage_state] = ode45( ...
    @(time, state) second_stage_dynamics(time, state, p), ...
    second_stage_time_span, ...
    fairing_separation_state, ...
    second_stage_options);

second_stage_cutoff_state = ...
    second_stage_state(end, :)';

second_stage_cutoff_orbit = ...
    calculate_orbital_elements(second_stage_cutoff_state, p);

fprintf('============================================================\n');
fprintf('SECOND-STAGE ASCENT SUMMARY\n');
fprintf('============================================================\n');
fprintf('Engine cutoff time:          %.3f s\n', ...
    second_stage_time(end));
fprintf('Altitude at cutoff:          %.3f km\n', ...
    second_stage_cutoff_orbit.altitude / 1000);
fprintf('Predicted apogee:            %.3f km\n', ...
    second_stage_cutoff_orbit.apogee_altitude / 1000);
fprintf('Predicted perigee:           %.3f km\n', ...
    second_stage_cutoff_orbit.perigee_altitude / 1000);
fprintf('Mass at cutoff:              %.3f kg\n', ...
    second_stage_cutoff_state(7));
fprintf('Orbital eccentricity:        %.8f\n\n', ...
    second_stage_cutoff_orbit.eccentricity);


%% 13. COAST TO APOGEE

% The coast phase stops at either:
% 1. Apogee, or
% 2. Ground impact.
coast_options = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-9, ...
    'Events', @(time, state) ...
        apogee_or_impact_events(time, state, p));

coast_time_span = [
    second_stage_time(end), ...
    second_stage_time(end) + 6000
];

[coast_time, coast_state] = ode45( ...
    @(time, state) coast_dynamics(time, state, p), ...
    coast_time_span, ...
    second_stage_cutoff_state, ...
    coast_options);

coast_final_state = ...
    coast_state(end, :)';

coast_final_orbit = ...
    calculate_orbital_elements(coast_final_state, p);

if coast_final_orbit.altitude <= 0

    error(['The vehicle impacted the Earth before reaching ', ...
           'apogee. Increase ascent performance or modify ', ...
           'the ascent guidance.']);
end

fprintf('============================================================\n');
fprintf('APOGEE BEFORE CIRCULARIZATION\n');
fprintf('============================================================\n');
fprintf('Apogee time:                 %.3f s\n', coast_time(end));
fprintf('Apogee altitude:             %.3f km\n', ...
    coast_final_orbit.altitude / 1000);
fprintf('Perigee altitude:            %.3f km\n', ...
    coast_final_orbit.perigee_altitude / 1000);
fprintf('Radial velocity:             %.6f m/s\n', ...
    coast_final_orbit.radial_velocity);
fprintf('Tangential velocity:         %.3f m/s\n', ...
    coast_final_orbit.tangential_velocity);
fprintf('Local circular velocity:     %.3f m/s\n\n', ...
    coast_final_orbit.circular_velocity);


%% 14. CIRCULARIZATION BURN

circularization_options = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-9, ...
    'Events', @(time, state) ...
        circularization_events(time, state, p));

circularization_time_span = [
    coast_time(end), ...
    coast_time(end) + 1500
];

[circularization_time, circularization_state] = ode45( ...
    @(time, state) ...
        circularization_dynamics(time, state, p), ...
    circularization_time_span, ...
    coast_final_state, ...
    circularization_options);

final_state = ...
    circularization_state(end, :)';

final_orbit = ...
    calculate_orbital_elements(final_state, p);


%% 15. FINAL ORBIT REPORT

fprintf('============================================================\n');
fprintf('FINAL ORBIT AFTER CIRCULARIZATION\n');
fprintf('============================================================\n');
fprintf('Final time:                  %.3f s\n', ...
    circularization_time(end));
fprintf('Final altitude:              %.3f km\n', ...
    final_orbit.altitude / 1000);
fprintf('Final apogee altitude:       %.3f km\n', ...
    final_orbit.apogee_altitude / 1000);
fprintf('Final perigee altitude:      %.3f km\n', ...
    final_orbit.perigee_altitude / 1000);
fprintf('Final radial velocity:       %.6f m/s\n', ...
    final_orbit.radial_velocity);
fprintf('Final tangential velocity:   %.3f m/s\n', ...
    final_orbit.tangential_velocity);
fprintf('Final circular velocity:     %.3f m/s\n', ...
    final_orbit.circular_velocity);
fprintf('Final eccentricity:           %.8f\n', ...
    final_orbit.eccentricity);
fprintf('Final vehicle mass:           %.3f kg\n', ...
    final_state(7));


%% 16. COMBINE OUTPUT DATA

time_output = [
    first_stage_time;
    second_stage_initial_time(2:end);
    second_stage_time(2:end);
    coast_time(2:end);
    circularization_time(2:end)
];

state_output = [
    first_stage_state;
    second_stage_initial_state_history(2:end, :);
    second_stage_state(2:end, :);
    coast_state(2:end, :);
    circularization_state(2:end, :)
];

position_output = state_output(:, 1:3);
velocity_output = state_output(:, 4:6);
mass_output = state_output(:, 7);

altitude_output = ...
    vecnorm(position_output, 2, 2) - p.R_E;

speed_output = ...
    vecnorm(velocity_output, 2, 2);


%% 17. DYNAMIC PRESSURE

dynamic_pressure_output = ...
    zeros(length(time_output), 1);

for index = 1:length(time_output)

    current_position = ...
        position_output(index, :)';

    current_velocity = ...
        velocity_output(index, :)';

    current_altitude = ...
        altitude_output(index);

    if current_altitude >= 0 && current_altitude < 84000

        [~, ~, ~, air_density] = ...
            atmosisa(current_altitude);

    else

        air_density = 0;
    end

    earth_rotation_vector = [0; 0; p.wE];

    atmospheric_velocity = ...
        cross(earth_rotation_vector, current_position);

    relative_velocity = ...
        current_velocity - atmospheric_velocity;

    relative_speed = norm(relative_velocity);

    dynamic_pressure_output(index) = ...
        0.5 * air_density * relative_speed^2;
end

[max_dynamic_pressure, max_q_index] = ...
    max(dynamic_pressure_output);

max_q_time = ...
    time_output(max_q_index);

max_q_altitude = ...
    altitude_output(max_q_index);

fprintf('\n============================================================\n');
fprintf('MAXIMUM DYNAMIC PRESSURE\n');
fprintf('============================================================\n');
fprintf('Maximum dynamic pressure:    %.3f kPa\n', ...
    max_dynamic_pressure / 1000);
fprintf('Time of Max-Q:               %.3f s\n', ...
    max_q_time);
fprintf('Altitude at Max-Q:           %.3f km\n\n', ...
    max_q_altitude / 1000);


%% 18. GROUND-RELATIVE DOWNRANGE DISTANCE

initial_position_vector = ...
    position_output(1, :)';

downrange_distance = ...
    zeros(length(time_output), 1);

for index = 1:length(time_output)

    current_time = ...
        time_output(index);

    current_position = ...
        position_output(index, :)';

    earth_rotation_angle = ...
        p.wE * current_time;

    eci_to_ecef_matrix = [
         cos(earth_rotation_angle),  sin(earth_rotation_angle), 0;
        -sin(earth_rotation_angle),  cos(earth_rotation_angle), 0;
         0,                          0,                         1
    ];

    current_ecef_position = ...
        eci_to_ecef_matrix * current_position;

    cosine_central_angle = ...
        dot(initial_position_vector, current_ecef_position) / ...
        (norm(initial_position_vector) * ...
         norm(current_ecef_position));

    cosine_central_angle = ...
        max(-1, min(1, cosine_central_angle));

    central_angle = ...
        acos(cosine_central_angle);

    downrange_distance(index) = ...
        p.R_E * central_angle;
end


%% 19. ALTITUDE, VELOCITY, MASS, AND DOWNRANGE PLOTS

figure( ...
    'Name', 'Falcon 9 Trajectory', ...
    'Position', [100 100 1100 750]);

subplot(2, 2, 1);

plot( ...
    time_output, ...
    altitude_output / 1000, ...
    'b', ...
    'LineWidth', 1.8);

grid on;
xlabel('Time [s]');
ylabel('Altitude [km]');
title('Altitude Profile');

subplot(2, 2, 2);

plot( ...
    time_output, ...
    speed_output / 1000, ...
    'r', ...
    'LineWidth', 1.8);

grid on;
xlabel('Time [s]');
ylabel('Inertial Speed [km/s]');
title('Velocity Profile');

subplot(2, 2, 3);

plot( ...
    time_output, ...
    mass_output / 1000, ...
    'k', ...
    'LineWidth', 1.8);

grid on;
xlabel('Time [s]');
ylabel('Mass [metric tons]');
title('Mass Depletion');

subplot(2, 2, 4);

plot( ...
    downrange_distance / 1000, ...
    altitude_output / 1000, ...
    'g', ...
    'LineWidth', 1.8);

grid on;
xlabel('Ground Distance [km]');
ylabel('Altitude [km]');
title('Ground-Relative Flight Path');


%% 20. DYNAMIC-PRESSURE PLOT

figure( ...
    'Name', 'Dynamic Pressure', ...
    'Position', [200 200 900 450]);

plot( ...
    time_output, ...
    dynamic_pressure_output / 1000, ...
    'm', ...
    'LineWidth', 1.8);

hold on;

plot( ...
    max_q_time, ...
    max_dynamic_pressure / 1000, ...
    'ro', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');

grid on;

xlim([0 min(150, time_output(end))]);

xlabel('Time [s]');
ylabel('Dynamic Pressure [kPa]');

title(sprintf( ...
    'Max-Q = %.2f kPa at t = %.2f s, altitude = %.2f km', ...
    max_dynamic_pressure / 1000, ...
    max_q_time, ...
    max_q_altitude / 1000));

hold off;


%% 21. THREE-DIMENSIONAL ECI TRAJECTORY

figure( ...
    'Name', 'ECI Orbital Trajectory', ...
    'Position', [250 150 850 850]);

earth_radius_km = ...
    p.R_E / 1000;

[sphere_x, sphere_y, sphere_z] = ...
    sphere(60);

surf( ...
    sphere_x * earth_radius_km, ...
    sphere_y * earth_radius_km, ...
    sphere_z * earth_radius_km, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.2 0.5 0.8], ...
    'FaceAlpha', 0.6);

hold on;

plot3( ...
    position_output(:, 1) / 1000, ...
    position_output(:, 2) / 1000, ...
    position_output(:, 3) / 1000, ...
    'r', ...
    'LineWidth', 2);

plot3( ...
    position_output(1, 1) / 1000, ...
    position_output(1, 2) / 1000, ...
    position_output(1, 3) / 1000, ...
    'go', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'g');

plot3( ...
    position_output(end, 1) / 1000, ...
    position_output(end, 2) / 1000, ...
    position_output(end, 3) / 1000, ...
    'ko', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'k');

axis equal;
grid on;

xlabel('ECI X [km]');
ylabel('ECI Y [km]');
zlabel('ECI Z [km]');

title('Falcon 9 Earth-Centered Inertial Trajectory');

legend( ...
    'Earth', ...
    'Rocket trajectory', ...
    'Launch position', ...
    'Final position', ...
    'Location', 'best');

view(3);
hold off;


%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================


function derivative = ...
    first_stage_dynamics(time, state, p)

    position = state(1:3);
    velocity = state(4:6);
    mass = state(7);

    radius = norm(position);
    altitude = radius - p.R_E;

    if altitude >= 0 && altitude < 84000

        [~, ~, atmospheric_pressure, air_density] = ...
            atmosisa(altitude);

    else

        atmospheric_pressure = 0;
        air_density = 0;
    end

    earth_rotation_vector = [0; 0; p.wE];

    atmospheric_velocity = ...
        cross(earth_rotation_vector, position);

    relative_velocity = ...
        velocity - atmospheric_velocity;

    relative_speed = ...
        norm(relative_velocity);

    dynamic_pressure = ...
        0.5 * air_density * relative_speed^2;

    % The previous value of 355400 Pa was 355.4 kPa.
    % The limit below is set to 30 kPa.
    max_q_limit = 30000;
    max_q_warning = 0.8 * max_q_limit;

    throttle_setting = 1.0;

    if dynamic_pressure > max_q_warning && altitude < 84000

        throttle_setting = ...
            1.0 - 0.4 * ...
            (dynamic_pressure - max_q_warning) / ...
            (max_q_limit - max_q_warning);

        throttle_setting = ...
            max(0.60, min(1.0, throttle_setting));
    end

    mass_flow_rate = ...
        p.first_stage_vacuum_thrust / ...
        (p.g0 * p.first_stage_isp_vacuum) * ...
        throttle_setting;

    effective_isp = ...
        p.first_stage_isp_vacuum - ...
        (p.first_stage_isp_vacuum - ...
         p.first_stage_isp_sea_level) * ...
        atmospheric_pressure / p.P_SL;

    thrust_magnitude = ...
        mass_flow_rate * p.g0 * effective_isp;

    radial_direction = ...
        position / radius;

    earth_axis = [0; 0; 1];

    east_direction = ...
        cross(earth_axis, radial_direction);

    if norm(east_direction) < 1e-8

        east_direction = [0; 1; 0];

    else

        east_direction = ...
            east_direction / norm(east_direction);
    end

    north_direction = ...
        cross(radial_direction, east_direction);

    horizontal_direction = ...
        cos(p.launch_azimuth) * north_direction + ...
        sin(p.launch_azimuth) * east_direction;

    if time < 10

        thrust_direction = radial_direction;

    else

        target_pitch_from_vertical_deg = 70;

        normalized_time = ...
            (time - 10) / (164 - 10);

        normalized_time = ...
            max(0, min(1, normalized_time));

        pitch_angle = deg2rad( ...
            target_pitch_from_vertical_deg * normalized_time^2);

        thrust_direction = ...
            cos(pitch_angle) * radial_direction + ...
            sin(pitch_angle) * horizontal_direction;

        thrust_direction = ...
            thrust_direction / norm(thrust_direction);
    end

    gravity_force = ...
        -p.mu * mass / radius^3 * position;

    if relative_speed > 1e-8

        drag_force = ...
            -0.5 * air_density * relative_speed^2 * ...
            p.reference_area * p.drag_coefficient * ...
            relative_velocity / relative_speed;

    else

        drag_force = [0; 0; 0];
    end

    if mass <= p.first_stage_dry_mass

        thrust_magnitude = 0;
        mass_flow_rate = 0;
    end

    thrust_force = ...
        thrust_magnitude * thrust_direction;

    position_derivative = velocity;

    velocity_derivative = ...
        (gravity_force + drag_force + thrust_force) / mass;

    mass_derivative = ...
        -mass_flow_rate;

    derivative = [
        position_derivative;
        velocity_derivative;
        mass_derivative
    ];
end


function derivative = ...
    second_stage_dynamics(~, state, p)

    position = state(1:3);
    velocity = state(4:6);
    mass = state(7);

    radius = norm(position);

    radial_direction = ...
        position / radius;

    angular_momentum = ...
        cross(position, velocity);

    if norm(angular_momentum) < 1e-8

        tangential_direction = ...
            velocity / max(norm(velocity), 1e-8);

    else

        tangential_direction = ...
            cross(angular_momentum, position);

        tangential_direction = ...
            tangential_direction / norm(tangential_direction);
    end

    radial_velocity = ...
        dot(velocity, radial_direction);

    tangential_velocity = ...
        dot(velocity, tangential_direction);

    % The vehicle burns mainly in the prograde direction.
    % A radial correction reduces excessive climb or descent velocity.
    radial_correction_gain = 0.0015;

    radial_velocity_correction = ...
        -radial_correction_gain * radial_velocity;

    thrust_direction = ...
        tangential_direction + ...
        radial_velocity_correction * radial_direction;

    thrust_direction = ...
        thrust_direction / norm(thrust_direction);

    % Add a small radial component during the early ascent.
    % This helps avoid an overly horizontal trajectory immediately
    % after stage separation.
    if radius < p.R_E + 150000

        vertical_component = 0.08;

        thrust_direction = ...
            thrust_direction + ...
            vertical_component * radial_direction;

        thrust_direction = ...
            thrust_direction / norm(thrust_direction);
    end

    thrust_magnitude = ...
        p.second_stage_vacuum_thrust;

    mass_flow_rate = ...
        p.second_stage_vacuum_thrust / ...
        (p.g0 * p.second_stage_isp_vacuum);

    % Stop the ascent burn while retaining propellant for
    % the circularization burn.
    if mass <= p.second_stage_ascent_cutoff_mass

        thrust_magnitude = 0;
        mass_flow_rate = 0;
    end

    gravity_force = ...
        -p.mu * mass / radius^3 * position;

    thrust_force = ...
        thrust_magnitude * thrust_direction;

    position_derivative = velocity;

    velocity_derivative = ...
        (gravity_force + thrust_force) / mass;

    mass_derivative = ...
        -mass_flow_rate;

    derivative = [
        position_derivative;
        velocity_derivative;
        mass_derivative
    ];
end


function derivative = ...
    coast_dynamics(~, state, p)

    position = state(1:3);
    velocity = state(4:6);
    mass = state(7);

    radius = norm(position);

    gravity_acceleration = ...
        -p.mu / radius^3 * position;

    position_derivative = velocity;
    velocity_derivative = gravity_acceleration;
    mass_derivative = 0;

    derivative = [
        position_derivative;
        velocity_derivative;
        mass_derivative
    ];
end


function derivative = ...
    circularization_dynamics(~, state, p)

    position = state(1:3);
    velocity = state(4:6);
    mass = state(7);

    radius = norm(position);

    angular_momentum = ...
        cross(position, velocity);

    if norm(angular_momentum) < 1e-8

        tangential_direction = ...
            velocity / max(norm(velocity), 1e-8);

    else

        tangential_direction = ...
            cross(angular_momentum, position);

        tangential_direction = ...
            tangential_direction / norm(tangential_direction);
    end

    tangential_velocity = ...
        dot(velocity, tangential_direction);

    local_circular_velocity = ...
        sqrt(p.mu / radius);

    if tangential_velocity < local_circular_velocity && ...
            mass > p.second_stage_dry_mass

        thrust_magnitude = ...
            p.second_stage_vacuum_thrust;

        mass_flow_rate = ...
            p.second_stage_vacuum_thrust / ...
            (p.g0 * p.second_stage_isp_vacuum);

        thrust_direction = ...
            tangential_direction;

    else

        thrust_magnitude = 0;
        mass_flow_rate = 0;
        thrust_direction = [0; 0; 0];
    end

    gravity_force = ...
        -p.mu * mass / radius^3 * position;

    thrust_force = ...
        thrust_magnitude * thrust_direction;

    position_derivative = velocity;

    velocity_derivative = ...
        (gravity_force + thrust_force) / mass;

    mass_derivative = ...
        -mass_flow_rate;

    derivative = [
        position_derivative;
        velocity_derivative;
        mass_derivative
    ];
end


function [value, is_terminal, direction] = ...
    second_stage_cutoff_events(~, state, p)

    orbital_data = ...
        calculate_orbital_elements(state, p);

    current_mass = ...
        state(7);

    target_apogee_error = ...
        orbital_data.apogee_altitude - p.target_altitude;

    propellant_reserve_error = ...
        current_mass - p.second_stage_ascent_cutoff_mass;

    % If the current orbit is unbound, disable the apogee event.
    if ~isfinite(orbital_data.apogee_altitude)

        target_apogee_error = 1e12;
    end

    value = [
        target_apogee_error;
        propellant_reserve_error
    ];

    is_terminal = [
        1;
        1
    ];

    direction = [
        1;
        -1
    ];
end


function [value, is_terminal, direction] = ...
    apogee_or_impact_events(~, state, p)

    position = state(1:3);
    velocity = state(4:6);

    radius = norm(position);

    radial_velocity = ...
        dot(position, velocity) / radius;

    altitude = ...
        radius - p.R_E;

    value = [
        radial_velocity;
        altitude
    ];

    is_terminal = [
        1;
        1
    ];

    direction = [
        -1;
        -1
    ];
end


function [value, is_terminal, direction] = ...
    circularization_events(~, state, p)

    position = state(1:3);
    velocity = state(4:6);

    radius = norm(position);

    angular_momentum = ...
        cross(position, velocity);

    tangential_velocity = ...
        norm(angular_momentum) / radius;

    local_circular_velocity = ...
        sqrt(p.mu / radius);

    radial_velocity = ...
        dot(position, velocity) / radius;

    velocity_error = ...
        tangential_velocity - local_circular_velocity;

    % Circularization is only allowed when the radial velocity is
    % already close to zero.
    if abs(radial_velocity) > 5

        value = -abs(velocity_error);

    else

        value = velocity_error;
    end

    is_terminal = 1;
    direction = 1;
end


function orbital_data = ...
    calculate_orbital_elements(state, p)

    position = state(1:3);
    velocity = state(4:6);

    radius = norm(position);
    speed = norm(velocity);

    radial_velocity = ...
        dot(position, velocity) / radius;

    angular_momentum_vector = ...
        cross(position, velocity);

    angular_momentum = ...
        norm(angular_momentum_vector);

    tangential_velocity = ...
        angular_momentum / radius;

    circular_velocity = ...
        sqrt(p.mu / radius);

    specific_orbital_energy = ...
        speed^2 / 2 - p.mu / radius;

    if specific_orbital_energy < 0

        semi_major_axis = ...
            -p.mu / (2 * specific_orbital_energy);

        eccentricity_vector = ...
            cross(velocity, angular_momentum_vector) / p.mu - ...
            position / radius;

        eccentricity = ...
            norm(eccentricity_vector);

        perigee_radius = ...
            semi_major_axis * (1 - eccentricity);

        apogee_radius = ...
            semi_major_axis * (1 + eccentricity);

        perigee_altitude = ...
            perigee_radius - p.R_E;

        apogee_altitude = ...
            apogee_radius - p.R_E;

    else

        semi_major_axis = Inf;
        eccentricity = Inf;
        perigee_altitude = -Inf;
        apogee_altitude = Inf;
    end

    orbital_data.radius = radius;
    orbital_data.altitude = radius - p.R_E;
    orbital_data.speed = speed;
    orbital_data.radial_velocity = radial_velocity;
    orbital_data.tangential_velocity = tangential_velocity;
    orbital_data.circular_velocity = circular_velocity;
    orbital_data.specific_orbital_energy = ...
        specific_orbital_energy;
    orbital_data.semi_major_axis = ...
        semi_major_axis;
    orbital_data.eccentricity = ...
        eccentricity;
    orbital_data.perigee_altitude = ...
        perigee_altitude;
    orbital_data.apogee_altitude = ...
        apogee_altitude;
end