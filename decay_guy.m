% US Standard Atmosphere 1976 (Table 4.2)
% Altitude (km)
altitude_km = [0, 150, 200, 250, 300, 350, 400, 450, 500, 550, ...
               600, 650, 700, 750, 800, 850, 900, 950, 1000];

% Density (kg/m^3)
density = [1.225, 2.076e-9, 2.541e-10, 6.073e-11, 1.916e-11, ...
           7.014e-12, 2.803e-12, 1.184e-12, 5.215e-13, 2.384e-13, ...
           1.137e-13, 5.712e-14, 3.070e-14, 1.788e-14, 1.136e-14, ...
           7.824e-15, 5.759e-15, 4.453e-15, 3.561e-15];

% Pressure (Pa)
pressure = [1.01e5, 4.54e-4, 8.47e-5, 2.48e-5, 8.77e-6, ...
            3.45e-6, 1.45e-6, 6.45e-7, 3.02e-7, 1.51e-7, ...
            8.21e-8, 4.89e-8, 3.19e-8, 2.26e-8, 1.70e-8, ...
            1.34e-8, 1.09e-8, 8.98e-9, 7.51e-9];

% Scale height (km)
scale_height_km = [8.4345, 23.380, 36.183, 44.924, 51.193, ...
                    55.832, 59.678, 63.644, 68.785, 76.427, ...
                    88.244, 105.992, 130.630, 161.074, 193.862, ...
                    224.737, 250.894, 271.754, 288.203];

% Molecular weight (kg/kmol)
mol_weight = [29.0, 24.1, 21.3, 19.2, 17.7, 16.7, 16.0, 15.3, ...
              14.3, 13.1, 11.5, 9.72, 8.00, 6.58, 5.54, 4.85, ...
              4.40, 4.12, 3.94];

% Optional: combine into a table for easy inspection/export
atm_table = table(altitude_km', density', pressure', scale_height_km', mol_weight', ...
    'VariableNames', {'Altitude_km', 'Density_kgm3', 'Pressure_Pa', 'ScaleHeight_km', 'MolWeight_kgkmol'});

mu_E = 3.986004418E14;
R_E = 6371E3;

C_D = 2.2;
A = 4.5; % Frontal area
m = 1435 + 900; % Combined debris + vehicle mass
B = C_D*A/m;
N = length(atm_table.Altitude_km);
%B = 2.252378E-03;
%atm_table.Properties.VariableNames

deltaT = zeros(N,2);

dt = @(H,Hinit,Hscale,rho_0) - Hscale / (sqrt(mu_E*R_E)*B*rho_0)*(exp(H/Hscale)-exp(Hinit/Hscale));

target_time = 25*365*24*60*60; % Target time in seconds

rho_0 = atm_table.Density_kgm3./exp(-atm_table.Altitude_km./atm_table.ScaleHeight_km);

for i = N-1:-1:1
    func_return = dt(atm_table.Altitude_km(i)*1000,atm_table.Altitude_km(i+1)*1000,atm_table.ScaleHeight_km(i)*1000,rho_0(i));
    deltaT(i,:) = [atm_table.Altitude_km(i) func_return];
end

time_elapse = deltaT(1,2);
for i = 2:N
    time_elapse = time_elapse + deltaT(i,2);

    if time_elapse >= target_time
        disp("Maximum orbit alt: " + atm_table.Altitude_km(i) + ", layer: " + i);
        break
    end
end
