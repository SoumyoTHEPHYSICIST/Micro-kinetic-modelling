function [T,Y] = solve_compass
%% kf (forward Rate Constant Values)
kf=[1.5*10^-3 9.8*10^-4 4*10^-20 4.2*10^-13 6.4*10^-18 3.0*10^-11 9.5*10^-23 7.5*10^-17 1.1*10^-19 2.3*10^-20 9.3*10^-12 4.2*10^-19 3.9*10^-14 5.7*10^-11 2.2*10^-10 4.5*10^-10 6.5*10^-10 2.1*10^-10 7.8*10^-11 2.8*10^-11 2.0*10^-10 3.8*10^-11];
npert = 22;

%% Time span for the simulation (0 to 10^5 seconds)
tspan_fixed = linspace(0, 0.03, 10000);
tspan = [0 2];

%moles of CH4 = 4.308*10^-3
%moles of ZnCl2 = 9.5736*10^-4
N=6.022*10^23*10^-6;
PCH4=0.45*101325; %Pressure of CH4 in pascals
v0=1*10^-3; % in m^3
y0_CH4=PCH4*v0/(8.314*1273);
PAr=0.45*101325; %Pressure of Ar in pascals
y0_Ar=PAr*v0/(8.314*1273);
PZncl2=0.1*101325; %Pressure of Na in pascals
y0_ZnCl2=PZncl2*v0/(8.314*1273);

%% Initialize concentration vector for all 16 species
yi = zeros(16, 1);  % 1x16 matrix for initial concentrations

% Set initial concentrations for CH4, ZnCl2, and Ar
yi(1) = y0_CH4;     % Methane (CH4)
yi(4) = y0_ZnCl2;   % ZnCl2
yi(16) = y0_Ar;     % Argon (inert gas)

yi = yi*N; % Converting to molecules/cm^3

%% Solving Differential Equations using ode15s
[T Y] = ode15s(@compass, tspan, yi);
y_fixed = deval(ode15s(@compass, tspan, yi), tspan_fixed);
assignin('base', 'y_fixed', y_fixed);

%Perturbation vector
pert = zeros(22,1);
n = length(tspan_fixed);
ro = zeros(1, n);

for i = 1:n
    ro(i) = r_CH4(y_fixed(:,i));
end
assignin('base', 'ro', ro);

n = length(tspan_fixed); % Number of evaluations
xi_all = zeros(22, n);  % To store xi values for each perturbation

for pert_index = 1:22
    pert(pert_index) = 1 * kf(pert_index);  % Set the perturbation value
    
    % Solve the ODE
    a_fixed = deval(ode15s(@(f, a) pertubations(f, a, pert), tspan, yi), tspan_fixed);  % Use the already solved 'a'
    
    r = zeros(1, n);
    for i = 1:n
        r(i) = r_CH4(a_fixed(:, i));  % Calculate r values
    end
    
    % Store r in the base workspace for the current perturbation
    assignin('base', ['r', num2str(pert_index)], r);
    
    % Calculate xi for the current perturbation
    xi = abs((r - ro) * 1 ./ ro); 
    xi_all(pert_index, :) = xi;  % Store xi in the array
    assignin('base', ['xi', num2str(pert_index)], xi);  % Store xi in the base workspace
    
    % Reset perturbation to 0 after this iteration
    pert(pert_index) = 0;  % Set the current perturbation to zero
    
    % Plot xi values; use hold on to overlay plots
    hold on;
    plot(tspan_fixed, xi, 'LineWidth',1.2);
    text(tspan_fixed(end), xi(end), ['Rxn.' num2str(pert_index)], ...
          'VerticalAlignment', 'bottom', ...
         'HorizontalAlignment', 'left', 'FontSize', 8);
    pert(:) = 0;
end
% Add labels and title for the plot
% Add labels and title for the plot
% Add labels and title for the plot
xlabel('Time');
ylabel('xi values');
title('xi values for Perturbations 1 to 22');

% Create the legend entries
legendEntries = arrayfun(@(x) ['Rxn.' num2str(x)], 1:22, 'UniformOutput', false);

% Create the legend with a specified number of columns
lgd = legend(legendEntries, 'Location', 'best', 'Box', 'off', 'Color', 'none', 'NumColumns', 4); % Adjust the number of columns as needed

% Set the legend to be horizontal (using NumColumns automatically makes it horizontal)
set(lgd, 'Orientation', 'vertical'); % This can be omitted as NumColumns handles the layout
grid on



 %% Plot methane concentration over time (mol/m3)
figure;
plot(T, Y(:,1)/N,'LineWidth',2);
xlabel('Time (s)');
ylabel('Concentration (mol/m^3)');
title('Temporal Variation of Methane Concentration');
grid on

%% Plot Other species concentration over time (mol/m3)
figure;
plot(T, Y(:,10)/N,T, Y(:,7)/N, T,Y(:,8)/N,T, Y(:,14)/N);  % CH4 and p
xlabel('Time (s)');
ylabel('Concentration (mol/m^3)');
title('Temporal Variation of Concentrations of other species');
legend('C_2H_6', 'CH_3ZnCl','HCl','H_2');
grid on;
end




