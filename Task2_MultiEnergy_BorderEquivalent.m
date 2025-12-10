function Task2_MultiEnergy_BorderEquivalent()
% ================================================================
% TASK 2: Border Equivalent Method for Multi-Energy Systems
% --------------------------------------------------------
% - Models a single energy hub coupled to the ADN.
% - Includes: electricity import/export, gas boiler, CHP, thermal storage.
% - Derives a linear feasible region (border equivalent) and
%   aggregated bidding curve via vertex-search style LP.

% ================================================================

clc; close all;
fprintf("=============================================\n");
fprintf("  TASK 2 – BORDER EQUIVALENT MULTI-ENERGY HUB\n");
fprintf("=============================================\n\n");

%% 1) System parameters (physically reasonable, no tuning to plots)
T       = 1;          % single representative hour at border
Pel_max = 2.5;        % [MW] max electric export to ADN
Pel_min = -1.0;       % [MW] min (i.e. can import up to 1 MW from ADN)

% Local heat demand in the building / district (service to be met)
Qdemand = 3.0;        % [MW_th]

% Device capacities
P_CHP_el_max = 2.0;   % [MW_el] CHP electric capacity
Q_CHP_th_max = 3.5;   % [MW_th] CHP thermal capacity
Q_boiler_max = 4.0;   % [MW_th] gas boiler capacity
Q_storage_max = 2.0;  % [MW_th] max discharge from thermal storage

% Efficiencies (simple static values)
eta_CHP_el = 0.35;    % electric efficiency of CHP (el / gas_in)
eta_CHP_th = 0.45;    % thermal efficiency of CHP (th / gas_in)
eta_boiler = 0.9;     % boiler efficiency (th / gas_in)

% Cost coefficients (per MWh)
c_el_buy  = 80;       % [£/MWh] buying electricity from upstream grid
c_el_sell = 60;       % [£/MWh] selling electricity to ADN (bid price)
c_gas     = 30;       % [£/MWh_gas] gas price
c_storage = 5;        % [£/MWh_th] cycling cost of thermal storage

%% 2) Decision variables for the hub (one hour)
% x = [Pel_export; Pel_grid; F_gas_CHP; F_gas_boiler; Q_storage];
% Pel_export : net electric power exported to ADN (>0 export, <0 import)
% Pel_grid   : power bought from upstream external grid (always >= 0)
% F_gas_*    : gas energy input [MW_gas]
% Q_storage  : positive = discharge thermal storage

nvar = 5;

% Bounds
lb = [Pel_min;  0;           0;                0;              0           ];
ub = [Pel_max;  inf;         inf;              inf;            Q_storage_max];

%% 3) Linear constraints A*x <= b and Aeq*x = beq

% --- Electric power balance at the hub bus ---
% Pel_export + Pel_load_local = Pel_CHP + Pel_grid
% Assume local *electric* demand in the building is fixed:
Pel_local = 1.0;  % [MW_el]

% CHP electric output:
% P_CHP_el = eta_CHP_el * F_gas_CHP
% So we enforce capacity and relation via inequalities:
%    0 <= eta_CHP_el * F_gas_CHP <= P_CHP_el_max

% For thermal side:
% Q_th = eta_CHP_th * F_gas_CHP + eta_boiler * F_gas_boiler + Q_storage
% must satisfy Q_th >= Qdemand.

% We'll enforce relations as linear inequalities by introducing them
% directly in A, Aeq.

A   = [];
b   = [];
Aeq = [];
beq = [];

% ---- 3a) Electric capacity of CHP: eta_CHP_el * F_gas_CHP <= P_CHP_el_max
a = zeros(1,nvar);
a(3) =  eta_CHP_el;
A = [A; a];
b = [b; P_CHP_el_max];

% ---- 3b) Thermal capacity of CHP: eta_CHP_th * F_gas_CHP <= Q_CHP_th_max
a = zeros(1,nvar);
a(3) =  eta_CHP_th;
A = [A; a];
b = [b; Q_CHP_th_max];

% ---- 3c) Electric balance at the border:
% Pel_export + Pel_local = eta_CHP_el*F_gas_CHP + Pel_grid
a = zeros(1,nvar);
% left-hand side moved to RHS -> equality
% Pel_export - P_CHP_el - Pel_grid + Pel_local = 0
a(1) =  1;              % Pel_export
a(2) = -1;              % -Pel_grid
a(3) = -eta_CHP_el;     % -eta_CHP_el * F_gas_CHP
Aeq = [Aeq; a];
beq = [beq; -Pel_local];

% ---- 3d) Thermal demand: Q_th >= Qdemand
% Q_th = eta_CHP_th*F_gas_CHP + eta_boiler*F_gas_boiler + Q_storage
% => -Q_th <= -Qdemand  for A*x <= b form:
a = zeros(1,nvar);
a(3) = -eta_CHP_th;
a(4) = -eta_boiler;
a(5) = -1;
A = [A; a];
b = [b; -Qdemand];

%% 4) Define cost function (used later for bidding curve)
% Operating cost in £/h (assuming 1-hour block)
% cost = c_el_buy*Pel_grid + c_gas*(F_gas_CHP + F_gas_boiler) + c_storage*Q_storage

c = zeros(nvar,1);
c(2) = c_el_buy;
c(3) = c_gas;
c(4) = c_gas;
c(5) = c_storage;

%% 5) Compute the feasible region via vertex search
fprintf("Computing feasible region vertices via LP...\n");

% We are mainly interested in the projection on Pel_export vs. total heat
% Actually the thermal output is completely tied to constraints, so we can
% just record Pel_export and Q_th for each extreme point.

directions = eye(nvar);                % basis directions
directions = [directions, -directions]; % positive & negative

V_points = [];

optsLP = optimoptions('linprog','Display','none','Algorithm','dual-simplex');

for k = 1:size(directions,2)
    f_dir = directions(:,k);
    [x_ext,~,exitflag] = linprog(-f_dir,A,b,Aeq,beq,lb,ub,optsLP); % max in direction
    if exitflag == 1
        V_points = [V_points, x_ext];
    end
end

% Remove duplicates (within tolerance)
V_points = unique(round(V_points.',6),'rows').';

Pel_vert = V_points(1,:);  % Pel_export
Qth_vert = eta_CHP_th*V_points(3,:) + eta_boiler*V_points(4,:) + V_points(5,:);

fprintf("  Found %d distinct vertices in projected space.\n\n", numel(Pel_vert));

%% 6) Aggregated bidding curve: min operating cost for each Pel_export

fprintf("Computing aggregated bidding curve...\n");

Pel_grid_vec = linspace(Pel_min,Pel_max,25);  % candidate export levels
cost_vec     = nan(size(Pel_grid_vec));

for i = 1:numel(Pel_grid_vec)
    Pel_set = Pel_grid_vec(i);

    % Fix Pel_export = Pel_set with an extra equality:
    Aeq2  = [Aeq; zeros(1,nvar)];
    Aeq2(end,1) = 1;
    beq2  = [beq; Pel_set];

    [x_opt,f_opt,exitflag] = linprog(c,A,b,Aeq2,beq2,lb,ub,optsLP);
    if exitflag == 1
        cost_vec(i) = f_opt;
    end
end

%% 7) Plot feasible region and bidding curve

figure('Color','w','Position',[80 80 1050 480]);

subplot(1,2,1);
hold on; grid on; box on;
scatter(Pel_vert,Qth_vert,50,'filled');
xlabel('P_{el, export} to ADN [MW]');
ylabel('Q_{th} supplied [MW]');
title('Task 2: Multi-Energy Feasible Region (Border Equivalent)');
set(gca,'FontSize',11);

subplot(1,2,2);
plot(Pel_grid_vec,cost_vec,'-o','LineWidth',1.8,'MarkerSize',5);
grid on; box on;
xlabel('P_{el, export} to ADN [MW]');
ylabel('Operating Cost [£/h]');
title('Task 2: Aggregated Bidding Curve of Multi-Energy Hub');
set(gca,'FontSize',11);

fprintf("Task 2 completed – feasible region and bidding curve plotted.\n");
fprintf("Use these figures in your report as the 'border equivalent'\n");
fprintf("and 'aggregated bidding curve' for the multi-energy system.\n\n");
end
