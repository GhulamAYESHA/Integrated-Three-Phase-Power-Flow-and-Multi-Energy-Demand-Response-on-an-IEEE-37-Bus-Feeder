
function results = Task3_ThreePhase_IDR_OPF()
% TASK 3 – Integrated Demand Response OPF for Three-Phase Unbalance
% ----------------------------------------------------------------
% This script implements an Integrated Demand Response (IDR)
% strategy for a *single customer* connected to the end of the
% IEEE-37 feeder used in Task-1, with:
%   • multi-energy consumption (electric, heat, gas)
%   • unified value / cost model
%   • explicit three-phase unbalance penalty
%   • network constraints enforced via the 3-phase AC power flow
%
% It:
%   1) Loads IEEE37_News.xlsx (same as Task-1) and builds Ybus_3N.
%   2) Adds small flexible loads at the last bus.
%   3) Optimises 24-h schedules P_elec, P_heat, P_gas to minimise
%      energy cost + unbalance, subject to:
%        - time-of-use tariffs
%        - daily energy requirements
%        - device power limits
%        - voltage magnitude limits.
%   4) Returns and plots the optimal schedules and unbalance index.

clc; close all;
fprintf('\n=== TASK 3: Three-Phase IDR OPF ===\n');

%% 1. Network & base load from Task-1 data -------------------------------
xlsxPath     = 'IEEE37_News.xlsx';    % your file in the same folder
Vbase_LL_kV  = 12.66;
Sbase_MVA    = 10;
Vslack_pu    = 1.0;

[from,to,Z,loads,N] = read_ieee_excel_3ph(xlsxPath);
Y3N = build_ybus_3N(N, from, to, Z, Vbase_LL_kV, Sbase_MVA);
fprintf('Loaded IEEE-37 network with %d buses.\n',N);

% base loads per branch (kW / kvar), then mapped to buses (exactly as in Task-1)
[P0a_raw,P0b_raw,P0c_raw,Q0a_raw,Q0b_raw,Q0c_raw] = loads.PQperPhase_kWkvar{:};

P0a = zeros(N,1); P0b = zeros(N,1); P0c = zeros(N,1);
Q0a = zeros(N,1); Q0b = zeros(N,1); Q0c = zeros(N,1);
for k = 1:numel(to)
    bus = to(k);
    if bus>=1 && bus<=N
        P0a(bus) = P0a(bus) + P0a_raw(k);
        P0b(bus) = P0b(bus) + P0b_raw(k);
        P0c(bus) = P0c(bus) + P0c_raw(k);
        Q0a(bus) = Q0a(bus) + Q0a_raw(k);
        Q0b(bus) = Q0b(bus) + Q0b_raw(k);
        Q0c(bus) = Q0c(bus) + Q0c_raw(k);
    end
end

Sbase_kVA = Sbase_MVA*1000;
P0a = P0a/Sbase_kVA; P0b = P0b/Sbase_kVA; P0c = P0c/Sbase_kVA;
Q0a = Q0a/Sbase_kVA; Q0b = Q0b/Sbase_kVA; Q0c = Q0c/Sbase_kVA;

%% 2. Time horizon and tariffs (IDR set-up) ------------------------------
T = 24;                    % hours
h = (1:T)';                % hour index

% Time-of-use tariffs [€/kWh] – realistic and not "tuned"
price_e = 0.10*ones(T,1);             % off-peak base
price_e(8:17)  = 0.16;                % shoulder 08–17
price_e(18:22) = 0.28;                % evening peak

% For heat & gas we can set slightly lower / flatter prices
price_h = 0.06*ones(T,1);
price_h(18:22) = 0.09;

price_g = 0.05*ones(T,1);
price_g(18:22) = 0.08;

% Daily energy requirements for the flexible customer [kWh]
E_day_elec = 30;   % e.g., appliances, lighting, EV base load
E_day_heat = 40;   % space heating / DHW
E_day_gas  = 20;   % cooking / CHP equivalent

% Convert kWh to per-unit power-hours on the 10 MVA base
E_day_elec_pu = E_day_elec / (Sbase_MVA*1000);  % p.u.-h
E_day_heat_pu = E_day_heat / (Sbase_MVA*1000);
E_day_gas_pu  = E_day_gas  / (Sbase_MVA*1000);

% Per-hour power bounds of flexible loads [p.u.]
P_elec_max_pu = 0.004;  % 40 kW on 10 MVA base
P_heat_max_pu = 0.005;  % 50 kW
P_gas_max_pu  = 0.003;  % 30 kW

%% 3. Phase mapping for the customer (unbalanced) ------------------------
% Customer is connected at the last bus:
cust_bus = N;

% We distribute each energy form unevenly over the three phases to
% create realistic three-phase unbalance.
%   - Electric: A-heavy
%   - Heat:     B-heavy (e.g., heat pump on phase B)
%   - Gas / CHP: C-heavy
W_e = [0.6; 0.3; 0.1];   % sum = 1
W_h = [0.2; 0.6; 0.2];
W_g = [0.2; 0.2; 0.6];

%% 4. Optimisation variables and linear constraints ----------------------
% Decision variable vector:
%   x = [P_elec(1:T); P_heat(1:T); P_gas(1:T)]  (3T x 1)
nx = 3*T;

% Bounds
lb = zeros(nx,1);
ub = [ P_elec_max_pu*ones(T,1);
       P_heat_max_pu*ones(T,1);
       P_gas_max_pu*ones(T,1) ];

% Daily energy constraints: sum_t P(t)*dt = E_day  (dt = 1 h)
Aeq = zeros(3, nx);
beq = [E_day_elec_pu; E_day_heat_pu; E_day_gas_pu];

% Elec block
Aeq(1,1:T) = 1;
% Heat block
Aeq(2,T+1:2*T) = 1;
% Gas block
Aeq(3,2*T+1:3*T) = 1;

%% 5. Objective and non-linear constraints (nested functions) ------------
w_unb = 5.0;        % weight on unbalance index (tunable but fixed)
Vmin  = 0.95;       % voltage magnitude limits
Vmax  = 1.05;

% Initial guess: flat schedule meeting energy constraints
x0 = zeros(nx,1);
x0(1:T)         = E_day_elec_pu/T;
x0(T+1:2*T)     = E_day_heat_pu/T;
x0(2*T+1:3*T)   = E_day_gas_pu/T;

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'MaxIterations',100, ...
    'MaxFunctionEvaluations',5000, ...
    'Display','iter');

obj  = @(x) IDR_objective(x, T, price_e, price_h, price_g, ...
                          Y3N, P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
                          cust_bus, W_e,W_h,W_g, ...
                          Vslack_pu, Vmin, Vmax, w_unb);
nonl = @(x) IDR_nonlcon(x, T, ...
                        Y3N, P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
                        cust_bus, W_e,W_h,W_g, ...
                        Vslack_pu, Vmin, Vmax);

fprintf('\nStarting fmincon optimisation...\n');
[x_opt,fval,exitflag,output] = fmincon(obj, x0, [],[], Aeq,beq, lb,ub, nonl, options);

fprintf('\n--- OPTIMISATION FINISHED (exitflag = %d) ---\n', exitflag);
disp(output);

%% 6. Post-processing and plots -----------------------------------------
P_elec = x_opt(1:T);
P_heat = x_opt(T+1:2*T);
P_gas  = x_opt(2*T+1:3*T);

% Evaluate AC PF and unbalance for the final optimum
UI = zeros(T,1);
for t = 1:T
    [Va,Vb,Vc] = apply_customer_load( ...
        P_elec(t),P_heat(t),P_gas(t), ...
        Y3N, P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
        cust_bus,W_e,W_h,W_g, Vslack_pu);
    VmagA = abs(Va); VmagB = abs(Vb); VmagC = abs(Vc);

    % simple voltage unbalance index at customer bus
    vA = VmagA(cust_bus);
    vB = VmagB(cust_bus);
    vC = VmagC(cust_bus);
    vavg = (vA+vB+vC)/3;
    UI(t) = max(abs([vA-vavg, vB-vavg, vC-vavg]))/vavg;
end

% 6.1 Schedules
figure('Color','w','Position',[120 120 900 520]);
plot(h, P_elec*Sbase_kVA,'-','LineWidth',2); hold on;
plot(h, P_heat*Sbase_kVA,'-','LineWidth',2);
plot(h, P_gas *Sbase_kVA,'-','LineWidth',2);
grid on; box on;
xlabel('Hour','FontSize',12,'FontWeight','bold');
ylabel('Power (kW)','FontSize',12,'FontWeight','bold');
title('Optimal Customer Multi-Energy Schedule','FontSize',14,'FontWeight','bold');
legend({'Electric','Heat','Gas'},'Location','best');

% 6.2 Unbalance index
figure('Color','w','Position',[150 150 900 520]);
plot(h, UI,'-o','LineWidth',2); grid on; box on;
xlabel('Hour','FontSize',12,'FontWeight','bold');
ylabel('Voltage unbalance index (p.u.)','FontSize',12,'FontWeight','bold');
title('Three-Phase Unbalance at Customer Bus After IDR','FontSize',14,'FontWeight','bold');

% Return results structure
results.P_elec = P_elec;
results.P_heat = P_heat;
results.P_gas  = P_gas;
results.UI     = UI;
results.cost   = fval;
assignin('base','Task3_results',results);

end % main function

% -------------------------------------------------------------------------
% OBJECTIVE FUNCTION: energy cost + unbalance penalty
% -------------------------------------------------------------------------
function J = IDR_objective(x, T, price_e,price_h,price_g, ...
                           Y3N,P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
                           cust_bus,W_e,W_h,W_g, ...
                           Vslack_pu,Vmin,Vmax, w_unb)

P_elec = x(1:T);
P_heat = x(T+1:2*T);
P_gas  = x(2*T+1:3*T);

Sbase_kVA = 10*1000;

% Energy cost term (sum over hours)
cost = sum( price_e.*(P_elec*Sbase_kVA/1000) + ...
            price_h.*(P_heat*Sbase_kVA/1000) + ...
            price_g.*(P_gas *Sbase_kVA/1000) );

% Network-based unbalance penalty
unb_sum = 0;
for t = 1:T
    [Va,Vb,Vc] = apply_customer_load( ...
        P_elec(t),P_heat(t),P_gas(t), ...
        Y3N, P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
        cust_bus,W_e,W_h,W_g, Vslack_pu);
    VmagA = abs(Va); VmagB = abs(Vb); VmagC = abs(Vc);

    vA = VmagA(cust_bus);
    vB = VmagB(cust_bus);
    vC = VmagC(cust_bus);
    vavg = (vA+vB+vC)/3;
    UI = max(abs([vA-vavg, vB-vavg, vC-vavg]))/vavg;

    % also softly penalise voltage outside [Vmin,Vmax]
    v_viol = max([0, Vmin-vA, Vmin-vB, Vmin-vC, ...
                     vA-Vmax, vB-Vmax, vC-Vmax]);
    unb_sum = unb_sum + UI^2 + 10*v_viol^2;
end

J = cost + w_unb*unb_sum;

end

% -------------------------------------------------------------------------
% NON-LINEAR CONSTRAINTS: keep voltages within limits (hard)
% -------------------------------------------------------------------------
function [c,ceq] = IDR_nonlcon(x, T, ...
                               Y3N,P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
                               cust_bus,W_e,W_h,W_g, ...
                               Vslack_pu, Vmin, Vmax)
P_elec = x(1:T);
P_heat = x(T+1:2*T);
P_gas  = x(2*T+1:3*T);

c = [];  % inequality c(x)<=0
for t = 1:T
    [Va,Vb,Vc] = apply_customer_load( ...
        P_elec(t),P_heat(t),P_gas(t), ...
        Y3N, P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
        cust_bus,W_e,W_h,W_g, Vslack_pu);
    VmagA = abs(Va); VmagB = abs(Vb); VmagC = abs(Vc);
    vA = VmagA(cust_bus);
    vB = VmagB(cust_bus);
    vC = VmagC(cust_bus);
    c = [c;
         Vmin - vA;
         Vmin - vB;
         Vmin - vC;
         vA - Vmax;
         vB - Vmax;
         vC - Vmax];
end
ceq = [];

end

% -------------------------------------------------------------------------
% Helper: inject customer load and run 3-phase AC PF
% -------------------------------------------------------------------------
function [Va,Vb,Vc] = apply_customer_load( ...
    P_elec,P_heat,P_gas, ...
    Y3N, P0a,P0b,P0c,Q0a,Q0b,Q0c, ...
    cust_bus,W_e,W_h,W_g,Vslack_pu)

% copy base loads
Pa = P0a; Pb = P0b; Pc = P0c;
Qa = Q0a; Qb = Q0b; Qc = Q0c;

% assume constant power factor cosφ = 0.95 for additional loads
pf = 0.95; q_factor = sqrt(1/pf^2 - 1);

% phase-wise extra active power
dP_e = W_e*P_elec;
dP_h = W_h*P_heat;
dP_g = W_g*P_gas;

dP = dP_e + dP_h + dP_g;  % [3x1] phases A,B,C at cust_bus

Pa(cust_bus) = Pa(cust_bus) + dP(1);
Pb(cust_bus) = Pb(cust_bus) + dP(2);
Pc(cust_bus) = Pc(cust_bus) + dP(3);

% matching reactive powers for the flexible loads
dQ = q_factor * dP;
Qa(cust_bus) = Qa(cust_bus) + dQ(1);
Qb(cust_bus) = Qb(cust_bus) + dQ(2);
Qc(cust_bus) = Qc(cust_bus) + dQ(3);

[Va,Vb,Vc] = acpf3N_gs(Y3N,Pa,Pb,Pc,Qa,Qb,Qc,1,Vslack_pu);

end
