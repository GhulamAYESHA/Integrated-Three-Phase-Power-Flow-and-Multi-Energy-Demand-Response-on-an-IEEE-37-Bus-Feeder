function [Va,Vb,Vc] = wang2017_lpf_task1(Y3N,Pa,Pb,Pc,Qa,Qb,Qc,Vslack,lambda)
% Linear three-phase power flow in the spirit of Wang et al. (2017).
% We build J_PV and J_QV from the 3N x 3N Y-bus and solve
%   (J_PV + λ J_QV) (|V| - 1) = - (P + λ Q)

N = numel(Pa);

G = real(Y3N);
B = imag(Y3N);

ia = 1:N;
ib = N+1:2*N;
ic = 2*N+1:3*N;

% Block partitions
GA = G(ia,ia); GB = G(ib,ib); GC = G(ic,ic);
BA = B(ia,ia); BB = B(ib,ib); BC = B(ic,ic);

GAB = G(ia,ib); GAC = G(ia,ic);
GBA = G(ib,ia); GBC = G(ib,ic);
GCA = G(ic,ia); GCB = G(ic,ib);

BAB = B(ia,ib); BAC = B(ia,ic);
BBA = B(ib,ia); BBC = B(ib,ic);
BCA = B(ic,ia); BCB = B(ic,ib);

rt3 = sqrt(3);

%% J_PV blocks (active power sensitivity to |V|)
JAAp = GA;
JBBp = GB;
JCCp = GC;

JABp = -0.5*GAB + (rt3/2)*BAB;
JACp = -0.5*GAC - (rt3/2)*BAC;
JBAp = -0.5*GBA - (rt3/2)*BBA;
JBCp = -0.5*GBC + (rt3/2)*BBC;
JCAp = -0.5*GCA + (rt3/2)*BCA;
JCBp = -0.5*GCB - (rt3/2)*BCB;

%% J_QV blocks (reactive power sensitivity to |V|)
JAAq = -BA;
JBBq = -BB;
JCCq = -BC;

JABq = -0.5*BAB - (rt3/2)*GAB;
JACq = -0.5*BAC + (rt3/2)*GAC;
JBAq = -0.5*BBA + (rt3/2)*GBA;
JBCq = -0.5*BBC - (rt3/2)*GBC;
JCAq = -0.5*BCA - (rt3/2)*GCA;
JCBq = -0.5*BCB + (rt3/2)*GCB;

% Assemble full matrices
JPV = [JAAp JABp JACp;
       JBAp JBBp JBCp;
       JCAp JCBp JCCp];

JQV = [JAAq JABq JACq;
       JBAq JBBq JBCq;
       JCAq JCBq JCCq];

M = JPV + lambda*JQV;

P_all = [Pa;Pb;Pc];
Q_all = [Qa;Qb;Qc];

% Linearized equation around |V| = 1 p.u.
rhs = - (P_all + lambda*Q_all) + M*ones(3*N,1);

%% Impose slack magnitude on phases A1, B1, C1
slack_idx = [1, N+1, 2*N+1];

for k = slack_idx
    M(k,:) = 0;
    M(k,k) = 1;
    rhs(k) = Vslack;
end

Vmag = M \ rhs;

Va = Vmag(ia);
Vb = Vmag(ib);
Vc = Vmag(ic);
end
