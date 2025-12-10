function [Va,Vb,Vc] = acpf3N_gs(Y3N, P0a,P0b,P0c, Q0a,Q0b,Q0c, slack_bus, Vslack)
% 3-phase Gauss–Seidel AC power flow (PQ everywhere except slack A/B/C).
N = numel(P0a); ia = 1:N; ib = N+1:2*N; ic = 2*N+1:3*N;

V = ones(3*N,1);
V(ib) = exp(-1j*2*pi/3);   % -120°
V(ic) = exp( 1j*2*pi/3);   % +120°
V( ia(1) ) = Vslack*exp(1j*0);
V( ib(1) ) = Vslack*exp(-1j*2*pi/3);
V( ic(1) ) = Vslack*exp( 1j*2*pi/3);

Sload = -[P0a+1j*Q0a; P0b+1j*Q0b; P0c+1j*Q0c]; % negative injections
Sload([ia(1) ib(1) ic(1)]) = 0;

maxit=200; tol=1e-8; PQ=setdiff(1:3*N,[ia(1) ib(1) ic(1)]);
for it=1:maxit
    Vprev=V;
    for k=PQ
        Ykk = Y3N(k,k);
        sumYV = Y3N(k,:)*V - Ykk*V(k);
        V(k) = (1/Ykk) * ( (conj(Sload(k))/conj(V(k))) - sumYV );
    end
    V( ia(1) ) = Vslack*exp(1j*0);
    V( ib(1) ) = Vslack*exp(-1j*2*pi/3);
    V( ic(1) ) = Vslack*exp( 1j*2*pi/3);
    if max(abs(V-Vprev))<tol, break; end
end
Va=V(ia); Vb=V(ib); Vc=V(ic);
end