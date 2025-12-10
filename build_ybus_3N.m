function Y3N = build_ybus_3N(N, from, to, Z, Vbase_LL_kV, Sbase_MVA)
% Builds the 3N×3N Ybus (ordering: [A1..AN, B1..BN, C1..CN]).
% Series-only stamping (shunts can be added if you have them).

% Per-unit bases
Vbase_ph_kV = Vbase_LL_kV/sqrt(3);
Zbase_ohm   = (Vbase_ph_kV*1e3)^2 / (Sbase_MVA*1e6/3); % per phase
toPU = @(x) x ./ Zbase_ohm;

Yaa = sparse(N,N); Ybb = sparse(N,N); Ycc = sparse(N,N);
Yab = sparse(N,N); Yac = sparse(N,N); Ybc = sparse(N,N);

E = numel(from);
for e = 1:E
    i = from(e); j = to(e);
    
    % --- 1. Build 3x3 Impedance Matrix Ze in PU ---
    Raa_pu = toPU(Z.Raa(e)); Xaa_pu = toPU(Z.Xaa(e));
    Rbb_pu = toPU(Z.Rbb(e)); Xbb_pu = toPU(Z.Xbb(e));
    Rcc_pu = toPU(Z.Rcc(e)); Xcc_pu = toPU(Z.Xcc(e));
    Rab_pu = toPU(Z.Rab(e)); Xab_pu = toPU(Z.Xab(e));
    Rac_pu = toPU(Z.Rac(e)); Xac_pu = toPU(Z.Xac(e));
    Rbc_pu = toPU(Z.Rbc(e)); Xbc_pu = toPU(Z.Xbc(e));
    
    Z_mat_e = [ (Raa_pu + 1j*Xaa_pu), (Rab_pu + 1j*Xab_pu), (Rac_pu + 1j*Xac_pu);
                (Rab_pu + 1j*Xab_pu), (Rbb_pu + 1j*Xbb_pu), (Rbc_pu + 1j*Xbc_pu);
                (Rac_pu + 1j*Xac_pu), (Rbc_pu + 1j*Xbc_pu), (Rcc_pu + 1j*Xcc_pu) ];
    
    % --- 2. Calculate 3x3 Admittance Matrix Ye (Matrix Inverse) ---
    Y_e = inv(Z_mat_e);

    % Extract the nine elements of the branch admittance matrix Ye
    Yaa_e = Y_e(1,1); Yab_e = Y_e(1,2); Yac_e = Y_e(1,3);
    Yba_e = Y_e(2,1); Ybb_e = Y_e(2,2); Ybc_e = Y_e(2,3);
    Yca_e = Y_e(3,1); Ycb_e = Y_e(3,2); Ycc_e = Y_e(3,3);

    % --- 3. Stamping the Ybus blocks (3x3 blocks) ---
    % Note: The stamping logic for Ybus remains the standard:
    % Y_ii += Y_e, Y_jj += Y_e, Y_ij -= Y_e, Y_ji -= Y_e
    
    % Diagonal Blocks (Self-Admittances)
    % A phase stamp (Yaa block) - uses Y_e(1,1)
    Yaa(i,i)=Yaa(i,i)+Yaa_e; Yaa(j,j)=Yaa(j,j)+Yaa_e;
    Yaa(i,j)=Yaa(i,j)-Yaa_e; Yaa(j,i)=Yaa(i,j);
    
    % B phase stamp (Ybb block) - uses Y_e(2,2)
    Ybb(i,i)=Ybb(i,i)+Ybb_e; Ybb(j,j)=Ybb(j,j)+Ybb_e;
    Ybb(i,j)=Ybb(i,j)-Ybb_e; Ybb(j,i)=Ybb(i,j);
    
    % C phase stamp (Ycc block) - uses Y_e(3,3)
    Ycc(i,i)=Ycc(i,i)+Ycc_e; Ycc(j,j)=Ycc(j,j)+Ycc_e;
    Ycc(i,j)=Ycc(i,j)-Ycc_e; Ycc(j,i)=Ycc(i,j);

    % Off-Diagonal Blocks (Mutual-Admittances between phases at bus level)
    % Yab block (A-B coupling) - uses Y_e(1,2)
    Yab(i,i)=Yab(i,i)+Yab_e; Yab(j,j)=Yab(j,j)+Yab_e;
    Yab(i,j)=Yab(i,j)-Yab_e; Yab(j,i)=Yab(i,j);

    % Yac block (A-C coupling) - uses Y_e(1,3)
    Yac(i,i)=Yac(i,i)+Yac_e; Yac(j,j)=Yac(j,j)+Yac_e;
    Yac(i,j)=Yac(i,j)-Yac_e; Yac(j,i)=Yac(i,j);

    % Ybc block (B-C coupling) - uses Y_e(2,3)
    Ybc(i,i)=Ybc(i,i)+Ybc_e; Ybc(j,j)=Ybc(j,j)+Ybc_e;
    Ybc(i,j)=Ybc(i,j)-Ybc_e; Ybc(j,i)=Ybc(i,j);

    % Note: We do not need to explicitly stamp Yba, Yca, Ycb, as the final 
    % assembly handles the Hermitian property: Yba=Yab', Yca=Yac', Ycb=Ybc'.
end

% Assemble 3N×3N, relying on the Hermitian property (Yba = Yab' etc.)
Y3N = [Yaa  Yab  Yac ;
       Yab' Ybb  Ybc ;
       Yac' Ybc' Ycc];
end