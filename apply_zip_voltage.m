function [Pz,Qz] = apply_zip_voltage(P0,Q0,Vmag,FZ,FI,FP)
% ZIP load: P(V)= (FZ V^2 + FI V + FP) P0; same for Q.
V = Vmag(:);
Pz = (FZ*(V.^2) + FI*V + FP).*P0(:);
Qz = (FZ*(V.^2) + FI*V + FP).*Q0(:);
end