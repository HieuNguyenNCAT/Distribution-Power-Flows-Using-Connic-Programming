%Implement SCED in CVX 
%this is for MATPOWR casefiles in https://github.com/power-grid-lib/pglib-opf
% run this one [Pd, Pgmax, Pgmin, Gen2Bus, GenVarCost0, GenVarCost, Bij, MVAlim, A]=Read_Matpower(pglib_opf_case5_pjm)
function [opt_obj, dcflow, Pg, lmp]= SCED(Pd, Pgmax, Pgmin, Gen2Bus, GenVarCost0, GenVarCost, Bij, MVAlim, A)
%if case_file ==''
B= (Bij.*A);
[L, I] =size(A);
Gen2Bus=double(Gen2Bus);
G= length(Pgmax);
 cvx_begin 
    cvx_precision best
    cvx_solver gurobi_5
    variable obj 
    variable Powerflow(L)
    variable theta(I)
    variable VOLL(I) nonnegative
    variable Pg(G) nonnegative
    dual variable lmp
    minimize obj
    subject to
%        flow_eq(l,i,j)$((ord(i) = from_bus(l)) and (ord(j)=to_bus(l))).. Flow(l) =e= 1/Xl(l)*(theta(i)-theta(j));
        Powerflow == B*theta
        Powerflow <= MVAlim
        Powerflow >= -MVAlim
%-sum(g$(Gen2Bus(g) = ord(i)), Pg(g)) +sum(l$(to_bus(l)= ord(i)), Flow(l)) -sum(l$(from_bus(l)= ord(i)), Flow(l))+VOLL(i) =e= Pd(i) ;
        Gen2Bus'*Pg - A'*Powerflow + VOLL  == Pd : lmp
        Pg <= Pgmax
        Pg >= Pgmin
        obj == sum(GenVarCost0)+ GenVarCost'*Pg +1000*sum(VOLL)
 cvx_end
 opt_obj = obj;
 dcflow = Powerflow;
 Pg = Pg;
 lmp = lmp;
end 

 
 
 
 
 
 