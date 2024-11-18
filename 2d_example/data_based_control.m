function u  = data_based_control(x, Fx, SS, u_SS)
lmbd = sdpvar(size(SS,2),1);
y    = sdpvar(size(SS,1),1);

constr = [lmbd >=0; ones(1, size(SS,2))*lmbd == 1; x==Fx(y);SS*lmbd == y];
cost   = norm(lmbd,1);

options = sdpsettings('verbose',0,'solver','ipopt');
options.ipopt.tol=1e-3;

Problem = optimize(constr,cost,options);

u = u_SS*double(lmbd);

end