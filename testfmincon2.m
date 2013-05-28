%a0=1.0e-06*[0.0291; -0.0485; -0.1299; -0.0384; 0.0844; 0.0309; 0.0622; -0.1391];

%options = optimset('Algorithm','active-set','Display','off');
options = optimset('Algorithm','interior-point','Display','off');
%func = @(x) exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
func = @(a) -(a'*N1-0.5*a'*M2*a);  % log likelihood ratio

% C is the constraint matrix
C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];
beq=[0.0; 0.0; 0.0; 0.0];

% a1, 8 by 1 vector
%[~,fval] = fmincon(func,a0,[],[],[],[],[],[],@nonlconst,options);
[a1,fval] = fmincon(func,a0,[],[],C,beq,[],[]); % ,@nonlconst,options);