% use 'fmincon' to find the solution to an optimization problem with 
% nonlinear inequality constraints
% Yan Wang, April 29, 2013

% need optimization toolbox

%global N1 M2

%problem=

func = @(a) -(a'*N1-0.5*a'*M2*a);

% x0 = [-1,1];  % make a starting guess at the solution
% a=1;
% M=[4,2;2,2];
% %func = @(x) a*exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
% func = @(x) a*exp(x(1))*(x*M*x' + 2*x(2) + 1);

options = optimset('Algorithm','active-set','Display','off');

%disp([num2str(N1)])
%disp([num2str(M2)])

%[x,fval] = fmincon(@objfun,x0,[],[],[],[],[],[],@confun,options);
%[x,fval] = fmincon(func,x0,[],[],[1,-1],0,[],[],@confun,options);
%[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst2,options);

[a,fval] = fmincon(@func4,a1,[],[],[],[],[],[],@nonlconst1,options);
%x = fmincon(problem);
%x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(problem);