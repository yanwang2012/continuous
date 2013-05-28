function [LLR,a]=nonlinearConst(N1,M2)

B=diag([1,-1,0,0,1,-1,0,0]);
% C is the constraint matrix
C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];

% inverse square matrix, warning message is printed if X is badly scaled or nearly singular.
a1=M2\N1;  % tmp variable, \ backslash or matrix left division, more efficient
a2=M2\C';
% a: Lagrangian multiplier solution with only 4 linear equality constraints, 
% ref eq.5 of my note on 4/23/13
a=a1-a2*((C*a2)\(C*a1));
abest=a;
fbest = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );

% LLR with 4 linear constraints, and additional nonlinear constraints are 
% subject if certain condition is violated. Here accually LLR = - LLR

% setup for local optimization 'fmincon'
%a0 = a;  % make a starting guess at the solution
%beq=[0.0; 0.0; 0.0; 0.0];  % right hand side of 4 linear constraints
%options = optimset('Algorithm','active-set','Display','off');
options = optimset('Algorithm','interior-point','Display','off'); %, ...
    %'Hessian','user-supplied','HessFcn',@hessianfcn);
%func = @(x) exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
%func = @(a) -(a'*N1-0.5*a'*M2*a);  % min log likelihood ratio
Ampfunc = @ (a) 0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));
a0=zeros(4,1); % initial value for nonlinear constraint quadratic problem
penalty=0.9;  % penalty factor when every trial is failed


if a'*B*a > 0
    Amp=Ampfunc(a);
    
    if a(2)^2-4*Amp^2*a(5)^2/(a(1)^2+a(5)^2) <= 0
        % CASE 1: linear constraits and both nonlinear constraints are satisfied
        
        %LLR=-0.5*dot(N1,a1);  % LLR without constraints
        %LLR = -0.5 * ( N1'*M2inv*N1 - (C*M2inv*N1)'*(C*M2inv*C')*(C*M2inv*N1) );
        %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*M2inv*C')*(C*a1) );
        LLR = fbest; % -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );  % maybe most efficient and accurate
        %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*a2)*(C*a1) );
        %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*(C*a2)\(C*a1) );
    else
        % CASE 2: second nonlinear constraint is violated
        %a0=a; %+a.*(0.1*rand(8,1)-0.05);  % make a starting guess at the solution
        a0(1)=a(1);
        a0(2)=a(2);
        a0(3)=a(5);
        a0(4)=a(6);
        
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst2,options);
        [y,fval] = fmincon(@(x) func4(x,N1,M2),a0,[],[],[],[],[],[],@nonlconst2,options);
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst1,options);
        a(1)=y(1);
        a(2)=y(2);
        a(3)=-y(1);
        a(4)=y(2);
        a(5)=y(3);
        a(6)=y(4);
        a(7)=-y(3);
        a(8)=y(4);
        
        if a'*B*a >= 0 && a(2)^2-4*Ampfunc(a)^2*a(5)^2/(a(1)^2+a(5)^2) <= 0
            LLR=fval;
        else
            %LLR=Inf;  % may be not able to get a decent solution
            %LLR=0.0;
            a=abest;
            LLR=fbest*penalty; %fval;
            %disp('In LoglikelihoodRatioNonc.m: Inf occurs at 1');
        end
        
    end
    
else
    
    %disp('Warning: In LogLikelihoodRatio.m: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
    % CASE 3:
    %a0=a+a.*(0.1*rand(8,1)-0.05);
    a0(1)=a(1);
    a0(2)=a(2);
    a0(3)=a(5);
    a0(4)=a(6);

    %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst1,options);
    [y,fval] = fmincon(@(x) func4(x,N1,M2),a0,[],[],[],[],[],[],@nonlconst1,options);
    a(1)=y(1);
    a(2)=y(2);
    a(3)=-y(1);
    a(4)=y(2);
    a(5)=y(3);
    a(6)=y(4);
    a(7)=-y(3);
    a(8)=y(4);
    
    % fmincon may be not accurate, ...
    
    if a'*B*a >= 0 && a(2)^2-4*Ampfunc(a)^2*a(5)^2/(a(1)^2+a(5)^2) <= 0 % >
        LLR=fval;
    else
        %a0=a+a.*(0.1*rand(8,1)-0.05);
        a0(1)=a(1);
        a0(2)=a(2);
        a0(3)=a(5);
        a0(4)=a(6);
        
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst2,options);
        %[a,fval] = fmincon(func,a0,[],[],C,beq,[],[],@nonlconst1,options);
        [y,fval] = fmincon(@(x) func4(x,N1,M2),a0,[],[],[],[],[],[],@nonlconst2,options);
        a(1)=y(1);
        a(2)=y(2);
        a(3)=-y(1);
        a(4)=y(2);
        a(5)=y(3);
        a(6)=y(4);
        a(7)=-y(3);
        a(8)=y(4);
        
        if a'*B*a >= 0 && a(2)^2-4*Ampfunc(a)^2*a(5)^2/(a(1)^2+a(5)^2) <= 0
            LLR=fval;
        else
            %LLR=Inf;  % may be not able to get a decent solution
            %LLR=0.0;
            a=abest;
            LLR=fbest*penalty;  %fval;
            %disp('In LoglikelihoodRatioNonc.m: Inf occurs at 2');

        end
        
    end
    
end

% END of function