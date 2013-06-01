% nonlinear constraint multiplier finder function
% Yan Wang, May 27, 2013

function fval=MultiFunc1(mu,M,N,norm1,norm2)

B=diag([1,-1,0,0,1,-1,0,0]);
% C is the constraint matrix
C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];

% D=diag([0,1,0,0,-1,0,0,0]);

%norm1=norm(N,Inf);
%norm2=norm(M,Inf);

% temporary variables, use normalized unit
tmp0=M/norm2-2.0*mu*B;
tmp1=tmp0\C'; 
tmp2=tmp0\(N/norm1);
tmp3=tmp0\B;

lambda=-(C*tmp1)\(C*tmp2);
fval=(N/norm1+C'*lambda)'*tmp3*(tmp2+tmp1*lambda);
%lambda=lambda*norm1;
%fval=fval*(norm1^2)*(norm2^2);

%fval=a'*B*a;

% a2=M2\C';
% 
% a=(M/norm2-mu*B)\(N/norm1);

%mu=mu*norm2;

% fval=fval*(norm1*norm2)^2; we are normalizing the fval, so that the possible
% value of mu is around one 

% END of function