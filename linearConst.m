function [LLR,a]=linearConst(N1,M2)

B=diag([1,-1,0,0,1,-1,0,0]);
% C is the constraint matrix
C=[1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0; ...
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0];

% inverse square matrix, warning message is printed if X is badly scaled or nearly singular.
a1=M2\N1;  % \ backslash or matrix left division, more efficient
a2=M2\C';
a=a1-a2*((C*a2)\(C*a1));
fbest = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );
%M2inv=zeros(8,8);
%M2inv=inv(M2);

% LLR with 4 constraints, here is accually LLR = - LLR
if a'*B*a > 0
   %LLR=-0.5*dot(N1,a1);  % LLR without constraints
   %LLR = -0.5 * ( N1'*M2inv*N1 - (C*M2inv*N1)'*(C*M2inv*C')*(C*M2inv*N1) );
   %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*M2inv*C')*(C*a1) );
   %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*(M2\C'))\(C*a1)) );  % maybe most efficient and accurate
   LLR = fbest;
   %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*inv(C*a2)*(C*a1) );
   %LLR = -0.5 * ( dot(N1,a1) - (C*a1)'*(C*a2)\(C*a1) );

else
    
    %disp('Warning: In LogLikelihoodRatio.m: a(1)^2+a(5)^2-a(2)^2-a(6)^2 < 0');
    LLR=fbest;

end

% END of function