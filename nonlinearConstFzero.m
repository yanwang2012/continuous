% Semi-analytic calculation of the nonlinear constraints
% Yan Wang, May 29, 2013

function [LLR,a]=nonlinearConstFzero(N1,M2)

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
%abest=a;
fbest = -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );
Ampfunc = @ (a) 0.5*(sqrt(a(1)^2+a(5)^2)+sqrt(a(1)^2+a(5)^2-a(2)^2-a(6)^2));

if a'*B*a > 0
    Amp=Ampfunc(a);
    
    if a(2)^2-4*Amp^2*a(5)^2/(a(1)^2+a(5)^2) <= 0
        
        LLR = fbest; % -0.5 * ( dot(N1,a1) - (C*a1)'*((C*a2)\(C*a1)) );  % maybe most efficient and accurate
        
    else
        
        
    end
    
else
    
    norm1=norm(N1,Inf);  % maximum element in the vector
    norm2=norm(M2,Inf);  % maximum element in the matrix
    lb=-5.0;  % lower boundary in normalized unit
    ub=5.0;  % upper boundary
    xzero=[];  % variable contains the roots of function in normalized unit
    ftmp=[];  % -log likelihood ratio
    atmp=[];
    j=0;  % number of roots
    mul=0;  % multiplier in true unit

    for i=lb:0.02:ub  % here 'i' is not an integer
        if MultiFunc1(i-0.01,M2,N1,norm1,norm2) * MultiFunc1(i+0.01,M2,N1,norm1,norm2) < 0
            % here is one root in the interval
            j=j+1;
            xzero(j)=fzero(@(x)MultiFunc1(x,M2,N1,norm1,norm2),[i-0.01 i+0.01]);
        end
    end
    
    if j>0
        for k=1:1:j
            mul=xzero(k);  % still in normalized unit
            lam=-(C*(M/norm2-2.0*mul*B)\C')\(C*(M/norm2-2.0*mul*B)\(N/norm1));
            atmp(8,k)=(M/norm2-2.0*mul*B)\(N/norm1+C'*lam)*(norm1*norm2);  % true unit
            ftmp(k)=-(atmp(k)'*N-0.5*atmp(k)'*M*atmp(k));  % true unit
        end
    end
    
    [LLR,I]=min(ftmp);  % minimum value
    a=atmp(I);

end

% END of function