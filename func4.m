% function for the maximum log likelihood ratio, with 4 linear equality
% constraints included
% Yan Wang, May 6, 2013

function func=func4(a,N1,M2)
%function func=func4(a,N1,M2,norm1,norm2)

% as is 8 by 1 vector, a is 4 by 1 vector
as(1,1)=a(1,1);
as(2,1)=a(2,1);
as(3,1)=-a(1,1);
as(4,1)=a(2,1);
as(5,1)=a(3,1);
as(6,1)=a(4,1);
as(7,1)=-a(3,1);
as(8,1)=a(4,1);

func = -(as'*N1-0.5*as'*M2*as);

% END OF FUNCTION