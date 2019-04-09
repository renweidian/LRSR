function  A   =  Nonnegative_SC( D, X, par )

A        =    zeros( size(D,2), size(X, 2) );
V        =    zeros( size(D,2), size(X, 2) );
T        =    100;
fun      =    zeros(T,1);
DTD       =   D'*D;
DTX       =   D'*X;
Ek        =   eye(par.K);
mu        =   0.01;
ro        =   1.07;
%  aaa=(DTD + mu*Ek)\eye(par.K);
for  i  =  1:T
    S         =   (DTD + mu*Ek)\(DTX + mu*(A-V/(2*mu)) );
     A         =   max( soft(S+V/(2*mu), par.lambda/(2*mu)), 0);
% A         =    soft(S+V/(2*mu), par.lambda/(2*mu));
    V         =   V + 2*mu*( S - A );
    mu        =   mu*ro;

    fun(i)    =   0.5*sum(sum((X-D*A).^2)) + par.lambda*sum(sum(abs(A)));
end
