function    [D,mc]    =  Nonnegative_DL( X, par )
Q            =    randperm( size(X,2) );
 D            =    X(:, Q(1:par.K));
% D            =    X(:, 1:par.K);
D            =    D./repmat((sqrt(sum(D.^2))+eps), size(D,1), 1);

Iter         =    10;
fun          =    zeros(Iter, 1);
X_s          =    X;
 par.lambda   =    0.001;

for  t   =  1 : Iter

    % Nonnegative Sparse coding        
    B     =    Nonnegative_SC( D, X_s, par );        
    b     =    sum(B.^2, 2);
    R     =    X - D*B;
    
    % update dictionary
    for k = 1 : par.K       
        d_k_pre   =   D(:,k);
         d_k       =   max( D(:,k) + R*B(k,:)'/b(k), 0 );
%  d_k       = D(:,k) + R*B(k,:)'/b(k);
        d_k       =   d_k./max(norm(d_k), 1);
        D(:,k)    =   d_k;
       
        R0        =   R;
        R         =   R0 - (d_k - d_k_pre)*B(k,:);
    end
    
    fun(t)   =   0.5*sum(sum((X_s-D*B).^2)) + par.lambda*sum( sum( abs(B) ) );    
end
