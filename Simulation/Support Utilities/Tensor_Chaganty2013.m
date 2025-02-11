function Solution = Tensor_Chaganty2013(y,X,Option)

% global N p K p2 p3 pK
[N,p]   = size(X);
p2      = p*p;
p3      = p*p2;

K       = Option.NumMix;
pK      = p*K;

%%
if isfield(Option,'RegressWeight')              % Outlier weight-out method to achieve "robust" regression
    RegressWeight	= Option.RegressWeight;
else
    RegressWeight.Name= '';
end
if isfield(Option,'RegressMethod')                  % Regression method
    RegressMethod   = Option.RegressMethod;
else
    RegressMethod.Name= 'ls';
end
if isfield(Option,'FactorizationMethod')            % Tensor factorization method
    FactorizationMethod = Option.FactorizationMethod;
else
    FactorizationMethod.Name= 'rpsf';
end

RegressWeight.Name      = lower(RegressWeight.Name);
RegressMethod.Name      = lower(RegressMethod.Name);
FactorizationMethod.Name= lower(FactorizationMethod.Name);

%% Precalculate tensors

y2 = y.^2;
y3 = y2.*y;
x2 = zeros(N,p2);
x3 = zeros(N,p3);
parfor i = 1 : N
    x2(i,:)	= reshape(X(i,:)'*X(i,:),[1,p2]); %#ok<PFGV>
    x3(i,:)	= reshape(x2(i,:)'*X(i,:),[1,p3]); %#ok<PFGV>
end

%%

% Estimate beta using the tensor method
switch RegressMethod.Name
    case 'ls'   % weighted least square solution
        
        M1 = (X'*X) \ (X'*y);
        
        fun2    = @(param) f_M2(param,x2,y2,1/sqrt(N));
        prini   = repmat(1/K,K,1);
        xini    = [randn(pK,1); prini; 0];%
        UB      = [  inf(pK,1); ones(K,1); inf];
        LB      = [ -inf(pK,1);zeros(K,1);   0];
        Aeq     = [zeros(1,pK), ones(1,K),   0];
        beq     = 1;
        xout    = fmincon(fun2,xini,[],[],Aeq,beq,LB,UB);
        %xout    = fmincon(fun2,xini);
        
        beta    = reshape(xout(1:pK),[p,K]);
        pr      = xout(pK+1:pK+K)';
        M2      = (pr.*beta)*beta';
        VarErr  = xout(end);
        
        [U,S,~] = svd(M2);
        S   = diag(S);
        W   = U(:,1:K)*diag( 1./sqrt(S(1:K)) );
        
        a       = 3*VarErr*(X*M1);
        fun3	= @(param) f_M3(param,x3,y3,a,1/sqrt(N));
        xini    = [beta(:); pr'];%
        UB      = [  inf(pK,1); ones(K,1)];
        LB      = [ -inf(pK,1);zeros(K,1)];
        Aeq     = [zeros(1,pK), ones(1,K)];
        beq     = 1;
        xout    = fmincon(fun3,xini,[],[],Aeq,beq,LB,UB);
        
        beta    = reshape(xout(1:pK),[p,K]);
        pr      = xout(pK+1:pK+K)';
        vecM3   = zeros(p3,1);
        for k = 1 : K
            b   = beta(:,k);
            vecM3	= vecM3 + reshape(reshape((pr(k)*b)*b',[p2,1])*b',[p3,1]);
        end
        M3	= reshape(vecM3,[p,p,p]);
        M3WWW = zeros(K,K,K);
        for k = 1 : K
            TensorMap_k	= sum(M3 .* reshape(W(:,k),[1,1,p]),3);
            for j = 1 : K
                TensorMap_j	= sum(TensorMap_k .* W(:,j)', 2);
                M3WWW(:,j,k)= W' * TensorMap_j;
            end
        end
        
        switch FactorizationMethod.Name
            case 'tpm'
                [eigenv,eigenval]   = TPM(M3WWW);
        case 'rpsf'
                [eigenv,eigenval]   = RPSF(M3WWW);
        end
        pr0     = 1./(eigenval.^2);
        pr      = pr0/(sum(pr0));
        beta    = (pinv(W')*(diag(eigenval)*eigenv)) .* (pr0./pr);
        Solution.pr     = pr;
        Solution.beta   = beta;
        
end

end

%% Private function

function f = f_M2(vecParam,x2,y2,lambda)
    global N p K p2 pK
    
    beta    = reshape(vecParam(1:pK),[p,K]);
    pr      = vecParam(pK+1:pK+K)';
    varerr	= vecParam(end);
    
    M2      = (pr.*beta)*beta';
    vecM2   = reshape(M2,[1,p2]);
    nucnorm = sum(svd(M2),1);
       
    f = 0.5*sum((x2*vecM2' - y2 + varerr).^2)/N + lambda*nucnorm;
end

function f = f_M3(vecParam,x3,y3,a,lambda)

    global N p K p2 p3 pK
    
    beta    = reshape(vecParam(1:pK),[p,K]);
    pr      = vecParam(pK+1:pK+K)';
    
    % C = tensorprod(A,B,dimA,dimB) in 2022
    vecM3 = zeros(p3,1);
    for k = 1 : K
        b   = beta(:,k);
        vecM3	= vecM3 + reshape(reshape((pr(k)*b)*b',[p2,1])*b',[p3,1]);
    end
    % Due to supersymmetry of bⓍbⓍb, 1,2,3 unfolding are identical
    M3_unfold = reshape(vecM3,[p,p2]);
    
    nucnorm = sum(svd(M3_unfold));
    f = 0.5*sum((x3*vecM3 + a - y3).^2)/N + lambda*nucnorm;
end