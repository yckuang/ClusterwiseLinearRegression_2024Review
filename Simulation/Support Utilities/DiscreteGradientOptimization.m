function [x,info] = DiscreteGradientOptimization(fun,x,lambda,optim)

evalcount = 0;

bExtrapolate= true;
bStationery = false;
p       = length(x);
% function value at the current location
funx    = fun(x);                                                           evalcount = evalcount + 1;

Drop    = [];
bxLoop  = true;
xcount	= 0;
while bxLoop
    xcount	= xcount + 1;
    
    % Find descend direction
    D       = [];
    g       = randn(1,p); g = g/norm(g);
    bgLoop  = true;
    gcount	= 0;
    while bgLoop
        gcount	= gcount + 1;
        
        v = DiscreteGradient(fun,x,g,lambda,optim.DGopt);
        D = [D,v']; %#ok<AGROW>
        
        [w,mgw] = MinNorm(D);
        if mgw < optim.delta    % Found stationery point
            bgLoop = false;
            bStationery = true;
        else
            g = -w/mgw;
            if fun(x+lambda*g) < funx - optim.c*lambda*mgw                  evalcount = evalcount + 1;
                bgLoop = false; % Found a viable descend direction
            end
        end        
        if gcount > 2*p
            bgLoop = false;
        end        
    end

    % Update x
    if bStationery
        bxLoop = false;
    else    % Find the biggest viable step in the new direction
        L = lambda*[1:0.5:3]';
        x_test  = x + L*g;
        bValid  = fun(x_test) < funx + optim.c*L*mgw;                       evalcount = evalcount + length(L);
        sig     = L(find(bValid,1,'last')); % maximum L that satisfies the descend condition
        if isempty(sig) % none satisfy the condition
            sig = lambda*0.5;
        end
        x = x + sig*g;
    end
    funxold = funx;
    funx    = fun(x);                                                       evalcount = evalcount + 1;
    
    abstol = abs(funxold-funx);
    reltol = abstol/abs(funxold);
    Drop = [Drop;log([abstol,lambda])]; %#ok<AGROW>
    
    if size(Drop,1) >= 3 
        if bExtrapolate 
            poly = polyfit(Drop(:,1),Drop(:,2),2);
            u = Drop(end,1)-0.7; % -3dB reduction in error
            v = polyval(poly,u);
            if v > Drop(end,1) || Drop(end,1) > Drop(end-1,1)
                % Overshoot detected, stop extrapolation use geometrically shrinking step instead
                bExtrapolate = false;
                lambda  = optim.beta*lambda;
            else
                lambda = exp(v);
            end
        else    % use geometrically shrinking step
            lambda  = optim.beta*lambda;
        end
    end
    
    if xcount > optim.maxiter
        bxLoop = false;
        warning('Optimization terminated because the maximum number of iteration exceeded')
    end
    if reltol < optim.reltol || abstol < optim.abstol
        bxLoop = false;
    end
    
end
% Drop
info.evalcount = evalcount;

end


%% Private function
function [w,mgw] = MinNorm(D)
    % This implementation uses constrained quadratic programming instead of
    % Wolfe's combinatorial algorithm. The algorithm normally converges in
    % a few steps.
    [~,n] = size(D);
    
    alpha0  = repmat(1/n,[n,1]);    
    options = optimoptions('quadprog','Display','off');
    alpha   = quadprog(0.5*(D'*D),[],[],[],...
                        ones(1,n),1,zeros(n,1),[],alpha0,options);

    w   = (D*alpha)';
    mgw = norm(w);
end



function [D,evalcount] = DiscreteGradient(fun,x,g,lambda,DGopt)
    evalcount = 0;
    
    if ~isempty(DGopt)
        a = DGopt;
    else
        a = 0.5;
    end
    z = lambda^(1.44224957031);

    p = length(x);
    [~,i] = max(abs(g));

    %%
    D = zeros(1,p);
    u = zeros(p+1,p);

    ea = repmat(a.^[1:p],p,1);
    for index = 1 : p-1
        ea(index,index+1:p) = 0;
    end
    e  = sign(randn(1,p));
    ea = e.*ea;

    u(1,:) = x + lambda*g;
    for j = 1 : p
        u(j+1,:) = u(1,:) + z*ea(j,:);
        if j ~= i
            D(j) = (fun(u(j+1,:))-fun(u(j,:)))/(z*ea(j,j));                 evalcount = evalcount + 2;
        end
    end
    D(i) = (fun(u(1,:)) - fun(x) - lambda*(D*g'))/(lambda*g(i));            evalcount = evalcount + 2;

end