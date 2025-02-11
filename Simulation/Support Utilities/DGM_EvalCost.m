function f = DGM_EvalCost(b,y,X1,pp1,k)

[NumEval,~] = size(b);
f = zeros(NumEval,1);
for index = 1 : NumEval
    beta	= reshape(b(index,:),[pp1,k]);
    aerr    = min(abs(y-X1*beta),[],2);
    f(index)= rms(aerr);
end
