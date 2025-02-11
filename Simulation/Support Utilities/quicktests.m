clc

Input.K = 4;
beta = randn(5,4);
Solution.beta = beta(:,[4,3,2,1]) + [0.05*randn([5,3]),1*randn([5,1])];
S = CorrectPermutation(Solution, beta, Input);




