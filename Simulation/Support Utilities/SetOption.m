function Option = SetOption(Option, Method, K)
% Sets default paramaters of options
% 25/11/2022
%%
Option.NumMix           = K;
Option.Initialization   = Method.init;
Option.DistType         = Method.distr;
Option.RegressWeight    = Method.regrweight;
Option.RegressMethod    = Method.regr;
Option.ScaleEstWeight   = Method.scaleweight;
Option.Initialization.FactorizationMethod   = Method.init.Name;



end