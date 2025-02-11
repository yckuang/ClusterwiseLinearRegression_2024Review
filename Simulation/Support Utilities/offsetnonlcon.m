function [c,ceq] = offsetnonlcon(u,set,offset)

c   = -sum((u-set).^2) + offset^2;
c   = c';
ceq = [];