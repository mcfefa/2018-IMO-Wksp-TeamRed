function [ xdot ] = bladder_ODE( t, x, params, BCG, tspan )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

HI = x(1);
HN = x(2);
LI = x(3);
LN = x(4);
A = x(5);
S = x(6);

rHI = params(1);
rHN = params(2);
rLI = params(3);
rLN = params(4);
K = params(5);
delta0 = params(6);
rhoS = params(7);
rhoI = params(8);
lambdaA = params(9);
lambdaS = params(10);
beta = params(11);
eta = params(12);
lambdaAS = params(13);

% [rPI, rPL, rNI, rNL, K, delta0, rhoS, rhoI, lambdaS, lambdaI] = params

BCG0 = BCG((tspan >= t) + (tspan <= t))

delta = delta0 * max((A/(S+A) - 0.5), 0);

HIdot = rHI * HI * (1 - (HI + HN + LI + LN) / K) - delta * HI;

HNdot = rHN * HN * (1 - (HI + HN + LI + LN) / K);

LIdot = rLI * LI * (1 - (HI + HN + LI + LN) / K);

LNdot = rLN * LN * (1 - (HI + HN + LI + LN) / K);

Adot = rhoS * (HI + LI) + beta * BCG / (eta + BCG) - lambdaA * A; % additional -lambdaAS * A * S term?

Sdot = rhoI * (HN + LN) - lambdaS * S;

xdot = [HIdot HNdot LIdot LNdot Adot Sdot]';


end

