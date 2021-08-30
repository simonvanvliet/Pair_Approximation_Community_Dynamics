function [rangeP, rangeT, varargout] = GR_analytical_fromSettings(settings)
%calculates analytical growth range from settings file

Deff=(settings.D0/(settings.cellSpacing)^2)...
    *(1-settings.rho)/(1+settings.rho/2)....
    *(1-settings.rho)/settings.rho;

rlP=settings.rl/sqrt(settings.delta_l);
rlT=settings.rl*sqrt(settings.delta_l);

ruP=settings.ru/sqrt(settings.delta_u);
ruT=settings.ru*sqrt(settings.delta_u);

DP=Deff/sqrt(settings.delta_D);
DT=Deff*sqrt(settings.delta_D);

r0P = sqrt(DP ./ ((ruP + rlP)));
r0T = sqrt(DT ./ ((ruT + rlT)));

ru = settings.ru;
rl = settings.rl;
r0 = sqrt(Deff ./ ((ru + rl)));

[GRsimpleP, GRfullP] = GR_analytical(r0P, settings.Iconst2, rlP);
[GRsimpleT, GRfullT] = GR_analytical(r0T, settings.Iconst1, rlT);
[GRsimple, ~] = GR_analytical(r0, settings.Iconst1, rl);


%SET: SELECT FULL OR SIMPLE RANGE
rangeP = GRsimpleP;
rangeT = GRsimpleT;

varargout{1} = GRsimple;