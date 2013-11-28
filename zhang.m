fs = 1024;
fc = [8, 9.75, 11.89, 14.49, 17.67, 21.53, 26.25, 32]; %center frequencies
qf = 0.33;
forder = 4;

f1 = @(f0,Q) f0*(sqrt(1+1/(4*Q^2))- 1/(2*Q));
f2 = @(f0,Q) f0*(sqrt(1+1/(4*Q^2))+ 1/(2*Q));

Wn = [f1(fc(1), qf), f2(fc(1), qf)]/(fs/2);
[b,a]  =cheby2(forder,30, Wn); %need to now r value, for now setting to 20
%s = dfilt.df2(b,a);
%fvtool(s);

