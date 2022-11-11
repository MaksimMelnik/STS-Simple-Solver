function out=Evibr(T, num_vibr_levels, f)
if nargin==2
   f=1;
end
e_i=levels_e(num_vibr_levels);
n=density_f(T, f, num_vibr_levels);
out=e_i*n';
end