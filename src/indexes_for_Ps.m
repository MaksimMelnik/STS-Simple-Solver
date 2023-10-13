function index = indexes_for_Ps(Ps)
% Function for calculation indexes of molecules for kinetics structure.
% Ps is the cell array of particles in the mixture compasition.
% 30.09.2023 by Maksim Melnik.

index = cell(1, length(Ps));
first = 1;
for ind = 1:length(Ps)
    num_states = sum(Ps{ind}.num_vibr_levels(1:Ps{ind}.num_elex_levels));
    index{ind} = first : (first+num_states-1);
    first = index{ind}(end) + 1;
end
end