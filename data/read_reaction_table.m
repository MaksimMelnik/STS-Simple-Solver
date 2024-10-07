function Reactions = read_reaction_table
% Reading table with reactions and converting to the code's format.
% 04.10.2024 by Maksim Melnik

% load('particles.mat') %#ok<LOAD>
load("reactions.mat", "Reactions")
% k = 1.380649e-23;           % Boltzmann constant, J/K
% eV_to_J = 1.6022e-19;       % eV to J

    % template
template.name = NaN;
template.source = NaN;
template.particles = NaN;
template.type = NaN;
template.neq_model = "equal";
template.direction_forward = true;
template.reverse = false;
template.index = NaN;
template.A   = 0;
template.n   = 0;
template.E   = 0;
template.d_T = 1;
template.E_th = 0;

     % zero reaction
% zero_r.name = 'zero';
% react1 = template;
% react1.type = "const";
% react1.A   = 0;
% react1.n   = 0;
% zero_r.data = react1;

% table_t = table2cell(readtable("reactions_data_table.csv"));
table_t = readtable("reactions_data_table.csv");
ind = 1;
temp_react = template;
temp_react.name     = table_t.ReactionName{ind};
temp_react.A        = table_t.A(ind);
temp_react.source   = table_t.Source{ind};
if isKey(Reactions, temp_react.name)
    error("reaction already exist")
else
    new_reaction.name = temp_react.name;
    new_reaction.particles = [convertCharsToStrings(table_t.Reagent1{ind}), ...
        convertCharsToStrings(table_t.Reagent2{ind}), ...
        convertCharsToStrings(table_t.Product1{ind}), ...
        convertCharsToStrings(table_t.Product2{ind})];
    if ~iscell(table_t.Product3(ind))
        if ~isnan(table_t.Product3(ind))
            error("5 components in the reaction")
        end
    end
    temp_react.particles    = new_reaction.particles;
    temp_react.index = {{table_t.Reagent1_e(ind), ...
                    convertCharsToStrings(table_t.Reagent1_v{ind})}, ...
        {table_t.Reagent2_e(ind), ...
                    convertCharsToStrings(table_t.Reagent2_v{ind})}, ...
        {table_t.Product1_e(ind), ...
                    convertCharsToStrings(table_t.Product1_v{ind})}, ...
        {table_t.Product2_e(ind), ...
                    convertCharsToStrings(table_t.Product2_v{ind})}};
    temp_react.type = table_t.type{ind};
    temp_react_Starik = temp_react;
    temp_react_Starik.source = strcat(temp_react_Starik.source,'_Starik');
    temp_react_Starik.neq_model = "Starik_test";
    keySet              = {temp_react.source, temp_react_Starik.source};
    valueSet            = {temp_react, temp_react_Starik};
    new_reaction.data   = containers.Map(keySet, valueSet);
    Reactions(new_reaction.name) = new_reaction.data;
end
end