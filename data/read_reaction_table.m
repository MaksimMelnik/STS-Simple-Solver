function Reactions = read_reaction_table
% Reading table with reactions and converting to the code's format.
% 04.10.2024 by Maksim Melnik

load("reactions.mat", "Reactions")
    % template
template.name = NaN;
template.source = NaN;
template.particles = NaN;
template.type = NaN;
template.neq_model = "equilibrium";
template.direction_forward = true;
template.reverse = false;
template.index = NaN;
template.A   = 0;
template.n   = 0;
template.E   = 0;
template.d_T = 1;
template.E_th = 0;
    % reading the table with reactions
table_t = readtable("reactions_data_table.csv");
for ind = 1:height(table_t)
 temp_react = template;
 temp_react.name     = table_t.ReactionName{ind};
 temp_react.A        = table_t.A(ind);
 if ~isnan(table_t.n(ind))
      temp_react.n = table_t.n(ind);
 end
 if ~isnan(table_t.E(ind))
      temp_react.E = table_t.E(ind);
 end
 if ~isnan(table_t.T_reaction(ind))
      temp_react.T_reaction = table_t.T_reaction(ind);
 end
 if ~isnan(table_t.E_threshold(ind))
      temp_react.E_threshold = table_t.E_threshold(ind);
 end
 if ~isempty(table_t.neq_model(ind))
      temp_react.neq_model = table_t.neq_model{ind};
 end
 temp_react.neq_model = table_t.neq_model{ind};
 temp_react.source   = table_t.Source{ind};
 if isKey(Reactions, temp_react.name)
    error("reaction already exist")
 else
    temp_react.particles = [convertCharsToStrings(table_t.Reagent1{ind}), ...
        convertCharsToStrings(table_t.Reagent2{ind}), ...
        convertCharsToStrings(table_t.Product1{ind}), ...
        convertCharsToStrings(table_t.Product2{ind})];
    R1v = index_convert_fun(table_t.Reagent1_v(ind));
    R2v = index_convert_fun(table_t.Reagent2_v(ind));
    P1v = index_convert_fun(table_t.Product1_v(ind));
    P2v = index_convert_fun(table_t.Product2_v(ind));
    temp_react.index = {{table_t.Reagent1_e(ind), R1v}, ...
       {table_t.Reagent2_e(ind), R2v}, {table_t.Product1_e(ind), P1v}, ...
                                        {table_t.Product2_e(ind), P2v}};
    if ~isempty(table_t.Product3{ind})
        temp_react.particles = [temp_react.particles, ...
                    convertCharsToStrings(table_t.Product3{ind})];
        P3v = index_convert_fun(table_t.Product3_v(ind));
        temp_react.index{end + 1} = {table_t.Product3_e(ind), P3v};
    end
    temp_react.type = table_t.type{ind};
    if ~isnan(table_t.direction_forward(ind))   % reading the direction
        if ~table_t.direction_forward(ind)      %   of the reaction
            temp_react.direction_forward = false;
        end
    end
    temp_react_Starik = temp_react;
    temp_react_Starik.source = strcat(temp_react_Starik.source,'_Starik');
    temp_react_Starik.neq_model = "Starik_test";
    keySet      = {temp_react.source, temp_react_Starik.source};
    valueSet    = {temp_react, temp_react_Starik};
    data        = containers.Map(keySet, valueSet);
    Reactions(temp_react.name) = data;
 end
end
end

function Pv = index_convert_fun(Pv)
if iscell(Pv)
 switch class(Pv{1})
     case 'char'
        Pv = convertCharsToStrings(Pv{1});
     case {'string', 'double'}
        Pv = Pv{1};
     otherwise
         error("unsupported class")
 end
 if Pv ~= "all"
    error("not all indexes")
 end
else
    if isnan(Pv)
        Pv = "all";
    else
        error("unsupported format")
    end
end
end