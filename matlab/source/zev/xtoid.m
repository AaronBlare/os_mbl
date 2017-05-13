function [id] = xtoid (x)
global states_id
id = states_id.xtoid(x + 1);
end