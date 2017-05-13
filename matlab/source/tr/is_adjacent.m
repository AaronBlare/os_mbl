function [adj] = is_adjacent (x, y)
global states_id
adj = states_id.adjacent(bitxor(x, y) + 1);
end