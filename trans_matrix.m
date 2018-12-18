function [matrix] = trans_matrix(fkin_array, i)
    if i == 0
        matrix = [0, 0, 0];
    else
        matrix = fkin_array(1:3, 4, i);
    end
end