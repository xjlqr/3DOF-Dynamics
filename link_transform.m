function [matrix] = link_transform(DH, i)
    syms a

    twist_transform(a) = [
        1,       0,       0, 0;
        0,  cos(a), -sin(a), 0;
        0,  sin(a),  cos(a), 0;
        0,       0,       0, 1;
    ];

    length_transform(a) = [
        1, 0, 0, a;
        0, 1, 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1;
    ];

    angle_transform(a) = [
        cos(a), -sin(a), 0, 0;
        sin(a),  cos(a), 0, 0;
            0,       0, 1, 0;
         0,       0, 0, 1;
    ];

    offset_transform(a) = [
        1, 0, 0, 0;
        0, 1, 0, 0;
        0, 0, 1, a;
        0, 0, 0, 1;
    ];

    matrix = twist_transform(DH(i, 1))*length_transform(DH(i, 2))*angle_transform(DH(i, 3))*offset_transform(DH(i, 4));
end






