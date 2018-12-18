function [fkin_array] = forward_kinematics(DH)
    links = size(DH,1);
    fkin_array = sym('o', [4, 4, links]);
    for i=1:links
        fkin_array(:,:,i) = sym(eye(4));
    end
    
    for i=1:links
        for j=i:links
            fkin_array(:,:,j)=fkin_array(:,:,j)*link_transform(DH, i);
        end
        fkin_array(:,:,i)=simplify(fkin_array(:,:,i), 10);
    end
end