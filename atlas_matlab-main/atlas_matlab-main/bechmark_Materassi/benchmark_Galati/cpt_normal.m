function N = cpt_normal(n)
% Compute the transformation matrix
    if length(n) ~= 3
        error('Input n must be a 3*1 vector');
    end

    n1 = n(1);
    n2 = n(2);
    n3 = n(3);

    N = [ n1   0    0   n2  n3   0;
          0    n2   0   n1   0  n3;
          0    0   n3   0   n1  n2 ];
%     N = repmat(N,4,1);
end