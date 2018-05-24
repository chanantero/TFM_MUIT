function n = quatrotate_custom( q, r )
% q. 1x4 or Mx4 matrix. Each row is a quaternion
% r. 1x3 or Mx3 matrix. Each row is a 3D vector

numQuat = size(q, 1);
numVec = size(r, 1);

if numQuat > 1 && numVec == 1
    n = zeros(numQuat, 3);
    for k = 1:numQuat

        M = quat2LinearTransformation(q(k, :));
        
        n(k, :) = (M * r')';
        
    end

elseif numQuat == 1 && numVec > 1
   M = quat2LinearTransformation(q);        
   n = (M * r')';
    
elseif numQuat == numVec
        n = zeros(numQuat, 3);

    for k = 1:numQuat

        M = quat2LinearTransformation(q(k, :));
        
        n(k, :) = (M * r(k, :)')';
        
    end
    
end
    
end

function M = quat2LinearTransformation(q)

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

M = [(1 - 2*q2^2 - 2*q3^2), 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2);
    2*(q1*q2 - q0*q3), (1 - 2*q1^2 - 2*q3^2), 2*(q2*q3 + q0*q1);
    2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), (1 - 2*q1^2 - 2*q2^2)];

end

