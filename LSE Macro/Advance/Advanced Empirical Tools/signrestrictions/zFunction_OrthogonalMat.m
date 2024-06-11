function Q_all = zFunction_OrthogonalMat(method, M)

% -----------------------------------------------------%

% DESCRIPTION 
% generates orthogonal matrices, i.e. matrices Q such that Q*Q' = I. Uses
% three methods, rotation, reflections or QR decomposition. Rotation and
% reflection required the specification of a radiant. Qr instead can be
% applied to any matrix. Here it is applied to matrices generate from a
% normal distribution. The code generates M orthogonal matrices
  
% -----------------------------------------------------%
       
Q_all = NaN*ones(2,2,M);

% used for rotations/reflection
theta_vec = linspace(0, 2*pi, M); 

% used for QR decomposition
reset(RandStream.getDefaultStream);
step = randn(2,2,M); 
        
for i = 1:M

    % ROTATIONS
    if method == 1
        Q = [cos(theta_vec(i)), -sin(theta_vec(i)); sin(theta_vec(i)), cos(theta_vec(i))];

    % REFLECTIONS
    elseif method == 2
        Q = [cos(theta_vec(i)), sin(theta_vec(i)); sin(theta_vec(i)), -cos(theta_vec(i))];

    % QR DECOMPOSITIONS
    elseif method == 3
        [Q_bis, R_bis] = qr(step(:,:,i));  % so that step1 = Q_bis*R_bis    
        Q = Q_bis';
        
    end

    assert(max(max(Q'*Q - eye(2))) < 0.001)

Q_all(:,:,i) = Q;

end

end

    