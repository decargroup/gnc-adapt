function phi = rotm2rotvec(R)
%ROTM2ROTVEC Converts Rotation matrix 'R' to rotation vector 'phi'

if ~all(ismembertol(mtimesx(R, 'T', R), eye(3), 1e-6), 'all')
    error('Input must be a rotation matrix');
end

phi = rotm2axang(R);
phi = phi(:, 1:3)' .* phi(:, 4)';

end

