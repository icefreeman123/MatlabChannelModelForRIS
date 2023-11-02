function [ AllPtAndCenter_new ] = RodriguesRotVec( OrgVec, RotVec, Pt )
%% Rodrigues' rotation formula
RotAng = acos( dot(OrgVec, RotVec) / (norm(OrgVec)*norm(RotVec)) ); % rotation angle
RotAxis = cross(OrgVec, RotVec); % rotation axis
RotAxis = RotAxis / norm(RotAxis);
Kmat = [0, -1*RotAxis(3), RotAxis(2); ...
        RotAxis(3), 0, -1*RotAxis(1); ...
        -1*RotAxis(2), RotAxis(1), 0]; % cross-product matrix
Rmat = eye(3) + sin(RotAng)*Kmat + (1-cos(RotAng))*Kmat*Kmat; % rotation matrix
AllPtAndCenter_new = Rmat * Pt; % rotated vector 

end

