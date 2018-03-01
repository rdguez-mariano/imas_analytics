function [ T,R_phi1,R_phi2,lambda ] = affine_decomp( A, opts)
%AFFINE_DECOMP Summary of this function goes here
%   Detailed explanation goes here

%Default Parameter
epsilon = 10^(-10);
HUMAN_OUTPUT = false;

% Input Parameters
switch nargin
    case 2
        if (~isstruct(opts))
            error('Optional parameter is not of type struct!');
        end
        names = fieldnames(opts);
        for i=1:length(names)
            name = names{i};
            switch name
                case 'epsilon'
                    epsilon = opts.epsilon;
                case 'human_output'
                    HUMAN_OUTPUT = opts.human_output;                                        
            end
        end
    otherwise
        if (nargin~=1)
            error('Wrong number of input parameters!');
        end
end


%   [U,S,V] = svd(A);
%   A = U*S*V'
[U,T,V] = svd(A);
K = [-1 0; 0 1];

% K*D*K = D
if ((norm(det(U)+1)<=epsilon)&&(norm(det(V)+1)<=epsilon))
    U = U*K;
    V = V*K;
end

% Computing First Rotation
phi1 = atan2( V(2,1), V(1,1));
R_phi1 = [cos(phi1) sin(phi1); -sin(phi1) cos(phi1)];
assert(norm(V'-R_phi1,'fro')<=epsilon,'V not equal to our rotation');

%Computing Second Rotation
phi2 = atan2( U(1,2),U(1,1));
R_phi2 = [cos(phi2) sin(phi2); -sin(phi2) cos(phi2)];
assert(norm(U-R_phi2,'fro')<=epsilon,'U not equal to our rotation');

% Computing Tilt and Lambda
lambda = T(2,2);
T(1,1)=T(1,1)/T(2,2);
T(2,2)=1;

assert(norm(A - lambda*R_phi2*T*R_phi1,'fro')<=epsilon,'Couldnt decompose A');

if (HUMAN_OUTPUT)
    T = T(1,1);
    R_phi1 = phi1;%*(180/pi);
    R_phi2 = phi2;%*(180/pi);    
end

end

