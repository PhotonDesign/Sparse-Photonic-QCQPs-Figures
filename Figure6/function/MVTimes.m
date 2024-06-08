
% A: Nx x Ny x 4 matrix
% v: Nx x Ny x 2 vector
function[C] = MVTimes(A,v)
    C = zeros(size(v));
    C(:,:,1,:) = A(:,:,1,:) * v(:,:,1,:) + A(:,:,2,:) * v(:,:,2,:);
    C(:,:,2,:) = A(:,:,3,:) * v(:,:,1,:) + A(:,:,4,:) * v(:,:,2,:);
end