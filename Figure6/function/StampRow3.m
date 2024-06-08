% Helper function to fill out 3rd column of 3x3 matrix
function[X] = StampRow3(X,A,B,C,D)
X(:,:,1) = A;
X(:,:,2) = B;
X(:,:,3) = C;
X(:,:,4) = D;
end