function[X] = Times3(A,B,C)
X = MTimes(A, MTimes(B,C));
end