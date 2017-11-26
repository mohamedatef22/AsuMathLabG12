# AsuMathLabG12
# Phase 1:Core Operation
Mathematical library software.
Include:
 -Class for matrix.
 -Support dynamic creation and destruction of matrices of any size.
 -Support addition, subtraction, multiplication, transpose and division.
 -Process input user commands and show results directly.
 -Process input file.

# Phase 2:Advanced operations and Tuning
-Support mathematical functions (Trigonometric,Logarethmic,Roots and Power)
-Support mathematical expressions (with dot or not)
-Support flexible matrix parser (accept matrix in matrix ,expressions and variables)
-Support error handling

#  Some Examples
- A = 5.5 + 12 * sin(0.4) + 2.2^4;
- A should equal 33.5986
- B=[1.2 2.3 A;[1.3 2.4;4.6 1.3],[3.2;7.8]];
- B should be a matrix 3x3 and its values 1.2 2.3 33.5986
                                          1.3 2.4 3.2
                                          4.6 1.3 7.8
- C=rand(4,4)

C =

    0.8147    0.6324    0.9575    0.9572
    0.9058    0.0975    0.9649    0.4854
    0.1270    0.2785    0.1576    0.8003
    0.9134    0.5469    0.9706    0.1419
    

- D=eye(4,4)

D =    1 0 0 0
       0 1 0 0
       0 0 1 0
       0 0 0 1

- E=zeros(2,3);


E =    0  0  0
       0  0  0
       
