function res = actionSE3(B,input)

    res = [B(1:3,1:3)*input(1:3); B(1:3,1:3)*input(4:6) + cross(B(:,end),(B(1:3,1:3)*input(1:3)))];


end