function A = hat(v)

    %It returns the skew symmetric matrix A associated to the vector v

    A = zeros(3);
    A(1,2) = -v(3);
    A(1,3) = v(2);
    A(2,3) = -v(1);
    
    A = A - A';
    
end