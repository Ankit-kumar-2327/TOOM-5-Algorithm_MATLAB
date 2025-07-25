function [mul_value] = TOOM5_MY_PROPO(a, b)

    a_sym = sym(a);
    b_sym = sym(b);

    a_str=char(a_sym);
    b_str=char(b_sym);

    if strlength(a_str) > strlength(b_str)
        b_str = pad(b_str, strlength(a_str), 'left', '0');
    elseif strlength(b_str) > strlength(a_str)
        a_str = pad(a_str, strlength(b_str), 'left', '0');
    end

    n = strlength(a_str);
    if mod(n, 5) ~= 0
        pad_len = 5 - mod(n, 5);
        a_str = pad(a_str, n + pad_len, 'left', '0');
        b_str = pad(b_str, n + pad_len, 'left', '0');
        n = strlength(a_str);
    end
    m = n / 5;
    x= sym(10)^m;

    a4 = sym(a_str(1:m));
    a3 = sym(a_str(m+1:2*m));
    a2 = sym(a_str(2*m+1:3*m));
    a1 = sym(a_str(3*m+1:4*m));
    a0 = sym(a_str(4*m+1:end));

    b4 = sym(b_str(1:m));
    b3 = sym(b_str(m+1:2*m));
    b2 = sym(b_str(2*m+1:3*m));
    b1 = sym(b_str(3*m+1:4*m));
    b0 = sym(b_str(4*m+1:end));

    evalA = @(v) a0 + a1*v + a2*v^2 + a3*v^3 + a4*v^4;
    evalB = @(v) b0 + b1*v + b2*v^2 + b3*v^3 +b4*v^4;

    points = sym([0,1,-1,2,-2,3,-3,4]);

    y = vpa(sym(zeros(9,1)));

    for i = 1:8
        y(i) = evalA(points(i)) * evalB(points(i));
    end

    y(9) = a4 * b4;

V = [ 1   0   0   0   0    0     0     0      0;
      1   1   1   1   1    1     1     1      1;
      1  -1   1  -1   1   -1     1    -1      1;
      1   2   4   8   16   32    64    128    256;
      1  -2   4  -8   16  -32    64   -128    256;
      1   3   9   27  81   243   729   2187   6561;
      1  -3   9  -27  81  -243   729  -2187   6561;
      1   4   16  64  256  1024  4096  16384  65536;
      0   0   0   0   0    0    0      0      1   ];

   syms c0 c1 c2 c3 c4 c5 c6 c7 c8
   C = [c0 ; c1 ; c2 ; c3 ; c4 ; c5 ; c6 ; c7 ; c8];

   equation = V*C == y;
   display(equation);

   solutions = solve(equation,[c0,c1,c2,c3,c4,c5,c6,c7,c8]);
   display(solutions);

   mul_value = solutions.c0 + solutions.c1*x + solutions.c2*x^2 + ...
             solutions.c3*x^3 + solutions.c4*x^4 + solutions.c5*x^5 + ...
             solutions.c6*x^6 + solutions.c7*x^7 + solutions.c8*x^8;

   syms p0 p1 p2 p3 p4 p5 p6 p7 p8  % here p0 to p8 has same values as y(1) to y(9)
   P = [p0 ; p1 ; p2 ; p3 ; p4 ; p5 ; p6 ; p7 ; p8];

   equation = V*C == P;

   format long g
   solutions = solve(equation,[c0,c1,c2,c3,c4,c5,c6,c7,c8]);
   fields = fieldnames(solutions);

   for i = 1:numel(fields)
       fprintf("%s : ",fields{i});
       disp(solutions.(fields{i}));
   end
   
   disp(mul_value);
end
