function F = get_sb_matrix_1(N, w1, w2, tau )

F = zeros(N, N) ;

d_zero_value = w2 - w1;

for n=1:N
    for k=1:N
        d = k - n + tau;
        if d == 0
            F(n, k) = d_zero_value;
        else
            F(n, k) = -1j/ d * (exp(1i*w2*d) - exp(1i*w1*d));
        end
    end;
end;