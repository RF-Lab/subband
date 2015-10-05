function F = get_sb_matrix_2(N, w, tau)

F = zeros(N, N, N, N) ;

d_zero_value = 2*w;

for n=1:N
    for k=1:N
        for m = 1:N
            for q = 1:N
                d = k - n + m - q + tau;
                if d == 0
                    F(n, k, m, q) = d_zero_value;
                else
                    F(n, k, m, q) = -1j/ d * (exp(1i*w*d) - exp(1i*(-w)*d));
                end;
            end;
        end;
    end;
end;