function rxx_sb2 = get_acf_from_sb_matrix_2(N, F0, F1, F2, x)

% we need 3 pointes to use Yule-Walker equation
rxx_sb2 = zeros(3, 1);

for n=1:N
    for k=1:N
        for m = 1:N
            for q = 1:N
                x_val =  x(n) * conj(x(k)) * x(m) * conj(x(q));
                rxx_sb2(1) = rxx_sb2(1) + x_val * F0(n, k, m, q);
                rxx_sb2(2) = rxx_sb2(2) + x_val * F1(n, k, m, q);
                rxx_sb2(3) = rxx_sb2(3) + x_val * F2(n, k, m, q);                
            end;
        end;
    end;
end;