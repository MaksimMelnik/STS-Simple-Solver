function k_VT=kvt_fho_fr(t, M1, M2, coefs_for_poly) 
    %display(M1.name, M2.name);

    if M1.name == 'O2'
        vibr_lvls = 35;
        if M2.name == 'O2'
            coefs_num = 1;
        end
        if M2.name == 'O'
            coefs_num = 2;
        end
        if M2.name == 'Ar'
            coefs_num = 3;
        end
        if M2.name == 'N'
            coefs_num = 4;
        end
        if M2.name == 'N2'
            coefs_num = 5;
        end
        if M2.name == 'NO'
            coefs_num = 6;
        end
    end

    if M1.name == 'N2'
        vibr_lvls = 46;
        if M2.name == 'O2'
            coefs_num = 7;
        end
        if M2.name == 'O'
            coefs_num = 8;
        end
        if M2.name == 'Ar'
            coefs_num = 9;
        end
        if M2.name == 'N'
            coefs_num = 10;
        end
        if M2.name == 'N2'
            coefs_num = 11;
        end
        if M2.name == 'NO'
            coefs_num = 12;
        end
    end

    if M1.name == 'NO'
        vibr_lvls = 37;
        if M2.name == 'O2'
            coefs_num = 13;
        end
        if M2.name == 'O'
            coefs_num = 14;
        end
        if M2.name == 'Ar'
            coefs_num = 15;
        end
        if M2.name == 'N'
            coefs_num = 16;
        end
        if M2.name == 'N2'
            coefs_num = 17;
        end
        if M2.name == 'NO'
            coefs_num = 18;
        end
    end

    coefs = coefs_for_poly(coefs_num);
    coefs = coefs{1};
    poly_coefs = coefs(1:vibr_lvls, 2:end);
    val_for_normalization = coefs(1:vibr_lvls, 1);
    
    poly_len = size(poly_coefs);
    poly_len = poly_len(1,2);
    poly = ones(1, poly_len); % a0*t^0 + a1*t^(-1/3) + a2*t^(-2/3)...
    for i = 1:poly_len-1
       poly(i+1) = t.^(-i/3);
    end
    
    k_VT = zeros(1, M1.num_vibr_levels(1)-1);
    for i = 1:M1.num_vibr_levels-1
       k_VT(i) = exp(dot(poly_coefs(i,:), poly)*val_for_normalization(i));
    end
end