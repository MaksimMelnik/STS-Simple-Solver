function k_VT=get_k_VT_FHO_FR(t, M1, M2, poly_coefs, val_for_normalization)
   poly_len = size(poly_coefs);
   poly_len = poly_len(1,2);
   poly = ones(1, poly_len); % a0*t^0 + a1*t^(-1/3) + a2*t^(-2/3)...
   for i = 1:poly_len-1
       poly(i+1) = t.^(-i/3);
   end

   k_VT = zeros(1, M1.num_vibr_levels(1)-1);
   for i = 1:M1.num_vibr_levels-1
       k_VT(i) = exp(dot(poly_coefs(i,:),poly)*val_for_normalization(i));
   end
end