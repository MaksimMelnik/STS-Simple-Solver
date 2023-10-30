function k_VT=kvt_fho_fr(t, M1, M2, coefs_for_poly)
   if M2.name == 'O2'
      %display('O2');
      coefs = coefs_for_poly(1);
      coefs = coefs{1};
      poly_coefs = coefs(1:35, 2:end);
      val_for_normalization = coefs(1:35, 1);

   else
       if M2.name == 'O'
          %display('O');
          coefs = coefs_for_poly(2);
          coefs = coefs{1};
          poly_coefs = coefs(1:35, 2:end);
          val_for_normalization = coefs(1:35, 1);      
       end
       if M2.name == "Ar"
          %display('Ar');
          coefs = coefs_for_poly(3);
          coefs = coefs{1};
          poly_coefs = coefs(1:35, 2:end);
          val_for_normalization = coefs(1:35, 1);
       end
   end

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