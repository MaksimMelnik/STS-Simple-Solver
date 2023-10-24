function k_VT=kvt_fho_fr(t, M1, M2, coefs_O2)
   if M2.name == 'O2'
      display('O2');
      poly_coefs = coefs_O2(1:35, 2:end);
      val_for_normalization = coefs_O2(1:35, 1);

   else
       if M2.name == 'O'
          display('O');
          poly_coefs = coefs_O2(35:70, 2:end);
          val_for_normalization = coefs_O2(35:70, 1);      
       end
       if M2.name == "Ar"
          display('Ar');
          poly_coefs = coefs_O2(70:105, 2:end);
          val_for_normalization = coefs_O2(70:105, 1);
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