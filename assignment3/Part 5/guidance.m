function chi_d = guidance(y_e_p, y_e_p_int, last_wp, next_wp, lookahead_delta, kappa)

    path_vector = next_wp - last_wp;
    pi_p = atan2(path_vector(2), path_vector(1));
    
    Kp = 1/lookahead_delta;
    Ki = kappa*Kp;
    chi_d = pi_p - atan(Kp*y_e_p + Ki*y_e_p_int);
 
end