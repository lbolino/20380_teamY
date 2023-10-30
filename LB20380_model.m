%% 20.380 modeling

f1 = figure;
k = [4.478e-5, 1.153e-8, 0.405, 1e6, 0.27, 25.2, 1.14, 1e5, 2.306e-9];
t_span = [0,10000];
y0 = zeros(9,1);

[t,y] = ode15s(@model, t_span, y0, [], k);
P = y(:,4);
plot(t, P)

function dydt = model(~,y,k)
    dydt = zeros(9,1);
    kg = k(1);
    kep = k(2);
    koff_i = k(3);
    kon_i = k(4);
    kcat_s = k(5);
    kcat_i = k(6);
    koff_s = k(7);
    kon_s = k(8);
    K = k(9);

    B = y(1);
    E = y(2);
    PQ = y(3);
    P = y(4);
    PQE = y(5);
    Q = y(6);
    S = y(7);
    SE = y(8);
    S_c = y(9);

    dydt(1,:) = kg*B;
    dydt(2,:) = kep*B+kcat_i*PQE+kcat_s*SE+koff_i*PQE-kon_i*PQ*E-kon_s*S*E;
    dydt(3,:) = koff_i*PQE-kon_i*PQ*E;
    dydt(4,:) = kcat_i*PQE;
    dydt(5,:) = kon_i*PQ*E-kcat_i*PQE;
    dydt(6,:) = kcat_i*PQE;
    dydt(7,:) = K+koff_s*SE-kon_s*S*E;
    dydt(8,:) = kon_s*S*E-koff_s*SE-kcat_s*SE;
    dydt(9,:) = kcat_s*SE;
end




