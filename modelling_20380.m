%% 20.380 modeling

k = [2.614e-4, 1.153e-8, 1e6, .27, 1.14, 1e5, 2.306e-9];
ode_options = odeset('RelTol',1e-12,'AbsTol',[1e-12]);
dosage=1;
t_span=linspace(0,12*60*60,1000000);
y0 = [1.43e-22, 0, dosage, 0, 0, 0.0606, 0, 0];
[t,y]=ode15s(@model,t_span, y0, ode_options, k);

figure
subplot(2,2,1)
plot(t,y(:,4))
xlabel("Time (s)")
ylabel("Fluorescent Probe-Enzyme Complex (M)")
subplot(2,2,2)
plot(t,y(:,3))
xlabel("Time (s)")
ylabel("Quenched Probe (M)")
subplot(2,2,3)
plot(t,y(:,1))
xlabel("Time (s)")
ylabel("P aeruginosa (M)")
subplot(2,2,4)
plot(t,y(:,2))
xlabel("Time (s)")
ylabel("Elastase (M)")

function dydt = model(t,y,k)
    dydt = zeros(8,1);
    kg = k(1);
    kep = k(2);
    kon_i = k(3);
    kcat_s = k(4);
    koff_s = k(5);
    kon_s = k(6);
    K = k(7);
    tau=3600;

    % Probe delivery

    B = y(1);
    E = y(2);
    PQ = y(3);
    PE = y(4);
    Q = y(5);
    S = y(6);
    SE = y(7);
    S_c = y(8);

    dydt(1,:) = kg*B;
    dydt(2,:) = kep*B+kcat_s*SE+koff_s*SE-kon_i*PQ*E-kon_s*S*E;
    if t>tau
        dydt(3,:) = -kon_i*PQ*E;
    else
        dydt(3,:)=0;
    end
    dydt(4,:) = kon_i*PQ*E;
    dydt(5,:) = kon_i*PQ*E;
    dydt(6,:) = K+koff_s*SE-kon_s*S*E;
    dydt(7,:) = kon_s*S*E-koff_s*SE-kcat_s*SE;
    dydt(8,:) = kcat_s*SE;
end