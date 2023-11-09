%% 20.380 modeling
k =[2.614e-4, 1.01e14, 1e6, .27, 9.33, 1e5, 2.306e-9,3600, 5.35e-6];
ode_options = odeset('RelTol',1e-12,'AbsTol',[1e-12]);
dosage=1e-4;
t_span=linspace(0,2.*60.*60,100000);
y0 = [1.43e-22, 0, dosage, 0, 0, 0.0606, 0, 0];
[t,y]=ode15s(@model,t_span, y0, ode_options, k);

figure
hold all
kon_i_vals=[100,1000,10000,100000,1000000];
C=['-k','-b','-r','-g','-y'];
for i=1:5
    k_iter=[2.614e-4, 1.01e14, kon_i_vals(i), .27, 9.33, 1e5, 2.306e-9,3600,5.35e-6];
    [t,y1]=ode15s(@model,t_span, y0, ode_options, k_iter);
    subplot(2,1,1)
    hold on
    plot(t,y1(:,4),C(i))
    xlim([3600 7200])
    xlabel('Time (s)')
    ylabel("Fluorescent Probe-Enzyme Complex (M)")
    subplot(2,1,2)
    hold on
    plot(t,y1(:,8),C(i))
    xlabel('Time (s)')
    ylabel("Cleaved Elastin (M)")
end
subplot(2,1,1)
legend('100 M^{-1} s^{-1}','1000 M^{-1} s^{-1}','10000 M^{-1} s^{-1}','100000 M^{-1} s^{-1}','1000000 M^{-1} s^{-1}')
subplot(2,1,2)
legend('100 M^{-1} s^{-1}','1000 M^{-1} s^{-1}','10000 M^{-1} s^{-1}','100000 M^{-1} s^{-1}','1000000 M^{-1} s^{-1}')
sgtitle('k_{on}^i Sensitivity Analysis')
hold off

figure
dosages=[1e-12,1e-9,1e-6,1e-4];
C=['-k','-b','-r','-g'];
t_span_2=linspace(0,2.*60.*60,1000000);
ode_options_2=odeset('RelTol',1e-8,'AbsTol',[1e-8]);
for i=1:4
    y0_iter=[1.43e-22, 0, dosages(i), 0, 0, 0.0606, 0, 0];
    [t,y2]=ode15s(@model,t_span_2, y0_iter, ode_options, k);
    semilogy(t,y2(:,4),C(i))
    hold on
    xlim([3500 7200])
    xlabel('Time (s)')
    ylabel("Fluorescent Probe-Enzyme Complex (M)")
end
legend('1 pM','1 nM','1 uM','100 uM')
title('Dosage Analysis ([PQ]_0)')
hold off

figure
dosages=[0,1e-6,1e-5,1e-4,1e-3];
C=['-k','-b','-r','-g','-y'];
for i=1:5
    y0_iter=[1.43e-22, 0, dosages(i), 0, 0, 0.0606, 0, 0];
    [t,y2]=ode15s(@model,t_span, y0_iter, ode_options, k);
    semilogy(t,y2(:,8),C(i))
    xlim([3500 7200])
    xlabel('Time (s)')
    ylabel("Cleaved Elastin (M)")
    hold on
end
legend('None','1 uM','10 uM','100 uM','1 mM')
title('Dosage Effect on Elastin Cleavage')
hold off

delivery_times=[1.*3600,3.*3600,5.*3600,10.*3600];
t_span_3=linspace(0,12.*60.*60,1000000);
ode_options_3 = odeset('RelTol',1e-8,'AbsTol',[1e-8]);
C=['-k','-b','-r','-g'];
figure
hold all
for i=1:4
    k_iter=[2.614e-4, 1.01e14, 1e6, .27, 9.33, 1e5, 2.306e-9, delivery_times(i),5.35e-6];
    [t3,y3]=ode15s(@model,t_span_3, y0, ode_options_3, k_iter);
    subplot(2,1,1)
    hold on
    plot(t3,y3(:,4),C(i))
    xlabel('Time (s)')
    ylabel("Fluorescent Probe-Enzyme Complex (M)")
    subplot(2,1,2)
    hold on
    plot(t3,y3(:,8),C(i))
    xlabel('Time (s)')
    ylabel("Cleaved Elastin (M)")
end
subplot(2,1,1)
legend('1 hr','3 hr','5 hr', '10 hr')
subplot(2,1,2)
legend('1 hr','3 hr','5 hr', '10 hr')
sgtitle('Delivery Time Analysis')
hold off

function dydt = model(t,y,k)
    dydt = zeros(8,1);
    kg = k(1);
    kep = k(2);
    kon_i = k(3);
    kcat_s = k(4);
    koff_s = k(5);
    kon_s = k(6);
    K = k(7);
    tau=k(8);
    k_probe = k(9);

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
    if t>tau
        dydt(2,:) = kep*B+kcat_s*SE+koff_s*SE-kon_i*PQ*E-kon_s*S*E;
        dydt(3,:) = k_probe*PQ-kon_i*PQ*E;
        dydt(4,:) = kon_i*PQ*E;
        dydt(5,:) = kon_i*PQ*E;
    else
        dydt(2,:) = kep*B+kcat_s*SE+koff_s*SE-kon_s*S*E;
        dydt(3,:)=0;
        dydt(4,:)=0;
        dydt(5,:)=0;
    end
    dydt(6,:) = K+koff_s*SE-kon_s*S*E;
    dydt(7,:) = kon_s*S*E-koff_s*SE-kcat_s*SE;
    dydt(8,:) = kcat_s*SE;
end