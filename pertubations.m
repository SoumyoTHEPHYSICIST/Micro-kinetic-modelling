function dydt = pertubations(t,y,pert)

kf=[1.5*10^-3 9.8*10^-4 4*10^-20 4.2*10^-13 6.4*10^-18 3.0*10^-11 9.5*10^-23 7.5*10^-17 1.1*10^-19 2.3*10^-20 9.3*10^-12 4.2*10^-19 3.9*10^-14 5.7*10^-11 2.2*10^-10 4.5*10^-10 6.5*10^-10 2.1*10^-10 7.8*10^-11 2.8*10^-11 2.0*10^-10 3.8*10^-11];

kb=[1.1*10^-9 5.0*10^-11 5.2*10^-15 2.6*10^-14 3.6*10^-14 6.1*10^-13 3.2*10^-17 2.9*10^-13 3.5*10^-13 4.4*10^-14 3.9*10^-13 2.4*10^-14 1.2*10^-17 5.0*10^-13 2.5*10^-14 1.4*10^-4 1.4*10^-4 12 9500 3.0*10^-5 80 1500];
% Calculate equilibrium constants (Keq) from unperturbed ki and ki_
    Keq1 = kf(1) / kb(1);
    Keq2 = kf(2) / kb(2);
    Keq3 = kf(3) / kb(3);
    Keq4 = kf(4) / kb(4);
    Keq5 = kf(5) / kb(5);
    Keq6 = kf(6) / kb(6);
    Keq7 = kf(7) / kb(7);
    Keq8 = kf(8) / kb(8);
    Keq9 = kf(9) / kb(9);
    Keq10 = kf(10) / kb(10);
    Keq11 = kf(11) / kb(11);
    Keq12 = kf(12) / kb(12);
    Keq13 = kf(13) / kb(13);
    Keq14 = kf(14) / kb(14);
    Keq15 = kf(15) / kb(15);
    Keq16 = kf(16) / kb(16);
    Keq17 = kf(17) / kb(17);
    Keq18 = kf(18) / kb(18);
    Keq19 = kf(19) / kb(19);
    Keq20 = kf(20) / kb(20);
    Keq21 = kf(21) / kb(21);
    Keq22 = kf(22) / kb(22);

    % Apply perturbations to the forward rate constants (ki)
    kf(1) = kf(1) + pert(1);
    kf(2) = kf(2) + pert(2);
    kf(3) = kf(3) + pert(3);
    kf(4) = kf(4) + pert(4);
    kf(5) = kf(5) + pert(5);
    kf(6) = kf(6) + pert(6);
    kf(7) = kf(7) + pert(7);
    kf(8) = kf(8) + pert(8);
    kf(9) = kf(9) + pert(9);
    kf(10) = kf(10) + pert(10);
    kf(11) = kf(11) + pert(11);
    kf(12) = kf(12) + pert(12);
    kf(13) = kf(13) + pert(13);
    kf(14) = kf(14) + pert(14);
    kf(15) = kf(15) + pert(15);
    kf(16) = kf(16) + pert(16);
    kf(17) = kf(17) + pert(17);
    kf(18) = kf(18) + pert(18);
    kf(19) = kf(19) + pert(19);
    kf(20) = kf(20) + pert(20);
    kf(21) = kf(21) + pert(21);
    kf(22) = kf(22) + pert(22);

    % Recalculate reverse rate constants (ki_) using Keq and perturbed ki
    kb(1) = kf(1) / Keq1;
    kb(2) = kf(2) / Keq2;
    kb(3) = kf(3) / Keq3;
    kb(4) = kf(4) / Keq4;
    kb(5) = kf(5) / Keq5;
    kb(6) = kf(6) / Keq6;
    kb(7) = kf(7) / Keq7;
    kb(8) = kf(8) / Keq8;
    kb(9) = kf(9) / Keq9;
    kb(10) = kf(10) / Keq10;
    kb(11) = kf(11) / Keq11;
    kb(12) = kf(12) / Keq12;
    kb(13) = kf(13) / Keq13;
    kb(14) = kf(14) / Keq14;
    kb(15) = kf(15) / Keq15;
    kb(16) = kf(16) / Keq16;
    kb(17) = kf(17) / Keq17;
    kb(18) = kf(18) / Keq18;
    kb(19) = kf(19) / Keq19;
    kb(20) = kf(20) / Keq20;
    kb(21) = kf(21) / Keq21;
    kb(22) = kf(22) / Keq22;

    dydt = zeros(16,1);
% Let Zi be the rates of the corresponding reactions 

Z1= kf(1)*y(1) - kb(1)*y(2)*y(3);

Z2= kf(2)*y(4) - kb(2)*y(5)*y(6);

Z3= kf(3)*y(4)*y(1) - kb(3)*y(7)*y(8);

Z4= kf(4)*y(1)*y(3) - kb(4)*y(2)*y(9);

Z5= kf(5)*y(1)*y(6) - kb(5)*y(10)*y(3);

Z6= kf(6)*y(1)*y(6) - kb(6)*y(2)*y(8);

Z7= kf(7)*y(1)*y(2) - kb(7)*y(11)*y(3);

Z8= kf(8)*y(5)*y(1) - kb(8)*y(12)*y(2);

Z9= kf(9)*y(5)*y(1) - kb(9)*y(7)*y(3);

Z10= kf(10)*y(5)*y(1) - kb(10)*y(13)*y(8);

Z11= kf(11)*y(4)*y(3) - kb(11)*y(5)*y(8);

Z12= kf(12)*y(4)*y(2) - kb(12)*y(14)*y(8);

Z13= kf(13)*y(12)*y(8) - kb(13)*y(4)*y(4);

Z14= kf(14)*y(7)*y(6) - kb(14)*y(14)*y(8);

Z15= kf(15)*y(7)*y(3) - kb(15)*y(14)*y(9);

Z16= kf(16)*y(3)*y(3) - kb(16)*y(9);

Z17= kf(17)*y(3)*y(6) - kb(17)*y(8);

Z18= kf(18)*y(3)*y(5) - kb(18)*y(12);

Z19= kf(19)*y(6)*y(6) - kb(19)*y(15);

Z20= kf(20)*y(6)*y(13) - kb(20)*y(7);

Z21= kf(21)*y(2)*y(2) - kb(21)*y(11);

Z22= kf(22)*y(2)*y(5) - kb(22)*y(7);

% Rates dydt(i) is taken as the rate of formation

dydt(1) = -(Z1+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10);

dydt(2) = (Z1+Z4-Z7+Z6+Z8-Z12-2*Z21-Z22);

dydt(3) = (Z1-Z4+Z5+Z7+Z9-Z11-Z15-2*Z16-Z17-Z18);

dydt(4) = (-Z2-Z3-Z11-Z12+Z13);

dydt(5) = (Z2-Z8-Z9-Z10+Z11-Z18-Z22);

dydt(6) = (Z2-Z5-Z6-Z14-Z17-2*Z19-Z20);

dydt(7) = (Z3+Z9+Z20+Z22);

dydt(8) = (Z3+Z6+Z10+Z11+Z12-Z13+Z14+Z17);

dydt(9) = (Z4+Z13+Z15+Z16);

dydt(10) = (Z5);

dydt(11) = (Z7+Z21);

dydt(12) = (Z8-Z13+Z18);

dydt(13) = (Z10-Z20);

dydt(14) = (Z12+Z14+Z15);

dydt(15) = (Z19);

dydt(16) = 0;

end