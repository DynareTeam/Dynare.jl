var Capital , Output, Labour, Consumption, Efficiency, efficiency ;

varexo EfficiencyInnovation;

parameters beta, theta, tau, alpha, Epsilon, delta, rho, effstar, sigma;

beta    =  0.990;
theta   =  0.357;
tau     = 30.000;
alpha   =  0.450;
delta   =  0.020;
rho     =  0.950;
effstar =  1.500;
sigma   =  0.010;
Epsilon =  0.500;

model;

#Psi = (Epsilon-1)/Epsilon;

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + sigma*EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency);

  // Eq. n°3:
  Output = Efficiency*(alpha*Capital(-1)^Psi+(1-alpha)*Labour^Psi)^(1/Psi);

  // Eq. n°4:
  Consumption + Capital - Output - (1-delta)*Capital(-1);

  // Eq. n°5:
  ((1-theta)/theta)*(Consumption/(1-Labour)) - (1-alpha)*Efficiency^((1-Psi))*(alpha*(Capital(-1)/Labour)^Psi+1-alpha)^((1-Psi)/Psi);

  // Eq. n°6:
  (((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption
  - beta*(Consumption(1)^theta*(1-Labour(1))^(1-theta))^(1-tau)/Consumption(1)*(alpha*Efficiency(1)^Psi*(Output(1)/Capital)^(1-Psi)+1-delta);

end;



steady_state_model;

  efficiency = 0;
  Efficiency = effstar;

  psi = (Epsilon-1)/Epsilon;

  Output_per_unit_of_Capital = ( (1/beta-1+delta) / (alpha*effstar^psi) )^(1/(1-psi));

  Consumption_per_unit_of_Capital = Output_per_unit_of_Capital-delta;

  Labour_per_unit_of_Capital =  ((Output_per_unit_of_Capital/Efficiency)^psi-alpha)^(1/psi)/(1-alpha)^(1/psi);

  gamma_1 = theta*(1-alpha)/(1-theta)*(Output_per_unit_of_Capital/Labour_per_unit_of_Capital)^(1-psi);
  gamma_2 = (Output_per_unit_of_Capital-delta)/Labour_per_unit_of_Capital;

  Labour = 1/(1+gamma_2/gamma_1);

  Output_per_unit_of_Labour=Output_per_unit_of_Capital/Labour_per_unit_of_Capital;

  Consumption_per_unit_of_Labour=Consumption_per_unit_of_Capital/Labour_per_unit_of_Capital;

  ShareOfCapital= alpha^(1/(1-psi))*effstar^psi/(1/beta-1+delta)^(psi/(1-psi));
  
  Consumption = Consumption_per_unit_of_Labour*Labour;
  Capital = Labour/Labour_per_unit_of_Capital;
  Output = Output_per_unit_of_Capital*Capital;

end;

shocks;
var EfficiencyInnovation = 1;
end;
