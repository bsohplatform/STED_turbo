import CoolProp.CoolProp as CP
import numpy as np
from dataclasses import dataclass
from math import *
from sys import *


class Radial_comp:
  def Turbo_design(self, Condition, Design, Stage):
    # Turbomachinery Boundary Condition
    Condition.Ho_in = CP.PropsSI("H","T",Condition.To_in,"P",Condition.Po_in, Condition.gas)
    Condition.So_in = CP.PropsSI("S","T",Condition.To_in,"P",Condition.Po_in, Condition.gas)
    Condition.Ho_out_ideal = CP.PropsSI("H","P",Condition.Po_out, "S",Condition.So_in, Condition.gas)
    
    # Calculation of Stage pressure ratio
    Design.total_P_ratio = Condition.Po_out / Condition.Po_in
    Design.stage_P_ratio = [Design.total_P_ratio**Design.ratio_fraction for i in range(Design.n_stage)]
    
    # Check if Ca brings the static condition below saturation state
    (Check_Ts, Check_Ps) = Aux.stagnation_to_static(Condition.To_in, Condition.Po_in, Design.Ca, Condition.gas)
    Tsat_at_Check_Ps = CP.PropsSI("T","P",Check_Ps,"Q",1.0,Condition.gas)
    if Check_Ts < Tsat_at_Check_Ps:
      print("The inlet velocity too high to make static temperature under the saturation temperature")
    
    Design.Ca_max = 
  
  def radial_compressor_stage_design(self, Inputs, Stage):
    Stage.Ts1 = Stage.Ts0
    Stage.Ps1 = Stage.Ps0
    Stage.Po10 = Stage.Po10_target
    Stage.mdot = Stage.stage_in_mdot
    
    # Start: Inlet velocity calculation
    Stage.impeller_in_mdot = Stage.mdot
    Stage = self.inlet_velocity_triangle_unshrouded(Inputs, Stage)
    
    Stage.leak_mdot = Stage.impeller_in_mdot*Inputs.leak_mdot_fraction
    Stage.beta2_guess = Inputs.beta2_method_value
    
    Stage.beta2i = Stage.beta2_guess
    beta2_iter_quit = 0
    while beta2_iter_quit =

  def alpha2i_search_for_Po10(Inputs, Stage):
  
  def unshrouded_impeller(Inputs, Stage):
    
  
  def inlet_velocity_triangle_unshrouded(self, Inputs, Stage):
    Stage.rho1 = CP.PropsSI("D","T",Stage.Ts1,"P",Stage.Ps1,Inputs.gas)
    temp_var = 
  
  def inlet_velocity_triangle_unshrouded_given_c(self, Inputs, Stage):
    Stage.r1_hub = Inputs.impeller_inlet_radius_calculation_method_value
    
    # finite thickness blade design inputs
    hth = Inputs.impeller_inlet_blade_hub_thickness
    tth = Inputs.impeller_inlet_blade_tip_thickness
    nfb = Inputs.n_full_blades
    
    # blade root capa check
    root_circum = Stage.r1_hub*2*pi
    if root_circum < nfb*hth:
      exit('inlet hub radius is not enough to contain blades')
      
    Stage.Cm1 = Inputs.impeller_input_velocity_design_value
    Stage.A1 = Stage.impeller_in_mdot / Stage.rho1 / Stage.Cm1
    Stage.r1_tip = (hth*nfb + tth*nfb + 2*((hth**2*nfb**2)/4+(hth*nfb**2*tth)/2 - 2*pi*hth*nfb*Stage.r1_hub+(nfb**2*tth**2)/4 - 2*pi*nfb*Stage.r1_hub*tth + 4*pi**2*Stage.r1_hub**2+4*pi*Stage.A1)**0.5)/(4*pi)
    
    # Impeller tip inlet velocity triangle
    Stage.U1_tip = Stage.omega*Stage.r1_tip
    Stage.C1_tip = Stage.Cm1/cos(Stage.alpha1*pi/180)
    Stage.Cw1_tip = Stage.C1_tip*sin(Stage.alpha1*pi/180)
    Stage.Ww1_tip = Stage.Cw1_tip - Stage.U1_tip
    Stage.W1_tip = sqrt(Stage.Cm1**2+Stage.Ww1_tip**2)
    Stage.beta1_tip = atan(Stage.Ww1_tip/Stage.Cm1)*180/pi
    
    Stage.r1 = sqrt((Stage.r1_tip**2+Stage.r1_hub**2)/2)
    Stage.U1 = Stage.omega*Stage.r1
    Stage.C1 = Stage.Cm1/cos(Stage.alpha1*pi/180)
    Stage.Cw1 = Stage.C1*sin(Stage.alpha1*pi/180)
    Stage.Ww1 = Stage.Cw1 - Stage.U1
    Stage.W1 = sqrt(Stage.Cm1**2+Stage.Ww1**2)
    Stage.beta1 = atan(Stage.Ww1/Stage.Cm1)*180/pi
    
    Stage.U1_hub = Stage.omega*Stage.r1_hub
    Stage.C1_hub = Stage.Cm1/cos(Stage.alpha1*pi/180)
    Stage.Cw1_hub = Stage.C1_hub*sin(Stage.alpha1*pi/180)
    Stage.Ww1_hub = Stage.Cw1_hub - Stage.U1_hub
    Stage.W1_hub = sqrt(Stage.Cm1**2+Stage.Ww1_hub**2)
    Stage.beta1_hub = atan(Stage.Ww1_hub/Stage.Cm1)*180/pi
    
    (Stage.To, Stage.Po) = Aux.static_to_stagnation(Stage.Ts1, Stage.Ps1, Stage.C1, Inputs.gas)
    Stage.Ho1 = CP.PropsSI("H","T",Stage.To1,"P",Stage.Po1,Inputs.gas)
    Stage.Hs1 = CP.PropsSI("H","T",Stage.Ts1,"P",Stage.Ps1,Inputs.gas)
    Stage.s1 = CP.PropsSI("S","T",Stage.To1,"P",Stage.Po1,Inputs.gas)
    
    return Stage

class Aux:
  @staticmethod
  def static_to_stagnation(Ts, Ps, v, gas):
      Z = CP.PropsSI("Z","T",Ts, "P", Ps, gas)
      gamma = CP.PropsSI("Cpmass","T",Ts, "P", Ps, gas)/CP.PropsSI("Cvmass","T",Ts, "P", Ps, gas)
      SS = CP.PropsSI("A","T",Ts, "P", Ps, gas)
      mach = v/SS
      dzdp = (CP.PropsSI("Z","T",Ts, "P", Ps*1.01, gas)-CP.PropsSI("Z","T",Ts, "P", Ps*0.99, gas))/(Ps*0.02)
      dzdt = (CP.PropsSI("Z","T",Ts*1.01, "P", Ps, gas)-CP.PropsSI("Z","T",Ts*0.99, "P", Ps, gas))/(Ts*0.02)
      beta_T = 1/Ps-1/Z*(dzdp)
      beta_P = 1/Ts-1/Z*(dzdt)
      ns = gamma/beta_T/Ps
      ms = (gamma-1)/gamma*beta_T/beta_P*Ps/Ts
      Po = Ps*((1+(ns-1)/2*mach**2)**(ns/(ns-1)))
      To = Ts*((1+(ns-1)/2*mach**2)**(ms*ns/(ns-1)))
      
      return (To, Po)

  @staticmethod
  def stagnation_to_static(To, Po, v, gas):
    Ts = To
    Ps = Po
    a = 1
    f = 0.01
    while a:
      
      (To_0, Po_0) = Aux.static_to_stagnation(Ts, Ps, v, gas)
      (To_1, Po_1) = Aux.static_to_stagnation(Ts*(1+f), Ps, v, gas)
      (To_2, Po_2) = Aux.static_to_stagnation(Ts, Ps*(1+f), v, gas)
      
      dTodTs_Ps = (To_1-To_0)/(Ts*f)
      dTodPs_Ts = (To_2-To_0)/(Ps*f)
      dPodTs_Ps = (Po_1-Po_0)/(Ts*f)
      dPodPs_Ts = (Po_2-Po_0)/(Ps*f)
      
      Terr = To_0-To
      Perr = Po_0-Po
      F = [Terr, Perr]
      X = [Ts, Ps]
      
      
      J = np.array([[dTodTs_Ps, dTodPs_Ts],[dPodTs_Ps, dPodPs_Ts]])
      P = np.dot(np.linalg.inv(J),F)
      
      X_new = X-P
      
      err = max(abs(Terr)/To, abs(Perr)/Po)
      if err < 1.0e-6:
        a = 0
      else:
        Ts = X_new[0]
        Ps = X_new[1]
    
    return(Ts, Ps)
    
if __name__ == '__main__':
  Ts = 300
  Ps = 101000
  v = 100
  gas = 'air'
  
  (To, Po) = Aux.static_to_stagnation(Ts, Ps, v, gas)
  print(To, Po)
  
  (Ts_val, Ps_val) = Aux.stagnation_to_static(To, Po, v, gas)
  print(Ts_val, Ps_val)