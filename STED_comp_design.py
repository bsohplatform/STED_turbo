import CoolProp.CoolProp as CP
import numpy as np
from dataclasses import dataclass
from math import *
from sys import *
from Turbo_dataclass import *

class Radial_comp:
  def Turbo_design(self, Condition, Design):
    # Turbomachinery Boundary Condition
    Condition.Ho_in = CP.PropsSI("H","T",Condition.To_in,"P",Condition.Po_in, Condition.gas)
    Condition.So_in = CP.PropsSI("S","T",Condition.To_in,"P",Condition.Po_in, Condition.gas)
    Condition.Ho_out_ideal = CP.PropsSI("H","P",Condition.Po_out, "S",Condition.So_in, Condition.gas)
    
    # Calculation of Stage pressure ratio
    Design.total_P_ratio = Condition.Po_out / Condition.Po_in
    Design.stage_P_ratio = [Design.total_P_ratio**Design.ratio_fraction for i in range(Design.n_stage)]
    
    # Check if Ca brings the static condition below saturation Design
    (Check_Ts, Check_Ps) = Aux.stagnation_to_static(Condition.To_in, Condition.Po_in, Design.Ca, Condition.gas)
    Tsat_at_Check_Ps = CP.PropsSI("T","P",Check_Ps,"Q",1.0,Condition.gas)
    if Check_Ts < Tsat_at_Check_Ps:
      print("The inlet velocity too high to make static temperature under the saturation temperature")
      exit()
    else:
      SS = CP.PropsSI("A","P",Check_Ps,"T",Check_Ts,Condition.gas)
      if Design.Ca > SS*0.9:
        print("The inlet velocity is too close to the speed of sound")
        exit()
      
  
  def radial_compressor_stage_design(self, Condition, Design, P_ratio, gas):
    Stage = Def_Stage()
    Stage.Ho1 = Condition.Ho_in
    Stage.To1 = Condition.To_in
    Stage.Po1 = Condition.Po_in
    
    Stage.So1 = CP.PropsSI("S","T",Stage.To1,"P",Stage.Po1,gas)
    Stage.Cr1 = Design.Ca 
    Stage.C1 = Design.Cr1/cos(Design.alpha1*pi/180) # axial composition
    Stage.Cw1 = Design.C1*sin(Design.alpha1*pi/180) # Whirl composition
    
    
    (Stage.Ts1, Stage.Ps1) = Aux.stagnation_to_static(Stage.To1, Stage.Po1, gas)
    Stage.rho1 = CP.PropsSI("D","T",Stage.To1,"P",Stage.Po1, gas)
    
    Stage.D1_tip = sqrt(Design.D1_root**2+4*Condition.mdot/(pi*Stage.rho1*Stage.Cr1))
    Stage.D1 = (Design.D1_root+Design.D1_tip)/2
    
    Stage.U1 = Stage.D1/2*Condition.rpm/60*2*pi
    Stage.U1_hub = Design.D1_root/2*Condition.rpm/60*2*pi
    Stage.U1_tip = Stage.D1_tip/2*Condition.rpm/60*2*pi
    
    Stage.Ww1 = Stage.U1 - Stage.Cw1
    Stage.W1_tip = sqrt(Stage.Cr1**2+(Stage.U1_tip-Stage.Cw1)**2)
    Stage.W1_hub = sqrt(Stage.Cr1**2+(Stage.U1_hub-Stage.Cw1)**2)
    
    Stage.beta1 = 180/pi*atan((Stage.U1-Stage.Cw1)/Stage.Cr1)
    Stage.beta1_hub = 180/pi*atan((Stage.U1_hub-Stage.Cw1)/Stage.Cr1)
    Stage.beta1_tip = 180/pi*atan((Stage.U1_tip-Stage.Cw1)/Stage.Cr1)
    
    Stage.Slip_Factor = 1-sqrt(cos(Design.BackSwept_beta*pi/180))/Design.n_vane**0.7
    
    # Initial Assumption of Impeller outlet
    Stage.Po2 = Stage.Po1*P_ratio
    Stage.Ho2_ideal = CP.PropsSI("H","P",Stage.Po2,"S",Stage.So1,gas)
    Stage.Ho2 = Stage.Ho2_ideal  
    Stage.To2 = CP.PropsSI("T","H",Stage.Ho2, "P", Stage.Po2, gas)
    
    Stage.Cr2 = Design.Ca # Initial Assumption of outlet axial velocity
    Stage.del_H = Stage.Ho2 - Stage.Ho1
    
    a = Stage.Slip_Factor
    b = -Stage.Cr2*tan(Design.BackSwept_beta*pi/180)
    c = -Stage.del_H-Stage.U1*Stage.Cw1
    
    Stage.U2 = (-b+sqrt(b**2-4*a*c))/2/a
    Stage.Cw2 = Stage.U2 + Stage.Cw2*tan(Design.BackSwept_beta*pi/180) - Stage.U2*(1-Stage.Slip_Factor)
    Stage.alpha2 = 180/pi*atan(Stage.Cw2/Stage.Cr2)
    
    Stage.Wr2 = Stage.Cr2
    Stage.Ww2 = Stage.U2-Stage.Cw2
    Stage.W2 = sqrt(Stage.Wr2**2+Stage.Ww2**2)
    
    (Stage.Ts2, Stage.Ps2) = Aux.stagnation_to_static(Stage.To2,Stage.Po2,Stage.C2,gas)
    Stage.rho2 = CP.PropsSI("D","T",Stage.Ts2,"P",Stage.Ps2,gas)
    Stage.mu2 = CP.PropsSI("V","T",Stage.Ts2,"P",Stage.Ps2,gas)
    Stage.D2 = 2*Stage.U2/(2*pi*Condition.rpm/60)
    Stage.Vane_depth = Condition.mdot/(Stage.rho2*Stage.Cr2)/(pi*Stage.D2) 
    # Flow area = impeller exit circumference*Vane depth
    
    # Initial Assumption
    Stage.To3 = Stage.To2
    Stage.Po3 = Stage.Po2
    
    Stage.Diffuser_depth = Stage.Vane_depth*Design.bstar 
    # bstar: ratio of inlet depth of vaneless diffuser to depth of impeller exit
    Stage.D2i = Stage.D2*Design.imp_to_dif 
    # imp_to_dif: ratio of diffuser inlet diameter to impeller exit diameter
    # Stage.D2i: Diffuser inlet diameter
    Stage.D3 = Stage.D2i*Design.dif_in_to_out
    # dif_in_to_out: ratio of diffuer outlet diamter to inlet diameter
    
    # Initial assumption of diffuser velocity
    Stage.Cr3 = Stage.Cr2/Design.ratio_dif_in_out
    Stage.Cw3 = Stage.Cw2/Design.ratio_dif_in_out
    Stage.C3 = sqrt(Stage.Cr3**2+Stage.Cw3**2)
    (Stage.Ts3, Stage.Ps3) = Aux.stagnation_to_static(Stage.To3, Stage.Po3, Stage.C3, gas)
    Stage.rho3 = CP.PropsSI("D","T",Stage.Ts3,"P",Stage.Ps3,gas)
    
    a = 1
    while a:
      Stage = radial_compressor_loss(Condition, Design, Stage)
      
  def radial_compressor_loss(self, Condition, Design, Stage):
    f_inc = 0.7
    
    

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