import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
from dataclasses import dataclass
from math import *
from sys import *
from Turbo_dataclass import *

class Radial_comp:
  def __init__(self, Condition, Design):
    self.Condition = Condition
    self.Design = Design
    
    self.com_diagram = pd.DataFrame({'ns': [0.163, 0.222, 0.318, 0.39, 0.475, 0.645],
                                     'ds': [17, 12.1, 8.1, 6.55, 5.5, 4.05],
                                    'eff': [50, 60, 70, 75, 80, 85]})
    
    self.tur_diagram = pd.DataFrame({'ns': [0.132, 0.205, 0.48, 0.545, 0.61, 0.81, 1],
                                     'ds': [10.4, 7.8, 3.7, 3.45, 3.2, 2.7, 2.35],
                                    'eff': [70, 80, 90, 95, 90, 80, 70]})
    
  def __call__(self):
    
    gas = 'R245FA'
    

    # Compressor inlet, outlet conditions
    in_t = 55+273.15 
    in_p = 0.34*1.0e6
    out_p = 1.0*1.0e6
    m_dot = 5
    
    Fluid_in = Def_Condition()
    Fluid_out = Def_Condition()
    
    Fluid_in.T = in_t
    Fluid_in.p = in_p
    Fluid_in.fluid = gas
    Fluid_in.m = m_dot
    
    Fluid_out.p = out_p
    Fluid_out.fluid = gas
    Fluid_out.m = m_dot
    
    ns = 0.645
    
    (rpm_tune, dia_tune, tip_tune, ns_tune, ds_tune, H_ad) = self.nsds_conversion(Fluid_in, Fluid_out, ns)

    a = 1
    
    '''for n in range(self.Design.n_stage):
      if n == 0:
        self.Preprocess(self.Condition, self.Design)
        Stage = self.radial_compressor_design_initialization(self.Condition, self.Design, self.Design.P_ratio_vec[n])
        Stage = self.radial_compressor_convergence(self.Condition, self.Design, Stage)'''
        
        
  #def nsds_finding_com(Condition, Design):
  #  Condition.        
  def nsds_conversion(self, Fluid_1, Fluid_2, ns):
    m_to_fit = 3.28084
    kg_to_lb = 2.205
    rho_conv = kg_to_lb/m_to_fit**3
    g_en = 9.80665 * m_to_fit
    
    mdot_en = Fluid_1.m * kg_to_lb # kg to lb
    
    h_in = CP.PropsSI("H","T",Fluid_1.T,"P",Fluid_1.p,Fluid_1.fluid)
    s_in = CP.PropsSI("S","T",Fluid_1.T,"P",Fluid_1.p,Fluid_1.fluid)
    h_out = CP.PropsSI("H","P",Fluid_2.p,"S",s_in,Fluid_2.fluid)
    H_ad = (h_out-h_in) / (g_en/m_to_fit) * m_to_fit
    vol_en = mdot_en/(CP.PropsSI("D","T",Fluid_1.T, "P", Fluid_1.p, Fluid_1.fluid)*rho_conv)
    
    omega = ns/sqrt(vol_en)*(H_ad*g_en)**(3/4)
    rpm = omega*60/2/pi
    
    rpm_tune = rpm//60*60
    omega_tune = rpm_tune/60*2*pi
    ns_tune = omega_tune*sqrt(vol_en)/(H_ad*g_en)**(3/4)
    (ds_tune, eff_tune, err_code) = self.nsds_diagram(ns_tune)
    dia_en = ds_tune * sqrt(vol_en) / ((H_ad*g_en)**0.25)
    dia_tune = dia_en / m_to_fit
    tip_tune = omega_tune*dia_tune/2
    
    return rpm_tune, dia_tune, tip_tune, ns_tune, ds_tune, H_ad
  
  def nsds_diagram(self, ns):
    ns_data = self.com_diagram.ns
    ds_data = self.com_diagram.ds
    eff_data = self.com_diagram.eff
    if ns < ns_data.iloc[0]:
      ds = 3.8
      eff = 80
      err_code = 1
    elif ns > ns_data.iloc[-1]:
      ds = 3.8
      eff = 80
      err_code = 2
    else:
      err_code = 0
      ns_lb_index = ns_data.loc[ns > ns_data].index[-1]
      ns_ub_index = ns_lb_index+1
      
      interp_x = (ns - ns_data[ns_lb_index])/(ns_data[ns_ub_index]-ns_data[ns_lb_index])
      ds = ds_data[ns_lb_index] + interp_x*(ds_data[ns_ub_index]-ds_data[ns_lb_index])
      eff = eff_data[ns_lb_index] + interp_x*(eff_data[ns_ub_index]-eff_data[ns_lb_index])
      
    return ds, eff, err_code
      
  def Preprocess(self, Condition, Design):
    # Turbomachinery Boundary Condition
    Condition.Ho_in = CP.PropsSI("H","T",Condition.To_in,"P",Condition.Po_in, Condition.gas)
    Condition.So_in = CP.PropsSI("S","T",Condition.To_in,"P",Condition.Po_in, Condition.gas)
    Condition.Ho_out_ideal = CP.PropsSI("H","P",Condition.Po_out, "S",Condition.So_in, Condition.gas)
    
    
    # Check if Ca brings the static condition below saturation Design
    try:
      (Check_Ts, Check_Ps) = Aux.stagnation_to_static(Condition.To_in, Condition.Po_in, Design.Ca, Condition.gas)
    except:
      print("The inlet velocity too high to make static temperature under the saturation temperature")
      exit()
      
    Tsat_at_Check_Ps = CP.PropsSI("T","P",Check_Ps,"Q",1.0,Condition.gas)
    if Check_Ts < Tsat_at_Check_Ps:
      print("The inlet velocity too high to make static temperature under the saturation temperature")
      exit()
    else:
      SS = CP.PropsSI("A","P",Check_Ps,"T",Check_Ts,Condition.gas)
      if Design.Ca > SS*0.9:
        print("The inlet velocity is too close to the speed of sound")
        exit()
      
  
  def radial_compressor_design_initialization(self, Condition, Design, P_ratio):
    Stage = Def_Stage()
    Stage.Ho1 = Condition.Ho_in
    Stage.To1 = Condition.To_in
    Stage.Po1 = Condition.Po_in
    
    
    Stage.So1 = CP.PropsSI("S","T",Stage.To1,"P",Stage.Po1,Condition.gas)
    Stage.Cr1 = Design.Ca 
    Stage.C1 = Stage.Cr1/cos(Design.alpha1*pi/180) # axial composition
    Stage.Cw1 = Stage.C1*sin(Design.alpha1*pi/180) # Whirl composition
    
    
    (Stage.Ts1, Stage.Ps1) = Aux.stagnation_to_static(Stage.To1, Stage.Po1, Stage.C1, Condition.gas)
    Stage.rho1 = CP.PropsSI("D","T",Stage.To1,"P",Stage.Po1, Condition.gas)
    
    Stage.D1_tip = sqrt(Design.D1_root**2+4*Condition.mdot/(pi*Stage.rho1*Stage.Cr1))
    Stage.D1 = (Design.D1_root+Stage.D1_tip)/2
    
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
    Stage.Ho2_ideal = CP.PropsSI("H","P",Stage.Po2,"S",Stage.So1,Condition.gas)
    Stage.Ho2 = Stage.Ho2_ideal  
    Stage.To2 = CP.PropsSI("T","H",Stage.Ho2, "P", Stage.Po2, Condition.gas)
    
    Stage = self.impeller_outlet_state_calculation(Condition, Design, Stage, Design.Ca)
    
    # Initial Assumption
    Stage.To3 = Stage.To2
    Stage.Po3 = Stage.Po2
    
    Stage.Diffuser_depth = Stage.Vane_depth*Design.bstar 
    # bstar: ratio of inlet depth of vaneless diffuser to depth of impeller exit
    Stage.D2i = Stage.D2*Design.ratio_imp_dif 
    # imp_to_dif: ratio of diffuser inlet diameter to impeller exit diameter
    # Stage.D2i: Diffuser inlet diameter
    Stage.D3 = Stage.D2i*Design.ratio_dif_in_out
    # dif_in_to_out: ratio of diffuer outlet diamter to inlet diameter
    
    # Initial assumption of diffuser velocity
    Stage.Cr3 = Stage.Cr2/Design.ratio_dif_in_out
    Stage.Cw3 = Stage.Cw2/Design.ratio_dif_in_out
    Stage.C3 = sqrt(Stage.Cr3**2+Stage.Cw3**2)
    (Stage.Ts3, Stage.Ps3) = Aux.stagnation_to_static(Stage.To3, Stage.Po3, Stage.C3, Condition.gas)
    Stage.rho3 = CP.PropsSI("D","T",Stage.Ts3,"P",Stage.Ps3,Condition.gas)
    
    return Stage
  
  def radial_compressor_convergence(self, Condition, Design, Stage):
    a = 1
    Cr2_iter = Stage.Cr2
    while a:
      Stage = self.radial_compressor_loss(Condition, Design, Stage)
      Stage.Ho2 = Stage.Ho2_ideal+Stage.loss_external+Stage.loss_internal
      # Only pressure is affected by internal loss but the enthalpy is conserved
      # So only the internal loss is reflected when calculating pressure
      Stage.Po2 = CP.PropsSI("P","H",Stage.Ho2 - Stage.loss_internal,"S",Stage.So1,Condition.gas)
      Stage.So2 = CP.PropsSI("S","H",Stage.Ho2, "P", Stage.Po2, Condition.gas)
      Stage.To2 = CP.PropsSI("T","H",Stage.Ho2, "P", Stage.Po2, Condition.gas)
      
      (Stage.Ts2, Stage.Ps2) = Aux.stagnation_to_static(Stage.To2, Stage.Po2, Stage.C2, Condition.gas)
      Stage.rho2 = CP.PropsSI("D","T",Stage.Ts2,"P",Stage.Ps2, Condition.gas)
      
      # Station 2 update
      Cr2_update = Condition.mdot/(pi*Stage.D2*Stage.Vane_depth*Stage.rho2)
      Stage = self.impeller_outlet_state_calculation(Condition, Design, Stage, Cr2_update)
      
      Stage.Ho3 = Stage.Ho2
      Stage.Po3 = CP.PropsSI("P","H",Stage.Ho3-Stage.loss_total,"S",Stage.So1,Condition.gas)
      Stage.So3 = CP.PropsSI("S","H",Stage.Ho3,"P",Stage.Po3,Condition.gas)
      Stage.To3 = CP.PropsSI("T","H",Stage.Ho3,"P",Stage.Po3,Condition.gas)
      
      Stage.Diffuser_depth = Stage.Vane_depth*Design.bstar 
      Stage.D2i = Stage.D2*Design.ratio_imp_dif 
      Stage.D3 = Stage.D2i*Design.ratio_dif_in_out
      
      b = 1
      
      C3_iter = Stage.C3
      while b:
        Stage.Cr3 = Condition.mdot/(pi*Stage.D3*Stage.Diffuser_depth*Design.ratio_dif_unblocked*Stage.rho3)
        Stage.Cw3 = Stage.Cw2/Design.ratio_dif_in_out
        Stage.C3 = sqrt(Stage.Cr3**2+Stage.Cw3**2)
        (Stage.Ts3, Stage.Ps3) = Aux.stagnation_to_static(Stage.To3, Stage.Po3, Stage.C3, Condition.gas)
        Stage.rho3 = CP.PropsSI("D","T",Stage.Ts3, "P", Stage.Ps3, Condition.gas)
        
        err_C3 = abs((C3_iter-Stage.C3)/Stage.C3)
        C3_iter = Stage.C3
        
        if err_C3 < 1.0e-4:
          b = 0
      
      err_C2 = abs((Cr2_iter-Stage.Cr2)/Stage.Cr2)
      Cr2_iter = Stage.Cr2
      
      if err_C2 < 1.0e-4:
        a = 0
      
  def impeller_outlet_state_calculation(self, Condition, Design, Stage, Cr2):
    Stage.Cr2 = Cr2 # Initial Assumption of outlet axial velocity
    Stage.del_H = Stage.Ho2 - Stage.Ho1
    
    a = Stage.Slip_Factor
    b = -Stage.Cr2*tan(Design.BackSwept_beta*pi/180)
    c = -Stage.del_H-Stage.U1*Stage.Cw1
    
    Stage.U2 = (-b+sqrt(b**2-4*a*c))/2/a
    Stage.Cw2 = Stage.U2 + Stage.Cr2*tan(Design.BackSwept_beta*pi/180) - Stage.U2*(1-Stage.Slip_Factor)
    Stage.C2 = sqrt(Stage.Cw2**2+Stage.Cr2**2)
    Stage.alpha2 = 180/pi*atan(Stage.Cw2/Stage.Cr2)
    
    Stage.Wr2 = Stage.Cr2
    Stage.Ww2 = Stage.U2-Stage.Cw2
    Stage.W2 = sqrt(Stage.Wr2**2+Stage.Ww2**2)
    
    (Stage.Ts2, Stage.Ps2) = Aux.stagnation_to_static(Stage.To2,Stage.Po2,Stage.C2,Condition.gas)
    Stage.rho2 = CP.PropsSI("D","T",Stage.Ts2,"P",Stage.Ps2,Condition.gas)
    Stage.mu2 = CP.PropsSI("V","T",Stage.Ts2,"P",Stage.Ps2,Condition.gas)
    Stage.D2 = 2*Stage.U2/(2*pi*Condition.rpm/60)
    Stage.Vane_depth = Condition.mdot/(Stage.rho2*Stage.Cr2)/(pi*Stage.D2) 
    # Flow area = impeller exit circumference*Vane depth
    
    return Stage
      
  def radial_compressor_loss(self, Condition, Design, Stage):
    #---------- Internal loss models ----------#
    
    if Design.Loss_incidence == 1:
      # Conrad
      f_inc = 0.7
      Stage.loss_incidence = Design.Incidence_tune*(f_inc*Stage.Ww1**2/2)
    else:
      Stage.loss_incidence = 0
    
    if Design.Loss_blade_loading == 1:
      # Coppage
      ratio_W = Stage.W2/Stage.W1_tip
      ratio_D = Stage.D1_tip/Stage.D2
      Stage.Df = 1-ratio_W+ratio_W*0.75*Stage.del_H/Stage.U2**2/(Design.n_vane/pi*(1-ratio_D)+2*ratio_D)
      Stage.loss_bladeloading = Design.Blade_loading_tune*(0.05*Stage.Df**2*Stage.U2**2)
    else:
      Stage.loss_bladeloading = 0
      
    if Design.Loss_skin_friction == 1:
      # Jansen
      Cf = 0.005
      Stage.Lb = Stage.D1_tip/2+Stage.D2-Design.D1_root
      Stage.Dh = pi*(Stage.D1_tip**2-Design.D1_root**2)/(pi*(Design.D1_root+Stage.D1_tip)+2*Design.n_vane*(Stage.D1_tip-Design.D1_root))
      W_bar = (Stage.Cw1+Stage.C2+Stage.W1_tip+2*Stage.W1_hub+3*Stage.W2)/8
      Stage.loss_skinfriction = Design.Skin_friction_tune*(2*Cf*Stage.Lb/Stage.Dh*W_bar**2)
    else:
      Stage.loss_skinfriction = 0
    
    if Design.Loss_clearance == 1:
      # Jansen
      C1c = 0.6*Design.Clearance/Stage.Vane_depth*Stage.Cw2
      C2c = (Stage.D1_tip**2-Design.D1_root**2)/(Stage.D2-Stage.D1_tip)/(1+Stage.rho2/Stage.rho1)
      C3c = sqrt(2*pi/Stage.Vane_depth/Design.n_vane*Stage.Cw2*Stage.Cr1*C2c)
      Stage.loss_clearance = Design.Clearance_tune * (C1c*C3c)
    else:
      Stage.loss_clearance = 0
      
    if Design.Loss_mixing == 1:
      # Johnston & Dean
      C1m = 1/(1+(tan(Stage.alpha2*pi/180))**2)*Stage.C2**2/2
      C2m = ((1-Design.eps_wake-Design.bstar)/(1-Design.eps_wake))**2
      Stage.loss_mixing = Design.Mixing_tune*(C1m*C2m)
    else:
      Stage.loss_mixing = 0
    #------------------------------------------#
    
    #---------- External loss models ----------#
    if Design.Loss_disk_friction == 1:
      Re_df = Stage.rho2*Stage.U2*Stage.D2/2/Stage.mu2
      if Re_df > 3.0e5:
        f_df = 0.0622/Re_df**0.2
      else:
        f_df = 2.67/Re_df**0.5
      
      Stage.loss_diskfriction = Design.Disk_friction_tune*(f_df*0.5*(Stage.rho1+Stage.rho2)*Stage.D2**2*Stage.U2**3/16/Condition.mdot)
    else:
      Stage.loss_diskfriction = 0
      
    # Recirculation loss <Oh et al.>    
    if Design.Loss_recirculation == 1:
      ratio_W = Stage.W2/Stage.W1_tip
      ratio_D = Stage.D1_tip/Stage.D2
      Stage.Df = 1-ratio_W+ratio_W*0.75*Stage.del_H/Stage.U2**2/(Design.n_vane/pi*(1-ratio_D)+2*ratio_D)
      Stage.loss_recirculation = Design.Recirculation_tune * (2e-5*sinh(3.5*(Stage.alpha2/180*pi)**3)*Stage.Df**2*Stage.U2**2)
    else:
      Stage.loss_recirculation = 0
    
    if Design.Loss_leakage == 1:
      Stage.Lb = Stage.D1_tip/2+Stage.D2-Design.D1_root
      r_bar = (Stage.D2+Stage.D1_tip)/4
      b_bar = (Stage.D1_tip-Design.D1_root+Stage.Vane_depth)/2
      del_Pcl = Condition.mdot/2*(Stage.D2*Stage.Cw2-Stage.D1_tip*Stage.Cw1)/Design.n_vane/r_bar/b_bar/Stage.Lb
      Ucl = 0.816*sqrt(2*del_Pcl/Stage.rho2)
      mdot_cl = Stage.rho2*Design.n_vane*Design.Clearance*Stage.Lb*Ucl
      Stage.loss_leakage = Design.Leakage_tune*mdot_cl*Ucl*Stage.U2/2/Condition.mdot
    else:
      Stage.loss_leakage = 0
      
    if Design.Loss_windage == 1:
      P_windage = 1.379e+6
      rho_max = CP.PropsSI("D","T", 320, "P", P_windage, Condition.gas)
      Re_r = Design.rotor_r*(Design.rotor_gap*rho_max*Condition.rpm/60)/Design.mu
      Cd = 0.0030818
      Stage.loss_windage = Design.Windage_tune*(pi*Cd*rho_max*Design.rotor_r**4*(Condition.rpm*2*pi/60)**3*Design.rotor_length)
    else:
      Stage.loss_windage = 0
      
    # Total internal loss
    Stage.loss_internal = Stage.loss_incidence + Stage.loss_bladeloading + Stage.loss_skinfriction + Stage.loss_clearance + Stage.loss_mixing
    
    # Total external loss
    Stage.loss_external = Stage.loss_diskfriction + Stage.loss_recirculation + Stage.loss_leakage + Stage.loss_windage
    
    # Total loss
    Stage.loss_total = Stage.loss_internal + Stage.loss_external
    
    return Stage
    
class Aux:
  @staticmethod
  def static_to_stagnation(Ts, Ps, v, gas, method="r"):
    if method == "r":
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
    elif method == "i":
      gamma = CP.PropsSI("Cpmass","T",Ts, "P", Ps, gas)/CP.PropsSI("Cvmass","T",Ts, "P", Ps, gas)
      SS = CP.PropsSI("A","T",Ts, "P", Ps, gas)
      mach = v/SS
      Po = Ps*((1+(gamma-1)/2*mach**2)**(gamma/(gamma-1)))
      To = Ts*(1+(gamma-1)/2*mach**2)
    
    return (To, Po)

  @staticmethod
  def stagnation_to_static(To, Po, v, gas):
    Ts = To*0.8
    Ps = Po*0.8
    a = 1
    f = 0.01
    n = 0
    while a:
      n = n+1
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
        if n < 100:
          Ts = X_new[0]
          Ps = X_new[1]
        else:
          a = 0
        
    
    return(Ts, Ps)
    
if __name__ == '__main__':
  Condition = Def_Condition()
  Design = Def_Design()
  
  Condition.gas = "air"
  Condition.mdot = 10  # kg/sec
  Condition.rpm = 50000  # rev/min
  Condition.To_in = 30 + 273.15  # Inlet stagnation temperature (K)
  Condition.Po_in = 0.1 * 1E6  # Inlet stagnation pressure (Pa)
  Condition.Po_out = 0.3 * 1E6 # Outlet stagnation pressure (Pa)
  
  # Specified Design values
  Design.P_ratio_vec = [Condition.Po_out/Condition.Po_in]
  Design.n_stage = 1  # number of stages
  Design.Ca = 100 # Axial velocity(m/s)
  Design.n_vane = 20  # number of vanes
  
  (Ts, Ps) = Aux.stagnation_to_static(Condition.To_in, Condition.Po_in, Design.Ca, Condition.gas )
  Condition.rho_in = CP.PropsSI("D","T",Ts,"P",Ps,Condition.gas)
  
  Design.D1_root = sqrt(Condition.mdot/(Condition.rho_in*Design.Ca*pi/4*2.24)) # Impeller eye root diameter (m)
  Design.Clearance = 0.000254*1.0 # Clearance (m)

  # Loss model
  Design.Loss_incidence = 1
  Design.Loss_blade_loading = 1
  Design.Loss_skin_friction = 1
  Design.Loss_clearance = 1
  Design.Loss_mixing = 1
  Design.Loss_disk_friction = 1
  Design.Loss_recirculation = 1
  Design.Loss_leakage = 1
  Design.Loss_windage = 0

  # Loss model tunning factors
  Design.Incidence_tune = 1
  Design.Blade_loading_tune = 1
  Design.Skin_friction_tune = 1
  Design.Clearance_tune = 1
  Design.Mixing_tune = 1
  Design.Disk_friction_tune = 1
  Design.Recirculation_tune = 1
  Design.Leakage_tune = 1

  # Assumed design parameters
  Design.eps_wake = 0.45 # Wake fraction: 
  Design.bstar = 1.2 # ratio of inlet depth of vaneless diffuser to depth of impeller exit
  Design.ratio_imp_dif = 1.03 # ratio of diameter of impeller exit to inlet diameter of diffuser: vaneless diffuser area
  Design.ratio_dif_in_out = 1.8 # ratio of diameter of diffuser inlet to outlet
  
  Design.ratio_dif_unblocked = 0.8 # ratio of blocked surface area to flow area at the outlet of diffuser: actual flow area


  # Design options for turbomachineries
  Design.alpha1 = 30 # Inlet guide vane angle (degree)
  Design.BackSwept_beta = 50 # Backswept angle (degree): meridional corodinator
  
  Test_comp = Radial_comp(Condition, Design)
  Test_comp()
  