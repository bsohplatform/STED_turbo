import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
from dataclasses import dataclass
from math import *
from sys import *
from Turbo_dataclass import *
from copy import *
import matplotlib.pyplot as plt

CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP')

class Radial_comp:
  def __init__(self, Condition, Design, Plot):
    self.Condition = Condition
    self.Off_Condition = deepcopy(Condition)
    self.Design = Design
    self.Plot = Plot
    
    self.com_diagram = pd.DataFrame({'ns': [0.163, 0.222, 0.318, 0.39, 0.475, 0.645],
                                     'ds': [17, 12.1, 8.1, 6.55, 5.5, 4.05],
                                    'eff': [50, 60, 70, 75, 80, 85]})
    
    self.tur_diagram = pd.DataFrame({'ns': [0.132, 0.205, 0.48, 0.545, 0.61, 0.81, 1],
                                     'ds': [10.4, 7.8, 3.7, 3.45, 3.2, 2.7, 2.35],
                                    'eff': [70, 80, 90, 95, 90, 80, 70]})
    
  def __call__(self):
    Fluid_in = Def_Condition()
    Fluid_out = Def_Condition()
    self.Plot.pr_mat = []
    self.Plot.eff_mat = []
    self.Plot.mdot_mat = []
    
    Fluid_in.T = self.Condition.To_in
    Fluid_in.p = self.Condition.Po_in
    Fluid_in.fluid = self.Condition.gas
    Fluid_in.m = self.Condition.mdot
    
    Fluid_out.p = self.Condition.Po_out
    Fluid_out.fluid = Fluid_in.fluid
    Fluid_out.m = Fluid_in.m
    
    ns = 0.645
    
    (rpm_tune, dia_tune, tip_tune, ns_tune, ds_tune, H_ad) = self.nsds_conversion(Fluid_in, Fluid_out, ns)
    
    self.Condition.rpm = rpm_tune
     
    mdot_frac_list = [(self.Plot.mdot_plot_lb+(self.Plot.mdot_plot_ub-self.Plot.mdot_plot_lb)/(self.Plot.mdot_num_range-1)*n) for n in range(self.Plot.mdot_num_range)]
    rpm_frac_list = [(self.Plot.rpm_plot_lb+(self.Plot.rpm_plot_ub-self.Plot.rpm_plot_lb)/(self.Plot.rpm_num_range-1)*n) for n in range(self.Plot.rpm_num_range)]
     
    
    for n in range(self.Design.n_stage):
      if n == 0:
        self.Condition.Ho_in = CP.PropsSI("H","T",self.Condition.To_in,"P",self.Condition.Po_in, self.Condition.gas)
        self.Condition.So_in = CP.PropsSI("S","T",self.Condition.To_in,"P",self.Condition.Po_in, self.Condition.gas)
        self.Condition.Ho_out_ideal = CP.PropsSI("H","P",self.Condition.Po_out, "S",self.Condition.So_in, self.Condition.gas)

        self.Preprocess(self.Condition, self.Design)
        Stage, Geo = self.radial_compressor_design_initialization(self.Condition, self.Design, self.Design.P_ratio_vec[n])
        Stage, Geo = self.radial_compressor_convergence(self.Condition, self.Design, Stage, Geo)
        
        self.Plot.eff_on_design = (Condition.Ho_out_ideal-Stage.Ho1)/(Stage.Ho3-Stage.Ho1)*100
        self.Plot.pr_on_design = Stage.Po3/Stage.Po1
        
        Off_Stage = deepcopy(Stage)
        for rpm in rpm_frac_list:
          pr_list = []
          eff_list = []
          mdot_list = []
          for mdot in mdot_frac_list:
            self.Off_Condition.rpm = rpm*self.Condition.rpm
            self.Off_Condition.mdot = mdot*rpm*self.Condition.mdot
            Off_Stage = self.radial_compressor_off_design(self.Off_Condition, self.Design, Off_Stage, Geo)
            self.Off_Condition.Ho_out_ideal = CP.PropsSI("H","P",Off_Stage.Po3,"S",Off_Stage.So1, self.Off_Condition.gas)
            
            pr_list.append(Off_Stage.Po3/Off_Stage.Po1)
            eff_list.append((self.Off_Condition.Ho_out_ideal-Off_Stage.Ho1)/(Off_Stage.Ho3-Off_Stage.Ho1)*100)
            mdot_list.append(self.Off_Condition.mdot)
          
          self.Plot.pr_mat.append(pr_list)
          self.Plot.eff_mat.append(eff_list)
          self.Plot.mdot_mat.append(mdot_list)
          self.Plot.rpm_list = rpm_frac_list
      
      self.plot_performance_map(self.Plot, self.Condition)
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
    Geo = Def_Geo()
    
    Stage.Ho1 = Condition.Ho_in
    Stage.To1 = Condition.To_in
    Stage.Po1 = Condition.Po_in
    
    
    Stage.So1 = CP.PropsSI("S","T",Stage.To1,"P",Stage.Po1,Condition.gas)
    Stage.Cr1 = Design.Ca 
    
    Stage, Geo = self.impeller_inlet_state_calculation(Condition, Design, Stage, Geo, Stage.Cr1, "design")
    
    Geo.beta1 = 180/pi*atan((Stage.U1-Stage.Cw1)/Stage.Cr1)
    Geo.beta1_hub = 180/pi*atan((Stage.U1_hub-Stage.Cw1)/Stage.Cr1)
    Geo.beta1_tip = 180/pi*atan((Stage.U1_tip-Stage.Cw1)/Stage.Cr1)
    
    Geo.Slip_Factor = 1-sqrt(cos(Design.BackSwept_beta*pi/180))/Design.n_vane**0.7
    
    # Initial Assumption of Impeller outlet
    Stage.Po2 = Stage.Po1*P_ratio
    Stage.Ho2_ideal = CP.PropsSI("H","P",Stage.Po2,"S",Stage.So1,Condition.gas)
    Stage.Ho2 = Stage.Ho2_ideal  
    Stage.To2 = CP.PropsSI("T","H",Stage.Ho2, "P", Stage.Po2, Condition.gas)
    
    Stage.Cr2 = Design.Ca
    Stage, Geo = self.impeller_outlet_state_calculation(Condition, Design, Stage, Geo, Stage.Cr2, "design")
    
    # Initial Assumption
    Stage.To3 = Stage.To2
    Stage.Po3 = Stage.Po2
    
    Geo.Diffuser_depth = Geo.Vane_depth*Design.bstar 
    # bstar: ratio of inlet depth of vaneless diffuser to depth of impeller exit
    Geo.D2i = Geo.D2*Design.ratio_imp_dif 
    # imp_to_dif: ratio of diffuser inlet diameter to impeller exit diameter
    # Stage.D2i: Diffuser inlet diameter
    Geo.D3 = Geo.D2i*Design.ratio_dif_in_out
    # dif_in_to_out: ratio of diffuer outlet diamter to inlet diameter
    
    # Initial assumption of diffuser velocity
    Stage.Cr3 = Stage.Cr2/Design.ratio_dif_in_out
    Stage = self.diffuser_inlet_state_calculation(Condition, Design, Stage, Stage.Cr3)
    
    return Stage, Geo
  
  def radial_compressor_convergence(self, Condition, Design, Stage, Geo):
    a = 1
    Cr2_iter = Stage.Cr2
    while a:
      Stage, Geo = self.radial_compressor_loss(Condition, Design, Stage, Geo, "design")
      Stage.Ho2 = Stage.Ho2_ideal+Stage.loss_external+Stage.loss_internal
      # Only pressure is affected by internal loss but the enthalpy is conserved
      # So only the internal loss is reflected when calculating pressure
      Stage.Po2 = CP.PropsSI("P","H",Stage.Ho2 - Stage.loss_internal,"S",Stage.So1,Condition.gas)
      Stage.To2 = CP.PropsSI("T","H",Stage.Ho2, "P", Stage.Po2, Condition.gas)
      Stage.So2 = CP.PropsSI("S","T",Stage.To2, "P", Stage.Po2, Condition.gas)
      
      
      (Stage.Ts2, Stage.Ps2) = Aux.stagnation_to_static(Stage.To2, Stage.Po2, Stage.C2, Condition.gas)
      Stage.rho2 = CP.PropsSI("D","T",Stage.Ts2,"P",Stage.Ps2, Condition.gas)
      
      # Station 2 update
      Stage.Cr2 = Condition.mdot/(pi*Geo.D2*Geo.Vane_depth*Stage.rho2)
      Stage, Geo = self.impeller_outlet_state_calculation(Condition, Design, Stage, Geo, Stage.Cr2, "design")
      
      Stage.Ho3 = Stage.Ho2
      Stage.Po3 = CP.PropsSI("P","H",Stage.Ho3-Stage.loss_total,"S",Stage.So1,Condition.gas)
      Stage.So3 = CP.PropsSI("S","H",Stage.Ho3,"P",Stage.Po3,Condition.gas)
      Stage.To3 = CP.PropsSI("T","H",Stage.Ho3,"P",Stage.Po3,Condition.gas)
      
      Geo.Diffuser_depth = Geo.Vane_depth*Design.bstar
      Geo.D2i = Geo.D2*Design.ratio_imp_dif 
      Geo.D3 = Geo.D2i*Design.ratio_dif_in_out
      
      b = 1
      
      C3_iter = Stage.C3
      while b:
        Stage.Cr3 = Condition.mdot/(pi*Geo.D3*Geo.Diffuser_depth*Design.ratio_dif_unblocked*Stage.rho3)
        Stage = self.diffuser_inlet_state_calculation(Condition, Design, Stage, Stage.Cr3)
        err_C3 = abs((C3_iter-Stage.C3)/Stage.C3)
        C3_iter = Stage.C3
  
        if err_C3 < Condition.tol:
          b = 0
      
      err_C2 = abs((Cr2_iter-Stage.Cr2)/Stage.Cr2)
      Cr2_iter = Stage.Cr2
      
      if err_C2 < Condition.tol:
        a = 0
        
    return Stage, Geo
  
  def radial_compressor_off_design(self, Off_Condition, Design, Off_Stage, Geo):    
    # Initial guess
    Off_Stage.Cr1 = 4*Off_Condition.mdot/pi/(Geo.D1_tip**2-Design.D1_root**2)/Off_Stage.rho1
    Off_Stage.Cw1 = Off_Stage.Cr1*tan(Design.alpha1*pi/180)
    Off_Stage.C1 = sqrt(Off_Stage.Cr1**2+Off_Stage.Cw1**2)
    
    Off_Stage.U1 = 2*pi*Off_Condition.rpm/60*Geo.D1/2
    Off_Stage.U1_hub = 2*pi*Off_Condition.rpm/60*Design.D1_root/2
    Off_Stage.U1_tip = 2*pi*Off_Condition.rpm/60*Geo.D1_tip/2
    
    print('****************************************')
    a = 1
    iter_Cr1 = Off_Stage.Cr1
    n_c1 = 0
    while a:
      [Off_Stage.Ts1, Off_Stage.Ps1] = Aux.stagnation_to_static(Off_Stage.To1,Off_Stage.Po1,Off_Stage.C1,Condition.gas)
      Off_Stage.rho1 = CP.PropsSI("D","T",Off_Stage.Ts1,"P", Off_Stage.Ps1, Off_Condition.gas)
      
      Off_Stage.Cr1 = 4*Off_Condition.mdot/pi/(Geo.D1_tip**2-Design.D1_root**2)/Off_Stage.rho1
      Off_Stage, Geo = self.impeller_inlet_state_calculation(Off_Condition, Design, Off_Stage, Geo, Off_Stage.Cr1, "off-design")
      err_Cr1 = abs(iter_Cr1-Off_Stage.Cr1)/Off_Stage.Cr1
      iter_Cr1 = Off_Stage.Cr1
      n_c1 = n_c1+1
      print('[Cr1 Calculation Stage]   No_iter: %.3f   err_iter: %.3f' %(n_c1, err_Cr1))
      if err_Cr1 < Off_Condition.tol:
        a = 0
    
    Off_Stage.U2 = (2*pi*Off_Condition.rpm/60)*Geo.D2/2
    n_c2 = 0
    iter_Cr2 = Off_Condition.mdot/pi/Geo.D2/Geo.Vane_depth/Off_Stage.rho2
    b = 1
    while b:
      Off_Stage, Geo = self.radial_compressor_loss(Off_Condition, Design, Off_Stage, Geo, "off-design")
      Off_Stage.Ho2 = Off_Stage.Ho1+Off_Stage.U2*Off_Stage.Cw2 - Off_Stage.U1*Off_Stage.Cw1
      Off_Stage.del_H = Off_Stage.Ho2 - Off_Stage.Ho1
      Off_Stage.Po2 = CP.PropsSI("P","H",Off_Stage.Ho2-Off_Stage.loss_internal,"S",Off_Stage.So1, Off_Condition.gas)
      Off_Stage.To2 = CP.PropsSI("T","P",Off_Stage.Po2,"H",Off_Stage.Ho2, Off_Condition.gas)
      Off_Stage.So2 = CP.PropsSI("S","T",Off_Stage.To2, "P", Off_Stage.Po2, Off_Condition.gas)
      [Off_Stage.Ts2, Off_Stage.Ps2] = Aux.stagnation_to_static(Off_Stage.To2, Off_Stage.Po2, Off_Stage.C2, Off_Condition.gas)
      Off_Stage.rho2 = CP.PropsSI("D","T",Off_Stage.Ts2, "P", Off_Stage.Ps2, Off_Condition.gas)
      Off_Stage.Cr2 = Off_Condition.mdot/pi/Geo.D2/Geo.Vane_depth/Off_Stage.rho2
      Off_Stage, Geo = self.impeller_outlet_state_calculation(Off_Condition, Design, Off_Stage, Geo, Off_Stage.Cr2, "off-design")
      
      
      Off_Stage.Ho3 = Off_Stage.Ho2
      Off_Stage.Po3 = CP.PropsSI("P","H",Off_Stage.Ho2-Off_Stage.loss_total,"S",Off_Stage.So1,Off_Condition.gas)
      Off_Stage.To3 = CP.PropsSI("T","H",Off_Stage.Ho3,"P",Off_Stage.Po3,Off_Condition.gas)
      Off_Stage.So3 = CP.PropsSI("S","T",Off_Stage.To3, "P", Off_Stage.Po3, Off_Condition.gas)
      
      c = 1
      n_c3 = 0
      iter_C3 = Off_Stage.C3
      while c:
        Off_Stage.Cr3 = Off_Condition.mdot/pi/Geo.D3/Geo.Diffuser_depth/Design.ratio_dif_unblocked/Off_Stage.rho3
        Off_Stage = self.diffuser_inlet_state_calculation(Off_Condition, Design, Off_Stage, Off_Stage.Cr3)
        
        err_C3 = abs(iter_C3-Off_Stage.C3)/Off_Stage.C3
        iter_C3 = Off_Stage.C3
        n_c3 = n_c3+1
        print('[C3 Calculation Stage]   No_iter: %.3f   err_iter: %.3f' %(n_c3, err_C3))
        if err_C3 < Off_Condition.tol:
          c = 0
      
      err_Cr2 = abs(iter_Cr2-Off_Stage.Cr2)/Off_Stage.Cr2
      iter_Cr2 = Off_Stage.Cr2
      n_c2 = n_c2+1
      print('[Cr2 Calculation Stage]   No_iter: %.3f   err_iter: %.3f' %(n_c2, err_Cr2))
      if err_Cr2 < Off_Condition.tol:
        b = 0
      
      
    return Off_Stage
      
  def impeller_inlet_state_calculation(self, Condition, Design, Stage, Geo, Cr1, mode):
    Stage.Cw1 = Cr1*tan(Design.alpha1*pi/180)
    Stage.C1 = sqrt(Cr1**2+Stage.Cw1**2)
    
    if mode == 'design':
      (Stage.Ts1, Stage.Ps1) = Aux.stagnation_to_static(Stage.To1, Stage.Po1, Stage.C1, Condition.gas)
      Stage.rho1 = CP.PropsSI("D","T",Stage.Ts1,"P",Stage.Ps1, Condition.gas)
      Geo.D1_tip = sqrt(Design.D1_root**2+4*Condition.mdot/(pi*Stage.rho1*Stage.Cr1))
      Geo.D1 = (Design.D1_root+Geo.D1_tip)/2  
      
      Stage.U1 = Geo.D1/2*Condition.rpm/60*2*pi
      Stage.U1_hub = Design.D1_root/2*Condition.rpm/60*2*pi
      Stage.U1_tip = Geo.D1_tip/2*Condition.rpm/60*2*pi
    
    Stage.Ww1 = Stage.U1 - Stage.Cw1
    Stage.W1_tip = sqrt(Stage.Cr1**2+(Stage.U1_tip-Stage.Cw1)**2)
    Stage.W1_hub = sqrt(Stage.Cr1**2+(Stage.U1_hub-Stage.Cw1)**2)
    
    return Stage, Geo
      
  def impeller_outlet_state_calculation(self, Condition, Design, Stage, Geo, Cr2, mode):
    if mode == "design":
      Stage.del_H = Stage.Ho2 - Stage.Ho1
      a = Geo.Slip_Factor
      b = -Cr2*tan(Design.BackSwept_beta*pi/180)
      c = -Stage.del_H-Stage.U1*Stage.Cw1
      Stage.U2 = (-b+sqrt(b**2-4*a*c))/2/a
      
    Stage.Cw2 = Stage.U2*Geo.Slip_Factor - Cr2*tan(Design.BackSwept_beta*pi/180)
    Stage.C2 = sqrt(Stage.Cw2**2+Cr2**2)
    Stage.alpha2 = 180/pi*atan(Stage.Cw2/Cr2)
    
    Stage.Wr2 = Cr2
    Stage.Ww2 = Stage.U2-Stage.Cw2
    Stage.W2 = sqrt(Stage.Wr2**2+Stage.Ww2**2)
        
    (Stage.Ts2, Stage.Ps2) = Aux.stagnation_to_static(Stage.To2,Stage.Po2,Stage.C2,Condition.gas)
    Stage.rho2 = CP.PropsSI("D","T",Stage.Ts2,"P",Stage.Ps2,Condition.gas)
    Stage.mu2 = CP.PropsSI("V","T",Stage.Ts2,"P",Stage.Ps2,Condition.gas)
    
    if mode == "design":  
      Geo.D2 = 2*Stage.U2/(2*pi*Condition.rpm/60)  
      Geo.Vane_depth = Condition.mdot/(Stage.rho2*Cr2)/(pi*Geo.D2) 
      # Flow area = impeller exit circumference*Vane depth
    
    return Stage, Geo
  
  def diffuser_inlet_state_calculation(self, Condition, Design, Stage, Cr3):
    Stage.Cw3 = Stage.Cw2/Design.ratio_dif_in_out
    Stage.C3 = sqrt(Cr3**2+Stage.Cw3**2)
    (Stage.Ts3, Stage.Ps3) = Aux.stagnation_to_static(Stage.To3, Stage.Po3, Stage.C3, Condition.gas)
    Stage.rho3 = CP.PropsSI("D","T",Stage.Ts3,"P",Stage.Ps3,Condition.gas)
    
    return Stage
  
  def radial_compressor_loss(self, Condition, Design, Stage, Geo, mode):
    #---------- Internal loss models ----------#
    if mode == "design":
      Geo.Lb = Geo.D1_tip/2+Geo.D2-Design.D1_root
      Geo.Dh = pi*(Geo.D1_tip**2-Design.D1_root**2)/(pi*(Design.D1_root+Geo.D1_tip)+2*Design.n_vane*(Geo.D1_tip-Design.D1_root))
      
    if Design.Loss_incidence == 1:
      # Conrad
      f_inc = 0.7
      Stage.loss_incidence = Design.Incidence_tune*(f_inc*Stage.Ww1**2/2)
    else:
      Stage.loss_incidence = 0
    
    if Design.Loss_blade_loading == 1:
      # Coppage
      ratio_W = Stage.W2/Stage.W1_tip
      ratio_D = Geo.D1_tip/Geo.D2
      Df = 1-ratio_W+ratio_W*0.75*Stage.del_H/Stage.U2**2/(Design.n_vane/pi*(1-ratio_D)+2*ratio_D)
      Stage.loss_bladeloading = Design.Blade_loading_tune*(0.05*Df**2*Stage.U2**2)
    else:
      Stage.loss_bladeloading = 0
      
    if Design.Loss_skin_friction == 1:
      # Jansen
      Cf = 0.005
      W_bar = (Stage.Cw1+Stage.C2+Stage.W1_tip+2*Stage.W1_hub+3*Stage.W2)/8
      Stage.loss_skinfriction = Design.Skin_friction_tune*(2*Cf*Geo.Lb/Geo.Dh*W_bar**2)
    else:
      Stage.loss_skinfriction = 0
    
    if Design.Loss_clearance == 1:
      # Jansen
      C1c = 0.6*Design.Clearance/Geo.Vane_depth*Stage.Cw2
      C2c = (Geo.D1_tip**2-Design.D1_root**2)/(Geo.D2-Geo.D1_tip)/(1+Stage.rho2/Stage.rho1)
      C3c = sqrt(2*pi/Geo.Vane_depth/Design.n_vane*Stage.Cw2*Stage.Cr1*C2c)
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
      Re_df = Stage.rho2*Stage.U2*Geo.D2/2/Stage.mu2
      if Re_df > 3.0e5:
        f_df = 0.0622/Re_df**0.2
      else:
        f_df = 2.67/Re_df**0.5
      
      Stage.loss_diskfriction = Design.Disk_friction_tune*(f_df*0.5*(Stage.rho1+Stage.rho2)*Geo.D2**2*Stage.U2**3/16/Condition.mdot)
    else:
      Stage.loss_diskfriction = 0
      
    # Recirculation loss <Oh et al.>    
    if Design.Loss_recirculation == 1:
      ratio_W = Stage.W2/Stage.W1_tip
      ratio_D = Geo.D1_tip/Geo.D2
      Df = 1-ratio_W+ratio_W*0.75*Stage.del_H/Stage.U2**2/(Design.n_vane/pi*(1-ratio_D)+2*ratio_D)
      Stage.loss_recirculation = Design.Recirculation_tune * (2e-5*sinh(3.5*(Stage.alpha2/180*pi)**3)*Df**2*Stage.U2**2)
    else:
      Stage.loss_recirculation = 0
    
    if Design.Loss_leakage == 1:
      r_bar = (Geo.D2+Geo.D1_tip)/4
      b_bar = (Geo.D1_tip-Design.D1_root+Geo.Vane_depth)/2
      del_Pcl = Condition.mdot/2*(Geo.D2*Stage.Cw2-Geo.D1_tip*Stage.Cw1)/Design.n_vane/r_bar/b_bar/Geo.Lb
      Ucl = 0.816*sqrt(2*del_Pcl/Stage.rho2)
      mdot_cl = Stage.rho2*Design.n_vane*Design.Clearance*Geo.Lb*Ucl
      Stage.loss_leakage = Design.Leakage_tune*mdot_cl*Ucl*Stage.U2/2/Condition.mdot
    else:
      Stage.loss_leakage = 0
      
    if Design.Loss_windage == 1:
      P_windage = 1.379e+6
      rho_max = CP.PropsSI("D","T", 320, "P", P_windage, Condition.gas)
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
    
    return Stage, Geo
  
  def plot_performance_map(self, Plot, Condition):
    fig, axs = plt.subplots(2,1, constrained_layout=True)
    axs[0].set_title('Pressure ratio map ('+str(round(Condition.rpm))+' RPM)')
    axs[0].set_xlabel('Mass flow rate [kg/s]')
    axs[0].set_ylabel('Pressure ratio')
    axs[1].set_title('Total to Total efficiency map ('+str(round(Condition.rpm))+' RPM)')
    axs[1].set_xlabel('Mass flow rate [kg/s]')
    axs[1].set_ylabel('Efficiency [%]')
    for r in range(Plot.rpm_num_range):
      axs[0].plot(Plot.mdot_mat[r][:],Plot.pr_mat[r][:],label=str(round(Plot.rpm_list[r]*100))+'[%] RPM')
      axs[1].plot(Plot.mdot_mat[r][:],Plot.eff_mat[r][:],label=str(round(Plot.rpm_list[r]*100))+'[%] RPM')
    
    axs[0].plot(Condition.mdot, Plot.pr_on_design, '*', markersize=12, label='Design Point')
    axs[1].plot(Condition.mdot, Plot.eff_on_design, '*', markersize=12, label='Design Point')
    axs[0].legend(fontsize=8)
    fig.savefig('./Performance_map')
  
class Aux:
  @staticmethod
  def static_to_stagnation(Ts, Ps, v, gas, method="h"):
    if method == "r":
      f = 0.001
      Z =  CP.PropsSI("Z","T",Ts, "P", Ps, gas)
      gamma = CP.PropsSI("Cpmass","T",Ts, "P", Ps, gas)/CP.PropsSI("Cvmass","T",Ts, "P", Ps, gas)
      SS = CP.PropsSI("A","T",Ts, "P", Ps, gas)
      mach = v/SS
      dzdp = (CP.PropsSI("Z","T",Ts, "P", Ps*(1+f), gas)-CP.PropsSI("Z","T",Ts, "P", Ps*(1-f), gas))/(Ps*(2*f))
      dzdt = (CP.PropsSI("Z","T",Ts*(1+f), "P", Ps, gas)-CP.PropsSI("Z","T",Ts*(1-f), "P", Ps, gas))/(Ts*(2*f))
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
    
    elif method == "h":
      Hs = CP.PropsSI("H","T",Ts,"P",Ps,gas)
      S = CP.PropsSI("S","T",Ts,"P",Ps,gas)
      Ho = Hs+0.5*v**2
      Po = CP.PropsSI("P","H",Ho,"S",S,gas)
      To = CP.PropsSI("T","H",Ho,"S",S,gas)
    
    return (To, Po)
  
  @staticmethod
  def stagnation_to_static(To, Po, v, gas, method = "h"):
    if method == "n":
      Ts = To*0.9
      Ps = Po*0.9
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
    elif method == "h":
      Ho = CP.PropsSI("H","T",To,"P",Po,gas)
      S = CP.PropsSI("S","T",To,"P",Po,gas)
      Hs = Ho-0.5*v**2
      Ps = CP.PropsSI("P","H",Hs,"S",S,gas)
      Ts = CP.PropsSI("T","H",Hs,"S",S,gas)
    
    return(Ts, Ps)
    
if __name__ == '__main__':
  Condition = Def_Condition()
  Design = Def_Design()
  Plot = Def_Plot()
  
  Condition.tol = 1.0e-6
  Condition.gas = "REFPROP::R245FA"
  Condition.mdot = 5  # kg/sec
  Condition.To_in = 65 + 273.15  # Inlet stagnation temperature (K)
  Condition.Po_in = 0.34 * 1E6  # Inlet stagnation pressure (Pa)
  Condition.Po_out = 1.0 * 1E6 # Outlet stagnation pressure (Pa)
  
  # Specified Design values
  Design.P_ratio_vec = [Condition.Po_out/Condition.Po_in]
  Design.n_stage = 1  # number of stages
  Design.Ca = 50 # Axial velocity(m/s)
  Design.n_vane = 10  # number of vanes
  
  (Ts, Ps) = Aux.stagnation_to_static(Condition.To_in, Condition.Po_in, Design.Ca, Condition.gas )
  Condition.rho_in = CP.PropsSI("D","T",Ts,"P",Ps,Condition.gas)
  
  Design.D1_root = sqrt(Condition.mdot/(Condition.rho_in*Design.Ca*pi/4*2.24)) # Impeller eye root diameter (m)
  Design.Clearance = 0.00354*1.0 # Clearance (m)

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
  
  Plot.mdot_plot_lb = 0.6
  Plot.mdot_plot_ub = 1.1
  Plot.rpm_plot_lb = 0.6
  Plot.rpm_plot_ub = 1.1
  Plot.mdot_num_range = 6
  Plot.rpm_num_range = 6
  
  Test_comp = Radial_comp(Condition, Design, Plot)
  Stage = Test_comp()
