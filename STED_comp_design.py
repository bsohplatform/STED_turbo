import CoolProp.CoolProp as CP
from dataclasses import dataclass
from math import *
from sys import *

class Radial_comp:
  def inlet_velocity_triangle_unshrouded(self, Inputs, Stage):
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
    
class Aux:
  @staticmethod
  def static_to_stagnation(self, Ts, Ps, v, gas, sts_model):
      Z = CP.PropsSI("Z","T",Ts, "P", Ps, gas)
      gamma = CP.PropsSI("Cpmass","T",Ts, "P", Ps, gas)/CP.PropsSI("Cvmass","T",Ts, "P", Ps, gas)
      SS = property("A","T",Ts, "P", Ps, gas)
      mach = v/SS
      dzdp = (CP.PropsSI("Z","T",Ts, "P", Ps*1.01, gas)-CP.PropsSI("Z","T",Ts, "P", Ps*0.99, gas))/(Ps*0.02)
      dzdt = (CP.PropsSI("Z","T",Ts*1.01, "P", Ps, gas)-CP.PropsSI("Z","T",Ts*0.99, "P", Ps, gas))/(Ts*0.02)
      beta_T = 1/Ps-1/Z*(dzdp)
      beta_P = 1/Ts-1/Z*(dzdt)
      ns = gamma/beta_T/Ps
      ms = (gamma-1)/gamma*beta_T/beta_P*Ps/Ts
      Po = Ps*((1+(ns-1)/2*mach**2)**(ns/(ns-1)))
      To = Ts*((1+(ns-1)/2*mach**2)**(ms*ns/(ns-1)))