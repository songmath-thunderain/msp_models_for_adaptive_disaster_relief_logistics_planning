#!/usr/bin/env python
# coding: utf-8

class solveParams:
  def __init__(self,max_iter,stall,cutviol_maxiter,time_limit,cutviol):
    self.max_iter = max_iter;
    self.stall = stall;
    self.cutviol_maxiter = cutviol_maxiter;
    self.time_limit = time_limit;
    self.cutviol = cutviol;

class hurricaneData:
  def __init__(self, P_intensity, P_location, P_landfall, T, Tmin, P_joint, absorbing_states, P_terminals):
    self.P_intensity = P_intensity;
    self.P_location = P_location;
    self.P_landfall = P_landfall;
    self.T = T;
    self.Tmin = Tmin;
    self.P_joint = P_joint;
    self.absorbing_states = absorbing_states;
    self.P_terminals = P_terminals;
    
class networkData:
  def __init__(self, Ni, Nj, SP, DP, fuel, cb, ca, ch, h, p, q, D_max, x_cap, x_0, f_cap, SCEN):
    self.N0 = Ni+1;
    self.Ni = Ni;
    self.Nj = Nj;
    self.SP = SP;
    self.DP = DP;
    self.fuel = fuel;
    self.cb = cb;
    self.ca = ca;
    self.ch = ch;
    self.h = h;
    self.p = p;
    self.q = q;
    self.D_max = D_max;
    self.x_cap = x_cap;
    self.x_0 = x_0;
    self.f_cap = f_cap;
    self.SCEN = SCEN;
    
#Ni: number of supply points (without MDC).
#Nj: number of demand points.
#N0: number of supply points including MDC.
#x_cap: the capacity of each supply point.
#cb: unit cost of transporting/rerouting relief items from the MDC/SP i to/between different SP i' 
#ch: unit cost for holding an item at SP i.
#h: unit cost for purchasing a relief item.
#-----------------------------------------------------------------------------------------
#ca: unit cost of transporting relief items from the MDC/SP i to/between a demand point j.
#p: penalty cost for failing to serve the demand.
#q: salvage cost for each unit of overstock.
#-----------------------------------------------------------------------------------------
#Na: number of states in the intensity MC
#Nb: number of states in the locations MC
#K: number of states in the joint distrubutions between intesities and locations
#P_intesity: transition probability matrix of intensity MC
#P_location: transition probability matrix of location MC
#P_joint: transition probability matrix of the joint distrubution beween intensity and location MC
