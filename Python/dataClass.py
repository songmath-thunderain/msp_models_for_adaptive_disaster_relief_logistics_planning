#!/usr/bin/env python
# coding: utf-8

class inputParams:
  def __init__(self,dissipate_option,absorbing_option,k_init,nbOS):
    self.dissipate_option = dissipate_option;
    self.absorbing_option = absorbing_option;
    self.k_init = k_init;
    self.nbOS = nbOS;

class solveParams:
  def __init__(self,max_iter,stall,cutviol_maxiter,time_limit,cutviol):
    self.max_iter = max_iter;
    self.stall = stall;
    self.cutviol_maxiter = cutviol_maxiter;
    self.time_limit = time_limit;
    self.cutviol = cutviol;

class hurricaneData:
  def __init__(self, Na, K, T, P_joint, states, absorbing_states, smallestTransProb, nodeLists, nodeScenList, nodeScenWeights):
    self.Na = Na;
    self.K = K;
    self.T = T;
    self.P_joint = P_joint;
    self.states = states;
    self.absorbing_states = absorbing_states;
    self.smallestTransProb = smallestTransProb;
    self.nodeLists = nodeLists;
    self.nodeScenList = nodeScenList;
    self.nodeScenWeights = nodeScenWeights;
    
class networkData:
  def __init__(self, Ni, Nj, fuel, cb, ca, ch, cp, p, q, x_cap, x_0, SCEN):
    self.N0 = Ni+1;
    self.Ni = Ni;
    self.Nj = Nj;
    self.fuel = fuel;
    self.cb = cb;
    self.ca = ca;
    self.ch = ch;
    self.cp = cp;
    self.p = p;
    self.q = q;
    self.x_cap = x_cap;
    self.x_0 = x_0;
    self.SCEN = SCEN;
    
#Ni: number of supply points (without MDC).
#Nj: number of demand points.
#N0: number of supply points including MDC.
#x_cap: the capacity of each supply point.
#cb: unit cost of transporting/rerouting relief items from the MDC/SP i to/between different SP i' 
#ch: unit cost for holding an item at SP i.
#cp: unit cost for relief item procurement.
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
