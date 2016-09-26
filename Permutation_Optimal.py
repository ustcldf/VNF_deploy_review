#!/usr/bin/python
from __future__ import division
import os
#from numpy import *
from scipy import *
import numpy as np
import operator
import time
import math
import copy
import json
import pickle
import random
from collections import defaultdict
from VNF_topo_generate import fat_tree
#from VNF_request_generate import REQUEST_NUM, get_request_seq, data_back
from VNF_DPT_source import x, VNF_SPE, brc_base, cpu_brc_per, mem_brc_per, CPU_TOTAL,MEM_TOTAL,h_band,band_alpha

ROOT_ID = -1
INFINITY = 1e10
EPSILON_F = 5.960465e-8

#x = np.arange(0,24,0.1)
class VNF_route(object):
    def __init__(self, ):
        self.switches = None
        self.adjacency = None #defaultdict(lambda:defaultdict(lambda:None))
        # [sw1][sw2] -> (distance, intermediate)
        self.path_map = defaultdict(lambda:defaultdict(lambda:(None,None)))
    
    def _calc_paths(self, ):
        def dump ():
          for i in sws.iteritems():
            for j in sws.iteritems():
              a = self.path_map[i][j][0]
              #a = adjacency[i][j]
              if a is None: a = "*"
              print a
            print
      
        sws = self.switches
        self.path_map.clear()
        for k,s in sws.iteritems():
          for j,value in self.adjacency[s.id].iteritems():
            if value is None: continue
            self.path_map[s.id][j] = (1,None)
          self.path_map[s.id][s.id] = (0,None) # distance, intermediate
        #dump()
      
        for k,s_1 in sws.iteritems():
          for i,s_2 in sws.iteritems():
            for j,s_3 in sws.iteritems():
              if self.path_map[s_2.id][s_1.id][0] is not None:
                if self.path_map[s_1.id][s_3.id][0] is not None:
                  # i -> k -> j exists
                  ikj_dist = self.path_map[s_2.id][s_1.id][0]+self.path_map[s_1.id][s_3.id][0]
                  if self.path_map[s_2.id][s_3.id][0] is None or ikj_dist < self.path_map[s_2.id][s_3.id][0]:
                    # i -> k -> j is better than existing
                    self.path_map[s_2.id][s_3.id] = (ikj_dist, s_1.id)
      
        #print "--------------------"
        #dump()
    
    def _get_raw_path (self, src, dst):
        """
        Get a raw path (just a list of nodes to traverse)
        """
        if len(self.path_map) == 0: self._calc_paths()
        if src is dst:
          # We're here!
          return []
        if self.path_map[src][dst][0] is None:
          return None
        intermediate = self.path_map[src][dst][1]
        if intermediate is None:
          # Directly connected
          return []
        return self._get_raw_path(src, intermediate) + [intermediate] + \
               self._get_raw_path(intermediate, dst)
    
    def _check_path (self, p):
        """
        Make sure that a path is actually a string of nodes
      
        returns True if path is valid
        """
        for a,b in zip(p[:-1],p[1:]):
          if self.adjacency[a][b] == None:
            return False
          if self.adjacency[b][a] == None:
            return False
        return True
    
    def _get_path (self, src, dst):
        """
        Gets a cooked path -- a list of (node)
        """
        # Start with a raw path...
        if src == dst:
          path = [src]
        else:
          path = self._get_raw_path(src, dst)
          if path is None: return None
          path = [src] + path + [dst]
        #print "what's src is:",src
        #print "what's dst is:",dst
        #print "what's in path is:",'\n',path
      
        assert self._check_path(path), "Illegal path!"
      
        return path
    
    def _get_switch(self,host_id,e_link):
        for id_1,iterm_1 in e_link.iteritems():
            for id_2,iterm_2 in iterm_1.iteritems():
                if id_2 == host_id:
                    return id_1
  
class VNF_Optimal(object):
    
    def __init__(self,):
        #self.container = []
        self.fattree = fat_tree()
        self.fattree.K = 16
        self.fattree.instant_topo()
        self.sws = self.fattree.switch
        self.adj = self.fattree.core_link
        self.core_link_nf = self.fattree.core_link
        self.e_link_nf = self.fattree.edge_link
        #self.request_lib_nf = data_back
        self.host_lib_nf = self.fattree.host_dict
        #self.NFR_lib_nf = NFCR_to_NFRs(request_lib_nf)
        self.host_info_ff_min = None
        self.host_num_min = 10000
        
        
    
    def verify_node_resource_left(self, path, link_map, request, host):
        #sw_id = route._get_switch(host.id)
        #path = route._get_path(ROOT_ID, sw_id)
        #modify host resource left
        req_len = len(request)               #here, request is a list
        req_cpu = 0
        req_mem = 0
        req_band = 0
        req_cpu += request[1].cpu
        req_mem += request[1].mem
        
        #modify the vnfs in each host
        vnf_in_req = []
        vnf_in_req.append(request[1].id)
        #print "vnf num in req", len(vnf_in_req), vnf_in_req
        vnf_add_new = []
        for vnf_i in vnf_in_req:
            if vnf_i in host.vnf_container:
                continue
            else:
                host.vnf_container.append(vnf_i)
                vnf_add_new.append(vnf_i)
                host.cpu_res -= cpu_brc_per
                host.mem_res -= mem_brc_per
        
        host.verify_cpu(req_cpu)
        host.verify_mem(req_mem)
    
    def check_value(self, list_ori, list_cmp):             #our data is a sequence, check value is to calc the break times of resource
        list_len = len(list_cmp)                #list_ori is the left resource sequence, list_cmp is the resource sequence that is to be setted to responding host
        break_t = 0
        for i in range(list_len):
            if list_ori[i] - list_cmp[i] < 0 :
                break_t += 1
        return break_t
        
    def check_positive(self, list_temp):
        list_len = len(list_temp)
        positive_num = 0
        for i in range(list_len):
            if list_temp[i] < 0:
                positive_num += 1
        return positive_num
        
    def calc_req_vnf(self, req_list):                #count the vnf species in req
        vnf_list = []
        req_len = len(req_list)
        for i in range(req_len):
            vnf_list.append(req_list[i].id)
        
        return vnf_list

    def check_end(self, path, link_map, request, host):          #judge if the request can be setted to host
        #host_backup = copy.deepcopy(host)
        req_len = len(request)             #request is a list
        req_cpu = 0
        req_mem = 0
        req_band = 0
        req_cpu = request[1].cpu
        req_mem = request[1].mem
        req_band = request[1].in_band + request[1].out_band                           #request band is determined by the first VNF
        vnf_species = request[1].id
        #print "vnf num in req", len(vnf_in_req), vnf_in_req
        vnf_add_new = []
        if vnf_species in host.vnf_container:
            pass
        else:
            host.vnf_container.append(vnf_species)
            vnf_add_new.append(vnf_species)
            host.cpu_res -= cpu_brc_per
            host.mem_res -= mem_brc_per
        host_vnf_container = copy.deepcopy(host.vnf_container)
        sw_band_bt_total = 0
        
        for a,b in zip(path[:-1],path[1:]):
            sw_band_bt_total  += self.check_value(link_map[a][b], req_band)
            #band_b_t_list.append(band_b_t)
        h_cpu_bt_total = self.check_value(host.cpu_res, req_cpu)                       #break times
        h_mem_bt_total = self.check_value(host.mem_res, req_mem)
        h_band_bt_total = self.check_value(host.band_res, req_band)
        
        
        #print "h_cpu_bt_total, h_mem_bt_total, h_band_bt_total",h_cpu_bt_total, h_mem_bt_total, h_band_bt_total
        if sw_band_bt_total == 0 and h_cpu_bt_total == 0 and h_mem_bt_total == 0 and h_band_bt_total == 0:                     #to be changed, a threshold can be defined
            #print "resource verified!"
            
            if len(vnf_add_new) != 0:
                for vnf_j in vnf_add_new:         #just check, no matter the host could hold the request or not, it'll should be restore back
                    host.vnf_container.remove(vnf_j)
                    host.cpu_res += cpu_brc_per
                    host.mem_res += mem_brc_per
                #host_vnf_container_2 = copy.deepcopy(host.vnf_container)
            
            return 1
        else:
            #host = copy.deepcopy(host_backup)
            if len(vnf_add_new) != 0:
                for vnf_k in vnf_add_new:         #if can not hold, restore the vnf_container
                    host.vnf_container.remove(vnf_k)
                    host.cpu_res += cpu_brc_per
                    host.mem_res += mem_brc_per
                #host_vnf_container_2 = copy.deepcopy(host.vnf_container)
                #print "after restore vnf in host_vnf_container",host_vnf_container_2
            
            return 0
        
    
    def verify_band_left_0(self, path, core_link, NFR, host): #the host does not contain the related NFCR
        
        if NFR[1].pre_vnf_id == None and NFR[1].next_vnf_id == None: #the NFCR contains only one VNF
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id == None and NFR[1].next_vnf_id != None:  #the VNF is the first
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id != None and NFR[1].next_vnf_id == None: #the VNF is the last
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id != None and NFR[1].next_vnf_id != None: #the VNF in the middle
            host.verify_band(NFR[1].in_band + NFR[1].out_band)
    
    def verify_band_left_1(self, path, core_link, NFR, host):  #the host contains the related NFCR
        vnf_in_NFCR = self.calc_req_vnf(host.req_container_dict[NFR[0]])
        if NFR[1].pre_vnf_id in vnf_in_NFCR and NFR[1].next_vnf_id not in vnf_in_NFCR: 
            host.band_res += NFR[1].in_band
            host.verify_band(NFR[1].out_band)
        if NFR[1].pre_vnf_id not in vnf_in_NFCR and NFR[1].next_vnf_id in vnf_in_NFCR:
            host.band_res += NFR[1].out_band
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id not in vnf_in_NFCR and NFR[1].next_vnf_id not in vnf_in_NFCR:
            host.verify_band(NFR[1].in_band + NFR[1].out_band)
        if NFR[1].pre_vnf_id in vnf_in_NFCR and NFR[1].next_vnf_id in vnf_in_NFCR:
            host.band_res += NFR[1].in_band + NFR[1].out_band
    
    def NF_map(self, route, core_link, e_link, NFR_lib, host_lib, NFR_id_list):     #map the request to host using Next Fit
        print "Next Fit"
        host_id = 0
        req_abandoned = []
        request_num = len(NFR_lib)
        #sw_id = route._get_switch(host_id,e_link)
        #path = route._get_path(ROOT_ID, sw_id)
        for i in NFR_id_list:
            sw_id = route._get_switch(host_id,e_link)
            path = route._get_path(ROOT_ID, sw_id)
            if self.check_end(path, core_link, NFR_lib[i], host_lib[host_id]):
                if NFR_lib[i][0] in host_lib[host_id].req_container:
                    host_lib[host_id].req_container_dict[NFR_lib[i][0]].append(NFR_lib[i][1])
                    self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_id])
                    self.verify_band_left_1(path,core_link,NFR_lib[i],host_lib[host_id])
                else:
                    host_lib[host_id].req_container.append(NFR_lib[i][0])
                    host_lib[host_id].req_container_dict[NFR_lib[i][0]] = [NFR_lib[i][1]]
                    self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_id])
                    self.verify_band_left_0(path,core_link,NFR_lib[i],host_lib[host_id])
            else:
                host_id += 1
                sw_id = route._get_switch(host_id,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_lib[i], host_lib[host_id]):
                    if NFR_lib[i][0] in host_lib[host_id].req_container:
                        host_lib[host_id].req_container_dict[NFR_lib[i][0]].append(NFR_lib[i][1])
                        self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_id])
                        self.verify_band_left_1(path,core_link,NFR_lib[i],host_lib[host_id])
                    else:
                        host_lib[host_id].req_container.append(NFR_lib[i][0])
                        host_lib[host_id].req_container_dict[NFR_lib[i][0]] = [NFR_lib[i][1]]
                        self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_id])
                        self.verify_band_left_0(path,core_link,NFR_lib[i],host_lib[host_id])
                else:
                    req_abandoned.append(i)
        print "final host id (NF):", host_id
        print "abandoned requeats are:",req_abandoned
        return (host_id,host_lib)     
        
    def search_min_dict(self, dict_exp):
        min_val = INFINITY
        for key, value in dict_exp.iteritems():
            if value <= min_val:
                min_val = value
                min_id = key
        return min_id
    

    def search_max_dict(self, dict_exp):
        max_val = 0
        for key, value in dict_exp.iteritems():
            if value >= max_val:
                max_val = value
                max_id = key
        return max_id
    
    def Find_optimal(self, NFR_lib, A, NFR_id_list_len, cur):  
        if cur==NFR_id_list_len:  
            Results = A
            print "A",A
            route = VNF_route()
            route.switches = self.sws
            route.adjacency = self.adj
            core_link = copy.deepcopy(self.core_link_nf)
            e_link = copy.deepcopy(self.e_link_nf)
            host_lib = copy.deepcopy(self.host_lib_nf)
            host_info_nf = self.NF_map(route, core_link, e_link, NFR_lib, host_lib, A)
            if host_info_nf[0] <= self.host_num_min:
                self.host_num_min = host_info_nf[0]
                self.host_info_nf_min = host_info_nf
        else:  
            #for i in xrange(1,n+1):
            for i in xrange(0,NFR_id_list_len):
                ok=True  
                for j in xrange(0,cur):  
                    if A[j]==i:  
                        ok=False  
                        break  
                if ok:  
                    A[cur]=i
                    #print "A[cur]", cur, A[cur]
                    self.Find_optimal(NFR_lib, A, NFR_id_list_len, cur+1)  
    
    
 


def NFCR_to_NFRs(request_lib):
    NFR_num = 0
    NFR_dict = {}
    for NFCR_id,NFCR_value in request_lib.iteritems():
        for NFR_i in NFCR_value:
            NFR_tuple = (NFCR_id,NFR_i)
            NFR_dict[NFR_num] = NFR_tuple
            NFR_num += 1
    return NFR_dict


    

    
def search_list_max(list_temp):
    max_value = 0
    for i in range(len(list_temp)):
        if list_temp[i] >= max_value:
            max_value = list_temp[i]
    return max_value

def search_list_min(list_temp):
    min_value = INFINITY
    for i in range(len(list_temp)):
        if list_temp[i] <= min_value:
            min_value = list_temp[i]
    return min_value

def reset_list(list_temp, new_value):
    for i in range(len(list_temp)):
        list_temp[i] = new_value


def pre_process_req_lib(request_lib):
    for req_id,req_value in request_lib.iteritems():
        max_cpu = 0
        max_band = 0
        max_mem = 0
        for i in range(len(req_value)):
            max_cpu = search_list_max(req_value[i].cpu)
            max_mem = search_list_max(req_value[i].mem)
            max_band = search_list_max(req_value[i].band)
            reset_list(req_value[i].cpu, max_cpu)
            reset_list(req_value[i].mem, max_mem)
            reset_list(req_value[i].band, max_band)
    
def check_resource_lefted_list(list_temp):
    verified_times = 0
    for i in range(len(list_temp)):
        if list_temp[i] <= 0:
            verified_times += 1
    return verified_times

def get_aver_list(list_temp):
    aver_value = 0
    total_value = 0
    for i in range(len(list_temp)):
        total_value += list_temp[i]
    aver_value = total_value/len(list_temp)
    return aver_value
    



def after_pod_adjust_check(host_lib, host_id_ori):
    #host_empty_list = []
    print "let's start to check the data in the host_lib"
    print "ID - host_id, A - len of req_container, B - keys of req_container_dict, C - len of vnf_container"
    print "D1 - max of cpu lefted + cpu_in_req, D2 - min of cpu lefted + cpu_in_req, D3 - aver of cpu lefted + cpu_in_req"
    print "E1 - max of mem lefted + mem_in_req, E2 - min of mem lefted + mem_in_req, E3 - aver of mem lefted + mem_in_req"
    print "F1 - max of band lefted + band_in_req, F2 - min of band lefted + band_in_req, F3 - aver of band lefted + band_in_req"
    print "ID**A**B**C**D1**D2**D3**E1**E2**E3**F1**F2**F3"
    for i in range(host_id_ori+1):
        #cpu_verified_t = check_resource_lefted_list(host_lib[i].cpu_res)
        #mem_verified_t = check_resource_lefted_list(host_lib[i].mem_res)
        #band_verified_t = check_resource_lefted_list(host_lib[i].band_res)
        A = len(host_lib[i].req_container)
        #B = len(host_lib[i].req_container_dict)
        B = host_lib[i].req_container_dict.keys()
        #C = len(host_lib[i].vnf_container)
        C = host_lib[i].vnf_container
        req_cpu = 0
        req_mem = 0
        req_band = 0
        if len(host_lib[i].req_container_dict.keys()):
            for req_id, req_value in host_lib[i].req_container_dict.iteritems():
                if len(req_value):
                    for vnf in req_value:
                        req_cpu += vnf.cpu
                        req_mem += vnf.mem
                    req_band += req_value[0].in_band
          
        
        
        D3 = cpu_aver_left = get_aver_list(host_lib[i].cpu_res + req_cpu + len(host_lib[i].vnf_container) * cpu_brc_per)
        E3 = mem_aver_left = get_aver_list(host_lib[i].mem_res + req_mem + len(host_lib[i].vnf_container) * mem_brc_per)
        F3 = band_aver_left = get_aver_list(host_lib[i].band_res + req_band)
        D1 = cpu_max_left = search_list_max(host_lib[i].cpu_res + req_cpu + len(host_lib[i].vnf_container) * cpu_brc_per)
        E1 = mem_max_left = search_list_max(host_lib[i].mem_res + req_mem + len(host_lib[i].vnf_container) * mem_brc_per)
        F1 = band_max_left = search_list_max(host_lib[i].band_res + req_band)
        D2 = cpu_min_left = search_list_min(host_lib[i].cpu_res + req_cpu + len(host_lib[i].vnf_container) * cpu_brc_per)
        E2 = mem_min_left = search_list_min(host_lib[i].mem_res + req_mem + len(host_lib[i].vnf_container) * mem_brc_per)
        F2 = band_min_left = search_list_min(host_lib[i].band_res + req_band)
        
        print i,'**',A,'**',B,'**',C,'**',D1,'**',D2,'**',D3,'**',E1,'**',E2,'**',E3,'**',F1,'**',F2,'**',F3


def host_resource_left_check(host_lib, host_id_ori, host_saved_list):
    host_empty_list = []
    print "let's start to check the data in the host_lib"
    print "ID - host_id, A - len of req_container, B - keys of req_container_dict, C - len of vnf_container"
    print "D1 - max of cpu lefted, D2 - min of cpu lefted, D3 - aver of cpu lefted"
    print "E1 - max of mem lefted, E2 - min of mem lefted, E3 - aver of mem lefted"
    print "F1 - max of band lefted, F2 - min of band lefted, F3 - aver of band lefted"
    print "ID**A**B**C**D1**D2**D3**E1**E2**E3**F1**F2**F3"
    for i in range(host_id_ori+1):
        #cpu_verified_t = check_resource_lefted_list(host_lib[i].cpu_res)
        #mem_verified_t = check_resource_lefted_list(host_lib[i].mem_res)
        #band_verified_t = check_resource_lefted_list(host_lib[i].band_res)
        A = len(host_lib[i].req_container)
        B = len(host_lib[i].req_container_dict)
        #B = host_lib[i].req_container_dict.keys()
        C = len(host_lib[i].vnf_container)
        #C = host_lib[i].vnf_container
        D3 = cpu_aver_left = get_aver_list(host_lib[i].cpu_res)
        E3 = mem_aver_left = get_aver_list(host_lib[i].mem_res)
        F3 = band_aver_left = get_aver_list(host_lib[i].band_res)
        D1 = cpu_max_left = search_list_max(host_lib[i].cpu_res)
        E1 = mem_max_left = search_list_max(host_lib[i].mem_res)
        F1 = band_max_left = search_list_max(host_lib[i].band_res)
        D2 = cpu_min_left = search_list_min(host_lib[i].cpu_res)
        E2 = mem_min_left = search_list_min(host_lib[i].mem_res)
        F2 = band_min_left = search_list_min(host_lib[i].band_res)
        
        print i,'**',A,'**',B,'**',C,'**',D1,'**',D2,'**',D3,'**',E1,'**',E2,'**',E3,'**',F1,'**',F2,'**',F3
        if len(host_lib[i].req_container) == 0:
            host_empty_list.append(i)
        #print "the host has been checked is:",i
        #print "the verified times of cpu", cpu_verified_t
        #print "the verified times of mem", mem_verified_t
        #print "the verified times of band", band_verified_t
        #print "the aver value of lefted cpu", cpu_aver_left
        #print "the aver value of lefted mem", mem_aver_left
        #print "the aver value of lefted band", band_aver_left
    """
    print "Now, we compare the empty list"
    for host_id in host_empty_list:
        if host_id in host_saved_list:
            print host_id, "is in host_saved_list"
        else:
            print "the host_id", host_id,"is not in host_saved_list"
            print "the host_empty_list is not the same with host_saved_list"
    """
    

def VNF_simulation_Optimal(data_back):
    fattree = fat_tree()
    fattree.K = 16
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    request_lib_nf = copy.deepcopy(data_back)
    
    NFR_lib_nf = NFCR_to_NFRs(request_lib_nf)
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_Optimal()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    print "Optimal Solution"
    A = []
    cur = 0
    for i in range(len(NFR_lib_nf)):
        A.append(len(NFR_lib_nf)+1)
    NFR_id_list_len = len(NFR_lib_nf)
    print "NFR_id_list_len", NFR_id_list_len
    req_2_host.__init__()
    req_2_host.Find_optimal(NFR_lib_nf, A, NFR_id_list_len, cur)
    print "*********start save info*********"
    
    host_id_nf = req_2_host.host_info_nf_min[0]
    host_lib_nf = req_2_host.host_info_nf_min[1]
    
    with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_NF_rw.pickle', 'w') as f:
        pickle.dump(host_lib_nf,f)
    saved_info_dict_nf = {}
    saved_info_dict_nf[0] = host_id_nf+1
    print "the used host num(nf) is:",host_id_nf+1
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_Optimal.json', 'w') as f:
        json.dump(saved_info_dict_nf,f)
    
    #print "check host_lib_nf"
    #saved_host_nf = []
    #host_resource_left_check(host_lib_nf,host_id_nf,saved_host_nf)

    
    print "***********Optimal Solution is over**********"
    



"""
if __name__ == '__main__':
    print "************Start the Process***************"
    time_start = time.time()
   
    #A=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  
    A = []
    Muta_num = 11
    for i in range(Muta_num):
        A.append(100)
    
    cur = 0  
    n = Muta_num
    print_permutation(n,A,cur)
    
    time_end = time.time()
    time_used = time_end - time_start
    
    time_second = time_used % 60
    time_ms = (time_used-time_second) / 60
    time_mininute = time_ms % 60
    time_hs = (time_ms - time_mininute) / 60
    time_hour = time_hs % 24
    time_day = (time_hs - time_hour) / 24
    print "****************************"
    print "Time of the process:"
    print "Used time (seconds)", time_used
    print time_day, "Days",time_hour, "Hours", time_mininute, "Mininutes", time_second, "Seconds"
    print "***********Over*************"
"""