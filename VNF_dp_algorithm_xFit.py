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

Cor_threshold = 1.0 #default is 0.3
Cor_threshold_shortest = 1.0

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

class VNF_PowerNets(object):
    
    def verify_node_resource_left(self, path, link_map, request, host):
        #sw_id = route._get_switch(host.id)
        #path = route._get_path(ROOT_ID, sw_id)
        #modify host resource left
        req_len = len(request)               #here, request is a list
        req_cpu = 0
        req_mem = 0
        req_band = 0
        req_cpu += copy.deepcopy(request[1].cpu)
        req_mem += copy.deepcopy(request[1].mem)
        
        #modify the vnfs in each host
        vnf_in_req = []
        vnf_in_req.append(request[1].id)
        #print "vnf num in req", len(vnf_in_req), vnf_in_req
        #vnf_add_new = []
        for vnf_i in vnf_in_req:
            if vnf_i in host.vnf_container:
                continue
            else:
                host.vnf_container.append(vnf_i)
                #vnf_add_new.append(vnf_i)
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
        req_cpu = copy.deepcopy(request[1].cpu)
        req_mem = copy.deepcopy(request[1].mem)
        req_band = copy.deepcopy(request[1].in_band + request[1].out_band)                           #request band is determined by the first VNF
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
        
        """
        #90_percentle threshold as the workload.
        Vector_len = len(req_cpu)
        cpu_90_percentile = self.get_90_percentile(req_cpu)
        mem_90_percentile = self.get_90_percentile(req_mem)
        band_90_percentile = self.get_90_percentile(req_band)
        
        req_cpu = np.array([float(cpu_90_percentile) for i in range(Vector_len)])
        req_mem = np.array([float(mem_90_percentile) for i in range(Vector_len)])
        req_band = np.array([float(band_90_percentile) for i in range(Vector_len)])
        """
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
        
    def get_90_percentile(self, value_list):
        value_list_copy = copy.deepcopy(value_list)
        value_list_copy.sort()
        list_len = len(value_list_copy)
        seq_90_percentile = math.ceil(0.9*list_len)
        return value_list_copy[seq_90_percentile]
        
    
    def calc_host_resource(self, host_candinate_list):
        host_resource = {}
        for i in range(len(host_candinate_list)):
            h_value = host_candinate_list[i]
            host_resource[h_value.id] = self.get_aver(h_value.cpu_res)/CPU_TOTAL + self.get_aver(h_value.mem_res)/MEM_TOTAL + self.get_aver(h_value.band_res)/h_band
        #min_h_id = self.search_min_dict(host_resource)
        return host_resource
    
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
    
    def sort_requests(self, NFR_lib):
        req_90_percentile_source_dict = {}
        for i in range(len(NFR_lib)):
            VNF_i_cpu = 0
            VNF_i_mem = 0
            VNF_i_band = 0
            VNF_i_cpu += self.get_90_percentile(NFR_lib[i][1].cpu)/CPU_TOTAL
            VNF_i_mem += self.get_90_percentile(NFR_lib[i][1].mem)/MEM_TOTAL
            VNF_i_band += self.get_90_percentile(NFR_lib[i][1].in_band + NFR_lib[i][1].out_band)/h_band  #the band resource deserves recosideration
            req_90_percentile_source_dict[i] = VNF_i_cpu + VNF_i_mem + VNF_i_band
        req_id_list = []
        while(len(req_90_percentile_source_dict)):
            max_id = self.search_max_dict(req_90_percentile_source_dict)
            req_id_list.append(max_id)
            req_90_percentile_source_dict.pop(max_id)
        
        return req_id_list
    
    def calc_d_val_pearson(self, mtr1,mtr2):              #calc d value
        #print "*****calc d value*****"
        [row,col] = mtr1.shape
        #print "row,col",row,col
        #mtr_temp = np.zeros([row,col])
        #normalization
        #normalization
        """
        normal_max_dict1 = {}
        normal_max_dict2 = {}
        for k in range(row):
            normal_max_dict1[k] = 0
            normal_max_dict2[k] = 0
            for l in range(col):
                if mtr1[k][l] >= normal_max_dict1[k]:
                    normal_max_dict1[k] = mtr1[k][l]
                if mtr2[k][l] >= normal_max_dict2[k]:
                    normal_max_dict2[k] = mtr2[k][l]
        
        for m in range(row):
            for n in range(col):
                mtr1[m][n] = mtr1[m][n]/(normal_max_dict1[m]+EPSILON_F)
                mtr2[m][n] = mtr2[m][n]/(normal_max_dict2[m]+EPSILON_F)
        """ 
        
        #pearson correlation coefficient
        pearson = {}
        for i in range(row):
            pearson[i] = 0
            factor_1 = 0
            factor_2 = 0
            factor_3 = 0
            factor_4 = 0
            factor_5 = 0
            j = 0
            for j in range(col):
                factor_1 += mtr1[i][j]*mtr2[i][j]
                factor_2 += mtr1[i][j]
                factor_3 += mtr2[i][j]
                factor_4 += mtr1[i][j]*mtr1[i][j]
                factor_5 += mtr2[i][j]*mtr2[i][j]
            #print "sqrt_1:",col*factor_4-math.pow(factor_2,2)
            Cross_x = col*factor_4-math.pow(factor_2,2)
            Cross_y = col*factor_5-math.pow(factor_3,2)
            #print "i, Row,Col", i, row, col
            #print "col*factor_4-math.pow(factor_2,2)", col*factor_4-math.pow(factor_2,2)
            #print "col*factor_5-math.pow(factor_3,2)", col*factor_5-math.pow(factor_3,2)
            if Cross_x < 0 and Cross_x > -1*e-8:
                Cross_x = 0
                #print "vector len", len(mtr2[i]), "Mistake vector:", mtr2[i]
                #print "col*factor_5",col*factor_5, "math.pow(factor_3,2)", math.pow(factor_3,2)
            if Cross_y < 0 and Cross_y > -1*e-8:
                Cross_y = 0
            if Cross_x*Cross_y == 0:
                pearson[i] = 1
            else:
                pearson[i] = (col * factor_1 - factor_2 * factor_3)/(math.sqrt(Cross_x*Cross_y))
            #print "i",i,"Pearson", pearson[i]
        #d = 0
        #for i in range(row):
            #d += pearson[i]
        
        return pearson
    
    def Find_shortest_host_for_vnf(self, route, core_link, e_link, host_id_cur, host_list_c, NFR_temp, host_lib):
        i = 1
        print "Find_shortest_host_for_vnf for:", NFR_temp[0], NFR_temp[1].id
        print "host_id_cur", host_id_cur
        host_id_list = []
        for k in range(len(host_list_c)):
            host_id_list.append(host_list_c[k].id)
        print "len(host_list_c)", len(host_list_c)
        bound = max(host_id_cur, len(host_list_c)-host_id_cur)
        print "bound is:",bound
        #sw_id = route._get_switch(host_id,e_link)
        #path = route._get_path(ROOT_ID, sw_id)
        while(1):
            if (host_id_cur - i) in host_id_list:
                sw_id = route._get_switch(host_id_cur - i,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_temp, host_lib[host_id_cur - i]) == 1:
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_id_cur - i].CPU - host_lib[host_id_cur - i].cpu_res
                    h_mtr[1] = host_lib[host_id_cur - i].Mem - host_lib[host_id_cur - i].mem_res
                    h_mtr[2] = host_lib[host_id_cur - i].Band - host_lib[host_id_cur - i].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] & h_NFR_cor[2] in first while(1) 1", h_NFR_cor[0],h_NFR_cor[2]
                    if (h_NFR_cor[0] <= Cor_threshold) and (h_NFR_cor[2] <= Cor_threshold): #and (h_NFR_cor[1] <= Cor_threshold)  #+h_NFR_cor[1])/2
                        print "host_id",host_id_cur - i
                        return host_id_cur - i
            if (host_id_cur + i) in host_id_list:
                sw_id = route._get_switch(host_id_cur + i,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_temp, host_lib[host_id_cur + i]) == 1:
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_id_cur + i].CPU - host_lib[host_id_cur + i].cpu_res
                    h_mtr[1] = host_lib[host_id_cur + i].Mem - host_lib[host_id_cur + i].mem_res
                    h_mtr[2] = host_lib[host_id_cur + i].Band - host_lib[host_id_cur + i].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] & h_NFR_cor[2] in first while(1) 2", h_NFR_cor[0],h_NFR_cor[2]
                    if (h_NFR_cor[0] <= Cor_threshold) and (h_NFR_cor[2] <= Cor_threshold):#and (h_NFR_cor[1] <= Cor_threshold) #+h_NFR_cor[1])/2
                        print "host_id",host_id_cur + i
                        return host_id_cur + i
            i += 1
            print "i in the first while(1)", i
            if i > bound:           # can not find the host satisfies all the cor constraints :(h_NFR_cor[0] <= Cor_threshold) and (h_NFR_cor[2] <= Cor_threshold) 
                break
        print "start to relax the cor constraints"
        i = 1 
        while(1):            #relax the cor constraints to:(h_NFR_cor[0] <= 0.3) and (h_NFR_cor[1] <= 0.3)
            if (host_id_cur - i) in host_id_list:
                sw_id = route._get_switch(host_id_cur - i,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_temp, host_lib[host_id_cur - i]) == 1:
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_id_cur - i].CPU - host_lib[host_id_cur - i].cpu_res
                    h_mtr[1] = host_lib[host_id_cur - i].Mem - host_lib[host_id_cur - i].mem_res
                    h_mtr[2] = host_lib[host_id_cur - i].Band - host_lib[host_id_cur - i].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in second while(1) 1", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
                        print "host_id",host_id_cur - i
                        return host_id_cur - i
            if (host_id_cur + i) in host_id_list:
                sw_id = route._get_switch(host_id_cur + i,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_temp, host_lib[host_id_cur + i]) == 1:
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_id_cur + i].CPU - host_lib[host_id_cur + i].cpu_res
                    h_mtr[1] = host_lib[host_id_cur + i].Mem - host_lib[host_id_cur + i].mem_res
                    h_mtr[2] = host_lib[host_id_cur + i].Band - host_lib[host_id_cur + i].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in second while(1) 2", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
                        print "host_id",host_id_cur + i
                        return host_id_cur + i
            i += 1
            print "i in the second while(1)", i
            if i > bound:           # can not find the host satisfies all the cor constraints :(h_NFR_cor[0] <= Cor_threshold)
                break        
        print "start to relax the cor constraints 2"
        i = 1
        while(1):            #relax the cor constraints to:(h_NFR_cor[0] <= 0.3) and (h_NFR_cor[1] <= 0.3)
            if (host_id_cur - i) in host_id_list:
                sw_id = route._get_switch(host_id_cur - i,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_temp, host_lib[host_id_cur - i]) == 1:
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_id_cur - i].CPU - host_lib[host_id_cur - i].cpu_res
                    h_mtr[1] = host_lib[host_id_cur - i].Mem - host_lib[host_id_cur - i].mem_res
                    h_mtr[2] = host_lib[host_id_cur - i].Band - host_lib[host_id_cur - i].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in third while(1) 1", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold_shortest):     #+h_NFR_cor[1])/2
                        print "host_id",host_id_cur - i
                        return host_id_cur - i
            if (host_id_cur + i) in host_id_list:
                sw_id = route._get_switch(host_id_cur + i,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_temp, host_lib[host_id_cur + i]) == 1:
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_id_cur + i].CPU - host_lib[host_id_cur + i].cpu_res
                    h_mtr[1] = host_lib[host_id_cur + i].Mem - host_lib[host_id_cur + i].mem_res
                    h_mtr[2] = host_lib[host_id_cur + i].Band - host_lib[host_id_cur + i].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in third while(1) 2", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold_shortest):       #+h_NFR_cor[1])/2
                        print "host_id",host_id_cur + i
                        return host_id_cur + i
            i += 1
            print "i in the third while(1)", i
            if i > bound:           # can not find the host satisfies all the cor constraints :(h_NFR_cor[0] <= 0.3) and (h_NFR_cor[1] <= 0.3),
                return None        # so, needs to start a new host
    
    def Find_vnf_id_in_chain(self, vnf_chain):
        vnf_id_list = []
        for vnf in vnf_chain:
            vnf_id_list.append(vnf.id)
        return vnf_id_list
        
    
    def Find_host_cor(self, NFR_temp, host_list, host_lib):
        j = 0
        for j in range(len(host_list)):
            if NFR_temp[0] in host_lib[host_list[j].id].req_container:   # VNF chain i_id is in the current host
                vnf_id_in_chain = self.Find_vnf_id_in_chain(host_lib[host_list[j].id].req_container_dict[NFR_temp[0]])
                if (NFR_temp[1].pre_vnf_id in vnf_id_in_chain) or (NFR_temp[1].next_vnf_id in vnf_id_in_chain):
                    return host_list[j].id                     #if vnf's pre_vnf or next_vnf in the host, return the host_id
        return None                #if there is no pre_vnf or next_vnf in the host, return None.
    
    def First_fit(self,route, core_link, e_link, req_abandoned, NFR_temp, host_id, host_lib, host_list):
        print "First Fit in PowerNets"
        inditor = 0
        j = 0
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        for j in range(len(host_list)):
            #print "request seq_num:", i_id, '\t', "host_id now", host_list[j].id
            sw_id = route._get_switch(host_list[j].id,e_link)
            path = route._get_path(ROOT_ID, sw_id)
            if self.check_end(path, core_link, NFR_temp, host_lib[host_list[j].id]) == 1: # the vnf can be embedded in the current host
                if NFR_temp[0] in host_lib[host_list[j].id].req_container:   # VNF chain i_id is in the current host
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_list[j].id].CPU - host_lib[host_list[j].id].cpu_res
                    h_mtr[1] = host_lib[host_list[j].id].Mem - host_lib[host_list[j].id].mem_res
                    h_mtr[2] = host_lib[host_list[j].id].Band - host_lib[host_list[j].id].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in FF:1", h_NFR_cor[0] 
                    if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
                        #print "VNF_chain", NFR_lib[i_id][0], "VNF_id", NFR_lib[i_id][1].id, "satisfies the cor constraints."
                        host_lib[host_list[j].id].req_container_dict[NFR_temp[0]].append(NFR_temp[1])
                        self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_list[j].id])
                        self.verify_band_left_1(path,core_link,NFR_temp,host_lib[host_list[j].id])
                        inditor = 1                  #an inditor, indite if there is a need to start a new host
                        break                       #find the first host to hold the request, and stop the process
                else:        #VNF chain i_id is not in the current host       
                    h_mtr = np.zeros([3,len(x)])
                    h_mtr[0] = host_lib[host_list[j].id].CPU - host_lib[host_list[j].id].cpu_res
                    h_mtr[1] = host_lib[host_list[j].id].Mem - host_lib[host_list[j].id].mem_res
                    h_mtr[2] = host_lib[host_list[j].id].Band - host_lib[host_list[j].id].band_res
                    NFR_mtr = np.zeros([3,len(x)])
                    NFR_mtr[0] = NFR_temp[1].cpu
                    NFR_mtr[1] = NFR_temp[1].mem
                    NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                    h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id,"h_NFR_cor[0] in FF:2", h_NFR_cor[0] 
                    if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold):
                        #print "VNF_chain", NFR_lib[i_id][0], "VNF_id", NFR_lib[i_id][1].id, "satisfies the cor constraints."
                        host_lib[host_list[j].id].req_container.append(NFR_temp[0])
                        host_lib[host_list[j].id].req_container_dict[NFR_temp[0]] = [NFR_temp[1]]
                        self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_list[j].id])
                        self.verify_band_left_0(path,core_link,NFR_temp,host_lib[host_list[j].id])
                        inditor = 1                  #an inditor, indite if there is a need to start a new host
                        break                       #find the first host to hold the request, and stop the process
        #again, First-Fit needs to relax the cor_threshold too.
        if inditor == 0:
            print "Start to relax the Cor_threshold in First-Fit"
            j = 0
            sw_id = route._get_switch(host_id,e_link)
            path = route._get_path(ROOT_ID, sw_id)
            for j in range(len(host_list)):
                #print "request seq_num:", i_id, '\t', "host_id now", host_list[j].id
                sw_id = route._get_switch(host_list[j].id,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_temp, host_lib[host_list[j].id]) == 1: # the vnf can be embedded in the current host
                    if NFR_temp[0] in host_lib[host_list[j].id].req_container:   # VNF chain i_id is in the current host
                        h_mtr = np.zeros([3,len(x)])
                        h_mtr[0] = host_lib[host_list[j].id].CPU - host_lib[host_list[j].id].cpu_res
                        h_mtr[1] = host_lib[host_list[j].id].Mem - host_lib[host_list[j].id].mem_res
                        h_mtr[2] = host_lib[host_list[j].id].Band - host_lib[host_list[j].id].band_res
                        NFR_mtr = np.zeros([3,len(x)])
                        NFR_mtr[0] = NFR_temp[1].cpu
                        NFR_mtr[1] = NFR_temp[1].mem
                        NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                        h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                        print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in FF:3", h_NFR_cor[0] 
                        if (h_NFR_cor[0] <= Cor_threshold_shortest):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
                            #print "VNF_chain", NFR_lib[i_id][0], "VNF_id", NFR_lib[i_id][1].id, "satisfies the cor constraints."
                            host_lib[host_list[j].id].req_container_dict[NFR_temp[0]].append(NFR_temp[1])
                            self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_list[j].id])
                            self.verify_band_left_1(path,core_link,NFR_temp,host_lib[host_list[j].id])
                            inditor = 1                  #an inditor, indite if there is a need to start a new host
                            break                       #find the first host to hold the request, and stop the process
                    else:        #VNF chain i_id is not in the current host       
                        h_mtr = np.zeros([3,len(x)])
                        h_mtr[0] = host_lib[host_list[j].id].CPU - host_lib[host_list[j].id].cpu_res
                        h_mtr[1] = host_lib[host_list[j].id].Mem - host_lib[host_list[j].id].mem_res
                        h_mtr[2] = host_lib[host_list[j].id].Band - host_lib[host_list[j].id].band_res
                        NFR_mtr = np.zeros([3,len(x)])
                        NFR_mtr[0] = NFR_temp[1].cpu
                        NFR_mtr[1] = NFR_temp[1].mem
                        NFR_mtr[2] = NFR_temp[1].in_band + NFR_temp[1].out_band
                        h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                        print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id,"h_NFR_cor[0] in FF:4", h_NFR_cor[0] 
                        if (h_NFR_cor[0] <= Cor_threshold_shortest):# and (h_NFR_cor[1] <= Cor_threshold):
                            #print "VNF_chain", NFR_lib[i_id][0], "VNF_id", NFR_lib[i_id][1].id, "satisfies the cor constraints."
                            host_lib[host_list[j].id].req_container.append(NFR_temp[0])
                            host_lib[host_list[j].id].req_container_dict[NFR_temp[0]] = [NFR_temp[1]]
                            self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_list[j].id])
                            self.verify_band_left_0(path,core_link,NFR_temp,host_lib[host_list[j].id])
                            inditor = 1                  #an inditor, indite if there is a need to start a new host
                            break                       #find the first host to hold the request, and stop the process
            
        
        if inditor == 0:
            print "start a new host"
            host_id += 1
            host_list.append(host_lib[host_id])
            #print "request seq_num:", i_id, '\t', "host_id now", host_id 
            sw_id = route._get_switch(host_id,e_link)
            path = route._get_path(ROOT_ID, sw_id)
            if self.check_end(path, core_link, NFR_temp, host_lib[host_id]) == 1:
                if NFR_temp[0] in host_lib[host_id].req_container:
                    host_lib[host_id].req_container_dict[NFR_temp[0]].append(NFR_temp[1])
                    self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_id])
                    self.verify_band_left_1(path,core_link,NFR_temp,host_lib[host_id])
                else:
                    host_lib[host_id].req_container.append(NFR_temp[0])
                    host_lib[host_id].req_container_dict[NFR_temp[0]] = [NFR_temp[1]]
                    self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_id])
                    self.verify_band_left_0(path,core_link,NFR_temp,host_lib[host_id])
            else:
                req_abandoned.append(NFR_temp[0])
        return host_id

    
    def PowerNets(self, route, core_link, e_link, key, vnf_chain_list, host_lib, host_id):
        print "PowerNets"
        print "sort the requests based on the 90-percentile resource"
        NFR_num = 0
        NFR_lib = {}
        for NFR_i in vnf_chain_list:
            NFR_tuple = (key,NFR_i)
            NFR_lib[NFR_num] = NFR_tuple
            NFR_num += 1
        VNF_move_id_list = self.sort_requests(NFR_lib)
        req_abandoned = []
        host_list = []
        for l in range(host_id+1):
            host_list.append(host_lib[l])
        NFR_num = len(NFR_lib)
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        for i in range(len(VNF_move_id_list)):
            i_id = VNF_move_id_list[i]
            if host_id == 0 and len(host_lib[host_id].req_container) == 0:   #the first vnf , put it in the current host directly
                #print "The first VNF"
                print "req_num",key,"vnf_id",NFR_lib[i_id][1].id
                host_lib[host_id].req_container.append(NFR_lib[i_id][0])
                host_lib[host_id].req_container_dict[NFR_lib[i_id][0]] = [NFR_lib[i_id][1]]
                self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id])
                self.verify_band_left_0(path,core_link,NFR_lib[i_id],host_lib[host_id])
            else:
                host_cor = self.Find_host_cor(NFR_lib[i_id],host_list,host_lib)
                if host_cor == None:
                    host_id = self.First_fit(route, core_link,e_link,req_abandoned, NFR_lib[i_id], host_id, host_lib, host_list)
                else:
                    sw_id = route._get_switch(host_cor,e_link)
                    path = route._get_path(ROOT_ID, sw_id)
                    if self.check_end(path, core_link, NFR_lib[i_id], host_lib[host_cor]) == 1: # the vnf can be embedded in the current host
                        h_mtr = np.zeros([3,len(x)])
                        h_mtr[0] = host_lib[host_cor].CPU - host_lib[host_cor].cpu_res
                        h_mtr[1] = host_lib[host_cor].Mem - host_lib[host_cor].mem_res
                        h_mtr[2] = host_lib[host_cor].Band - host_lib[host_cor].band_res
                        NFR_mtr = np.zeros([3,len(x)])
                        NFR_mtr[0] = NFR_lib[i_id][1].cpu
                        NFR_mtr[1] = NFR_lib[i_id][1].mem
                        NFR_mtr[2] = NFR_lib[i_id][1].in_band + NFR_lib[i_id][1].out_band
                        h_NFR_cor = self.calc_d_val_pearson(h_mtr, NFR_mtr)
                        print "VNF_CHAIN:",NFR_lib[i_id][0],"VNF:",NFR_lib[i_id][1].id, "find cor host", "Corefficient:", h_NFR_cor[0]
                        if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold):       #vnf satisfies the correlation threshold #+h_NFR_cor[1])/2
                            #print "the vnf",NFR_lib[i_id][1].id, "satisfies the cor constraints"
                            host_lib[host_cor].req_container_dict[NFR_lib[i_id][0]].append(NFR_lib[i_id][1])
                            self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_cor])
                            self.verify_band_left_1(path,core_link,NFR_lib[i_id],host_lib[host_cor])
                        else:      #vnf can not satisfy the correlation threshold, find the nearest proper host
                            #print "the vnf",NFR_lib[i_id][1].id, "can not satisfy the cor constraints"
                            host_list_copy = copy.deepcopy(host_list)
                            #host_list_copy.remove(host_list[j].id)
                            host_id_t = self.Find_shortest_host_for_vnf(route, core_link, e_link, host_cor, host_list_copy, NFR_lib[i_id], host_lib)
                            if host_id_t == None:
                                print "can not find proper near host for the vnf(start new host): CHAIN_ID,VNF_ID",NFR_lib[i_id][0], NFR_lib[i_id][1].id
                                #start new host
                                host_id += 1
                                host_list.append(host_lib[host_id])
                                print "request seq_num:", NFR_lib[i_id][0], '\t', "host_id now", host_id 
                                sw_id = route._get_switch(host_id,e_link)
                                path = route._get_path(ROOT_ID, sw_id)
                                if self.check_end(path, core_link, NFR_lib[i_id], host_lib[host_id]) == 1:
                                    if NFR_lib[i_id][0] in host_lib[host_id].req_container:
                                        host_lib[host_id].req_container_dict[NFR_lib[i_id][0]].append(NFR_lib[i_id][1])
                                        self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id])
                                        self.verify_band_left_1(path,core_link,NFR_lib[i_id],host_lib[host_id])
                                    else:
                                        host_lib[host_id].req_container.append(NFR_lib[i_id][0])
                                        host_lib[host_id].req_container_dict[NFR_lib[i_id][0]] = [NFR_lib[i_id][1]]
                                        self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id])
                                        self.verify_band_left_0(path,core_link,NFR_lib[i_id],host_lib[host_id])
                                else:
                                    req_abandoned.append(i_id)
                                #break
                            else:
                                #print "Find the proper near host:",host_id_t, "for the vnf",NFR_lib[i_id][1].id
                                sw_id = route._get_switch(host_id_t,e_link)
                                path = route._get_path(ROOT_ID, sw_id)
                                if NFR_lib[i_id][0] in host_lib[host_id_t].req_container:
                                    host_lib[host_id_t].req_container_dict[NFR_lib[i_id][0]].append(NFR_lib[i_id][1])
                                    self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id_t])
                                    self.verify_band_left_1(path,core_link,NFR_lib[i_id],host_lib[host_id_t])
                                else:
                                    host_lib[host_id_t].req_container.append(NFR_lib[i_id][0])
                                    host_lib[host_id_t].req_container_dict[NFR_lib[i_id][0]] = [NFR_lib[i_id][1]]
                                    self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id_t])
                                    self.verify_band_left_0(path,core_link,NFR_lib[i_id],host_lib[host_id_t])
                    else:
                        #print "the vnf",NFR_lib[i_id][1].id, "can not satisfy the cor constraints"
                        host_list_copy = copy.deepcopy(host_list)
                        #host_list_copy.remove(host_list[j].id)
                        host_id_t = self.Find_shortest_host_for_vnf(route, core_link, e_link, host_cor, host_list_copy, NFR_lib[i_id], host_lib)
                        if host_id_t == None:
                            print "can not find proper near host for the vnf(start new host): CHAIN_ID,VNF_ID",NFR_lib[i_id][0], NFR_lib[i_id][1].id
                            #start new host
                            host_id += 1
                            host_list.append(host_lib[host_id])
                            print "request seq_num:", NFR_lib[i_id][0], '\t', "host_id now", host_id 
                            sw_id = route._get_switch(host_id,e_link)
                            path = route._get_path(ROOT_ID, sw_id)
                            if self.check_end(path, core_link, NFR_lib[i_id], host_lib[host_id]) == 1:
                                if NFR_lib[i_id][0] in host_lib[host_id].req_container:
                                    host_lib[host_id].req_container_dict[NFR_lib[i_id][0]].append(NFR_lib[i_id][1])
                                    self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id])
                                    self.verify_band_left_1(path,core_link,NFR_lib[i_id],host_lib[host_id])
                                else:
                                    host_lib[host_id].req_container.append(NFR_lib[i_id][0])
                                    host_lib[host_id].req_container_dict[NFR_lib[i_id][0]] = [NFR_lib[i_id][1]]
                                    self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id])
                                    self.verify_band_left_0(path,core_link,NFR_lib[i_id],host_lib[host_id])
                            else:
                                req_abandoned.append(i_id)
                            #break
                        else:
                            #print "Find the proper near host:",host_id_t, "for the vnf",NFR_lib[i_id][1].id
                            sw_id = route._get_switch(host_id_t,e_link)
                            path = route._get_path(ROOT_ID, sw_id)
                            if NFR_lib[i_id][0] in host_lib[host_id_t].req_container:
                                    host_lib[host_id_t].req_container_dict[NFR_lib[i_id][0]].append(NFR_lib[i_id][1])
                                    self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id_t])
                                    self.verify_band_left_1(path,core_link,NFR_lib[i_id],host_lib[host_id_t])
                            else:
                                host_lib[host_id_t].req_container.append(NFR_lib[i_id][0])
                                host_lib[host_id_t].req_container_dict[NFR_lib[i_id][0]] = [NFR_lib[i_id][1]]
                                self.verify_node_resource_left(path, core_link, NFR_lib[i_id],host_lib[host_id_t])
                                self.verify_band_left_0(path,core_link,NFR_lib[i_id],host_lib[host_id_t])
                
        print "req_abandoned is:",req_abandoned,len(req_abandoned)
        #print "host_list",host_list
        #host_obsoleted = []
        #host_id_list = []
        #for l in range(len(host_list)):
            #host_id_list.append(host_list[l].id)
        #print "host_list length is:",len(host_list)
        #print "host_id in the host_list is:", host_id_list
        #for k in range(host_id+1):
            #if len(host_lib[k].req_container) == 0:
                #host_obsoleted.append(k)
        #host_num = host_id+1 - len(host_obsoleted)
        #print "host_num", host_num
        print "final host id (PN):", host_id
        return host_id
    
  
class VNF_xFit(object):
    
    def verify_resource_left(self, path, link_map, request, host):
        #sw_id = route._get_switch(host.id)
        #path = route._get_path(ROOT_ID, sw_id)
        #modify host resource left
        req_len = len(request)               #here, request is a list
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for i in range(req_len):
            #if request[i].id in host.vnf_container == False:
                #host.vnf_container.append(request[i].id)
            req_cpu += request[i].cpu
            req_mem += request[i].mem
        req_band = request[0].in_band                            #request band is determined by the first VNF
        
        #modify the vnfs in each host
        vnf_in_req = []
        vnf_in_req = self.calc_req_vnf(request)
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
        host.verify_band(req_band)
        #modify link
        for a,b in zip(path[:-1],path[1:]):
            link_map[a][b] = link_map[a][b] - req_band
            link_map[b][a] = link_map[b][a] - req_band
        
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
        for i in range(req_len):
            req_cpu += request[i].cpu
            req_mem += request[i].mem
        req_band = request[0].in_band                            #request band is determined by the first VNF
        vnf_in_req = []
        vnf_in_req = self.calc_req_vnf(request)
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
        
    
    def NF_map(self, route, core_link, e_link, request_lib, host_lib, req_id_list):     #map the request to host using Next Fit
        print "Next Fit"
        host_id = 0
        req_abandoned = []
        request_num = len(request_lib)
        PowerNets_C = VNF_PowerNets()
        #sw_id = route._get_switch(host_id,e_link)
        #path = route._get_path(ROOT_ID, sw_id)
        for i in req_id_list:
            sw_id = route._get_switch(host_id,e_link)
            path = route._get_path(ROOT_ID, sw_id)
            if self.check_end(path, core_link, request_lib[i], host_lib[host_id]):
                host_lib[host_id].req_container.append(i)
                host_lib[host_id].req_container_dict[i] = request_lib[i]
                self.verify_resource_left(path, core_link, request_lib[i],host_lib[host_id])
            else:
                host_id += 1
                sw_id = route._get_switch(host_id,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, request_lib[i], host_lib[host_id]):
                    host_lib[host_id].req_container.append(i)
                    host_lib[host_id].req_container_dict[i] = request_lib[i]
                    self.verify_resource_left(path, core_link, request_lib[i],host_lib[host_id])
                else:
                    req_abandoned.append(i)
                    host_id = PowerNets_C.PowerNets(route, core_link, e_link, i, request_lib[i], host_lib, host_id)
        print "final host id (NF):", host_id
        print "abandoned requeats are:",req_abandoned
        return host_id
    
    def FF_map(self, route, core_link, e_link, request_lib, host_lib, req_id_list):      #map the request to host using First Fit
        print "First Fit"
        req_abandoned = []
        host_id = 0
        host_list = []
        PowerNets_C = VNF_PowerNets()
        host_list.append(host_lib[host_id])
        request_num = len(request_lib)
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        #request_lib_backup = copy.deepcopy(request_lib)
        for i in req_id_list:
            j = 0
            inditor = 0
            for j in range(len(host_list)):
                #print "request seq_num:", i, '\t', "host_id now", host_list[j].id
                if self.check_end(path, core_link, request_lib[i], host_list[j]) == 1:
                    #print "map the request to the host",'\t',i, host_list[j].id
                    host_lib[host_list[j].id].req_container.append(i)
                    host_lib[host_list[j].id].req_container_dict[i] = request_lib[i]
                    self.verify_resource_left(path, core_link, request_lib[i],host_list[j])
                    inditor = 1                  #an inditor, indite if there is a need to start a new host
                    break                       #find the first host to hold the request, and stop the process
            if inditor == 0:
                host_id += 1
                host_list.append(host_lib[host_id])
                #print "request seq_num:", i, '\t', "host_id now", host_id 
                sw_id = route._get_switch(host_id,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, request_lib[i], host_lib[host_id]) == 1:
                    #print "map the request to the host",'\t',i, host_id
                    host_lib[host_id].req_container.append(i)
                    host_lib[host_id].req_container_dict[i] = request_lib[i]
                    self.verify_resource_left(path, core_link, request_lib[i],host_lib[host_id])
                else:
                    req_abandoned.append(i)
                    host_id = PowerNets_C.PowerNets(route, core_link, e_link, i, request_lib[i], host_lib, host_id)
        print "req_abandoned is:",req_abandoned,'\n',len(req_abandoned)
        host_obsoleted = []
        for k in range(host_id+1):
            if len(host_lib[k].req_container) == 0:
                host_obsoleted.append(k)
        host_id = host_id - len(host_obsoleted)
        print "final host id (FF):", host_id
        return host_id
                
            
        
    def get_aver(self, value_list):
        list_len = len(value_list)
        v_total = 0
        for i in range(list_len):
            v_total += value_list[i]
        aver = v_total/list_len
        return aver
    
    def calc_host_resource(self, host_candinate_list):
        host_resource = {}
        for i in range(len(host_candinate_list)):
            h_value = host_candinate_list[i]
            host_resource[h_value.id] = self.get_aver(h_value.cpu_res)/CPU_TOTAL + self.get_aver(h_value.mem_res)/MEM_TOTAL + self.get_aver(h_value.band_res)/h_band
        #min_h_id = self.search_min_dict(host_resource)
        return host_resource
    
    def search_min_dict(self, dict_exp):
        min_val = INFINITY
        for key, value in dict_exp.iteritems():
            if value <= min_val:
                min_val = value
                min_id = key
        return min_id
    
    def BF_map(self, route, core_link, e_link, request_lib, host_lib, req_id_list):      #map the request to host using Best Fit
        print "Best Fit"
        
        host_id = 0
        host_list = []
        req_abandoned = []
        PowerNets_C = VNF_PowerNets()
        host_list.append(host_lib[host_id])
        request_num = len(request_lib)
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        #request_lib_backup = copy.deepcopy(request_lib)
        for i in req_id_list:
            j = 0
            host_candinate_list = []
            
            for j in range(len(host_list)):
                if self.check_end(path, core_link, request_lib[i], host_list[j]):
                    host_candinate_list.append(host_list[j])
            if len(host_candinate_list):
                host_resource = self.calc_host_resource(host_candinate_list)
                host_chosen_id = self.search_min_dict(host_resource)
                host_lib[host_chosen_id].req_container.append(i)
                host_lib[host_chosen_id].req_container_dict[i] = request_lib[i]
                self.verify_resource_left(path, core_link, request_lib[i],host_lib[host_chosen_id])
            
            if len(host_candinate_list) == 0:
                host_id += 1
                host_list.append(host_lib[host_id])
                sw_id = route._get_switch(host_id,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, request_lib[i], host_lib[host_id]) == 1:
                    host_lib[host_id].req_container.append(i)
                    host_lib[host_id].req_container_dict[i] = request_lib[i]
                    self.verify_resource_left(path, core_link, request_lib[i],host_lib[host_id])
                else:
                    req_abandoned.append(i)
                    host_id = PowerNets_C.PowerNets(route, core_link, e_link, i, request_lib[i], host_lib, host_id)
        print "req_abandoned is:",req_abandoned,'\n',len(req_abandoned)
        host_obsoleted = []
        for k in range(host_id+1):
            if len(host_lib[k].req_container) == 0:
                host_obsoleted.append(k)
        host_id = host_id - len(host_obsoleted)
        print "final host id (BF):", host_id
        return host_id
    
    def search_max_dict(self, dict_exp):
        max_val = 0
        for key, value in dict_exp.iteritems():
            if value >= max_val:
                max_val = value
                max_id = key
        return max_id
    
    def sort_requests(self, request_lib):
        req_aver_source_dict = {}
        for i in range(len(request_lib)):
            VNF_i_cpu = 0
            VNF_i_mem = 0
            VNF_i_band = 0
            for VNF_i in request_lib[i]:
                VNF_i_cpu += self.get_aver(VNF_i.cpu)/CPU_TOTAL
                VNF_i_mem += self.get_aver(VNF_i.mem)/MEM_TOTAL
                VNF_i_band += self.get_aver(VNF_i.in_band)/h_band
            req_aver_source_dict[i] = VNF_i_cpu + VNF_i_mem + VNF_i_band
        req_id_list = []
        while(len(req_aver_source_dict)):
            max_id = self.search_max_dict(req_aver_source_dict)
            req_id_list.append(max_id)
            req_aver_source_dict.pop(max_id)
        
        return req_id_list
    
    def FFD_map(self, route, core_link, e_link, request_lib, host_lib):
        print "FFD"
        print "sort the requests based on the aver resource"
        req_id_list = self.sort_requests(request_lib)
        host_id = self.FF_map(route, core_link, e_link, request_lib, host_lib, req_id_list)
        return host_id
    
    def BFD_map(self, route, core_link, e_link, request_lib, host_lib):
        print "BFD"
        print "sort the requests based on the aver resource"
        req_id_list = self.sort_requests(request_lib)
        host_id = self.BF_map(route, core_link, e_link, request_lib, host_lib, req_id_list)
        return host_id
 

class VNF_adjust_pod(object):
    
    def count_species(self, host_id, host_lib):           #count the vnf species in each host
        VNF_dict = {}
        for i in range(host_id+1):
            VNF_list = []
            for req_id,req_value in host_lib[i].req_container_dict.iteritems():
                for VNF in req_value:
                    if VNF.id in VNF_list:
                        pass
                    else:
                        VNF_list.append(VNF.id)
                
            VNF_dict[i] = VNF_list
            host_lib[i].vnf_container = copy.deepcopy(VNF_list)          #each host has a vnf list to hold all the different vnf
        return VNF_dict
    
    def re_modify_host_res(self, host_id, host_lib):
        for i in range(host_id+1):
            vnf_num = len(host_lib[i].vnf_container)
            if vnf_num < VNF_SPE:
                host_lib[i].cpu_res = host_lib[i].cpu_res + (VNF_SPE - vnf_num)*cpu_brc_per      #if the vnf species < total num, then recover the responding resource
                host_lib[i].mem_res = host_lib[i].mem_res + (VNF_SPE - vnf_num)*mem_brc_per
    
    def get_aver(self, value_list):     #get the average of a list
        list_len = len(value_list)
        v_total = 0
        for i in range(list_len):
            v_total += value_list[i]
        aver = v_total/list_len
        return aver
        
    
    def count_resource_vnf(self, host_id, VNF_dict, host_lib): #count the resource of vnf,(aver traffic, traffic, mem, cpu, the out side band)
        
        VNF_t_dict = {}
        #print "VNF_t_dict in count_resource_vnf:"
        for i in range(host_id+1):
            #VNF_t_list = []
            VNF_t_dict2 = {}
            for VNF_id in VNF_dict[i]:
                #VNF_band_t = np.array([0 for j in range(len(x))])
                #VNF_band_0 = np.array([0 for j in range(len(x))])
                #VNF_mem = np.array([0 for j in range(len(x))])
                #VNF_cpu = np.array([0 for j in range(len(x))])
                
                VNF_band_t = 0
                VNF_band_0 = 0
                VNF_mem = 0
                VNF_cpu = 0
                for req_id,req_value in host_lib[i].req_container_dict.iteritems():
                    #req_id = host_lib[i].req_container[m]
                    n = 0
                    for n in range(len(req_value)):
                        VNF = req_value[n]
                        if VNF.id == VNF_id:
                            VNF_band_t += VNF.out_band                     #to calc cpu & mem
                            VNF_band_t += VNF.in_band
                            VNF_mem += VNF.mem
                            VNF_cpu += VNF.cpu
                            if n == 0:
                                VNF_band_0 += VNF.in_band
                VNF_t_dict2[VNF_id] = [self.get_aver(VNF_band_t),VNF_band_t,VNF_mem,VNF_cpu, VNF_band_0]         #use the average indict the traffic volume
                #VNF_t_list.append(VNF_t_dict2)
            #VNF_t_dict[i] = VNF_t_list
            VNF_t_dict[i] = VNF_t_dict2
            #print 'host_id:',i,'\t','VNF species:',VNF_t_dict[i].keys(), "VNF Number:", len(VNF_t_dict[i].keys())
        return VNF_t_dict
    
    
    def count_traffic_host(self, host):
        band_comsumed = host.Band - host.band_res
        return band_consumed
    
    def search_min_dict(self, dict_exp):    #find the min value in dict, and return the key
        min_val = INFINITY
        for key, value in dict_exp.iteritems():
            if isinstance(value,list):           #if value is a list
                if value[0] <= min_val:
                    min_val = value[0]
                    min_id = key
            else:                                # not list
                if value <= min_val:
                    min_val = value
                    min_id = key
                
        return min_id
        
    def search_max_dict(self, dict_exp):
        max_val = 0
        for key, value in dict_exp.iteritems():
            if isinstance(value,list):
                if value[0] >= max_val:
                    max_val = value[0]
                    max_id = key
            else:
                #print "value in factor:",value
                if value >= max_val:
                    max_val = value
                    max_id = key
        return max_id
 
    def sort_vnf(self, VNF_t_dict):
        for VNF_id, VNF_v in VNF_t_dict.iteritems():
            VNF_v_new = sorted(VNF_v.iteritems(), key=operator.itemgetter(1))    
            VNF_t_dict[VNF_id] = VNF_v_new                                                    #after sorted, the dict turn to list with tuple in it    
                                                                                            #x = {1:2, 3:4, 4:3, 2:1, 0:0}  
        return  VNF_t_dict                                                                  #sorted_x = sorted(x.iteritems(), key=operator.itemgetter(1))  
                                                                                            #print sorted_x  
                                                                                            #[(0, 0), (2, 1), (1, 2), (4, 3), (3, 4)]
    
    
    def verify_req_container_src(self, host_src):
        
        def search_vnf(vnf_list, vnf_id_s):
            for vnf_s in vnf_list:
                if vnf_s.id == vnf_id_s:
                    return vnf_s
            
        
        
        
        #container_len = len(host_src.req_container)
        #print "container len is", container_len
        #req_list_in_host = host_src.req_container
        req_to_removed = []
        for req_id,req_value in host_src.req_container_dict.iteritems():
            #req_id = req_list_in_host[m]
            vnf_list_req = []
            for vnf in req_value:
                vnf_id = vnf.id
                vnf_list_req.append(vnf_id)
            #if len(vnf_list_req) > len(host_src.vnf_pop_container):
                #pass
            #else:
            inditor = 0
            for k in vnf_list_req:
                #print k
                if k in host_src.vnf_pop_container:
                    vnf_removed = search_vnf(req_value, k)
                    host_src.req_container_dict[req_id].remove(vnf_removed) #the vnf has been removed, and the resource has been verified, so momdify the request in the request_lib
                else:                                                       #it is necessary to minus the same vnf in the requests, to in consistent with "verify_hvc_Vtd"
                    inditor += 1
                    continue
            if inditor == 0:                  #all vnfs in the request are included in the vnf_pop_container
                req_to_removed.append(req_id)
                #host_src.req_container.remove(req_id)
                #host_src.req_container_dict.pop(req_id)
        for r2r_id in req_to_removed:
            host_src.req_container.remove(r2r_id)
            host_src.req_container_dict.pop(r2r_id)
        
    
    def verify_band_and_req_dict(self, transformed_request, host_dst, host_src):
        
        def vnfID_in_req(vnf_list):
            vnfID_list = []
            for vnf_i in vnf_list:
                vnfID_list.append(vnf_i.id)
            return vnfID_list
            
        
        #the host_dst has no the same req
        def vnf_first_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id): #the req has not been moved clear
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list:       #next_vnf_id in the src_host
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                else:            #next_vnf_id in the third part, because the host_dst has no the same req_id,so it has no the next_vnf_id definitely.
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
            else:
                #the next_vnf_id is not none, and the next_vnf_id is in the third part,because the host_dst has no the same req_id,so it has no the next_vnf_id definitely.
                dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
        
        def vnf_last_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list:
                    dst_host.band_res -= trans_vnf.in_band
                    src_host.band_res -= trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
            else:
                src_host.band_res += trans_vnf.in_band
                dst_host.band_res -= trans_vnf.in_band
        
        def vnf_middle_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list and trans_vnf.next_vnf_id in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list and trans_vnf.next_vnf_id not in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id not in vnf_id_list and trans_vnf.next_vnf_id in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
            else:
                dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
        
        #the host_dst has the same req
        def vnf_first_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
        
        def vnf_last_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_src:
                    src_host.band_res -= trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res += trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_dst:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res += trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
        
        def vnf_middle_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id not in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.out_band - trans_vnf.in_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.out_band + trans_vnf.in_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_src and trans_vnf.pre_vnf_id not in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.out_band - trans_vnf.in_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
        
        for trans_id, trans_value in transformed_request.iteritems():
            if host_dst.req_container_dict.has_key(trans_id):
                #self.verify_req_container_src(host_src)
                host_dst.req_container_dict[trans_id].append(trans_value)
                if trans_value.pre_vnf_id == None:
                    vnf_first_1(host_dst, host_src, trans_id, trans_value)
                elif trans_value.next_vnf_id == None:
                    vnf_last_1(host_dst, host_src, trans_id, trans_value)
                else:
                    vnf_middle_1(host_dst, host_src, trans_id, trans_value)
                #host_dst.req_container_dict[trans_id] = self.right_sequence_vnf(host_dst.req_container_dict[trans_id], requests_lib[trans_id])
            else:
                #self.verify_req_container_src(host_src)
                host_dst.req_container_dict[trans_id] = [trans_value]
                host_dst.req_container.append(trans_id)
                if trans_value.pre_vnf_id == None and trans_value.next_vnf_id == None:
                    host_dst.band_res -= trans_value.in_band
                    host_src.band_res += trans_value.in_band
                elif trans_value.pre_vnf_id == None and trans_value.next_vnf_id != None:
                    vnf_first_0(host_dst, host_src, trans_id, trans_value)
                elif trans_value.next_vnf_id == None and trans_value.pre_vnf_id != None:
                    vnf_last_0(host_dst, host_src, trans_id, trans_value)
                else:
                    vnf_middle_0(host_dst, host_src, trans_id, trans_value)
    
    
    
    def mignate_vnf(self, src, dst, VNF_t_dict, host_src, host_dst, vnf_id, vnf, requests_lib):           #mignate the vnf and modify responding resource
        """
        src is the src host id, dst is the dst host id
        VNF_t_dict is a global pramater
        host_src and host_dst both are VNF_server class
        vnf_id is the id of vnf
        vnf is a list
        """
        
        VNF_t_dict[src].pop(vnf_id)
        host_src.vnf_container.remove(vnf_id)
        #print "what's in host_src.vnf_container are:", host_src.vnf_container
        host_src.vnf_pop_container.append(vnf_id)
        
        VNF_t_dict[dst][vnf_id][0] = VNF_t_dict[dst][vnf_id][0] + vnf[0]
        VNF_t_dict[dst][vnf_id][1] = VNF_t_dict[dst][vnf_id][1] + vnf[1]
        VNF_t_dict[dst][vnf_id][2] = VNF_t_dict[dst][vnf_id][2] + vnf[2]
        VNF_t_dict[dst][vnf_id][3] = VNF_t_dict[dst][vnf_id][3] + vnf[3]
        VNF_t_dict[dst][vnf_id][4] = VNF_t_dict[dst][vnf_id][4] + vnf[4]
        
        #verify the host_dst
        #calc the related requests in responding to the transformed vnf in host_src
        transformed_req = {}
        for req_src_id,req_src_value in host_src.req_container_dict.iteritems():
            for vnf_src in req_src_value:
                if vnf_src.id == vnf_id:
                    transformed_req[req_src_id] = vnf_src
        #print "transformed_req",transformed_req.keys()
        
        #modify the requests in host_dst after the vnf mignating
        cpu_resource_removed = 0
        mem_resource_removed = 0
        band_resource_removed = 0
        for trans_id, trans_value in transformed_req.iteritems():
            cpu_resource_removed += trans_value.cpu
            mem_resource_removed += trans_value.mem
        self.verify_req_container_src(host_src)
        self.verify_band_and_req_dict(transformed_req, host_dst, host_src)
        
        
        #modify the related resources lefted in host_src
        host_src.cpu_res = host_src.cpu_res + cpu_resource_removed + cpu_brc_per
        host_src.mem_res = host_src.mem_res + mem_resource_removed + mem_brc_per
        
        
        #self.verify_req_container_src(host_src)
        host_dst.cpu_res = host_dst.cpu_res - cpu_resource_removed
        host_dst.mem_res = host_dst.mem_res - mem_resource_removed
    
    
    def rich_or_poor(self, host):                                #if true: rich host, False:poor host
        gama = 0.0                                   #the threshold that determine rich or poor
        times_cpu = 0
        times_mem = 0
        times_band = 0
        for i in range(len(x)):
            if host.cpu_res[i] < 0:
                times_cpu += 1
            if host.mem_res[i] < 0:
                times_mem += 1
            if host.band_res[i] < 0:
                times_band += 1
        if times_cpu/len(x) > gama or times_mem/len(x) > gama or times_band/len(x) > gama:
            return False
        else:
            return True
    
    def calc_req_vnf(self, req_list):                #count the vnf species in req
        vnf_list = []
        req_len = len(req_list)
        for i in range(req_len):
            vnf_list.append(req_list[i].id)
        
        return vnf_list
    
    def sl_in_ll(self, s_list, l_list):               #judge if the req_vnf in the host_vnf_container
        sl_len = len(s_list)
        for i in range(sl_len):
            if (s_list[i] in l_list) == False:
                return False
        return True
        
    def check_value(self, list_ori, list_cmp): # calc the break times of resources constraints
        list_len = len(list_ori)
        break_t = 0
        for i in range(list_len):
            if list_cmp[i] >= list_ori[i]:
                break_t += 1
        return break_t
        
    
    
    def req_to_host(self, request, host):
        """
        judge if the request can be setted to host
        """
        gama = 0.0
        req_len = len(request)
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for i in range(req_len):
            req_cpu += request[i].cpu
            req_mem += request[i].mem
            req_band += request[i].in_band +  request[i].out_band                          #ensure the received host can host the request
        
        h_cpu_bt_total = self.check_value(host.cpu_res, req_cpu)/len(host.cpu_res)                       #break times
        h_mem_bt_total = self.check_value(host.mem_res, req_mem)/len(host.mem_res)
        h_band_bt_total = self.check_value(host.band_res, req_band)/len(host.band_res)
        
        if  h_cpu_bt_total > gama or h_mem_bt_total > gama or h_band_bt_total > gama:                     #a threshold be defined
            return False
        else:
            return True
       
    def verify_hvc_Vtd(self, host_start, host_end, host_lib, VNF_t_dict):           #refresh the host_vnf_container & VNF_t_dict
        for i in range(host_start, host_end):
            VNF_list = []
            for req_id, req_value in host_lib[i].req_container_dict.iteritems():
                for VNF in req_value:
                    if VNF.id in VNF_list:
                        pass
                    else:
                        VNF_list.append(VNF.id)
            host_lib[i].vnf_container = copy.deepcopy(VNF_list)
            j = 0
            for j in range(len(host_lib[i].vnf_pop_container)):
                if host_lib[i].vnf_pop_container[j] in host_lib[i].vnf_container:
                    host_lib[i].vnf_container.remove(host_lib[i].vnf_pop_container[j])
            #print "vnfs of host i", host_lib[i].vnf_container, i
        i = host_start
        for i in range(host_start, host_end):
            #VNF_t_list = []
            if len(host_lib[i].vnf_container) == 0:
                VNF_t_dict[i] = {}
            else:
                VNF_t_dict2 = {}
                for VNF_id in host_lib[i].vnf_container:
                    m = 0
                    n = 0
                    VNF_band_t = 0
                    VNF_band_0 = 0
                    VNF_mem = 0
                    VNF_cpu = 0
                    
                    #VNF_band_t = np.array([0 for j in range(len(x))])
                    #VNF_band_0 = np.array([0 for j in range(len(x))])
                    #VNF_mem = np.array([0 for j in range(len(x))])
                    #VNF_cpu = np.array([0 for j in range(len(x))])
                    for req_id_2,req_value_2 in host_lib[i].req_container_dict.iteritems():
                        n = 0
                        for n in range(len(req_value_2)):
                            VNF = req_value_2[n]
                            if VNF.id == VNF_id:
                                VNF_band_t += VNF.out_band                     #to calc cpu & mem
                                VNF_band_t += VNF.in_band
                                VNF_mem += VNF.mem
                                VNF_cpu += VNF.cpu
                                if n == 0:
                                    VNF_band_0 += VNF.in_band
                    VNF_t_dict2[VNF_id] = [self.get_aver(VNF_band_t),VNF_band_t,VNF_mem,VNF_cpu, VNF_band_0]         #use the average indict the traffic volume
                    
                VNF_t_dict[i] = VNF_t_dict2
                #print "keys of VNF_t_dict in verify_hvc_Vtd", i, VNF_t_dict2.keys()
    
    def calc_d_val_pearson(self, mtr1,mtr2):              #calc pearson value
        #print "*****calc d value*****"
        [row,col] = mtr1.shape
        #print "row,col",row,col
        #mtr_temp = np.zeros([row,col])
        #normalization
        #normalization
        """
        normal_max_dict1 = {}
        normal_max_dict2 = {}
        for k in range(row):
            normal_max_dict1[k] = 0
            normal_max_dict2[k] = 0
            for l in range(col):
                if mtr1[k][l] > normal_max_dict1[k]:
                    normal_max_dict1[k] = mtr1[k][l]
                if mtr2[k][l] > normal_max_dict2[k]:
                    normal_max_dict2[k] = mtr2[k][l]
        
        for m in range(row):
            for n in range(col):
                mtr1[m][n] = mtr1[m][n]/normal_max_dict1[m]
                mtr2[m][n] = mtr2[m][n]/normal_max_dict2[m]
            
        """
        #pearson correlation coefficient
        pearson = {}
        for i in range(row):
            pearson[i] = 0
            factor_1 = 0
            factor_2 = 0
            factor_3 = 0
            factor_4 = 0
            factor_5 = 0
            j = 0
            for j in range(col):
                factor_1 += mtr1[i][j]*mtr2[i][j]
                factor_2 += mtr1[i][j]
                factor_3 += mtr2[i][j]
                factor_4 += mtr1[i][j]*mtr1[i][j]
                factor_5 += mtr2[i][j]*mtr2[i][j]
            #print "sqrt_1:",col*factor_4-math.pow(factor_2,2)
            Cross_x = col*factor_4-math.pow(factor_2,2)
            Cross_y = col*factor_5-math.pow(factor_3,2)
            #print "i, Row,Col", i, row, col
            #print "col*factor_4-math.pow(factor_2,2)", col*factor_4-math.pow(factor_2,2)
            #print "col*factor_5-math.pow(factor_3,2)", col*factor_5-math.pow(factor_3,2)
            if Cross_x < 0 and Cross_x > -1*e-8:
                Cross_x = 0
                #print "vector len", len(mtr2[i]), "Mistake vector:", mtr2[i]
                #print "col*factor_5",col*factor_5, "math.pow(factor_3,2)", math.pow(factor_3,2)
            if Cross_y < 0 and Cross_y > -1*e-8:
                Cross_y = 0
            if Cross_x*Cross_y == 0:
                pearson[i] = 1
            else:
                pearson[i] = (col * factor_1 - factor_2 * factor_3)/(math.sqrt(Cross_x*Cross_y))
        d = 0
        for i in range(row):
            d += pearson[i]/3
        
        return d
            
    
    def calc_d_val(self, mtr1,mtr2):              #calc d value
        [row,col] = mtr1.shape
        #mtr_temp = np.zeros([row,col])
        #normalization
        #normalization
        normal_max_dict1 = {}
        normal_max_dict2 = {}
        for k in range(row):
            normal_max_dict1[k] = 0
            normal_max_dict2[k] = 0
            for l in range(col):
                if mtr1[k][l] > normal_max_dict1[k]:
                    normal_max_dict1[k] = mtr1[k][l]
                if mtr2[k][l] > normal_max_dict2[k]:
                    normal_max_dict2[k] = mtr2[k][l]
        
        for m in range(row):
            for n in range(col):
                mtr1[m][n] = mtr1[m][n]/(normal_max_dict1[m]+EPSILON_F)  #+EPSILON_F
                mtr2[m][n] = mtr2[m][n]/(normal_max_dict2[m]+EPSILON_F)
        
        mtr_abs = 0
        for i in range(row):
            for j in range(col):
                mtr_abs += math.fabs(mtr1[i][j]-mtr2[i][j])
        d = 1/(mtr_abs+EPSILON_F)
        return d
    
    
    def verify_brc(self, host_lib, host_start, host_end, VNF_t_dict_old, VNF_t_dict_new):
        for i in range(host_start,host_end):
            keys_len_old = len(VNF_t_dict_old[i].keys())
            keys_len_new = len(VNF_t_dict_new[i].keys())
            #print "keys_len_old",'\n',keys_len_old
            #print "keys_len_new",'\n',keys_len_new
            added_keys = []
            if keys_len_old > keys_len_new:
                host_lib[i].cpu_res += (keys_len_old-keys_len_new)*cpu_brc_per
                host_lib[i].mem_res += (keys_len_old-keys_len_new)*mem_brc_per
                """
                for key_id in VNF_t_dict_old[i].keys():
                    if key_id in VNF_t_dict_new[i].keys():
                        continue
                    else:
                        added_keys.append(key_id)
                
                for added_key_id in added_keys:
                    if added_key_id in host_lib[i].vnf_pop_container:
                        continue
                    else:
                        host_lib[i].vnf_container.remove(added_key_id)
                        host_lib[i].vnf_pop_container.append(added_key_id)
                """
            
        
    
    
    def eliminate_exceed(self, host_start,host_end, host_lib, VNF_t_dict, requests_lib):
        """
        host_rich , host that has residual resource
        host_poor , host that resource exceeded
        take some request in host_poor to host_rich to erase the exceeded resource
        """
        #print "start to check the exceeded resource in each host"
        host_rich = {}
        host_poor = {}
        for i in range(host_start,host_end):
            if self.rich_or_poor(host_lib[i]):
                host_rich[i] = host_lib[i]
            else:
                host_poor[i] = host_lib[i]
        if len(host_poor) == 0:
            print "there is no poor host"
            return True                           #if there is no host_poor, no more operations, return true
        else:
            """
            mignate the request from the last in the src host.req_container,
            to reduce the influence on the host usage in src host
            """
            #print "host_poor are:", host_poor.keys()
            for hp_id,hp_value in host_poor.iteritems():
                hp_c_backup = copy.deepcopy(hp_value.req_container_dict.keys())
                #print "poor host id:",hp_id
                while(self.rich_or_poor(hp_value) == False and len(hp_c_backup) != 0):     
                    req_id = hp_c_backup.pop()               #pick out the last req_request in host poor
                    req_vnf_list = self.calc_req_vnf(hp_value.req_container_dict[req_id])
                    req_matrix = np.zeros((3,len(x)))
                    req_len = len(hp_value.req_container_dict[req_id])
                    req_cpu = 0
                    req_mem = 0
                    i = 0
                    for i in range(req_len):
                        req_cpu += hp_value.req_container_dict[req_id][i].cpu
                        req_mem += hp_value.req_container_dict[req_id][i].mem
                    #req_band = hp_value.req_container_dict[req_id][0].in_band          #request band is determined by the first VNF
                    #req_band = np.zeros([1,len(x)])      #the modify of req_band is complicated and dubious, so just ignore it when making a choice
                    req_matrix[0] = req_cpu
                    req_matrix[1] = req_mem
                    #req_matrix[2] = req_band
                    condidate_d_dict = {}
                    for hr_id,hr_value in host_rich.iteritems():
                        if self.sl_in_ll(req_vnf_list, hr_value.vnf_container) and self.req_to_host(hp_value.req_container_dict[req_id], hr_value):
                            """
                            1.all the vnfs in request going to be removed must be included in the dst host
                            2.the resource left must could hold the request
                            """
                            host_matrix = np.zeros((3,len(x)))        #dst host matrix
                            host_matrix[0] = hr_value.CPU - hr_value.cpu_res
                            host_matrix[1] = hr_value.Mem - hr_value.mem_res
                            #host_matrix[2] = hr_value.Band - hr_value.band_res
                            
                            condidate_d_dict[hr_id] = self.calc_d_val_pearson(host_matrix, req_matrix)
                    if len(condidate_d_dict) != 0:
                        """
                        there is dst host can receive the responding request
                        Notice that:
                        here, we have not modify the band of links along the way to the dst host, so we assume that,
                        the bottleneck of resource restrict is in the link between the last sw and the dst host.
                        if the bottlenek can satisfy the contraints, then the link along the way will be satisfied too.
                        """
                        h_rece_id = self.search_min_dict(condidate_d_dict)
                        #print "received host_id is:", h_rece_id
                        #print "removed req_id",req_id,"vnf in the req", req_vnf_list
                        trans_request = []
                        trans_request = copy.deepcopy(hp_value.req_container_dict[req_id])
                        #print "trans_request", trans_request
                        #print "hp_value.req_container_dict[req_id]",hp_value.req_container_dict[req_id]
                        for trans_vnf in trans_request:
                            trans_req = {}
                            trans_req[req_id] = trans_vnf
                            
                            for vnf_index in hp_value.req_container_dict[req_id]:
                                if vnf_index.id == trans_vnf.id:
                                    vnf_removed = vnf_index
                            
                            hp_value.req_container_dict[req_id].remove(vnf_removed)
                            if len(hp_value.req_container_dict[req_id]) == 0:
                                hp_value.req_container_dict.pop(req_id)
                                hp_value.req_container.remove(req_id)
                            self.verify_band_and_req_dict(trans_req, host_lib[h_rece_id], hp_value)
                        
                        hp_value.cpu_res += req_cpu
                        hp_value.mem_res += req_mem
                        
                        host_lib[h_rece_id].cpu_res -= req_cpu
                        host_lib[h_rece_id].mem_res -= req_mem
                        #host_lib[h_rece_id].band_res -= req_band
                
                if self.rich_or_poor(hp_value) == True:
                    VNF_t_dict_old = copy.deepcopy(VNF_t_dict)
                    #print "verify the host_vnf_container & VNF_t_dict"
                    self.verify_hvc_Vtd(host_start, host_end, host_lib, VNF_t_dict)
                    #modify brc the vnf used on each host based on "VNF_t_dict_old & VNF_t_dict"
                    #print "verify the brc"
                    self.verify_brc(host_lib, host_start, host_end, VNF_t_dict_old, VNF_t_dict)
                    continue
                else:
                    return False                               #if there is one poor host can not be turned to rich, then the process of "eliminate_exceed"
            return True                                        #will not be successful
            
    
    def check_host_obsoleted(self, host_candinate, host_lib):
        pass
    
        
    def host_adjust(self, host_id, K, host_lib, requests_lib):
        alpha = 0.5                    #ratio of d
        beta = 0.5                     #ratio of aver traffic volume
        print "***count vnf species in each host***"
        VNF_t_dict_initial = self.count_species(host_id, host_lib)
        #print "what's in VNF_t_dict_initial",'\n', VNF_t_dict_initial
        print "count species over"
        #print "re_modify host resource lefted"
        #self.re_modify_host_res(host_id,host_lib)
        print "***count the traffic go through each vnf in each host***"
        VNF_t_dict = self.count_resource_vnf(host_id, VNF_t_dict_initial, host_lib)
        #print "what's in VNF_t_dict",'\n', VNF_t_dict.keys()
        print "count traffic over"
        
        host_start = 0
        host_end = host_id + 1
        host_candinate = []              #store all host that can be used to adjust the vnf
        host_i = host_start
        for host_i in range(host_start,host_end):
            host_candinate.append(host_i)
        #print "what's in host_candinate:",'\n',host_candinate
        #jump_out_while = 0
        #jump_out_while_host_id = []
        i = host_start
        host_lib_ori = copy.deepcopy(host_lib)
        while(len(host_candinate) > 1):
            jump_out_while = 0
            jump_out_while_host_id = []
            VNF_min_dict = {}
            for i in range(host_start,host_end):         #find the min vnf(based on aver traffic through it) in each host
                if i in host_candinate:
                    if len(VNF_t_dict[i]):               #i must not be empty
                        VNF_min_dict[i] = self.search_min_dict(VNF_t_dict[i])
                    else:             #all the vnfs have been moved
                        host_candinate.remove(i)
            #print "what's in VNF_min_dict:",'\n',VNF_min_dict
            #d & traffic volume
            #backup to restore
            if len(VNF_t_dict) == 0:
                break
            
            VNF_t_dict_backup = copy.deepcopy(VNF_t_dict)                        
            host_lib_backup = copy.deepcopy(host_lib)
            #print "what's in host_candinate (in while):",'\n',host_candinate
            i = host_start
            for i in range(host_start,host_end):
                #get the indictor to choose
                if i in host_candinate:    #the id must in host_candinate
                    VNF_to_match = {}        #store the same kind vnf in other's host
                    VNF_i_mtr = np.zeros([3,len(x)])
                    j = host_start
                    for j in range(host_start,host_end):
                        if j in host_candinate:    #the id must in host_candinate
                            if j == i: #establish the resource matrix
                                #continue
                                VNF_i_mtr[0] = VNF_t_dict[j][VNF_min_dict[i]][3]                 #cpu
                                VNF_i_mtr[1] = VNF_t_dict[j][VNF_min_dict[i]][2]                 #mem
                                VNF_i_mtr[2] = VNF_t_dict[j][VNF_min_dict[i]][1]                 #band
                            else:
                                if VNF_t_dict[j].has_key(VNF_min_dict[i]):                    #ensure that host_to_receive has the same vnf
                                    VNF_to_match[j] = VNF_t_dict[j][VNF_min_dict[i]]
                    #print "what's in VNF_to_math:",'\n',VNF_to_match.keys()
                    VNF_2m_bt = 0          #indicate the total aver trsffic in the VNF_to_match dict, it is used to normalize the traffic
                    for k,v_t in VNF_to_match.iteritems():
                        VNF_2m_bt += v_t[0]
                    host_2r_d = {}
                    k = host_start
                    for k in range(host_start,host_end):
                        if k in host_candinate:    #the id must in host_candinate
                            if k == i:
                                continue
                            else:
                                if VNF_t_dict[k].has_key(VNF_min_dict[i]):                    #ensure that host_to_receive has the same vnf
                                    h_2r_mtr = np.zeros([3,len(x)])
                                    h_2r_mtr[0] = host_lib[k].CPU - host_lib[k].cpu_res
                                    h_2r_mtr[1] = host_lib[k].Mem - host_lib[k].mem_res
                                    h_2r_mtr[2] = host_lib[k].Band - host_lib[k].band_res
                                    host_2r_d[k] = self.calc_d_val_pearson(VNF_i_mtr, h_2r_mtr)
                    #print "what's in host_2r_d:",'\n', host_2r_d.keys()
                    host_2r_dt = 0                        #like the VNF_2m_bt
                    for l,v_d in host_2r_d.iteritems():
                        host_2r_dt += v_d
                    factor_dict = {}
                    m = host_start
                    #for m in range(host_start,host_end):
                    for m in host_2r_d.keys():
                        if m in host_candinate:    #the id must in host_candinate
                            if m == i:
                                continue
                            else:
                                factor_dict[m] = alpha * host_2r_d[m]/host_2r_dt + beta * VNF_to_match[m][0]/(VNF_2m_bt+EPSILON_F)               #it is a balance between d and aver traffic volume
                    #find the host id that is going to receive the request
                    if len(factor_dict):
                        host_r_id = self.search_max_dict(factor_dict)       #set the host has max factor as the receive host
                        #mignate the vnf that has the min aver traffic
                        #print "mignate vnf:",'\t',"src",i,"dst",host_r_id,"vnf_id",VNF_min_dict[i]
                        self.mignate_vnf(i, host_r_id, VNF_t_dict, host_lib[i], host_lib[host_r_id], VNF_min_dict[i], VNF_t_dict[i][VNF_min_dict[i]],requests_lib)
                        #self.verify_req_container_src(host_lib[i])
                    else:
                        jump_out_while += 1
                        if i not in jump_out_while_host_id:
                            jump_out_while_host_id.append(i)
                        #print "jump_out_while_host_id(before eliminate_exceed)", jump_out_while_host_id
            #obsoleted_host = []    
            if (self.eliminate_exceed(host_start,host_end, host_lib, VNF_t_dict, requests_lib) == True):
                if len(jump_out_while_host_id) == len(host_candinate): #does it need further action?
                    #print "jump out the while loops(pod_i)",jump_out_while
                    #print "jump_out_while_host_id is:", jump_out_while_host_id
                    #for h_id,vnf_t_value in VNF_t_dict.iteritems():
                        #print "host id",h_id,'\t',"vnf species:",vnf_t_value.keys(),' ',"vnf number:", len(vnf_t_value.keys())
                    #print "host_candinates are:",host_candinate
                    break
                else:
                    #if all host that overuse the resource can be solved
                    print "pick out the host that can not be further used"
                    #obsoleted_host = self.check_host(host_lib_ori, host_candinate, host_lib)
                    #print "obsloleted_host",obsoleted_host
                    #k = 0
                    #for k in range(len(obsoleted_host)):                  #pick out the host that can not be further used
                    #    host_candinate.remove(obsoleted_host[k])
                    
                    #pick out the host that has 0 vnf that is to say the host is empty
                    obsoleted_host_2 = []
                    for host_id_2 in host_candinate:
                        if len(host_lib[host_id_2].vnf_container) == 0:
                            obsoleted_host_2.append(host_id_2)
                    print "obsoleted_host_2",obsoleted_host_2
                    k2 = 0
                    for k2 in range(len(obsoleted_host_2)):
                        host_candinate.remove(obsoleted_host_2[k2])
                    print "host_candinate", host_candinate
                    
            else:
                print "all host cannot be eliminated in a circle, restore the ori host_lib & VNF_t_dicts"
                host_lib = copy.deepcopy(host_lib_backup)
                VNF_t_dict = copy.deepcopy(VNF_t_dict_backup)
                for h_id,vnf_t_value in VNF_t_dict.iteritems():
                    print "host id",h_id,'\t',"vnf species:",vnf_t_value.keys(),' ',"vnf number:", len(vnf_t_value.keys())
                #if all host cannot be eliminated in a circle, what to do next?
                #the simplest way
                break
        
        print "the adjustment process in pods has been done!"
        
        for host_index in range(host_id+1):
            print "host_id", host_index, "vnf_vontainer", host_lib[host_index].vnf_container
        
        return host_lib
    
    def pod_adjust(self, host_id, K, host_lib, requests_lib):
        alpha = 0.5                    #ratio of d
        beta = 0.5                     #ratio of aver traffic volume
        host_rem = (host_id+1) % (K*K/4)
        full_use_pod_num = (host_id+1 - host_rem)/(K*K/4)
        pod_i = 0
        print "***count vnf species in each host***"
        VNF_t_dict_initial = self.count_species(host_id, host_lib)
        #print "what's in VNF_t_dict_initial",'\n', VNF_t_dict_initial
        print "count species over"
        #print "re_modify host resource lefted"
        #self.re_modify_host_res(host_id,host_lib)
        print "***count the traffic go through each vnf in each host***"
        VNF_t_dict = self.count_resource_vnf(host_id, VNF_t_dict_initial, host_lib)
        #print "what's in VNF_t_dict",'\n', VNF_t_dict.keys()
        print "count traffic over"
        """
        deal with full_use_pod_num
        """
        print "full use pod num is:", full_use_pod_num
        print "host_rem is:",host_rem
        print "*****start to adjust the full use pod*****"
        while(pod_i <= full_use_pod_num-1):
            #print "pod id is:",pod_i
            host_start = int(pod_i*K*K/4)
            host_end = int(host_start + K*K/4)
            host_candinate = []              #store all host that can be used to adjust the vnf
            host_i = host_start
            for host_i in range(host_start,host_end):
                host_candinate.append(host_i)
            #print "what's in host_candinate:",'\n',host_candinate
            #jump_out_while = 0
            #jump_out_while_host_id = []
            i = host_start
            host_lib_ori = copy.deepcopy(host_lib)
            while(len(host_candinate) > 1):
                jump_out_while = 0
                jump_out_while_host_id = []
                VNF_min_dict = {}
                for i in range(host_start,host_end):         #find the min vnf(based on aver traffic through it) in each host
                    if i in host_candinate:
                        if len(VNF_t_dict[i]):               #i must not be empty
                            VNF_min_dict[i] = self.search_min_dict(VNF_t_dict[i])
                        else:             #all the vnfs have been moved
                            host_candinate.remove(i)
                #print "what's in VNF_min_dict:",'\n',VNF_min_dict
                #d & traffic volume
                #backup to restore
                if len(VNF_t_dict) == 0:
                    break
                
                VNF_t_dict_backup = copy.deepcopy(VNF_t_dict)                        
                host_lib_backup = copy.deepcopy(host_lib)
                #print "what's in host_candinate (in while):",'\n',host_candinate
                i = host_start
                for i in range(host_start,host_end):
                    #get the indictor to choose
                    if i in host_candinate:    #the id must in host_candinate
                        VNF_to_match = {}        #store the same kind vnf in other's host
                        VNF_i_mtr = np.zeros([3,len(x)])
                        j = host_start
                        for j in range(host_start,host_end):
                            if j in host_candinate:    #the id must in host_candinate
                                if j == i: #establish the resource matrix
                                    #continue
                                    VNF_i_mtr[0] = VNF_t_dict[j][VNF_min_dict[i]][3]                 #cpu
                                    VNF_i_mtr[1] = VNF_t_dict[j][VNF_min_dict[i]][2]                 #mem
                                    VNF_i_mtr[2] = VNF_t_dict[j][VNF_min_dict[i]][1]                 #band
                                else:
                                    if VNF_t_dict[j].has_key(VNF_min_dict[i]):                    #ensure that host_to_receive has the same vnf
                                        VNF_to_match[j] = VNF_t_dict[j][VNF_min_dict[i]]
                        #print "what's in VNF_to_math:",'\n',VNF_to_match.keys()
                        VNF_2m_bt = 0          #indicate the total aver trsffic in the VNF_to_match dict, it is used to normalize the traffic
                        for k,v_t in VNF_to_match.iteritems():
                            VNF_2m_bt += v_t[0]
                        host_2r_d = {}
                        k = host_start
                        for k in range(host_start,host_end):
                            if k in host_candinate:    #the id must in host_candinate
                                if k == i:
                                    continue
                                else:
                                    if VNF_t_dict[k].has_key(VNF_min_dict[i]):                    #ensure that host_to_receive has the same vnf
                                        h_2r_mtr = np.zeros([3,len(x)])
                                        h_2r_mtr[0] = host_lib[k].CPU - host_lib[k].cpu_res
                                        h_2r_mtr[1] = host_lib[k].Mem - host_lib[k].mem_res
                                        h_2r_mtr[2] = host_lib[k].Band - host_lib[k].band_res
                                        host_2r_d[k] = self.calc_d_val_pearson(VNF_i_mtr, h_2r_mtr)
                        #print "what's in host_2r_d:",'\n', host_2r_d.keys()
                        host_2r_dt = 0                        #like the VNF_2m_bt
                        for l,v_d in host_2r_d.iteritems():
                            host_2r_dt += v_d
                        factor_dict = {}
                        m = host_start
                        #for m in range(host_start,host_end):
                        for m in host_2r_d.keys():
                            if m in host_candinate:    #the id must in host_candinate
                                if m == i:
                                    continue
                                else:
                                    factor_dict[m] = alpha * host_2r_d[m]/host_2r_dt + beta * VNF_to_match[m][0]/(VNF_2m_bt+EPSILON_F)               #it is a balance between d and aver traffic volume
                        #find the host id that is going to receive the request
                        if len(factor_dict):
                            host_r_id = self.search_max_dict(factor_dict)       #set the host has max factor as the receive host
                            #mignate the vnf that has the min aver traffic
                            #print "mignate vnf:",'\t',"src",i,"dst",host_r_id,"vnf_id",VNF_min_dict[i]
                            self.mignate_vnf(i, host_r_id, VNF_t_dict, host_lib[i], host_lib[host_r_id], VNF_min_dict[i], VNF_t_dict[i][VNF_min_dict[i]],requests_lib)
                            #self.verify_req_container_src(host_lib[i])
                        else:
                            jump_out_while += 1
                            if i not in jump_out_while_host_id:
                                jump_out_while_host_id.append(i)
                            #print "jump_out_while_host_id(before eliminate_exceed)", jump_out_while_host_id
                #obsoleted_host = []    
                if (self.eliminate_exceed(host_start,host_end, host_lib, VNF_t_dict, requests_lib) == True):
                    if len(jump_out_while_host_id) == len(host_candinate): #does it need further action?
                        #print "jump out the while loops(pod_i)",jump_out_while
                        #print "jump_out_while_host_id is:", jump_out_while_host_id
                        #for h_id,vnf_t_value in VNF_t_dict.iteritems():
                            #print "host id",h_id,'\t',"vnf species:",vnf_t_value.keys(),' ',"vnf number:", len(vnf_t_value.keys())
                        #print "host_candinates are:",host_candinate
                        break
                    else:
                        #if all host that overuse the resource can be solved
                        #print "pick out the host that can not be further used"
                        obsoleted_host = self.check_host(host_lib_ori, host_candinate, host_lib)
                        #print "obsloleted_host",obsoleted_host
                        k = 0
                        for k in range(len(obsoleted_host)):                  #pick out the host that can not be further used
                            host_candinate.remove(obsoleted_host[k])
                        
                        #pick out the host that has 0 vnf that is to say the host is empty
                        obsoleted_host_2 = []
                        for host_id_2 in host_candinate:
                            if len(host_lib[host_id_2].vnf_container) == 0:
                                obsoleted_host_2.append(host_id_2)
                        #print "obsoleted_host_2",obsoleted_host_2
                        k2 = 0
                        for k2 in range(len(obsoleted_host_2)):
                            host_candinate.remove(obsoleted_host_2[k2])
                        #print "host_candinate", host_candinate
                        """
                        for obsoleted_host_id in obsoleted_host:
                            if obsoleted_host_id in jump_out_while_host_id:
                                jump_out_while_host_id.remove(obsoleted_host_id)
                        if len(jump_out_while_host_id) == len(host_candinate):
                            break
                        """
                        
                else:
                    print "all host cannot be eliminated in a circle, restore the ori host_lib & VNF_t_dicts"
                    host_lib = copy.deepcopy(host_lib_backup)
                    VNF_t_dict = copy.deepcopy(VNF_t_dict_backup)
                    #for h_id,vnf_t_value in VNF_t_dict.iteritems():
                        #print "host id",h_id,'\t',"vnf species:",vnf_t_value.keys(),' ',"vnf number:", len(vnf_t_value.keys())
                    #if all host cannot be eliminated in a circle, what to do next?
                    #the simplest way
                    break  
            pod_i += 1
        
        """
        deal with the rest host,
        host_candinate store all the host id, if a host can not stiasfy the check_host condition,
        then the host will be removed from the host_candinate list, and repeat above process until the list is empty
        """
        
        if host_rem != 0:
            print "****start to deal with the rest host(host_rem)****"
            host_start = int(pod_i*K*K/4)
            host_end = int(host_start + host_rem)
            host_candinate = []                #store all the host id
            host_i = host_start
            for host_i in range(host_start,host_end):
                host_candinate.append(host_i)
            #print "what's in host_candinate:",'\n',host_candinate    
            i = host_start
            host_lib_ori = copy.deepcopy(host_lib)
            #jump_out_while = 0
            #jump_out_while_host_id = []
            while(len(host_candinate)):
                jump_out_while = 0
                jump_out_while_host_id = []
                VNF_min_dict = {}
                for i in range(host_start,host_end):         #find the min vnf(based on aver traffic going through it) in each host
                    if i in host_candinate:
                        if len(VNF_t_dict[i]):
                            VNF_min_dict[i] = self.search_min_dict(VNF_t_dict[i])
                        else:
                            host_candinate.remove(i)
                #print "what's in VNF_min_dict (host_rem):",'\n',VNF_min_dict
                #d & traffic volume
                #backup to restore
                if len(VNF_t_dict) == 0:
                    break
                VNF_t_dict_backup = copy.deepcopy(VNF_t_dict)                        
                host_lib_backup = copy.deepcopy(host_lib)
                
                #print "what's in host_candinate(in while host_rem):",'\n',host_candinate    
                i = host_start
                for i in range(host_start,host_end):
                    #get the indictor to choose
                    if i in host_candinate:    #the id must in host_candinate
                        VNF_to_match = {}        #store the same kind vnf in other's host
                        VNF_i_mtr = np.zeros([3,len(x)])
                        j = host_start
                        for j in range(host_start,host_end):
                            if j in host_candinate:    #the id must in host_candinate
                                if j == i: #establish the resource matrix
                                    #continue
                                    VNF_i_mtr[0] = VNF_t_dict[j][VNF_min_dict[i]][3]                 #cpu
                                    VNF_i_mtr[1] = VNF_t_dict[j][VNF_min_dict[i]][2]                 #mem
                                    VNF_i_mtr[2] = VNF_t_dict[j][VNF_min_dict[i]][1]                 #band
                                else:
                                    if VNF_t_dict[j].has_key(VNF_min_dict[i]):                    #ensure that host_to_receive has the same vnf
                                        VNF_to_match[j] = VNF_t_dict[j][VNF_min_dict[i]]
                        #print "what's in VNF_to_math:",'\n',VNF_to_match.keys()
                        VNF_2m_bt = 0          #indicate the total aver trsffic in the VNF_to_match dict, it is used to normalize the traffic
                        for k,v_t in VNF_to_match.iteritems():
                            VNF_2m_bt += v_t[0]
                        host_2r_d = {}
                        k = host_start
                        for k in range(host_start,host_end):
                            if k in host_candinate:    #the id must in host_candinate
                                if k == i:
                                    continue
                                else:
                                    if VNF_t_dict[k].has_key(VNF_min_dict[i]):                    #ensure that host_to_receive has the same vnf
                                        h_2r_mtr = np.zeros([3,len(x)])
                                        h_2r_mtr[0] = host_lib[k].CPU - host_lib[k].cpu_res
                                        h_2r_mtr[1] = host_lib[k].Mem - host_lib[k].mem_res
                                        h_2r_mtr[2] = host_lib[k].Band - host_lib[k].band_res
                                        host_2r_d[k] = self.calc_d_val_pearson(VNF_i_mtr, h_2r_mtr)
                        #print "what's in host_2r_d:",'\n', host_2r_d.keys()
                        host_2r_dt = 0                        #like the VNF_2m_bt
                        for l,v_d in host_2r_d.iteritems():
                            host_2r_dt += v_d
                        factor_dict = {}
                        m = host_start
                        #for m in range(host_start,host_end):
                        for m in host_2r_d.keys():
                            if m in host_candinate:    #the id must in host_candinate
                                if m == i:
                                    continue
                                else:
                                    factor_dict[m] = alpha * host_2r_d[m]/host_2r_dt + beta * VNF_to_match[m][0]/(VNF_2m_bt+EPSILON_F)               #it is a balance between d and aver traffic volume
                        #find the host id that is going to receive the request
                        if (len(factor_dict)):
                            host_r_id = self.search_max_dict(factor_dict)
                            #mignate the vnf that has the min aver traffic
                            #print "mignate vnf(host_rem != 0):",'\t',"src",i,"dst",host_r_id,"vnf_id",VNF_min_dict[i]
                            self.mignate_vnf(i, host_r_id, VNF_t_dict, host_lib[i], host_lib[host_r_id], VNF_min_dict[i], VNF_t_dict[i][VNF_min_dict[i]], requests_lib)
                            #self.verify_req_container_src(host_lib[i])
                        else:
                            #that is to say, there is no host has the same vnf with the host processing,
                            #but it is not means that the vnf in the processing host can not be mignate,
                            #the vnf still can be mignate to the host that doesn't has the same vnf
                            jump_out_while += 1
                            if i not in jump_out_while_host_id:
                                jump_out_while_host_id.append(i)
                            #print "jump_out_while_host_id(before eliminate_exceed in host_rem != 0)", jump_out_while_host_id
                
                if (self.eliminate_exceed(host_start,host_end, host_lib, VNF_t_dict, requests_lib) == True):
                    if len(jump_out_while_host_id) == len(host_candinate):
                        #print "jump out the while loops(host_rem != 0)",jump_out_while
                        #print "jump_out_while_host_id is:", jump_out_while_host_id
                        #for h_id,vnf_t_value in VNF_t_dict.iteritems():
                            #print "host id",h_id,'\t',"vnf species:",vnf_t_value.keys(),' ',"vnf number:", len(vnf_t_value.keys())
                        #print "host_candinates are:",host_candinate
                        break
                    else:
                        #if all host that overuse the resource can be solved
                        
                        obsoleted_host = self.check_host(host_lib_ori, host_candinate, host_lib)
                        k = 0
                        for k in range(len(obsoleted_host)):                  #pick out the host that can not be further used
                            host_candinate.remove(obsoleted_host[k])
                        
                        #pick out the host that has 0 vnf that is to say the host is empty
                        obsoleted_host_2 = []
                        for host_id_2 in host_candinate:
                            if len(host_lib[host_id_2].vnf_container) == 0:
                                obsoleted_host_2.append(host_id_2)
                        k2 = 0
                        for k2 in range(len(obsoleted_host_2)):
                            host_candinate.remove(obsoleted_host_2[k2])
                else:
                        #print "all host cannot be eliminated in a circle, restore the ori host_lib & VNF_t_dicts"
                        host_lib = copy.deepcopy(host_lib_backup)
                        VNF_t_dict = copy.deepcopy(VNF_t_dict_backup)
                        #for h_id,vnf_t_value in VNF_t_dict.iteritems():
                            #print "host id",h_id,'\t',"vnf species:",vnf_t_value.keys(),' ',"vnf number:", len(vnf_t_value.keys())
                        #if all host cannot be eliminated in a circle, what to do next?
                        #the simplest way
                        break    
        print "the adjustment process in pods has been done!"
        
        #for host_index in range(host_id+1):
            #print "host_id", host_index, "vnf_vontainer", host_lib[host_index].vnf_container
        
        return host_lib
        
    
    def pod_adjust_saved(self, host_id, host_lib):
        host_empty = []
        host_non_empty = []
        for i in range(host_id+1):
            if len(host_lib[i].vnf_container) == 0:                                 
                host_empty.append(i)
            else:
                host_non_empty.append(i)
        pod_adjust_saved_num = len(host_empty)
        print "pod_adjust_saved_num is:",'\t', pod_adjust_saved_num
        return pod_adjust_saved_num
    
    def check_host(self, host_lib_ori, host_candinate, host_lib):
        """
        check if the host can be further used,
        if the aver of res_band less than 0.5 times of the original,then the host can not be further used
        """
        #band_alpha = 0.5
        ob_h_list = []
        h_candinate = len(host_candinate)
        for i in range(h_candinate):
            if self.get_aver(host_lib[host_candinate[i]].band_res) < band_alpha*self.get_aver(host_lib_ori[host_candinate[i]].band_res):
                ob_h_list.append(host_candinate[i])
        return ob_h_list
        
    
class VNF_adjust_whole(object):
    
    def search_max_dict(self, dict_exp):
        max_val = 0
        for key, value in dict_exp.iteritems():
            if value >= max_val:
                max_val = value
                max_id = key
        return max_id
    
    
    def get_aver(self, value_list):
        list_len = len(value_list)
        v_total = 0
        for i in range(list_len):
            v_total += value_list[i]
        aver = v_total/list_len
        return aver
    
    def search_min_dict(self, dict_exp):
        min_val = INFINITY
        for key, value in dict_exp.iteritems():
            if value <= min_val:
                min_val = value
                min_id = key
        return min_id
    
    def calc_host_resource(self, host_id, host_lib):
        host_resource = {}
        for i in range(host_id+1):
            h_value = host_lib[i]
            host_resource[i] = self.get_aver(h_value.cpu_res)/CPU_TOTAL + self.get_aver(h_value.mem_res)/MEM_TOTAL + self.get_aver(h_value.band_res)/h_band
        #min_h_id = self.search_min_dict(host_resource)
        return host_resource
    
    def sort_host(self, host_resource):
        h_sorted_list = []
        while(len(host_resource)):
            h_max_id = self.search_max_dict(host_resource)
            h_sorted_list.append(h_max_id)
            host_resource.pop(h_max_id)
        return h_sorted_list    
    
    def calc_d_val_pearson(self, mtr1,mtr2):              #calc pearson value
        #print "*****calc d value*****"
        [row,col] = mtr1.shape
        #print "row,col",row,col
        #mtr_temp = np.zeros([row,col])
        #normalization
        #normalization
        """
        normal_max_dict1 = {}
        normal_max_dict2 = {}
        for k in range(row):
            normal_max_dict1[k] = 0
            normal_max_dict2[k] = 0
            for l in range(col):
                if mtr1[k][l] > normal_max_dict1[k]:
                    normal_max_dict1[k] = mtr1[k][l]
                if mtr2[k][l] > normal_max_dict2[k]:
                    normal_max_dict2[k] = mtr2[k][l]
        
        for m in range(row):
            for n in range(col):
                mtr1[m][n] = mtr1[m][n]/normal_max_dict1[m]
                mtr2[m][n] = mtr2[m][n]/normal_max_dict2[m]
            
        """
        #pearson correlation coefficient
        pearson = {}
        for i in range(row):
            pearson[i] = 0
            factor_1 = 0
            factor_2 = 0
            factor_3 = 0
            factor_4 = 0
            factor_5 = 0
            j = 0
            for j in range(col):
                factor_1 += mtr1[i][j]*mtr2[i][j]
                factor_2 += mtr1[i][j]
                factor_3 += mtr2[i][j]
                factor_4 += mtr1[i][j]*mtr1[i][j]
                factor_5 += mtr2[i][j]*mtr2[i][j]
            #print "sqrt_1:",col*factor_4-math.pow(factor_2,2)
            Cross_x = col*factor_4-math.pow(factor_2,2)
            Cross_y = col*factor_5-math.pow(factor_3,2)
            #print "i, Row,Col", i, row, col
            #print "col*factor_4-math.pow(factor_2,2)", col*factor_4-math.pow(factor_2,2)
            #print "col*factor_5-math.pow(factor_3,2)", col*factor_5-math.pow(factor_3,2)
            if Cross_x < 0 and Cross_x > -1*e-8:
                Cross_x = 0
                #print "vector len", len(mtr2[i]), "Mistake vector:", mtr2[i]
                #print "col*factor_5",col*factor_5, "math.pow(factor_3,2)", math.pow(factor_3,2)
            if Cross_y < 0 and Cross_y > -1*e-8:
                Cross_y = 0
            if Cross_x*Cross_y == 0:
                pearson[i] = 1
            else:
                pearson[i] = (col * factor_1 - factor_2 * factor_3)/(math.sqrt(Cross_x*Cross_y))
        d = 0
        for i in range(row):
            d += pearson[i]/3
        
        return d
     
    def calc_d_val(self, mtr1,mtr2):
        [row,col] = mtr1.shape
        #mtr_temp = np.zeros([row,col])
        #normalization
        #normalization
        normal_max_dict1 = {}
        normal_max_dict2 = {}
        for k in range(row):
            normal_max_dict1[k] = 0
            normal_max_dict2[k] = 0
            for l in range(col):
                if mtr1[k][l] > normal_max_dict1[k]:
                    normal_max_dict1[k] = mtr1[k][l]
                if mtr2[k][l] > normal_max_dict2[k]:
                    normal_max_dict2[k] = mtr2[k][l]
        
        for m in range(row):
            for n in range(col):
                mtr1[m][n] = mtr1[m][n]/normal_max_dict1[m]
                mtr2[m][n] = mtr2[m][n]/normal_max_dict2[m]
        
        mtr_abs = 0
        for i in range(row):
            for j in range(col):
                mtr_abs += math.fabs(mtr1[i][j]-mtr2[i][j])
        d = 1/(mtr_abs+EPSILON_F)
        return d
    
    def search_min_request(self, req_in_host, req_left_lib):
        min_d = INFINITY
        mtr_temp = np.zeros([3,len(x)])
        for req_id_h, req_value_h in req_in_host.iteritems():
            #print "req_value_h is:",'\n',req_value_h
            for req_j in req_value_h:
                mtr_temp[0] += req_j.cpu
                mtr_temp[1] += req_j.mem
            mtr_temp[2] += req_value_h[0].in_band                              #request band is determined by the first VNF
        for req_id, req_value in req_left_lib.iteritems():
            mtr_temp_1 = np.zeros([3,len(x)])
            for req_i in req_value:
                mtr_temp_1[0] += req_i.cpu
                mtr_temp_1[1] += req_i.mem
            mtr_temp_1[2] = req_value[0].in_band                           #request band is determined by the first VNF
            d = self.calc_d_val_pearson(mtr_temp,mtr_temp_1)
            if d <= min_d:
                min_d = d
                erase_id = req_id
        return erase_id
        
    
    def verify_resource_left_dst(self, request, host):
        req_len = len(request)
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for i in range(req_len):
            req_cpu += request[i].cpu
            req_mem += request[i].mem
        #req_band = request[0].band                            #request band is determined by the first VNF
        
        host.verify_cpu(req_cpu)
        host.verify_mem(req_mem)
        
        
    def verify_resource_left_src(self, request, host):
        req_len = len(request)
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for i in range(req_len):
            req_cpu += request[i].cpu
            req_mem += request[i].mem
        #req_band = request[0].band                            #request band is determined by the first VNF
        host.cpu_res += req_cpu
        host.mem_res += req_mem
        
        
    
        
    def check_value(self, list_ori, list_cmp):
        list_len = len(list_ori)
        break_t = 0
        for i in range(list_len):
            if list_cmp[i] >= list_ori[i]:
                break_t += 1
        return break_t
        
    
    
    def check_end(self, path, link_map, request, host):
        sw_band_bt_total = 0
        req_len = len(request)
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for i in range(req_len):
            req_cpu += request[i].cpu
            req_mem += request[i].mem
            req_band = request[i].in_band                            #request band is determined by the first VNF
        
        for a,b in zip(path[:-1],path[1:]):
            sw_band_bt_total  += self.check_value(link_map[a][b], req_band)
        h_cpu_bt_total = self.check_value(host.cpu_res, req_cpu)                       #break times
        h_mem_bt_total = self.check_value(host.mem_res, req_mem)
        h_band_bt_total = self.check_value(host.band_res, req_band)
        
        if sw_band_bt_total or h_cpu_bt_total or h_mem_bt_total or h_band_bt_total:                     #to be changed, a threshold can be defined
            return False
        else:
            return True
    
    def check_host(self, vnf_temp, host_dst):
        req_cpu = vnf_temp.cpu
        req_mem = vnf_temp.mem
        req_band = vnf_temp.in_band + vnf_temp.out_band
        #req_band = vnf_temp.in_band
        
        h_cpu_bt_total = self.check_value(host_dst.cpu_res, req_cpu)                       #break times
        h_mem_bt_total = self.check_value(host_dst.mem_res, req_mem)
        h_band_bt_total = self.check_value(host_dst.band_res, req_band)
        
        if h_cpu_bt_total or h_mem_bt_total or h_band_bt_total:                     #to be changed, a threshold can be defined
            return False
        else:
            host_dst.verify_band(req_band)
            host_dst.verify_cpu(req_cpu)
            host_dst.verify_mem(req_mem)
            return True
    
    def verify_band_and_req_dict(self, transformed_request, host_dst, host_src):
        
        def vnfID_in_req(vnf_list):
            vnfID_list = []
            for vnf_i in vnf_list:
                vnfID_list.append(vnf_i.id)
            return vnfID_list
            
        
        #the host_dst has no the same req
        def vnf_first_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id): #the req has not been moved clear
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list:       #next_vnf_id in the src_host
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                else:            #next_vnf_id in the third part, because the host_dst has no the same req_id,so it has no the next_vnf_id definitely.
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
            else:
                #the next_vnf_id is not none, and the next_vnf_id is in the third part,because the host_dst has no the same req_id,so it has no the next_vnf_id definitely.
                dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
        
        def vnf_last_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list:
                    dst_host.band_res -= trans_vnf.in_band
                    src_host.band_res -= trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
            else:
                src_host.band_res += trans_vnf.in_band
                dst_host.band_res -= trans_vnf.in_band
        
        def vnf_middle_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list and trans_vnf.next_vnf_id in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list and trans_vnf.next_vnf_id not in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id not in vnf_id_list and trans_vnf.next_vnf_id in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
            else:
                dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
        
        #the host_dst has the same req
        def vnf_first_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
        
        def vnf_last_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_src:
                    src_host.band_res -= trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res += trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_dst:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res += trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
        
        def vnf_middle_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id not in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.out_band - trans_vnf.in_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.out_band + trans_vnf.in_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_src and trans_vnf.pre_vnf_id not in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.out_band - trans_vnf.in_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
        
        for trans_id, trans_value in transformed_request.iteritems():
            if host_dst.req_container_dict.has_key(trans_id):
                #self.verify_req_container_src(host_src)
                host_dst.req_container_dict[trans_id].append(trans_value)
                if trans_value.pre_vnf_id == None:
                    vnf_first_1(host_dst, host_src, trans_id, trans_value)
                elif trans_value.next_vnf_id == None:
                    vnf_last_1(host_dst, host_src, trans_id, trans_value)
                else:
                    vnf_middle_1(host_dst, host_src, trans_id, trans_value)
                #host_dst.req_container_dict[trans_id] = self.right_sequence_vnf(host_dst.req_container_dict[trans_id], requests_lib[trans_id])
            else:
                #self.verify_req_container_src(host_src)
                host_dst.req_container_dict[trans_id] = [trans_value]
                host_dst.req_container.append(trans_id)
                if trans_value.pre_vnf_id == None and trans_value.next_vnf_id == None:
                    host_dst.band_res -= trans_value.in_band
                    host_src.band_res += trans_value.in_band
                elif trans_value.pre_vnf_id == None and trans_value.next_vnf_id != None:
                    vnf_first_0(host_dst, host_src, trans_id, trans_value)
                elif trans_value.next_vnf_id == None and trans_value.pre_vnf_id != None:
                    vnf_last_0(host_dst, host_src, trans_id, trans_value)
                else:
                    vnf_middle_0(host_dst, host_src, trans_id, trans_value)
    
    
    def mignate_request(self, route, core_link, e_link, request_lib, host_id, host_lib):
        """
        After the "pod_adjust", we can know that vnfs in each host tend to be different, that is to say, requests in the host can not be mignated freely.
        the host that intends to receive the mignated request must has all the vnfs in the requests, or it'll have to establish some new vnfs, which will
        consume more "brc". and this break our saving rules.
        """
        #first, we should process the empty host
        host_empty = []
        host_non_empty = []
        for i in range(host_id+1):
            if len(host_lib[i].vnf_container) == 0:                                 
                host_empty.append(i)
            else:
                host_non_empty.append(i)
        
        for j in host_non_empty:
            if len(host_lib[j].vnf_container) > 1:
                continue
            if len(host_lib[j].vnf_container) == 1:        #only mignate the requests in the host having only one vnf
                host_lib_backup = copy.deepcopy(host_lib)
                vnf_id_j = host_lib[j].vnf_container[0]
                req_removed_list = []
                for req_id,req_value in host_lib[j].req_container_dict.iteritems():
                    for k in host_non_empty:
                        if k in host_empty:
                            continue
                        else:
                            if k == j:
                                continue
                            if vnf_id_j not in host_lib[k].vnf_container:
                                continue
                            if vnf_id_j in host_lib[k].vnf_container:
                                sw_id = route._get_switch(k,e_link)
                                path = route._get_path(ROOT_ID, sw_id)
                                if self.check_end(path, core_link, req_value, host_lib[k]):
                                    req_removed_list.append((req_id,k))
                                    self.verify_resource_left_dst(req_value, host_lib[k])
                                    #verify resource left in src host
                                    self.verify_resource_left_src(req_value, host_lib[j])
                                    
                                    break
                
                for req_id_removed in req_removed_list:
                    trans_request = []
                    trans_request = copy.deepcopy(host_lib[j].req_container_dict[req_id_removed[0]])
                    #print "trans_request", trans_request
                    #print "hp_value.req_container_dict[req_id]",hp_value.req_container_dict[req_id]
                    for trans_vnf in trans_request:
                        trans_req = {}
                        trans_req[req_id_removed[0]] = trans_vnf
                        
                        for vnf_index in host_lib[j].req_container_dict[req_id_removed[0]]:
                            if vnf_index.id == trans_vnf.id:
                                vnf_removed = vnf_index
                        
                        host_lib[j].req_container_dict[req_id_removed[0]].remove(vnf_removed)
                        if len(host_lib[j].req_container_dict[req_id_removed[0]]) == 0:
                            host_lib[j].req_container_dict.pop(req_id_removed[0])
                            host_lib[j].req_container.remove(req_id_removed[0])
                        self.verify_band_and_req_dict(trans_req, host_lib[req_id_removed[1]], host_lib[j])
                    
                    
                
                if len(host_lib[j].req_container) == 0:
                    host_empty.append(j)
                    host_lib[j].vnf_container = []
                    host_lib[j].cpu_res += cpu_brc_per
                    host_lib[j].mem_res += mem_brc_per
                else:
                    host_lib = copy.deepcopy(host_lib_backup)           #if can not mignate all the requests, restore the "host_lib"
            
        print "saved host after VNF_adjust_whole is:",host_empty
        return host_empty    
        
    
    def count_species(self, host_lib):           #count the vnf species in the host
        VNF_list = []
        for req_id,req_value in host_lib.req_container_dict.iteritems():
            for VNF in req_value:
                if VNF.id in VNF_list:
                    pass
                else:
                    VNF_list.append(VNF.id)
            
        host_lib.vnf_container = copy.deepcopy(VNF_list)          #each host has a vnf list to hold all the different vnf
        return VNF_list


    def mignate_request_all(self, route, core_link, e_link, request_lib, host_id, host_lib):
        """
        first, sorted the host based on the aver resource
        
        h_sorted_list is a sorted list of host ids
        #start_p points the first position of the h_sorted_list, then it increses if all the requests in the responding host are moved
        end_p points the last position of the h_sorted_list, then it decreses if the responding host can not receive the removed requests
        
        secondly, pick the min host, then pick out the requests in the min_host,
        put them into the others' host, until the min host is empty , then we save one host.
        thirdly, repeat, until the start_p >= end_p 
        """
        host_empty = []
        host_non_empty = []
        for i in range(host_id+1):
            if len(host_lib[i].vnf_container) == 0:                                 
                host_empty.append(i)
            else:
                host_non_empty.append(i)
        print "****start to calc the rest resource in each host****"
        host_resource = self.calc_host_resource(host_id, host_lib)
        print "****start to sort the host list based on aver resource****"
        h_sorted_list = self.sort_host(host_resource)
        hsl_len = len(h_sorted_list)
        
        start_p = len(host_empty)
        end_p = hsl_len - 1
        #print "****the start_p and end_p are****",start_p, end_p
        shift_p = 0
        abandoned_host_id = []
        while(1):
            #print "****the start_p and end_p are(start)****",start_p, end_p
            host_src_id = h_sorted_list[start_p]
            host_src = host_lib[host_src_id]
            host_lib_backup = copy.deepcopy(host_lib)
            while(start_p < end_p):
                #print "****the start_p and end_p are****",start_p, end_p
                host_dst_id = h_sorted_list[end_p]
                host_dst = host_lib[host_dst_id]
                vnf_len_in_src_ori = self.count_species(host_src)
                
                #print "req_in_host_dst is:",req_in_host_dst
                #sw_id = route._get_switch(host_id,e_link)
                #path = route._get_path(ROOT_ID, sw_id)
                
                req_removed_list = []
                host_dst_for_check = copy.deepcopy(host_dst)
                for req_src_id,req_src_value in host_src.req_container_dict.iteritems():
                    for vnf_index in req_src_value:
                        if vnf_index.id in host_dst.vnf_container and self.check_host(vnf_index,host_dst_for_check):
                            req_removed_list.append((req_src_id,vnf_index))
                            req_list = [vnf_index]
                            self.verify_resource_left_src(req_list, host_src)
                            self.verify_resource_left_dst(req_list, host_dst)
                
                for req_id_removed in req_removed_list:
                    trans_req = {}
                    trans_req[req_id_removed[0]] = req_id_removed[1]
                        
                    for vnf_index in host_src.req_container_dict[req_id_removed[0]]:
                        if vnf_index.id == req_id_removed[1].id:
                            vnf_removed = vnf_index
                        
                    host_src.req_container_dict[req_id_removed[0]].remove(vnf_removed)
                    if len(host_src.req_container_dict[req_id_removed[0]]) == 0:
                        host_src.req_container_dict.pop(req_id_removed[0])
                        host_src.req_container.remove(req_id_removed[0])
                    self.verify_band_and_req_dict(trans_req, host_dst, host_src)
                vnf_len_in_src_new = self.count_species(host_src)
                host_src.cpu_res += (len(vnf_len_in_src_ori) - len(vnf_len_in_src_new))*cpu_brc_per
                host_src.mem_res += (len(vnf_len_in_src_ori) - len(vnf_len_in_src_new))*mem_brc_per
                for vnf_id_ori in vnf_len_in_src_ori:
                    if vnf_id_ori in vnf_len_in_src_new:
                        continue
                    else:
                        host_src.vnf_pop_container.append(vnf_id_ori)
                
                if len(host_src.req_container_dict.keys()) == 0:
                    #print "host_src", host_src_id, "has been moved clear"
                    host_empty.append(host_src_id)
                    #print "****the start_p and end_p are(len(host_src.req_container) == 0)****",start_p, end_p
                    break
                else:
                    end_p -= 1
            if start_p >= end_p:
                #print "****the start_p and end_p are(start_p >= end_p)****",start_p, end_p
                #host_lib = copy.deepcopy(host_lib_backup)
                #break
                #print "****the start_p and end_p are(if start_p >= end_p)****",start_p, end_p
                shift_p += 1
                if start_p + shift_p >= hsl_len - 1:
                    host_lib = copy.deepcopy(host_lib_backup)
                    #print "****the start_p and start_p + shift_p(before while) are****",start_p, start_p + shift_p
                    break
                abandoned_host_id.append(host_src_id)
                #host_src_id = h_sorted_list[start_p]
                while h_sorted_list[start_p + shift_p] in abandoned_host_id:
                    shift_p += 1
                    if start_p + shift_p >= hsl_len - 1:
                        break
                #abandoned_host_id.append(host_src_id)
                if start_p + shift_p >= hsl_len - 1:
                    host_lib = copy.deepcopy(host_lib_backup)
                    #print "****the start_p and start_p + shift_p(after while) are****",start_p, start_p + shift_p
                    break
                else:
                    host_lib = copy.deepcopy(host_lib_backup)
                    h_sorted_list[start_p] = h_sorted_list[start_p+shift_p]
                    h_sorted_list[start_p+shift_p] = host_src_id
                    end_p = hsl_len - 1
                
            else:
                start_p += 1 + shift_p
                end_p = hsl_len - 1
                shift_p = 0
                if start_p >= end_p:
                    #print "****the start_p and end_p are****",start_p, end_p
                    host_lib = copy.deepcopy(host_lib_backup)
                    break
                
        print "saved host after VNF_adjust_whole is:",len(host_empty)
        #for host_index in range(host_id+1):
            #print "host_id", host_index, "req_vontainer_dict", host_lib[host_index].req_container_dict.keys()
        return (host_empty,host_lib)
     
    
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
    

def VNF_simulation_xFit_randomized(data_back):
    fattree = fat_tree()
    fattree.K = 8
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    core_link_ff = copy.deepcopy(fattree.core_link)
    e_link_ff = copy.deepcopy(fattree.edge_link)
    core_link_ffd = copy.deepcopy(fattree.core_link)
    e_link_ffd = copy.deepcopy(fattree.edge_link)
    core_link_bf = copy.deepcopy(fattree.core_link)
    e_link_bf = copy.deepcopy(fattree.edge_link)
    core_link_bfd = copy.deepcopy(fattree.core_link)
    e_link_bfd = copy.deepcopy(fattree.edge_link)
    core_link_nf = copy.deepcopy(fattree.core_link)
    e_link_nf = copy.deepcopy(fattree.edge_link)
    
    request_lib_ff = copy.deepcopy(data_back)
    host_lib_ff = copy.deepcopy(fattree.host_dict)
    request_lib_ffd = copy.deepcopy(data_back)
    host_lib_ffd = copy.deepcopy(fattree.host_dict)
    request_lib_bf = copy.deepcopy(data_back)
    host_lib_bf = copy.deepcopy(fattree.host_dict)
    request_lib_bfd = copy.deepcopy(data_back)
    host_lib_bfd = copy.deepcopy(fattree.host_dict)
    request_lib_nf = copy.deepcopy(data_back)
    host_lib_nf = copy.deepcopy(fattree.host_dict)
    
    
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_xFit()
    adjust_pod = VNF_adjust_pod()
    adjust_whole = VNF_adjust_whole()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    print "Original xFit"
    """
    req_sequence_list = [0,1,2,3,4,5,6]  #for test
    for req_id in req_sequence_list:
        vnf_in_req_list = []
        for vnf in request_lib[req_id]:
            vnf_in_req_list.append(vnf.id)
        print "request",req_id,'\t',"vnf num:", vnf_in_req_list
    
    """
    req_sequence_list = []
    for i in range(len(request_lib_nf)):
        req_sequence_list.append(i)
    
    random.shuffle(req_sequence_list)
    #host_id = 128
    host_id_nf = req_2_host.NF_map(route, core_link_nf, e_link_nf, request_lib_nf, host_lib_nf, req_sequence_list)
    host_id_ff = req_2_host.FF_map(route, core_link_ff, e_link_ff, request_lib_ff, host_lib_ff, req_sequence_list)
    host_id_bf = req_2_host.BF_map(route, core_link_bf, e_link_bf, request_lib_bf, host_lib_bf, req_sequence_list)
    print "Descending xFit"
    host_id_ffd = req_2_host.FFD_map(route, core_link_ffd, e_link_ffd, request_lib_ffd, host_lib_ffd)
    host_id_bfd = req_2_host.BFD_map(route, core_link_bfd, e_link_bfd, request_lib_bfd, host_lib_bfd)
    #after_pod_adjust_check(host_lib,host_id)
    
    print "*****start pod adjust*****"
    host_lib_nf = adjust_pod.pod_adjust(host_id_nf, fattree.K, host_lib_nf, request_lib_nf)
    host_lib_ff = adjust_pod.pod_adjust(host_id_ff, fattree.K, host_lib_ff, request_lib_ff)
    host_lib_bf = adjust_pod.pod_adjust(host_id_bf, fattree.K, host_lib_bf, request_lib_bf)
    host_lib_ffd = adjust_pod.pod_adjust(host_id_ffd, fattree.K, host_lib_ffd, request_lib_ffd)
    host_lib_bfd = adjust_pod.pod_adjust(host_id_bfd, fattree.K, host_lib_bfd, request_lib_bfd)
    #host_lib = adjust_pod.host_adjust(host_id, fattree.K, host_lib, request_lib_xfit)
    pod_adjust_saved_num_nf = adjust_pod.pod_adjust_saved(host_id_nf, host_lib_nf)
    pod_adjust_saved_num_ff = adjust_pod.pod_adjust_saved(host_id_ff, host_lib_ff)
    pod_adjust_saved_num_bf = adjust_pod.pod_adjust_saved(host_id_bf, host_lib_bf)
    pod_adjust_saved_num_ffd = adjust_pod.pod_adjust_saved(host_id_ffd, host_lib_ffd)
    pod_adjust_saved_num_bfd = adjust_pod.pod_adjust_saved(host_id_bfd, host_lib_bfd)
    
    #after_pod_adjust_check(host_lib,host_id)
    
    print "*****start adjust whole*****"
    #saved_host = adjust_whole.mignate_request(route, core_link, e_link, request_lib_xfit, host_id, host_lib)
    host_info_nf = adjust_whole.mignate_request_all(route, core_link_nf, e_link_nf, request_lib_nf, host_id_nf, host_lib_nf)
    host_info_ff = adjust_whole.mignate_request_all(route, core_link_ff, e_link_ff, request_lib_ff, host_id_ff, host_lib_ff)
    host_info_bf = adjust_whole.mignate_request_all(route, core_link_bf, e_link_bf, request_lib_bf, host_id_bf, host_lib_bf)
    host_info_ffd = adjust_whole.mignate_request_all(route, core_link_ffd, e_link_ffd, request_lib_ffd, host_id_ffd, host_lib_ffd)
    host_info_bfd = adjust_whole.mignate_request_all(route, core_link_bfd, e_link_bfd, request_lib_bfd, host_id_bfd, host_lib_bfd)
    
    print "*********start save info*********"
    saved_host_nf = host_info_nf[0]
    host_lib_nf = host_info_nf[1]
    print "after adjust_whole, the saved host is:",saved_host_nf
    print "saved num is:", len(saved_host_nf)
    print "before adjust, the used host num is:",host_id_nf+1
    
    #with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_NF_r.pickle', 'w') as f:
        #pickle.dump(host_lib_nf,f)
    saved_info_dict_nf = {}
    saved_info_dict_nf[0] = host_id_nf+1
    saved_info_dict_nf[1] = pod_adjust_saved_num_nf
    saved_info_dict_nf[2] = len(saved_host_nf)
    saved_info_dict_nf[3] = saved_host_nf
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_NF_r.json', 'w') as f:
        json.dump(saved_info_dict_nf,f)
    
    saved_host_ff = host_info_ff[0]
    host_lib_ff = host_info_ff[1]
    print "after adjust_whole, the saved host is:",saved_host_ff
    print "saved num is:", len(saved_host_ff)
    print "before adjust, the used host num is:",host_id_ff+1
    
    #with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_FF_r.pickle', 'w') as f:
        #pickle.dump(host_lib_ff,f)
    saved_info_dict_ff = {}
    saved_info_dict_ff[0] = host_id_ff+1
    saved_info_dict_ff[1] = pod_adjust_saved_num_ff
    saved_info_dict_ff[2] = len(saved_host_ff)
    saved_info_dict_ff[3] = saved_host_ff
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_FF_r.json', 'w') as f:
        json.dump(saved_info_dict_ff,f)
    
    saved_host_ffd = host_info_ffd[0]
    host_lib_ffd = host_info_ffd[1]
    print "after adjust_whole, the saved host is:",saved_host_ffd
    print "saved num is:", len(saved_host_ffd)
    print "before adjust, the used host num is:",host_id_ffd+1
    
    #with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_FFD_r.pickle', 'w') as f:
        #pickle.dump(host_lib_ffd,f)
    saved_info_dict_ffd = {}
    saved_info_dict_ffd[0] = host_id_ffd+1
    saved_info_dict_ffd[1] = pod_adjust_saved_num_ffd
    saved_info_dict_ffd[2] = len(saved_host_ffd)
    saved_info_dict_ffd[3] = saved_host_ffd
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_FFD_r.json', 'w') as f:
        json.dump(saved_info_dict_ffd,f)
    
    saved_host_bf = host_info_bf[0]
    host_lib_bf = host_info_bf[1]
    print "after adjust_whole, the saved host is:",saved_host_bf
    print "saved num is:", len(saved_host_bf)
    print "before adjust, the used host num is:",host_id_bf+1
    
    #with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_BF_r.pickle', 'w') as f:
        #pickle.dump(host_lib_bf,f)
    saved_info_dict_bf = {}
    saved_info_dict_bf[0] = host_id_bf+1
    saved_info_dict_bf[1] = pod_adjust_saved_num_bf
    saved_info_dict_bf[2] = len(saved_host_bf)
    saved_info_dict_bf[3] = saved_host_bf
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_BF_r.json', 'w') as f:
        json.dump(saved_info_dict_bf,f)
    
    saved_host_bfd = host_info_bfd[0]
    host_lib_bfd = host_info_bfd[1]
    print "after adjust_whole, the saved host is:",saved_host_bfd
    print "saved num is:", len(saved_host_bfd)
    print "before adjust, the used host num is:",host_id_bfd+1
    
    #with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_BFD_r.pickle', 'w') as f:
        #pickle.dump(host_lib_bfd,f)
    saved_info_dict_bfd = {}
    saved_info_dict_bfd[0] = host_id_bfd+1
    saved_info_dict_bfd[1] = pod_adjust_saved_num_bfd
    saved_info_dict_bfd[2] = len(saved_host_bfd)
    saved_info_dict_bfd[3] = saved_host_bfd
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_BFD_r.json', 'w') as f:
        json.dump(saved_info_dict_bfd,f)
    
    #saved_host = []
    #host_resource_left_check(host_lib_nf,host_id_nf,saved_host_nf)
    print "***********xFit_random is over**********"
    
def VNF_simulation_xFit_non_randomized(data_back):
    fattree = fat_tree()
    fattree.K = 8
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    core_link_ff = copy.deepcopy(fattree.core_link)
    e_link_ff = copy.deepcopy(fattree.edge_link)
    core_link_ffd = copy.deepcopy(fattree.core_link)
    e_link_ffd = copy.deepcopy(fattree.edge_link)
    core_link_bf = copy.deepcopy(fattree.core_link)
    e_link_bf = copy.deepcopy(fattree.edge_link)
    core_link_bfd = copy.deepcopy(fattree.core_link)
    e_link_bfd = copy.deepcopy(fattree.edge_link)
    core_link_nf = copy.deepcopy(fattree.core_link)
    e_link_nf = copy.deepcopy(fattree.edge_link)
    
    request_lib_ff = copy.deepcopy(data_back)
    host_lib_ff = copy.deepcopy(fattree.host_dict)
    request_lib_ffd = copy.deepcopy(data_back)
    host_lib_ffd = copy.deepcopy(fattree.host_dict)
    request_lib_bf = copy.deepcopy(data_back)
    host_lib_bf = copy.deepcopy(fattree.host_dict)
    request_lib_bfd = copy.deepcopy(data_back)
    host_lib_bfd = copy.deepcopy(fattree.host_dict)
    request_lib_nf = copy.deepcopy(data_back)
    host_lib_nf = copy.deepcopy(fattree.host_dict)
    
    
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_xFit()
    adjust_pod = VNF_adjust_pod()
    adjust_whole = VNF_adjust_whole()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    print "Original xFit"
    """
    req_sequence_list = [0,1,2,3,4,5,6]  #for test
    for req_id in req_sequence_list:
        vnf_in_req_list = []
        for vnf in request_lib[req_id]:
            vnf_in_req_list.append(vnf.id)
        print "request",req_id,'\t',"vnf num:", vnf_in_req_list
    
    """
    req_sequence_list = []
    for i in range(len(request_lib_nf)):
        req_sequence_list.append(i)
    
    #random.shuffle(req_sequence_list)
    #host_id = 128
    host_id_nf = req_2_host.NF_map(route, core_link_nf, e_link_nf, request_lib_nf, host_lib_nf, req_sequence_list)
    host_id_ff = req_2_host.FF_map(route, core_link_ff, e_link_ff, request_lib_ff, host_lib_ff, req_sequence_list)
    host_id_bf = req_2_host.BF_map(route, core_link_bf, e_link_bf, request_lib_bf, host_lib_bf, req_sequence_list)
    print "Descending xFit"
    host_id_ffd = req_2_host.FFD_map(route, core_link_ffd, e_link_ffd, request_lib_ffd, host_lib_ffd)
    host_id_bfd = req_2_host.BFD_map(route, core_link_bfd, e_link_bfd, request_lib_bfd, host_lib_bfd)
    #after_pod_adjust_check(host_lib,host_id)
    
    print "*****start pod adjust*****"
    host_lib_nf = adjust_pod.pod_adjust(host_id_nf, fattree.K, host_lib_nf, request_lib_nf)
    host_lib_ff = adjust_pod.pod_adjust(host_id_ff, fattree.K, host_lib_ff, request_lib_ff)
    host_lib_bf = adjust_pod.pod_adjust(host_id_bf, fattree.K, host_lib_bf, request_lib_bf)
    host_lib_ffd = adjust_pod.pod_adjust(host_id_ffd, fattree.K, host_lib_ffd, request_lib_ffd)
    host_lib_bfd = adjust_pod.pod_adjust(host_id_bfd, fattree.K, host_lib_bfd, request_lib_bfd)
    #host_lib = adjust_pod.host_adjust(host_id, fattree.K, host_lib, request_lib_xfit)
    pod_adjust_saved_num_nf = adjust_pod.pod_adjust_saved(host_id_nf, host_lib_nf)
    pod_adjust_saved_num_ff = adjust_pod.pod_adjust_saved(host_id_ff, host_lib_ff)
    pod_adjust_saved_num_bf = adjust_pod.pod_adjust_saved(host_id_bf, host_lib_bf)
    pod_adjust_saved_num_ffd = adjust_pod.pod_adjust_saved(host_id_ffd, host_lib_ffd)
    pod_adjust_saved_num_bfd = adjust_pod.pod_adjust_saved(host_id_bfd, host_lib_bfd)
    
    #after_pod_adjust_check(host_lib,host_id)
    
    print "*****start adjust whole*****"
    #saved_host = adjust_whole.mignate_request(route, core_link, e_link, request_lib_xfit, host_id, host_lib)
    host_info_nf = adjust_whole.mignate_request_all(route, core_link_nf, e_link_nf, request_lib_nf, host_id_nf, host_lib_nf)
    host_info_ff = adjust_whole.mignate_request_all(route, core_link_ff, e_link_ff, request_lib_ff, host_id_ff, host_lib_ff)
    host_info_bf = adjust_whole.mignate_request_all(route, core_link_bf, e_link_bf, request_lib_bf, host_id_bf, host_lib_bf)
    host_info_ffd = adjust_whole.mignate_request_all(route, core_link_ffd, e_link_ffd, request_lib_ffd, host_id_ffd, host_lib_ffd)
    host_info_bfd = adjust_whole.mignate_request_all(route, core_link_bfd, e_link_bfd, request_lib_bfd, host_id_bfd, host_lib_bfd)
    
    print "*********start save info*********"
    saved_host_nf = host_info_nf[0]
    host_lib_nf = host_info_nf[1]
    print "after adjust_whole, the saved host is:",saved_host_nf
    print "saved num is:", len(saved_host_nf)
    print "before adjust, the used host num is:",host_id_nf+1
    
    with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_NF_nr.pickle', 'w') as f:
        pickle.dump(host_lib_nf,f)
    saved_info_dict_nf = {}
    saved_info_dict_nf[0] = host_id_nf+1
    saved_info_dict_nf[1] = pod_adjust_saved_num_nf
    saved_info_dict_nf[2] = len(saved_host_nf)
    saved_info_dict_nf[3] = saved_host_nf
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_NF_nr.json', 'w') as f:
        json.dump(saved_info_dict_nf,f)
    
    saved_host_ff = host_info_ff[0]
    host_lib_ff = host_info_ff[1]
    print "after adjust_whole, the saved host is:",saved_host_ff
    print "saved num is:", len(saved_host_ff)
    print "before adjust, the used host num is:",host_id_ff+1
    
    with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_FF_nr.pickle', 'w') as f:
        pickle.dump(host_lib_ff,f)
    saved_info_dict_ff = {}
    saved_info_dict_ff[0] = host_id_ff+1
    saved_info_dict_ff[1] = pod_adjust_saved_num_ff
    saved_info_dict_ff[2] = len(saved_host_ff)
    saved_info_dict_ff[3] = saved_host_ff
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_FF_nr.json', 'w') as f:
        json.dump(saved_info_dict_ff,f)
    
    saved_host_ffd = host_info_ffd[0]
    host_lib_ffd = host_info_ffd[1]
    print "after adjust_whole, the saved host is:",saved_host_ffd
    print "saved num is:", len(saved_host_ffd)
    print "before adjust, the used host num is:",host_id_ffd+1
    
    with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_FFD_nr.pickle', 'w') as f:
        pickle.dump(host_lib_ffd,f)
    saved_info_dict_ffd = {}
    saved_info_dict_ffd[0] = host_id_ffd+1
    saved_info_dict_ffd[1] = pod_adjust_saved_num_ffd
    saved_info_dict_ffd[2] = len(saved_host_ffd)
    saved_info_dict_ffd[3] = saved_host_ffd
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_FFD_nr.json', 'w') as f:
        json.dump(saved_info_dict_ffd,f)
    
    saved_host_bf = host_info_bf[0]
    host_lib_bf = host_info_bf[1]
    print "after adjust_whole, the saved host is:",saved_host_bf
    print "saved num is:", len(saved_host_bf)
    print "before adjust, the used host num is:",host_id_bf+1
    
    with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_BF_nr.pickle', 'w') as f:
        pickle.dump(host_lib_bf,f)
    saved_info_dict_bf = {}
    saved_info_dict_bf[0] = host_id_bf+1
    saved_info_dict_bf[1] = pod_adjust_saved_num_bf
    saved_info_dict_bf[2] = len(saved_host_bf)
    saved_info_dict_bf[3] = saved_host_bf
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_BF_nr.json', 'w') as f:
        json.dump(saved_info_dict_bf,f)
    
    saved_host_bfd = host_info_bfd[0]
    host_lib_bfd = host_info_bfd[1]
    print "after adjust_whole, the saved host is:",saved_host_bfd
    print "saved num is:", len(saved_host_bfd)
    print "before adjust, the used host num is:",host_id_bfd+1
    
    with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_BFD_nr.pickle', 'w') as f:
        pickle.dump(host_lib_bfd,f)
    saved_info_dict_bfd = {}
    saved_info_dict_bfd[0] = host_id_bfd+1
    saved_info_dict_bfd[1] = pod_adjust_saved_num_bfd
    saved_info_dict_bfd[2] = len(saved_host_bfd)
    saved_info_dict_bfd[3] = saved_host_bfd
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_BFD_nr.json', 'w') as f:
        json.dump(saved_info_dict_bfd,f)
    
    #saved_host = []
    #host_resource_left_check(host_lib_nf,host_id_nf,saved_host_nf)
    print "***********xFit_non_random is over**********"
  
"""    
if __name__ == '__main__':
    fattree = fat_tree()
    fattree.K = 4
    #fattree.__init__()
    fattree.instant_topo()
    #print "Topo Switches:",fattree.switch
    #print "Topo host:", fattree.host_dict
    #print "core_s num, agg_s num, edge_s num, host num", fattree.c_s_num, fattree.agg_s_num, fattree.edge_s_num, fattree.host_num
    #print "core switches:", fattree.core_switch_dict
    #print "agg switches:", fattree.agg_switch_dict
    #print "edge switches:", fattree.edge_switch_dict
    #print "core link:",fattree.core_link
    print "edge link:",fattree.edge_link
    print "pod:", fattree.pod_dict
    sws = fattree.switch
    adj = fattree.core_link
    #print "sws:", sws
    route_al = VNF_route()
    route_al.switches = sws
    route_al.adjacency = adj
    route_al._get_path(7,18)
    sw1 = route_al._get_switch(0,fattree.edge_link)
    sw2 = route_al._get_switch(6,fattree.edge_link)
    print "host 0 switch:", sw1
    print "host 6 switch:", sw2
"""
    
