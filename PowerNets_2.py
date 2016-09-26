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
Cor_threshold_shortest_1 = 1.0
Cor_threshold_shortest_2 = 1.0

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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] & h_NFR_cor[2] in first while(1)-1", h_NFR_cor[0], h_NFR_cor[2]
                    if (h_NFR_cor[0] <= Cor_threshold) and (h_NFR_cor[2] <= Cor_threshold): #and (h_NFR_cor[1] <= Cor_threshold) #+h_NFR_cor[1])/2
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] & h_NFR_cor[2] in first while(1)-2", h_NFR_cor[0], h_NFR_cor[2]
                    if (h_NFR_cor[0] <= Cor_threshold) and (h_NFR_cor[2] <= Cor_threshold):#and (h_NFR_cor[1] <= Cor_threshold) #+h_NFR_cor[1])/2
                        return host_id_cur + i
            i += 1
            print "i in the first while(1)", i
            if i > bound:           # can not find the host satisfies all the cor constraints :(h_NFR_cor[0] <= Cor_threshold) and (h_NFR_cor[2] <= Cor_threshold) 
                break
        print "start to relax the cor constraints 1"
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in second while(1)-1", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in second while(1)-2", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in third while(1)-1", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold_shortest_1):# and (h_NFR_cor[1] <= Cor_threshold):#+h_NFR_cor[1])/2
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in third while(1)-2", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold_shortest_1):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
                        return host_id_cur + i
            i += 1
            print "i in the third while(1)", i
            if i > bound:           # can not find the host satisfies all the cor constraints :(h_NFR_cor[0] <= Cor_threshold)
                break     
        print "start to relax the cor constraints 3"
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in forth while(1)-1", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold_shortest_2): #+h_NFR_cor[1])/2
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in forth while(1)-2", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold_shortest_2): #+h_NFR_cor[1])/2
                        return host_id_cur + i
            i += 1
            print "i in the forth while(1)", i
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
                    print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in FF:2", h_NFR_cor[0]
                    if (h_NFR_cor[0] <= Cor_threshold):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
                        #print "VNF_chain", NFR_lib[i_id][0], "VNF_id", NFR_lib[i_id][1].id, "satisfies the cor constraints."
                        host_lib[host_list[j].id].req_container.append(NFR_temp[0])
                        host_lib[host_list[j].id].req_container_dict[NFR_temp[0]] = [NFR_temp[1]]
                        self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_list[j].id])
                        self.verify_band_left_0(path,core_link,NFR_temp,host_lib[host_list[j].id])
                        inditor = 1                  #an inditor, indite if there is a need to start a new host
                        break                       #find the first host to hold the request, and stop the process
        #again, First-Fit needs to relax the cor_threshold too.
        
        if inditor == 0:
            print "Start to relax the Cor_threshold in First-Fit 1"
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
                        if (h_NFR_cor[0] <= Cor_threshold_shortest_1):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
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
                        if (h_NFR_cor[0] <= Cor_threshold_shortest_1):# and (h_NFR_cor[1] <= Cor_threshold):
                            #print "VNF_chain", NFR_lib[i_id][0], "VNF_id", NFR_lib[i_id][1].id, "satisfies the cor constraints."
                            host_lib[host_list[j].id].req_container.append(NFR_temp[0])
                            host_lib[host_list[j].id].req_container_dict[NFR_temp[0]] = [NFR_temp[1]]
                            self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_list[j].id])
                            self.verify_band_left_0(path,core_link,NFR_temp,host_lib[host_list[j].id])
                            inditor = 1                  #an inditor, indite if there is a need to start a new host
                            break                       #find the first host to hold the request, and stop the process
        #again, First-Fit needs to relax the cor_threshold too.
        if inditor == 0:
            print "Start to relax the Cor_threshold in First-Fit 2"
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
                        print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id, "h_NFR_cor[0] in FF:5", h_NFR_cor[0] 
                        if (h_NFR_cor[0] <= Cor_threshold_shortest_2):# and (h_NFR_cor[1] <= Cor_threshold): #+h_NFR_cor[1])/2
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
                        print "VNF_CHAIN:",NFR_temp[0],"VNF:",NFR_temp[1].id,"h_NFR_cor[0] in FF:6", h_NFR_cor[0] 
                        if (h_NFR_cor[0] <= Cor_threshold_shortest_2):# and (h_NFR_cor[1] <= Cor_threshold):
                            #print "VNF_chain", NFR_lib[i_id][0], "VNF_id", NFR_lib[i_id][1].id, "satisfies the cor constraints."
                            host_lib[host_list[j].id].req_container.append(NFR_temp[0])
                            host_lib[host_list[j].id].req_container_dict[NFR_temp[0]] = [NFR_temp[1]]
                            self.verify_node_resource_left(path, core_link, NFR_temp,host_lib[host_list[j].id])
                            self.verify_band_left_0(path,core_link,NFR_temp,host_lib[host_list[j].id])
                            inditor = 1                  #an inditor, indite if there is a need to start a new host
                            break                       #find the first host to hold the request, and stop the process
        
        #after 2 relax there still no proper host, so it needs to start a new one
        if inditor == 0:
            print "start a new host"
            host_id += 1
            host_list.append(host_lib[host_id])
            print "request seq_num:", NFR_temp[0], '\t', "host_id now", host_id 
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

    
    def PowerNets(self, route, core_link, e_link, NFR_lib, host_lib, VNF_move_id_list):
        print "PowerNets"
        print "sort the requests based on the 90-percentile resource"
        #VNF_move_id_list = self.sort_requests(NFR_lib)
        req_abandoned = []
        host_id = 0
        host_list = []
        host_list.append(host_lib[host_id])
        NFR_num = len(NFR_lib)
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        for i in range(len(VNF_move_id_list)):
            i_id = VNF_move_id_list[i]
            if i == 0:   #the first vnf , put it in the current host directly
                #print "The first VNF"
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
                
        print "req_abandoned is:",req_abandoned,'\n',len(req_abandoned)
        #print "host_list",host_list
        host_obsoleted = []
        host_id_list = []
        for l in range(len(host_list)):
            host_id_list.append(host_list[l].id)
        print "host_list length is:",len(host_list)
        print "host_id in the host_list is:", host_id_list
        for k in range(host_id+1):
            if len(host_lib[k].req_container) == 0:
                host_obsoleted.append(k)
        host_num = host_id+1 - len(host_obsoleted)
        print "host_num", host_num
        print "final host id (PN):", host_id
        host_info = (host_id,host_lib)
        return host_info

def search_max_dict(dict_exp):
    max_val = 0
    for key, value in dict_exp.iteritems():
        if value >= max_val:
            max_val = value
            max_id = key
    return max_id


    
def get_90_percentile(value_list):
    value_list_copy = copy.deepcopy(value_list)
    value_list_copy.sort()
    list_len = len(value_list_copy)
    seq_90_percentile = math.ceil(0.9*list_len)
    return value_list_copy[seq_90_percentile]

def sort_requests(NFR_lib):
    req_90_percentile_source_dict = {}
    for i in range(len(NFR_lib)):
        VNF_i_cpu = 0
        VNF_i_mem = 0
        VNF_i_band = 0
        VNF_i_cpu += get_90_percentile(NFR_lib[i][1].cpu)/CPU_TOTAL
        VNF_i_mem += get_90_percentile(NFR_lib[i][1].mem)/MEM_TOTAL
        VNF_i_band += get_90_percentile(NFR_lib[i][1].in_band + NFR_lib[i][1].out_band)/h_band  #the band resource deserves recosideration
        req_90_percentile_source_dict[i] = VNF_i_cpu + VNF_i_mem + VNF_i_band
    req_id_list = []
    while(len(req_90_percentile_source_dict)):
        max_id = search_max_dict(req_90_percentile_source_dict)
        req_id_list.append(max_id)
        req_90_percentile_source_dict.pop(max_id)
    
    return req_id_list    
    
 


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
    

def VNF_simulation_PowerNets(data_back):
    print "***************PowerNets is Start**********************"
    fattree = fat_tree()
    fattree.K = 16
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    core_link_PN = copy.deepcopy(fattree.core_link)
    e_link_PN = copy.deepcopy(fattree.edge_link)
    
    
    request_lib_PN = copy.deepcopy(data_back)
    host_lib_PN = copy.deepcopy(fattree.host_dict)
    
    
    NFR_lib_PN = NFCR_to_NFRs(request_lib_PN)
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_PowerNets()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    VNF_move_id_list = sort_requests(NFR_lib_PN)
    
    host_info_PN = req_2_host.PowerNets(route, core_link_PN, e_link_PN, NFR_lib_PN, host_lib_PN, VNF_move_id_list)
    
    print "*********start save info*********"
    
    host_id_PN = host_info_PN[0]
    host_lib_PN = host_info_PN[1]
    
    with open('F:/VNF_Review/ldf/simulation/data_store/host_lib_PN_rw.pickle', 'w') as f:
        pickle.dump(host_lib_PN,f)
    saved_info_dict_PN = {}
    saved_info_dict_PN[0] = host_id_PN+1
    print "the used host num(PN) is:",host_id_PN+1
    with open('F:/VNF_Review/ldf/simulation/data_store/saved_info_dict_PowerNets.json', 'w') as f:
        json.dump(saved_info_dict_PN,f)
    
    #print "check host_lib_PN"
    #saved_host_PN = []
    #host_resource_left_check(host_lib_PN,host_id_PN,saved_host_PN)
    
    print "***********PowerNets is over**********"
    

