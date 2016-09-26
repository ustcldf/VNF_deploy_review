# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:35:56 2015

@author: ustcldf
"""
#!/usr/bin/python
from __future__ import division
import os
#from numpy import *
from collections import defaultdict
from scipy import *
import numpy as np
import operator
import time
import math
import pickle

"""
#for optimal test
root_band = 10000
core_band = 10000
a_g_band = 10000
h_band = 10
x = np.arange(0,24.1,0.1)
#x = np.arange(0,144,1)
array_length = len(x)
#VNF species
VNF_SPE = 4
CPU_TOTAL = 10
MEM_TOTAL = 10
band_alpha = 0.5
brc_base = 0.5
cpu_brc_per = np.array([float(brc_base) for i in range(len(x))])
mem_brc_per = np.array([float(brc_base) for i in range(len(x))])
"""

#for Google_Cluster_data
root_band = 100
core_band = 100
a_g_band = 100
h_band = 1.6
x = np.arange(0,144,1)
array_length = len(x)
#VNF species
VNF_SPE = 20
CPU_TOTAL = 1
MEM_TOTAL = 1
band_alpha = 0.5
brc_base = 0.020
cpu_brc_per = np.array([float(brc_base) for i in range(len(x))])
mem_brc_per = np.array([float(brc_base) for i in range(len(x))])
"""
#for Gaussian data set
root_band = 100000
core_band = 100000
a_g_band = 100000
h_band = 150
x = np.arange(0,24.1,0.1)
#x = np.arange(0,144,1)
array_length = len(x)
#VNF species
VNF_SPE = 20
CPU_TOTAL = 100
MEM_TOTAL = 100
band_alpha = 0.5
brc_base = 2.75
cpu_brc_per = np.array([float(brc_base) for i in range(len(x))])
mem_brc_per = np.array([float(brc_base) for i in range(len(x))])
"""

class VNF_server(object):
    def __init__(self, ):
        #self.container = []
        self.req_container = []
        self.req_container_dict = {}
        self.vnf_container = []
        self.vnf_pop_container = []
        self.id = ''
        self.CPU = np.array([float(CPU_TOTAL) for i in range(len(x))])
        self.Mem = np.array([float(MEM_TOTAL) for i in range(len(x))])
        self.Band = np.array([float(h_band) for i in range(len(x))])
        self.cpu_res = np.array([float(CPU_TOTAL) for i in range(len(x))])
        self.mem_res = np.array([float(MEM_TOTAL) for i in range(len(x))])
        self.band_res = np.array([float(h_band) for i in range(len(x))])
    
    def verify_cpu(self, request_cpu):
        self.cpu_res = self.cpu_res - request_cpu           
    
    def verify_mem(self, request_mem):
        self.mem_res = self.mem_res - request_mem             
    
    def verify_band(self, request_band):
        self.band_res = self.band_res - request_band
               
class VNF_switch(object):
    def __init__(self, ):
        self.id = ''
    
    def verify_band(self, ):
        pass
    
class VNF_link(object):
    def __init__(self, ):
        self.adjacency = defaultdict(lambda:defaultdict(lambda:None))
        
    def verify_band(self, ):
        pass
    
class VNF_request(object):
    def __init__(self, ):
        self.id = ''
        self.cpu = None
        self.mem = None
        self.in_band = None
        self.out_band = None
        self.pre_vnf_id = None                #two inditor
        self.next_vnf_id = None
        
    


