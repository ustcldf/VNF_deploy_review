#!/usr/bin/python
from __future__ import division
import os
import re
import os
import shutil
import sys
from scipy import *
import numpy as np
import operator
import time
import math
import pickle
import copy
from collections import defaultdict

from VNF_dp_algorithm_abs import VNF_simulation_abs
from VNF_dp_algorithm_square import VNF_simulation_square
from VNF_dp_algorithm_cosin import VNF_simulation_cosin
from VNF_dp_algorithm_weighted_cosin import VNF_simulation_w_cosin
from VNF_dp_algorithm_pearson import VNF_simulation_pearson
from VNF_dp_algorithm_xFit import VNF_simulation_xFit_randomized, VNF_simulation_xFit_non_randomized
from VNF_dp_algorithm_xFit_based_vnfs import VNF_simulation_xFit_randomized_w, VNF_simulation_xFit_non_randomized_w
from Algorithm_H import Algorithm_H_r, Algorithm_H_nr
from PowerNets_2 import VNF_simulation_PowerNets
from Sher_wood_FF import VNF_simulation_SherWood_FF_BF
from Permutation_Optimal import VNF_simulation_Optimal
from T_SAT_PowerNets import VNF_simulation_T_SAT_PN
from VNF_request_generate import REQUEST_NUM



if __name__ == '__main__':
    time_start = time.time()
    for i in range(10):
    #with open('F:/VNF_Review/ldf\simulation\data\Task_2-20\Sample_num_144\Google_Cluster_task_usage_144_v_2x.pickle', 'r') as f:
        #data_back = pickle.load(f)
    #with open('F:/VNF_Review/ldf\simulation\data\Task_2-20\Sample_num_144\Seg_100\part_900.pickle', 'r') as f:
        #data_back = pickle.load(f)
    #with open('F:/VNF_Review/ldf\simulation\data\Task_2-20\Sample_num_144\Seg_100_v\part_900.pickle', 'r') as f:
        #data_back = pickle.load(f)
    #with open('F:/VNF_Review/ldf\simulation\data\Task_2-20\Sample_num_144\Seg_100_v_2x\part_900.pickle', 'r') as f:
        #data_back = pickle.load(f)
    #with open('F:/VNF_Review/ldf\simulation\data\Task_2-20\Sample_num_144\Seg_10\part_50.pickle', 'r') as f:
        #data_back = pickle.load(f)
    #with open('E:\Google_Cluster_Trace\data_collect\Optimal_test\'+'part_5_0.pickle', 'r') as f:
        #data_back = pickle.load(f)
        with open('F:/VNF_Review/ldf\simulation\data/'+'part_'+str(REQUEST_NUM)+'_'+str(i)+'.pickle', 'r') as f:
            data_back = pickle.load(f)
        """
        print "ABS"
        VNF_simulation_abs(data_back)
        
        print "square"
        VNF_simulation_square(data_back)
        """
        print "cosin"
        VNF_simulation_cosin(data_back)
        
        #print "weighted cosin"
        #VNF_simulation_w_cosin(data_back)
        
        #print "pearson"
        #VNF_simulation_pearson(data_back)
        
        print "xFit_randonm"
        VNF_simulation_xFit_randomized(data_back)
        
        #print "xFit_non_randonm"
        #VNF_simulation_xFit_non_randomized(data_back)
        
        print "xFit_randonm_w"
        VNF_simulation_xFit_randomized_w(data_back)
        
        #print "xFit_non_randonm_w"
        #VNF_simulation_xFit_non_randomized_w(data_back)
        
        print "Algorithm H r"
        Algorithm_H_r(data_back)
        
        #print "Algorithm H nr"
        #Algorithm_H_nr(data_back)
        
        print "PowerNets"
        VNF_simulation_PowerNets(data_back)
        
        print "T_SAT_PowerNets"
        VNF_simulation_T_SAT_PN(data_back)
        
        #print "Optimal"
        #VNF_simulation_Optimal(data_back)
        
        #print "Sher_Wood_FF_BF"
        #VNF_simulation_SherWood_FF_BF(data_back)
        
        
        store_path = 'F:/VNF_Review/ldf/simulation/data_store/'
        os.makedirs(store_path+'data_h_band_1.6_'+str(REQUEST_NUM)+'_'+str(i))
        for subj in os.listdir(store_path):
            #print subj
            #lists = subj.split('.')
            #file_ext = lists[-1]
            #ext_set = ['pickle','json']
            #if os.path.isdir(store_path+subjsubj):
            #if file_ext in ext_set:
            if os.path.isfile(store_path+subj):
                #print "the subj is a file"
                shutil.move(store_path+subj,store_path+'data_h_band_1.6_'+str(REQUEST_NUM)+'_'+str(i))
        
    time_end = time.time()
    time_used = time_end - time_start
    time_second = time_used % 60
    time_ms = (time_used-time_second) / 60
    time_mininute = time_ms % 60
    time_hs = (time_ms - time_mininute) / 60
    time_hour = time_hs % 24
    time_day = (time_hs - time_hour) / 24
    print "used time (seconds)", time_used
    print time_day, "Days",time_hour, "Hours", time_mininute, "Mininutes", time_second, "Seconds"
    print "over"



