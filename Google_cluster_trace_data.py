# -*- coding:utf-8 -*-
from __future__ import division
from scipy import stats
import re
import os
import shutil
import sys
import time
import csv
import json
import pickle
#import xlrd
#import xlwt
import datetime
import math
import numpy as np
from VNF_DPT_source import VNF_request, x, VNF_SPE, CPU_TOTAL, MEM_TOTAL , h_band, brc_base

"""
class VNF_request(object):
    def __init__(self, ):
        self.id = ''
        self.cpu = None
        self.mem = None
        self.in_band = None
        self.out_band = None
        self.pre_vnf_id = None                #two inditor
        self.next_vnf_id = None
"""

def Pre_process():
    Job_ID_dict = {} #jod_ID-->task_num
    #Google_cluster_file = file('e:/Google_Cluster_Trace/task_usage/part-00000-of-00500.csv','rb')
    Google_cluster_file = file('e:/Google_Cluster_Trace/test/part-00000-of-00500.csv','rb')
    Google_reader = csv.reader(Google_cluster_file)
    for G_line in Google_reader:
        start_time = int(G_line[0])
        end_time = int(G_line[1])
        if (start_time >= 600000000) and (start_time <= 610000000) and (end_time > 890000000) and (end_time <= 910000000): #and end_time == '900000000'
            print start_time, end_time
            #print "Job Counting"
            if Job_ID_dict.has_key(G_line[2]):
                Job_ID_dict[G_line[2]] += 1
            else:
                Job_ID_dict[G_line[2]] = 1
        else:
            if start_time > 890000000:
                break
    Task_num_list = []
    Job_ID_list = []
    for key, value in Job_ID_dict.iteritems():
        #if value>=2 and value <= 20:
        if value>=1 and value <= 10:
            Job_ID_list.append(key)
        Task_num_list.append(value)
    Task_num_list.sort()
    Task_1 = 0
    Task_2 = 0
    Task_3 = 0
    Task_4 = 0
    Task_5 = 0
    Task_1_5 = 0
    Task_5_10 = 0
    Task_10_20 = 0
    Task_20_30 = 0
    Task_30_40 = 0
    Task_40_50 = 0
    Task_50_100 = 0
    Task_supr_100 = 0
    for values in Task_num_list:
        if values <= 5:
            if values == 1:
                Task_1 += 1
            elif values == 2:
                Task_2 += 1
            elif values == 3:
                Task_3 += 1
            elif values == 4:
                Task_4 += 1
            else:
                Task_5 += 1
            Task_1_5 += 1
        elif values >5 and values <= 10:
            Task_5_10 += 1
        elif values > 10 and values <= 20:
            Task_10_20 += 1
        elif values > 20 and values <= 30:
            Task_20_30 += 1
        elif values > 30 and values <= 40:
            Task_30_40 += 1
        elif values > 40 and values <= 50:
            Task_40_50 += 1
        elif values > 50 and values <= 100:
            Task_50_100 += 1
        else:
            Task_supr_100 += 1
    
    print "***********************************"
    print "Task num equal 1:", Task_1
    print "Task num equal 2:", Task_2
    print "Task num equal 3:", Task_3
    print "Task num equal 4:", Task_4
    print "Task num equal 5:", Task_5
    print "Task num between 1 and 5:", Task_1_5
    print "Task num between 5 and 10:", Task_5_10
    print "Task num between 10 and 20:", Task_10_20
    print "Task num between 20 and 30:", Task_20_30
    print "Task num between 30 and 40:", Task_30_40
    print "Task num between 40 and 50:", Task_40_50
    print "Task num between 50 and 100:", Task_50_100
    print "Task num superior 100:", Task_supr_100
    print "***********************************"
    
    Job_ID_list.sort()
    print "Job_ID_list is:", Job_ID_list
    return  Job_ID_list  

def Find_ID_in_TL(Task_list, ID):
    for i in range(len(Task_list)):
        if ID == Task_list[i].id:
            return i

    
def Pick_up_max_value_in_dict(Sample_count):
    max_value = 0
    for key, value in Sample_count.iteritems():
        if value >= max_value:
            max_value = value
    return max_value
    

def Google_data_process(path, Job_ID_list):
    Job_ID_dict = {} #jod_ID-->task_ID-->task_usage
    Sample_num = 144 #144/288 144*5min/60 = 12 hours
    Job_num = len(Job_ID_list)
    Finish_flag = {}
    Task_ID_dict = {}
    Sample_count = {}
    Job_count_Sam = 0
    for temp_i in Job_ID_list:
        Finish_flag[temp_i] = 1
    
    
    for t_files in os.listdir(path):
        file_path = os.path.join(path, t_files)
        file_path_sp = os.path.split(file_path)
        lists = file_path_sp[1].split('.')
        file_ext = lists[-1]
        print "File name is:", lists[0]
        if file_ext == 'csv':
            Google_cluster_file = file(file_path,'rb')
            Google_reader = csv.reader(Google_cluster_file)
            for G_line in Google_reader:
                #if Finish_flag:
                    #print "Start Derive Data"
                if G_line[2] in Job_ID_list:
                    if Finish_flag[G_line[2]]:
                        if Job_ID_dict.has_key(G_line[2]):
                            print "Job_ID that has exists:", G_line[2]
                            if Task_ID_dict[G_line[2]].has_key(G_line[3]):
                                #max_count = Pick_up_max_value_in_dict(Sample_count)
                                print "Sample_count num in (if) is:", G_line[2], G_line[3], Sample_count[G_line[2]][G_line[3]]
                                Task_ID_dict[G_line[2]][G_line[3]].cpu[Sample_count[G_line[2]][G_line[3]]] = G_line[5]
                                Task_ID_dict[G_line[2]][G_line[3]].mem[Sample_count[G_line[2]][G_line[3]]] = G_line[7]
                            else:
                                print "Task ID is:", G_line[2], G_line[3]
                                Sample_count[G_line[2]][G_line[3]] = 0
                                VNF_temp = VNF_request()
                                Task_ID_dict[G_line[2]][G_line[3]] = VNF_temp
                                Task_ID_dict[G_line[2]][G_line[3]].id = G_line[3]
                                Task_ID_dict[G_line[2]][G_line[3]].cpu = np.array([float(0) for x in range(Sample_num)])
                                Task_ID_dict[G_line[2]][G_line[3]].mem = np.array([float(0) for x in range(Sample_num)])
                                Task_ID_dict[G_line[2]][G_line[3]].cpu[Sample_count[G_line[2]][G_line[3]]] = G_line[5]
                                Task_ID_dict[G_line[2]][G_line[3]].mem[Sample_count[G_line[2]][G_line[3]]] = G_line[7]
                                Job_ID_dict[G_line[2]].append(Task_ID_dict[G_line[2]][G_line[3]])
                        else:
                            Job_ID_dict[G_line[2]] = []
                            Task_ID_dict[G_line[2]] = {}
                            Sample_count[G_line[2]] = {}
                            print "Job_ID is:", G_line[2]
                            if Task_ID_dict[G_line[2]].has_key(G_line[3]):
                                #max_count = Pick_up_max_value_in_dict(Sample_count)
                                print "Sample_count num in (if) is:", Sample_count[G_line[2]][G_line[3]]
                                Task_ID_dict[G_line[2]][G_line[3]].cpu[Sample_count[G_line[2]][G_line[3]]] = G_line[5]
                                Task_ID_dict[G_line[2]][G_line[3]].mem[Sample_count[G_line[2]][G_line[3]]] = G_line[7]
                            else:
                                print "Task ID is:", G_line[2], G_line[3]
                                Sample_count[G_line[2]][G_line[3]] = 0
                                VNF_temp = VNF_request()
                                Task_ID_dict[G_line[2]][G_line[3]] = VNF_temp
                                #Task_ID_dict[G_line[2]][G_line[3]] = VNF_request()
                                Task_ID_dict[G_line[2]][G_line[3]].id = G_line[3]
                                Task_ID_dict[G_line[2]][G_line[3]].cpu = np.array([float(0) for x in range(Sample_num)])
                                Task_ID_dict[G_line[2]][G_line[3]].mem = np.array([float(0) for x in range(Sample_num)])
                                Task_ID_dict[G_line[2]][G_line[3]].cpu[Sample_count[G_line[2]][G_line[3]]] = G_line[5]
                                Task_ID_dict[G_line[2]][G_line[3]].mem[Sample_count[G_line[2]][G_line[3]]] = G_line[7]
                                Job_ID_dict[G_line[2]].append(Task_ID_dict[G_line[2]][G_line[3]])
                            
                        Sample_count[G_line[2]][G_line[3]] += 1
                        if Sample_count[G_line[2]][G_line[3]] == Sample_num:
                            print "Job", G_line[2], "Counts", Sample_num, "times."
                            Finish_flag[G_line[2]] = 0
                            Job_count_Sam += 1
                    else:
                        continue
        #if Finish_flag == 0:
        #    break
    
    Job_ID_dict_2 = {}
    i = 0
    j = 0
    for key_2, value_2 in Job_ID_dict.iteritems():
        delta_Sam = 1
        for m in range(len(value_2)):
            if (Sample_num - Sample_count[key_2][value_2[m].id]) >= 5:
                delta_Sam = 0
        if delta_Sam:
            if len(value_2) <= 20:
                Job_ID_dict_2[j] = value_2
                j += 1
        i += 1
    
    print "Job Num. ori is", i
    print "Job that has reach", Sample_num,  "num is:", Job_count_Sam
    print "Job Num. validated:", j
    for key_3, Task_list in Job_ID_dict_2.iteritems():
        l = 0
        if len(Task_list) == 1: #Job only has one task, so the pre_vnf_id & next_vnf_id are all none
            Task_list[l].in_band = 0.5*Task_list[l].cpu + 0.5*Task_list[l].mem
            Task_list[l].out_band = 0.5*Task_list[l].cpu + 0.5*Task_list[l].mem
        else:
            for l in range(len(Task_list)):
                if l == 0:
                    Task_list[l].next_vnf_id = Task_list[l+1].id
                    Task_list[l].in_band = 0.5*Task_list[l].cpu + 0.5*Task_list[l].mem
                    Task_list[l].out_band = 0.5*Task_list[l].cpu + 0.5*Task_list[l].mem
                elif l == (len(Task_list) - 1):
                    Task_list[l].pre_vnf_id = Task_list[l-1].id
                    Task_list[l].in_band = Task_list[l-1].out_band
                    Task_list[l].out_band = 0.5*Task_list[l].cpu + 0.5*Task_list[l].mem
                else:
                    Task_list[l].next_vnf_id = Task_list[l+1].id
                    Task_list[l].pre_vnf_id = Task_list[l-1].id
                    Task_list[l].in_band = Task_list[l-1].out_band
                    Task_list[l].out_band = 0.5*Task_list[l].cpu + 0.5*Task_list[l].mem
                        
    with open('E:/Google_Cluster_Trace/'+'Google_Cluster_task_usage_144.pickle', 'w') as f:
        pickle.dump(Job_ID_dict_2,f)
    #with open('E:/Google_Cluster_Trace/'+'Google_Cluster_task_usage_3.json', 'w') as f:
        #json.dump(Job_ID_dict,f)
            


#Google_cluster_file = file('e:/Google_Cluster_Trace/task_usage/part-00000-of-00500.csv','rb')

#Google_reader = csv.reader(Google_cluster_file)

#for G_line in Google_reader:
    #print G_line[0]

if __name__ == '__main__':
    time_start = time.time()
    Job_ID_list = Pre_process()
    print "******************Pre_process is Over****************"
    print "****************Start Collect Data*******************"
    #path = 'E:/Google_Cluster_Trace/task_usage'
    path = 'E:/Google_Cluster_Trace/test_2'
    Google_data_process(path, Job_ID_list)
    """
    with open('E:/Google_Cluster_Trace/'+'Google_Cluster_task_usage_4.pickle', 'r') as f:
        data_back = pickle.load(f)
    for i in range(1,2):
        for j in range(len(data_back[i])):
            print "Task ID:", data_back[i][j].id, "Pre_task_ID:", data_back[i][j].pre_vnf_id, "Next_task_ID:", data_back[i][j].next_vnf_id
            print "CPU Usage:", data_back[i][j].cpu
            print "Mem Usage:", data_back[i][j].mem
            print "In_band Usage:", data_back[i][j].in_band
            print "Out_band Usage:", data_back[i][j].out_band
    
    with open('E:/Google_Cluster_Trace/'+'Google_Cluster_task_usage_2.pickle', 'r') as f:
        data_back = pickle.load(f)       
    print "*********************Next is the 33 sample***************"
    for i in range(1,2):
        for j in range(len(data_back[i])):
            print "Task ID:", data_back[i][j].id, "Pre_task_ID:", data_back[i][j].pre_vnf_id, "Next_task_ID:", data_back[i][j].next_vnf_id
            print "CPU Usage:", data_back[i][j].cpu
            print "Mem Usage:", data_back[i][j].mem
            print "In_band Usage:", data_back[i][j].in_band
            print "Out_band Usage:", data_back[i][j].out_band    
    """
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
