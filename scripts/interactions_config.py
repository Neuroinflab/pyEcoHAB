# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 12:05:03 2017

@author: JanKMlaptop
"""

def load_experiments_info():
    experiments = []
    with open("experiments_desc.csv", "r") as cmp_file:
        for i,line in enumerate(cmp_file):
            if i==0:
                keys = line.split()[0].split(',')
            else:
                linedata = line.split("\n")[0].split(',')
                #print linedata
                experiments.append({})
                for j, key in enumerate(keys):
                    experiments[-1][key] = linedata[j]
    return experiments

def load_comparisons_info():
    comparisons = {}
    with open("comparisons.csv", "r") as cmp_file:
        for i,line in enumerate(cmp_file):
            if i==0:
                keys = line.split()[0].split(',')
            else:
                linedata = line.split()[0].split(',')
                #print linedata
                for j, key in enumerate(keys):
                    if "-" in linedata[j]:
                        diff = linedata[j]
                cmp_name = diff
                for data in linedata:
                    if data!="" and not "-" in data:
                        cmp_name+="|"+data
                comparisons[cmp_name] = {}
                comparisons[cmp_name]["division"] = {}
                comparisons[cmp_name]["conditions"] = {}
                for j, key in enumerate(keys):
                    if "-" in linedata[j]:
                       comparisons[cmp_name]["division"][key] = diff.split("-")
                    else:
                        if linedata[j]!="":
                            comparisons[cmp_name]["conditions"][key] = linedata[j]
    return comparisons
def group_data(comparison,comparisons,experiments):
    division_param = comparisons[comparison]['division'].keys()[0]
    division_values = comparisons[comparison]['division'].values()[0]
    names = {}
    for i in division_values:
        names[i] = []
    for exp in experiments:
        valid = True
        for key in comparisons[comparison]['conditions'].keys():
            if exp[key]!=comparisons[comparison]['conditions'][key]:
                print exp["path"], exp[key],comparisons[comparison]['conditions'][key]
                valid = False
                break
        if valid and exp[division_param] in division_values:
            names[exp[division_param]].append(exp["path"])
    print names        


                        
                    
                    
