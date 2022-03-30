##!/usr/bin/env python

###################################################
##      File name: mainQMP.py
##         Author: Cheng ZHANG
##          Email: cheng.zhang.abbott@vc.ibaraki.ac.jp
##           Date: 2021-12-01
###################################################

# import modules used here -- sys is a very standard one
import sys
#import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import utility
from operator import itemgetter, attrgetter
import random
import time
import os
import math
from scipy.special import lambertw

nCount = 0
nT = 0

def params_def_init():
    # This function read parameters from parameters from parameters.ini file
    # and define the corresponding global variables
    utility.log_info("Enter params_def_init()")
    parmfilename = os.path.join(utility.get_current_dir(),"parameters.ini")
    if os.path.isfile(parmfilename) == False:
        log_error("parameters.ini file is not exist!")
        return False
    utility.log_info("------Parameters from ini file------")
    parmfile = open(parmfilename,"r")
    line = parmfile.readline()
    while line:
        stripeline = line.strip()
        if ( stripeline[0:1] == "#" ):
        #skip the comment line with "#" in the head of the line
            line = parmfile.readline()
            continue
        key_value = stripeline.split("#")[0].split("=")
        if ( len(key_value) >= 2 ):
            globals()[key_value[0].strip()] = float(key_value[1].strip())
            #define the global variables, the value is asummed as float type
            utility.log_info("\t\t"+key_value[0].strip()+"\t=\t"+key_value[1].strip())
        line = parmfile.readline()
    parmfile.close()
    utility.log_info("------------End------------")


#Function name: isActionValid
#         [in]: l@int, location
#      return]: @bool, True if valid, False if not valid
#  Description: judge whether action is valid or not
#       Author: Cheng ZHANG

#Function name: def monetary_cost(l,b,a)
#         [in]: l@int, location
#         [in]: b@int, remaining data size
#         [in]: l@int, MU action
#      return]: @float
#  Description: return the monetary cost with location l and action a
#       Author: Cheng ZHANG

class BAP(object):
    """Class for budget allocation problem

    Attributes
    ----------
    B : int
        total budget
    L : int
        The number of locations
    """
    # The relation ship between quality (y) and budget (x) with different beta
    # beta = 1: y = 0.138 x + 18607
    # beta = 2: y = 0.1387 x + 50935
    # beta = 3: y = 0.1389 x + 89080
    # beta = 4: y = 0.1391 x + 131086
    # beta = 5: y = 0.1392 x + 175998
    # beta = 6: y = 0.1392 x + 223248
    # beta = 7: y = 0.1394 x + 272455
    # beta = 8: y = 0.1395 x + 323347
    # beta = 9: y = 0.1395 x + 375717
    # beta = 10: y = 0.1396 x + 429407
    '''
    budget: 10000.0
    i = 1
    added
    quality = 19987.0
    total
    quality = 19987.0
    i = 2
    added
    quality = 52322.0
    total
    quality = 72309.0
    i = 3
    added
    quality = 90469.0
    total
    quality = 162778.0
    i = 4
    added
    quality = 132477.0
    total
    quality = 295255.0
    i = 5
    added
    quality = 177390.0
    total
    quality = 472645.0
    i = 6
    added
    quality = 224640.0
    total
    quality = 697285.0
    i = 7
    added
    quality = 273849.0
    total
    quality = 971134.0
    i = 8
    added
    quality = 324742.0
    total
    quality = 1295876.0
    i = 9
    added
    quality = 377112.0
    total
    quality = 1672988.0
    i = 10
    added
    quality = 430803.0
    total
    quality = 2103791.0
    '''
    def __init__(self):
        print ("__init__():class BAP")

        self._TB                = int(g_TB) #total budget
        self._L                 = int(g_L_N) # number of locations
#        self._listBudget        = [100,120,140,160,180,200,220,240,260,280,300]#None
#        self._listBudget        = [100,100,100,100,100,100,100,100,100,100]#None
        #self._listQuality       = [100,120,139,158,177,196,215,234,253,272,391]#None
        self._listBudget        = [1818,3636,5454,7272,9090,10909,12727,14545,16363,18181]#None
        self._listQuality = [18857, 51439, 89837, 132097, 177263, 224766, 274229, 325376, 377999, 431945]  # None
    # beta = 1: y = 0.138 x + 18607
    # beta = 2: y = 0.1387 x + 50935
    # beta = 3: y = 0.1389 x + 89080
    # beta = 4: y = 0.1391 x + 131086
    # beta = 5: y = 0.1392 x + 175998
    # beta = 6: y = 0.1392 x + 223248
    # beta = 7: y = 0.1394 x + 272455
    # beta = 8: y = 0.1395 x + 323347
    # beta = 9: y = 0.1395 x + 375717
    # beta = 10: y = 0.1396 x + 429407
    def getQualityByLocationAndBudget(self,l,budget):
        quality = 0.0
        if l == 1:
            quality = 0.138 * budget + 18607
        elif l == 2:
            quality = 0.1387 * budget + 50935
        elif l == 3:
            quality = 0.1389 * budget + 89080
        elif l == 4:
            quality = 0.1391* budget + 131086
        elif l == 5:
            quality = 0.1392* budget + 175998
        elif l == 6:
            quality = 0.1392* budget + 223248
        elif l == 7:
            quality = 0.1394* budget + 272455
        elif l == 8:
            quality = 0.1395* budget + 323347
        elif l == 9:
            quality = 0.1395* budget + 375717
        elif l == 10:
            quality = 0.1396* budget + 429407

        return quality

    #set location weight for all the locations by uniform distribution
    def dynamicProgramingCore(self):
        dp_table = [[0 for j in range(self._TB + 1)] for i in range(len(self._listBudget) + 1)]
        for i in range(1, len(dp_table)):
            for j in range(1, len(dp_table[0])):
                if self._listBudget[i - 1] <= j:
                    dp_table[i][j] = max(self._listQuality[i - 1] + dp_table[i][j - self._listBudget[i - 1]], dp_table[i - 1][j])
                elif self._listBudget[i - 1] > j:
                    dp_table[i][j] = dp_table[i - 1][j]
        print("dp allocation:"+str(dp_table[-1][-1]))
#        print(dp_table[-1][-1])
        #print(dp_table)
        return dp_table[-1][-1]

    def dynamicProgramingCoreNew(self):
        dp_table = [[0 for j in range(self._TB + 1)] for i in range(len(self._listBudget) + 1)]
        for i in range(1, len(dp_table)):
            for j in range(1, len(dp_table[0])):
                if self._listBudget[i - 1] <= j:
                    dp_table[i][j] = max(self.getQualityByLocationAndBudget(i-1,self._listBudget[i-1]) + dp_table[i][j - self._listBudget[i - 1]], dp_table[i - 1][j])
                elif self._listBudget[i - 1] > j:
                    dp_table[i][j] = dp_table[i - 1][j]
        print(dp_table[-1][-1])
        for i in range(1, len(dp_table)):
            print ("i="+str(i))
            print(dp_table[i][-1])
        #print(len(dp_table))
        return dp_table[-1][-1]

    def averageAllocationNew(self):
        print("total budget:"+ str(self._TB))
        budget = self._TB / self._L
        #print ("budget:"+ str(budget))
        totalq=0
        for i in range(1,self._L+1,1):
            qualityadded = self.getQualityByLocationAndBudget(i,budget)
            totalq = totalq + qualityadded
            print ("i="+str(i)+" added quality=" + str(qualityadded)+ " total quality="+str(totalq))
        print("average allocation:"+str(totalq))
        print (self.getQualityByLocationAndBudget(10,self._TB))

    def weightBasedAllocation(self):
        #budget = self._TB / self._L
        weight_sum = 0
        for i in range(1, self._L + 1, 1):
            weight_sum = weight_sum + i
        #print ("weight sum:"+str(weight_sum))

        totalq=0
        for i in range(1,self._L+1,1):
            budget = self._TB * i / weight_sum
            #print("budget:" + str(budget))
            qualityadded = self.getQualityByLocationAndBudget(i,self._TB * i / weight_sum)
            totalq = totalq + qualityadded
            #print ("i="+str(i)+" added quality=" + str(qualityadded)+ " total quality="+str(totalq))
        print("weighted allocation:"+str(totalq))

    def averageAllocation(self):
        budget = self._TB / self._L
        totalq=0
        for i in range(self._L):
            totalq = totalq + budget * self._listQuality[i] / self._listBudget[i]
        print(totalq)

class QMP(object):

    """A Markov Decision Problem.
    Attributes
    ----------
    A : int
        The number of actions
    S : int
        The number of states
    """
    def __init__(self):
        print ("__init__():class QMP")

        self._N                 = int(g_N)
        #self._L                 = int(g_L_N)
        self._location_weight   = 8
        self._reputation        = None
        self._cost_parameters   = None

    #set location weight for all the locations by uniform distribution
    def setLocWeight(self, location_weight):
        self._location_weight = location_weight
        print (self._location_weight)

    #set reputation parameter for all users by uniform distribution
    def setRepuByUniform(self):
        self._reputation = np.random.randint(g_gamma_low, g_gamma_high, self._N)
        print (self._reputation)

    #set cost parameter for all users by uniform distribution
    def setCostByUniform(self):
        self._cost_parameters = np.random.randint(g_lambda_low, g_lambda_high, self._N)
        print (self._cost_parameters)

    def setRepuByNorm(self, mu, sigma):
        for i in range(self._N):
            self._reputation[i]= np.random.gauss(mu,sigma)
            print (reputaion)
        return

    def locationWeight(self):
        return self._location_weight

    def reputation(self):
        return self._reputation

    def reputationByIndex(self, i):
        if i < 0 and i >= self._N:
            utility.log_error("the index of reputation is out of range")
            return
        else:
            return self._reputation[i]

    def cost_parameters(self):
        return self._cost_parameters

    def cost_parametersByIndex(self, i):
        if i < 0 and i >= self._N:
            utility.log_error("the index of cost parameters is out of range")
            return
        else:
            return self._cost_parameters[i]

    def quality_of_data(self,i,effort):
        if i < 0 and i >= self._N:
            utility.log_error("the index is out of range")
            return
        else:
            return self._reputation[i]* self._location_weight * np.log(1+effort)

    def cost_of_data(self,i,effort):
        if i < 0 and i >= self._N:
            utility.log_error("the index is out of range")
            return
        else:
            return self._cost_parameters[i]* effort

    def effort_from_p(self,i,p):
        return self._reputation[i] * self._location_weight * p/self._cost_parameters[i] - 1

    def utility_DU(self,i, pi, effort):
        if i < 0 and i >= self._N:
            utility.log_error("the index is out of range")
            return
        else:
            return self._reputation[i]* self._location_weight * self.quality_of_data(i,effort) - self.cost_of_data(i,effort)

    def optimal_effort(self, i, pi):
        if i < 0 and i >= self._N:
            utility.log_error("the index is out of range")
            return
        else:
            return (self._reputation[i] * self._location_weight * pi) / self._cost_parameters[i] - 1

    def W(self, i, alpha):
        if i < 0 and i >= self._N:
            utility.log_error("the index is out of range")
            return
        else:
            return lambertw((self._reputation[i] * self._location_weight * np.e) / (self._cost_parameters[i] * alpha))

    def G(self, alpha):
        sumg = 0
        sumquality = 0
        sumpayment = 0
        opt_price_list = []
        opt_payment_list = []
        opt_effort_list = []
        opt_quality_list = []
        #utility.create_output_file("alpha_price_list.csv")
        #utility.save_to_output(str(alpha))
        if alpha == 0.0:
            return
        for i in range(self._N):
            w= self.W(i,alpha).real
            #first part of function G
            part1 = self._reputation[i]*self._location_weight *(w-1) / (alpha * (w + 1))
            #second part of function G
            part2 = (self._reputation[i]* self._location_weight / (alpha * w * (1+w))) * np.log(self._reputation[i]*self._location_weight/(self._cost_parameters[i]*alpha*w))
            sumg = sumg + part1 + part2
            price = 1/(alpha*w)
            effort = self.effort_from_p(i,price) # user i' effort
            quality = self.quality_of_data(i, effort) # user i's quality
            payment = price*quality # payment to user i
            sumquality = sumquality + quality # sum quality of all users
            sumpayment = sumpayment + payment # sum payment to all users

            opt_price_list.append(price)
            opt_payment_list.append(payment)
            opt_effort_list.append(effort)
            opt_quality_list.append(quality)
            #print ("price is:" + str(price))
            #print ("effort is:" + str(effort))
            DU=self.utility_DU(i,price,effort)
            if DU < 0:
                print ("User's utlity:"+ str(DU))
        #utility.save_to_output(",".join(str(r) for r in opt_price_list))
        #print ("When alpha is:" + str(alpha) + " G / B is:" + str(sumg))
        #print ("When alpha is:" + str(alpha) + " Total quality is:" + str(sumquality))
        #print("When alpha is:" + str(alpha) + " Total payment is:" + str(sumpayment))
        return sumpayment,sumquality, opt_price_list, opt_payment_list,opt_effort_list,opt_quality_list

    def G2(self, alpha):
        sumg = 0
        sumquality = 0
        sumpayment = 0
        opt_price_list = []
        opt_payment_list = []
        opt_effort_list = []
        opt_quality_list = []
        #utility.create_output_file("alpha_price_list.csv")
        #utility.save_to_output(str(alpha))
        if alpha == 0.0:
            return
        for i in range(self._N):
            w= self.W(i,alpha).real
            #first part of function G
            part1 = self._reputation[i]*self._location_weight *(w-1) / (alpha * (w + 1))
            #second part of function G
            part2 = (self._reputation[i]* self._location_weight / (alpha * w * (1+w))) * np.log(self._reputation[i]*self._location_weight/(self._cost_parameters[i]*alpha*w))
            sumg = sumg + part1 + part2
            price = 1/(alpha*w)
            effort = self.effort_from_p(i,price) # user i' effort
            quality = self.quality_of_data(i, effort) # user i's quality
            payment = price*quality # payment to user i
            sumquality = sumquality + quality # sum quality of all users
            sumpayment = sumpayment + payment # sum payment to all users

            opt_price_list.append(price)
            opt_payment_list.append(payment)
            opt_effort_list.append(effort)
            opt_quality_list.append(quality)
            #print ("price is:" + str(price))
            #print ("effort is:" + str(effort))
            DU=self.utility_DU(i,price,effort)
            if DU < 0:
                print ("User's utlity:"+ str(DU))
        #utility.save_to_output(",".join(str(r) for r in opt_price_list))
        #print ("When alpha is:" + str(alpha) + " G / B is:" + str(sumg))
        #print ("When alpha is:" + str(alpha) + " Total quality is:" + str(sumquality))
        #print("When alpha is:" + str(alpha) + " Total payment is:" + str(sumpayment))
        #return sumpayment,sumquality, opt_price_list, opt_payment_list,opt_effort_list,opt_quality_list
        return sumg,sumquality

    def price(self,i, alpha):
        return 1 / (alpha * self.W(i,alpha))

    def fix_price_sum(self, p):
        sumquality = 0
        sumbudget  = 0
        list_effort = []
        list_quality = []
        list_payment = []

        for i in range(self._N):
            effort = self.effort_from_p(i,p)
            list_effort.append(effort)
            if self.utility_DU(i,p,effort) >0:
                quality = self.quality_of_data(i,effort)
                list_quality.append(quality)
                budget = quality * p
                list_payment.append(budget)
                sumquality += quality
                sumbudget  += budget

        return sumbudget, sumquality, list_payment, list_effort,list_quality

    def run(self):
        print (self.reputation())
#        print (self.reputationByIndex(9))
        print (self._cost_parameters)
#        print (self._reputation)
 #       print ("sum reputation:", sum(self._reputation))
  #      print ("sum cost:", sum(self._cost_parameters))
        self.G(0.15)
        for p in np.arange(0,5.0,0.1):
            print (p)
            print (self.fix_price_sum(p))
        #for alpha in np.arange(0.1,0.26,0.01):
            #print (self.price(a))
        #    self.G(alpha)

    def optimal_price(self):
        filename = "optimal_price_alpha" + str(self._location_weight) + ".csv"
        utility.create_output_file(filename)
        utility.save_to_output("alpha, budget,total_quality")
        alpha_list = []

        for alpha in np.arange(0.05,0.46,0.005):
            alpha_list.append(alpha)

            result_G = self.G2(alpha)
            utility.save_to_output(str(alpha)+","+str(result_G[0])+","+str(result_G[1]))
        utility.close_output_file()

    def optimal_price_diff_location(self):
        for locweight in np.arange(1, 11, 1):
            self.setLocWeight(locweight)
            self.optimal_price()

    def sum_quality_diff_budget(self):

        utility.create_output_file("sum_quality_diff_budget.csv")

#--------for fixed price--------------------
        result_fix_p    = []
        result_sum_qulaity1    = []
        result_sum_budget1    = []

        #for p in np.arange(0,2.1,0.1):
        utility.save_to_output("price, payment, effort, quality")
        for p in np.arange(0.1,5.1,0.01):
            result_fix_p.append(p)
            result_sum = self.fix_price_sum(p)
            result_sum_budget1.append(result_sum[0])
            result_sum_qulaity1.append(result_sum[1])
            utility.save_to_output("price:,"+str(p))
            utility.save_to_output(",".join(str(r) for r in result_sum[2]))
            utility.save_to_output(",".join(str(r) for r in result_sum[3]))
            utility.save_to_output(",".join(str(r) for r in result_sum[4]))

        print (result_fix_p)
        print (result_sum_budget1)
        print (result_sum_qulaity1)

        utility.save_to_output("fixed price, sum budget, sum quality")
        utility.save_to_output(",".join(str(r) for r in result_fix_p))
        utility.save_to_output(",".join(str(r) for r in result_sum_budget1))
        utility.save_to_output(",".join(str(r) for r in result_sum_qulaity1))

        # --------for optimal price--------------------
        alpha_list = []
        result_sum_budget2  = []
        result_sum_qulaity2 = []

        for alpha in np.arange(0.05,0.46,0.005):
            alpha_list.append(alpha)

            result_G = self.G(alpha)
            result_sum_budget2.append(result_G[0])
            result_sum_qulaity2.append(result_G[1])
            utility.save_to_output("alpha="+str(alpha)+" and price, payment, effort, quality")
            utility.save_to_output(",".join(str(r) for r in result_G[2]))
            utility.save_to_output(",".join(str(r) for r in result_G[3]))
            utility.save_to_output(",".join(str(r) for r in result_G[4]))
            utility.save_to_output(",".join(str(r) for r in result_G[5]))

            #print (alpha)
           # self.G(alpha)

        print(alpha_list)
        print(result_sum_budget2)
        print(result_sum_qulaity2)

        utility.save_to_output("optimal price, sum budget, sum quality")

        utility.save_to_output(",".join(str(r) for r in alpha_list))
        utility.save_to_output(",".join(str(r) for r in result_sum_budget2))
        utility.save_to_output(",".join(str(r) for r in result_sum_qulaity2))

        utility.close_output_file()

        return

def main():
#    global g_switch_of_known_consumption, g_switch_of_unknown_consumption
    print ('\n=============================================================')
    print ('\n\tThis is simulation program for paper "Incentivizing Mobile Users for Crowdsensing on IoT  Platform"')
    print ('\n\tAuthor: Cheng ZHANG')
    print ('\n\t  Date: 2021/12/01')
    print ('\n=============================================================')
    utility.init_log_file("log.log")
    utility.log_info("Enter function: "+main.__name__)
##
    if ( params_def_init() is False ):
        utility.log_error("Fail to initialize parameters, abort the program")
        print ("[ERROR] Fail to initialize parameters, abort the program, see log.log file for detail information.")
        utility.close_log_file()
        return

    utility.create_output_file()
    utility.save_to_output("#flows,"+"baseline cost,"+"finished,"+"Heuristic cost,"+"finished")

    #qmp = QMP()
    #qmp.setCostByUniform()
    #qmp.setRepuByUniform()

    #qmp.optimal_price()
    #qmp.sum_quality_diff_budget()
    #qmp.optimal_price_diff_location()

    bap = BAP()
    #bap.dynamicProgramingCoreNew()
    bap.averageAllocationNew()
    bap.weightBasedAllocation()
    bap.dynamicProgramingCore()
    #

if __name__ == '__main__':
    main()

