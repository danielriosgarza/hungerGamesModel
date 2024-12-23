# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 12:28:47 2022

@author: danie
"""

from pathlib import Path
import os
import sys

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))

from updateParameters import *
from lmfit import Parameters



def getPramsFromFile(species, filePath):
    lmfit_params = Parameters()    
    with open(filePath) as f:
        f.readline()
        for line in f:
            a = line.strip().split('\t')
            lmfit_params.add(a[1], value = float(a[2]), min = float(a[3]), max = float(a[4]), vary=True)
    
    return lmfit_params
    
    


def assignBhParams(lmfit_params, conn):
    
    
    
    #z1
    
    num = str(lmfit_params['z1_l_s1'].value**lmfit_params['z1_h_s1'].value)
    
    denom = num + " + metObj.metD['trehalose'].concentration**" + str(lmfit_params['z1_h_s1'].value)
    
       
    zeta1 = "(" + num + "/(" + denom + "))"
    
    
    #z2
    
    num = "metObj.metD['trehalose'].concentration**50.0" #+ str(lmfit_params['z1_h_s1'].value)
    
    denom = num + " + 0.21**50.0" #str(lmfit_params['z1_l_s1'].value**lmfit_params['z1_h_s1'].value)
    
       
    zeta2 = "(" + num + ")/(" + denom + ")"
    
    
    
    #z3
    
    
    
    zeta3 = '""'# "(5.75**100)/((metObj.pH**100)+(5.75**100))" #roughly pH<5.5 
    
    #z4
    
    num = str(lmfit_params['z4_l_s3_s4'].value**lmfit_params['z4_h_s3_s4'].value)
    
    
    absD = "((metObj.metD['glucose'].concentration - metObj.metD['glutamate'].concentration)**2)**0.5"
    
    denom = num + " + ((metObj.metD['glucose'].concentration + metObj.metD['glutamate'].concentration - " + absD + ")/2)**"  + str(lmfit_params['z4_h_s3_s4'].value)
    
    
    zeta4 = "(" + num + "/(" + denom + "))"
    
    
    
    #z4
    
    zeta5 = '""'
    
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xa_mumax'].value, lmfit_params['xa_pHopt'].value, lmfit_params['xa_pHalpha'].value, 'xa'))
        
        update_subpopulations(conn, (lmfit_params['xb_mumax'].value, lmfit_params['xb_pHopt'].value, lmfit_params['xb_pHalpha'].value, 'xb'))
        
        
        
        update_subpopulations2subpopulations(conn, (zeta1, lmfit_params['z1_r'].value, 1))
        
        update_subpopulations2subpopulations(conn, (zeta2, 1.5, 2))
        
        update_subpopulations2subpopulations(conn, (zeta3, lmfit_params['z3_r'].value, 3))
        
        update_subpopulations2subpopulations(conn, (zeta4, lmfit_params['z4_r'].value, 4))
        
        update_subpopulations2subpopulations(conn, (zeta5, lmfit_params['z5_r'].value, 5))
        
        
        
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s1'].value, lmfit_params['xa_k_s1'].value, 1))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s5_s1'].value, 0, 2))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s6_s1'].value, 0, 3))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s2'].value, lmfit_params['xa_k_s2'].value, 4))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s6_s2'].value, 0, 5))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s3'].value, lmfit_params['xb_k_s3'].value, 6))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s4'].value, lmfit_params['xb_k_s4'].value, 7))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s6_s3_s4'].value, 0, 8))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s2'].value, lmfit_params['xb_k_s2'].value, 9))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s6_s2'].value, 0, 10))
        
        
        
        
        
        
    
def assignBtParams(lmfit_params, conn):
    
    #z6
    num1 = str(lmfit_params['z6_l_s3'].value**lmfit_params['z6_h_s3'].value)
    denom1 = num1 + " + (metObj.metD['glucose'].concentration**" + str(lmfit_params['z6_h_s3'].value) +")"
    
    num2 = "metObj.metD['mannose'].concentration**" + str(lmfit_params['z6_h_s7'].value)
    denom2 = num2 + " + " + str(lmfit_params['z6_l_s7'].value**lmfit_params['z6_h_s7'].value)
    
       
    zeta6 = "(" + num1 + "/(" + denom1 + ")) * (" + num2 + "/(" + denom2 + "))"
    
    #z7
    
    num = str(lmfit_params['z7_l_s3'].value**lmfit_params['z7_h_s3'].value)
    
    denom = num + " + (metObj.metD['glucose'].concentration**" + str(lmfit_params['z7_h_s3'].value) +")"
    
    zeta7 = "(" + num + "/(" + denom + "))" + "*(5.5**10)/((metObj.pH**10)+(5.5**10))"
    
    
    #z8
    
    num = str(lmfit_params['z8_l_s7'].value**lmfit_params['z8_h_s7'].value)
    
    denom = num + " + (metObj.metD['mannose'].concentration**" + str(lmfit_params['z8_h_s7'].value) +")"
    
    zeta8 = "(" + num + "/(" + denom + "))" + "*(5.5**10)/((metObj.pH**10)+(5.5**10))" #roughly pH<5.5 "*float(metObj.pH<5.5)"
    
    #z9

    zeta9 = '""'
    
    #z10
    num = "(metObj.metD['glucose'].concentration)**10" 
   
    denom = num + " + 0.5**10" 
    
    zeta10 = "(" + num + "/(" + denom + "))"
    
    
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xe_mumax'].value, lmfit_params['xe_pHopt'].value, lmfit_params['xe_pHalpha'].value, 'xe'))
        
        update_subpopulations(conn, (lmfit_params['xf_mumax'].value, lmfit_params['xf_pHopt'].value, lmfit_params['xf_pHalpha'].value, 'xf'))
    
        update_subpopulations2subpopulations(conn, (zeta6, lmfit_params['z6_r'].value, 6))
    
        update_subpopulations2subpopulations(conn, (zeta7, lmfit_params['z7_r'].value, 7))
    
        update_subpopulations2subpopulations(conn, (zeta8, lmfit_params['z8_r'].value, 8))
    
        update_subpopulations2subpopulations(conn, (zeta9, lmfit_params['z9_r'].value, 9))
        
        update_subpopulations2subpopulations(conn, (zeta10, 0.0001, 10))

        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s2'].value, lmfit_params['xe_k_s2'].value, 11))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s6_s2'].value, 0, 12))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s5_s2'].value, 0, 13))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s9_s2'].value, 0, 14))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s3'].value, lmfit_params['xe_k_s3'].value, 15))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s6_s3'].value, 0, 16))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s8_s3'].value, 0, 17))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xf_g_s7'].value, lmfit_params['xf_k_s7'].value, 18))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xf_g_s6_s7'].value, 0, 19))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xf_g_s8_s7'].value, 0, 20))

    
def assignRiParams(lmfit_params, conn):
    
    #z11
    num = "(metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**" + str(lmfit_params['z11_h_s5'].value)
    denom = num + " + " + str(lmfit_params['z11_l_s5'].value**lmfit_params['z11_h_s5'].value)
    
    zeta11 = "(" + num + "/(" + denom + "))"
    
    #z12
    
    num = str(lmfit_params['z12_l_s3_s2'].value**lmfit_params['z12_h_s3_s2'].value)
    denom = num + " + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**" + str(lmfit_params['z12_h_s3_s2'].value)
    
    zeta12 = "(" + num + "/(" + denom + "))"
    
    #z13
    zeta13 = '""'
    
    #z14
    
    zeta14 = '""'
    
    
    
    #z15
    
    num = "(metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**50" #+ str(lmfit_params['z11_h_s3'].value)
    
   
    denom = num + " + 8.0**50" #+ str(lmfit_params['z11_l_s3'].value**lmfit_params['z11_h_s3'].value)
    
    zeta15 = "(" + num + "/(" + denom + "))"
    
    
    #zeta14 = '""'
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xi_mumax'].value, lmfit_params['xi_pHopt'].value, lmfit_params['xi_pHalpha'].value, 'xi'))
        
        update_subpopulations(conn, (lmfit_params['xj_mumax'].value, lmfit_params['xj_pHopt'].value, lmfit_params['xj_pHalpha'].value, 'xj'))
    
        update_subpopulations2subpopulations(conn, (zeta11, lmfit_params['z11_r'].value, 11))
        
        update_subpopulations2subpopulations(conn, (zeta12, lmfit_params['z12_r'].value, 12))
        
        update_subpopulations2subpopulations(conn, (zeta13, lmfit_params['z13_r'].value, 13))
        
        update_subpopulations2subpopulations(conn, (zeta14, lmfit_params['z14_r'].value, 14))
        
        update_subpopulations2subpopulations(conn, (zeta15, 0.1, 15))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s2'].value, lmfit_params['xi_k_s2'].value, 21))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s6_s2'].value, 0, 22))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s10_s2'].value, 0, 23))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s3'].value, lmfit_params['xi_k_s3'].value, 24))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s6_s3'].value, 0, 25))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s5_s3'].value, 0, 26))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s10_s3'].value, 0, 27))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s5'].value, lmfit_params['xj_k_s5'].value, 28))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s10_s5'].value, 0, 29))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s6'].value, lmfit_params['xj_k_s6'].value, 30))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s10_s6'].value, 0, 31))

    

def assignParameters2db(species, lmfit_params, conn):
    if species == 'bh':
        assignBhParams(lmfit_params, conn)
    
    elif species == 'bt':
        assignBtParams(lmfit_params, conn)
        
    elif species == 'ri':
        assignRiParams(lmfit_params, conn)