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
            if a[0] == species:
                lmfit_params.add(a[1], value = float(a[2]), min = float(a[3]), max = float(a[4]), vary=True)
    
    return lmfit_params
    
    


def assignBhParams(lmfit_params, conn):
    
    
    
    #z1
    
    num1 = str(lmfit_params['z1_l_s1'].value**lmfit_params['z1_h_s1'].value)
    denom1 = num1 + " + metObj.metD['trehalose'].concentration**" + str(lmfit_params['z1_h_s1'].value)
    
    num2 = "metObj.metD['glucose'].concentration**" + str(lmfit_params['z1_h_s3'].value)
    denom2 = num2 + " + " + str(lmfit_params['z1_l_s3'].value**lmfit_params['z1_h_s3'].value)
    
       
    zeta1 = "(" + num1 + "/(" + denom1 + ")) * (" + num2 + "/(" + denom2 + "))"
    #zeta1 = "metObj.metD['trehalose'].concentration < " +  str(lmfit_params['z1_l_s1'].value)
    
    
    
    #z2
    zeta2 = '""'
    
    
    #z3
    
    num3 = str(lmfit_params['z3_l_s3_s4'].value**lmfit_params['z3_h_s3_s4'].value)
    
    
    absD = "((metObj.metD['glucose'].concentration - metObj.metD['glutamate'].concentration)**2)**0.5"
    denom3 = num3 + " + (metObj.metD['glucose'].concentration + metObj.metD['glutamate'].concentration + " + absD + ")/2"
    
    
    zeta3 = "(" + num3 + "/(" + denom3 + "))"
    
    
    
    #z4
    
    zeta4 = '""'
    
    num1 = str(lmfit_params['z14_l_s1'].value**lmfit_params['z14_h_s1'].value)
    denom1 = num1 + " + (metObj.metD['trehalose'].concentration + metObj.metD['trehalose'].concentration)**" + str(lmfit_params['z14_h_s1'].value)
    
    
    #z5
    #zeta14 = "(" + num1 + "/(" + denom1 + "))" 
    zeta5 = '""'
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xa_mumax'].value, lmfit_params['xa_pHopt'].value, lmfit_params['xa_pHalpha'].value, 'xa'))
        
        update_subpopulations(conn, (lmfit_params['xb_mumax'].value, lmfit_params['xb_pHopt'].value, lmfit_params['xb_pHalpha'].value, 'xb'))
        
        
        
        update_subpopulations2subpopulations(conn, (zeta1, lmfit_params['z1_r'].value, 1))
        
        update_subpopulations2subpopulations(conn, (zeta2, lmfit_params['z2_r'].value, 2))
        
        update_subpopulations2subpopulations(conn, (zeta3, lmfit_params['z3_r'].value, 3))
        
        update_subpopulations2subpopulations(conn, (zeta4, lmfit_params['z4_r'].value, 4))
        
        update_subpopulations2subpopulations(conn, (zeta5, lmfit_params['z5_r'].value, 14))
        
        
        
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s1'].value, lmfit_params['xa_k_s1'].value, 1))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s5_s1'].value, 0, 2))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s6_s1'].value, 0, 3))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s2'].value, lmfit_params['xa_k_s2'].value, 4))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xa_g_s6_s2'].value, 0, 5))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s3_s4'].value, lmfit_params['xb_k_s3'].value, 6))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s3_s4'].value, lmfit_params['xb_k_s4'].value, 7))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xb_g_s6_s3_s4'].value, 0, 8))
        
        
        
        
    
def assignBtParams(lmfit_params, conn):
    
    #z6
    num1 = str(lmfit_params['z6_l_s3'].value**lmfit_params['z6_h_s3'].value)
    denom1 = num1 + " + (metObj.metD['glucose'].concentration**" + str(lmfit_params['z6_h_s3'].value) +")"
    
    num2 = "metObj.metD['mannose'].concentration**" + str(lmfit_params['z6_h_s7'].value)
    denom2 = num2 + " + " + str(lmfit_params['z6_l_s7'].value**lmfit_params['z6_h_s7'].value)
    
       
    zeta6 = "(" + num1 + "/(" + denom1 + ")) * (" + num2 + "/(" + denom2 + "))"
    
    #z7
    
    # num1 = "metObj.metD['glucose'].concentration**" + str(lmfit_params['z6_h_s3'].value)
    # denom1 = num1 + " + " + str(lmfit_params['z6_l_s3'].value**lmfit_params['z6_h_s3'].value)
    
    # num2 = str(lmfit_params['z6_l_s7'].value**lmfit_params['z6_h_s7'].value)
    # denom2 = num2 + " + metObj.metD['mannose'].concentration**" + str(lmfit_params['z6_h_s7'].value)
    
    
    
       
    # zeta6 = "(" + num1 + "/(" + denom1 + ")) * (" + num2 + "/(" + denom2 + "))"
    zeta7 = 0
    
    
    #z8
    
    num1 = str(lmfit_params['z8_l_s3'].value**lmfit_params['z8_h_s3'].value)
    denom1 = num1 + " + metObj.metD['glucose'].concentration**" + str(lmfit_params['z8_h_s3'].value)
    
    zeta8 = "(" + num1 + "/(" + denom1 + "))"
    
    #z9
    
    num1 = str(lmfit_params['z9_l_s7'].value**lmfit_params['z9_h_s7'].value)
    denom1 = num1 + " + metObj.metD['mannose'].concentration**" + str(lmfit_params['z9_h_s7'].value)
    
    zeta9 = "(" + num1 + "/(" + denom1 + "))"
    
    
    #z10
    
    zeta10 = '""'
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xe_mumax'].value, lmfit_params['xe_pHopt'].value, lmfit_params['xe_pHalpha'].value, 'xe'))
        
        update_subpopulations(conn, (lmfit_params['xf_mumax'].value, lmfit_params['xf_pHopt'].value, lmfit_params['xf_pHalpha'].value, 'xf'))
    
        update_subpopulations2subpopulations(conn, (zeta5, lmfit_params['z6_r'].value, 6))
    
        update_subpopulations2subpopulations(conn, (zeta6, lmfit_params['z7_r'].value, 7))
    
        update_subpopulations2subpopulations(conn, (zeta7, lmfit_params['z8_r'].value, 8))
    
        update_subpopulations2subpopulations(conn, (zeta8, lmfit_params['z9_r'].value, 9))
    
        update_subpopulations2subpopulations(conn, (zeta9, lmfit_params['z10_r'].value, 10))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s2'].value, lmfit_params['xe_k_s2'].value, 9))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s6_s2'].value, 0, 10))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s5_s2'].value, 0, 11))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s9_s2'].value, 0, 12))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s3'].value, lmfit_params['xe_k_s3'].value, 13))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s6_s3'].value, 0, 14))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xe_g_s8_s3'].value, 0, 15))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xf_g_s7'].value, lmfit_params['xf_k_s7'].value, 16))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xf_g_s6_s7'].value, 0, 17))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xf_g_s8_s7'].value, 0, 18))

    
def assignRiParams(lmfit_params, conn):
    
    #z10
    num1 = "(metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**" + str(lmfit_params['z10_h_s5'].value)
    denom1 = num1 + " + " + str(lmfit_params['z10_l_s5'].value**lmfit_params['z10_h_s5'].value)
    
    zeta10 = "(" + num1 + "/(" + denom1 + "))"
    
    #z11
    
    num1 = str(lmfit_params['z11_l_s3'].value**lmfit_params['z11_h_s3'].value)
    denom1 = num1 + " + metObj.metD['glucose'].concentration**" + str(lmfit_params['z11_h_s3'].value)
    
    zeta11 = "(" + num1 + "/(" + denom1 + "))"
    
    #z12
    zeta12 = '""'
    
    #z13
    
    zeta13 = '""'
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xi_mumax'].value, lmfit_params['xi_pHopt'].value, lmfit_params['xi_pHalpha'].value, 'xi'))
        
        update_subpopulations(conn, (lmfit_params['xj_mumax'].value, lmfit_params['xj_pHopt'].value, lmfit_params['xj_pHalpha'].value, 'xj'))
    
        update_subpopulations2subpopulations(conn, (zeta10, lmfit_params['z10_r'].value, 10))
        
        update_subpopulations2subpopulations(conn, (zeta11, lmfit_params['z11_r'].value, 11))
        
        update_subpopulations2subpopulations(conn, (zeta12, lmfit_params['z12_r'].value, 12))
        
        update_subpopulations2subpopulations(conn, (zeta13, lmfit_params['z13_r'].value, 13))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s2'].value, lmfit_params['xi_k_s2'].value, 19))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s6_s2'].value, 0, 20))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s10_s2'].value, 0, 21))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s3'].value, lmfit_params['xi_k_s3'].value, 22))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s6_s3'].value, 0, 23))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s5_s3'].value, 0, 24))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xi_g_s10_s3'].value, 0, 25))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s5'].value, lmfit_params['xj_k_s5'].value, 26))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s10_s5'].value, 0, 27))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s6'].value, lmfit_params['xj_k_s6'].value, 28))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['xj_g_s10_s6'].value, 0, 29))

    

def assignParameters2db(species, lmfit_params, conn):
    if species == 'bh':
        assignBhParams(lmfit_params, conn)
    
    elif species == 'bt':
        assignBtParams(lmfit_params, conn)
        
    elif species == 'ri':
        assignRiParams(lmfit_params, conn)