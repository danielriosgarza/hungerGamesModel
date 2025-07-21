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
    '''
    file expected: #, species, parameter, value, min, max
    '''
    lmfit_params = Parameters()
    with open(filePath) as f:
        f.readline()
        for line in f:
            a = line.strip().split('\t')
            if a[0] == species:
              lmfit_params.add(a[1], value = float(a[2]), min = float(a[3]), max = float(a[4]), vary=True)

    return lmfit_params
    
    


def assignBhParams(lmfit_params, conn):


    #transition functions


    #z1

    z1_l_s1 = str(lmfit_params['z1_l_s1'].value)
    z1_h_s1 = " ** " + str(lmfit_params['z1_h_s1'].value)
    z1 = "(" + z1_l_s1 + z1_h_s1 + "/(" + z1_l_s1 + z1_h_s1 + " + metObj.metD['trehalose'].concentration" + z1_h_s1 + "))"



    #z2
    z2_l_s1 = str(lmfit_params['z2_l_s1'].value)
    z2_h_s1 = " ** " + str(lmfit_params['z2_h_s1'].value)
    z2 = "(metObj.metD['trehalose'].concentration" + z2_h_s1 + ")/(metObj.metD['trehalose'].concentration" + z2_h_s1 + " + " + z2_l_s1 + z2_h_s1 + ")"


    #z3

    z3 = '""'

    #z4

    z4_l_s3_s4 = str(lmfit_params['z4_l_s3_s4'].value)
    z4_h_s3_s4 = " ** " + str(lmfit_params['z4_h_s3_s4'].value)
    z4 = "(" + z4_l_s3_s4 + z4_h_s3_s4 + "/(" + z4_l_s3_s4 + z4_h_s3_s4 + " + ((metObj.metD['glucose'].concentration + metObj.metD['glutamate'].concentration - ((metObj.metD['glucose'].concentration - metObj.metD['glutamate'].concentration)**2)**0.5)/2)" + z4_h_s3_s4 + "))"

    z5 = '""'
    with conn:


        update_subpopulations(conn, (lmfit_params['xa_mumax'].value, lmfit_params['xa_pHopt'].value, lmfit_params['xa_pHalpha'].value, 'xa'))

        update_subpopulations(conn, (lmfit_params['xb_mumax'].value, lmfit_params['xb_pHopt'].value, lmfit_params['xb_pHalpha'].value, 'xb'))



        update_subpopulations2subpopulations(conn, (z1, lmfit_params['z1_r'].value, 1))

        update_subpopulations2subpopulations(conn, (z2, lmfit_params['z2_r'].value, 2))

        update_subpopulations2subpopulations(conn, (z3, lmfit_params['z3_r'].value, 3))

        update_subpopulations2subpopulations(conn, (z4, lmfit_params['z4_r'].value, 4))

        update_subpopulations2subpopulations(conn, (z5, lmfit_params['z5_r'].value, 5))




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
    
    z6_l_s3 = str(lmfit_params['z6_l_s3'].value)
    z6_h_s3 = " ** " + str(lmfit_params['z6_h_s3'].value)
    z6_l_s7 = str(lmfit_params['z6_l_s7'].value)
    z6_h_s7 = " ** " + str(lmfit_params['z6_h_s7'].value)
    z6 = "(" + z6_l_s3 + z6_h_s3 + "/(" + z6_l_s3 + z6_h_s3 + " + (metObj.metD['glucose'].concentration" + z6_h_s3 + "))) * (metObj.metD['mannose'].concentration" + z6_h_s7 + "/(" + z6_l_s7 + z6_h_s7 + " + metObj.metD['mannose'].concentration" + z6_h_s7 + "))"

    #z7
    
    z7_l_s3 = str(lmfit_params['z7_l_s3'].value)
    z7_h_s3 = " ** " + str(lmfit_params['z7_h_s3'].value)
    z7_l_pH = str(lmfit_params['z7_l_pH'].value)
    z7_h_pH = " ** " + str(lmfit_params['z7_h_pH'].value)
    
    z7 = "(" + z7_l_s3 + z7_h_s3 + "/(" + z7_l_s3 + z7_h_s3 + " + (metObj.metD['glucose'].concentration" + z7_h_s3+ "))) * (" + z7_l_pH + z7_h_pH + "/(" + z7_l_pH + z7_h_pH + " + metObj.pH" + z7_h_pH + "))"

    
    
    #z8
    
    z8_l_s7 = str(lmfit_params['z8_l_s7'].value)
    z8_h_s7 = " ** " + str(lmfit_params['z8_h_s7'].value)
    z8_l_pH = str(lmfit_params['z8_l_pH'].value)
    
    z8_h_pH = " ** " + str(lmfit_params['z8_h_pH'].value)
    
    z8 = "(" + z8_l_s7 + z8_h_s7 + "/(" + z8_l_s7 + z8_h_s7 + " + (metObj.metD['mannose'].concentration" + z8_h_s7+ "))) * (" + z8_l_pH + z8_h_pH + "/(" + z8_l_pH + z8_h_pH + " + metObj.pH" + z8_h_pH + "))"
    
    
    #z9
    
    z9 = '""'
    
    
    #z10
    
    z10_l_s3 = str(lmfit_params['z10_l_s3'].value)
    z10_h_s3 = " ** " + str(lmfit_params['z10_h_s3'].value)
    z10 = "((metObj.metD['glucose'].concentration)" + z10_h_s3 + "/((metObj.metD['glucose'].concentration)" + z10_h_s3 + " + " + z10_l_s3 + z10_h_s3 + "))"
    
    
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xe_mumax'].value, lmfit_params['xe_pHopt'].value, lmfit_params['xe_pHalpha'].value, 'xe'))
        
        update_subpopulations(conn, (lmfit_params['xf_mumax'].value, lmfit_params['xf_pHopt'].value, lmfit_params['xf_pHalpha'].value, 'xf'))
    
        update_subpopulations2subpopulations(conn, (z6, lmfit_params['z6_r'].value, 6))
    
        update_subpopulations2subpopulations(conn, (z7, lmfit_params['z7_r'].value, 7))
    
        update_subpopulations2subpopulations(conn, (z8, lmfit_params['z8_r'].value, 8))
    
        update_subpopulations2subpopulations(conn, (z9, lmfit_params['z9_r'].value, 9))
        
        update_subpopulations2subpopulations(conn, (z10, lmfit_params['z10_r'].value, 10))

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
    
    z11_l_s5_s6 = str(lmfit_params['z11_l_s5_s6'].value)
    z11_h_s5_s6 = " ** " + str(lmfit_params['z11_h_s5_s6'].value)
    z11 = "((metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)" + z11_h_s5_s6 + "/(" + z11_l_s5_s6 + z11_h_s5_s6 + " + (metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)" + z11_h_s5_s6 + "))"

    
    #z12
    z12_l_s3_s2 = str(lmfit_params['z12_l_s3_s2'].value)
    z12_h_s3_s2 = " ** " + str(lmfit_params['z12_h_s3_s2'].value)
    z12 = "(" + z12_l_s3_s2 + z12_h_s3_s2 + "/(" + z12_l_s3_s2 + z12_h_s3_s2 + " + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)" + z12_h_s3_s2 + "))"
    
    
    #z13
    
    z13 = '""'
    
    #z14
    
    z14 = '""'
    
    #z15
    
    z15_l_s3_s2 = str(lmfit_params['z15_l_s3_s2'].value)
    z15_h_s3_s2 = " ** " + str(lmfit_params['z15_h_s3_s2'].value)
    z15 = "((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)" + z15_h_s3_s2 + "/((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) " + z15_h_s3_s2 + " + " + z15_l_s3_s2 + z15_h_s3_s2 + "))"

    
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['xi_mumax'].value, lmfit_params['xi_pHopt'].value, lmfit_params['xi_pHalpha'].value, 'xi'))
        
        update_subpopulations(conn, (lmfit_params['xj_mumax'].value, lmfit_params['xj_pHopt'].value, lmfit_params['xj_pHalpha'].value, 'xj'))
    
        update_subpopulations2subpopulations(conn, (z11, lmfit_params['z11_r'].value, 11))
        
        update_subpopulations2subpopulations(conn, (z12, lmfit_params['z12_r'].value, 12))
        
        update_subpopulations2subpopulations(conn, (z13, lmfit_params['z13_r'].value, 13))
        
        update_subpopulations2subpopulations(conn, (z14, lmfit_params['z14_r'].value, 14))
        
        update_subpopulations2subpopulations(conn, (z15, lmfit_params['z15_r'].value, 15))
        
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