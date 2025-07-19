# -*- coding: utf-8 -*-
"""
Created on Thu May 29 11:34:11 2025

@author: drgarza
"""

def getTransitionFunction(statement):
    
    def tf(metObj):
        if statement == '""':
            return 1.0
        return eval(statement)
    return tf
    


def update_bh_Params(parameter, new_value, reactor):
    
    #1
    if parameter == 'xa_mumax':
        reactor.microbiome.subpopD['xa'].mumax = new_value
        return reactor
        

    #2
    elif parameter == 'xb_mumax':
        reactor.microbiome.subpopD['xb'].mumax = new_value
        return reactor
    
    #3
    elif parameter == 'xa_pHopt':
        reactor.microbiome.subpopD['xa'].pHopt = new_value
        return reactor
    
    #4
    elif parameter == 'xb_pHopt':
        reactor.microbiome.subpopD['xb'].pHopt = new_value
        return reactor
    
    #5
    elif parameter == 'xa_pHalpha':
        reactor.microbiome.subpopD['xa'].pHalpha = new_value
        return reactor
    
    #6
    elif parameter == 'xb_pHalpha':
        reactor.microbiome.subpopD['xb'].pHalpha = new_value
        return reactor
    
    #7
    elif parameter == 'xa_k_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('trehalose')] = new_value
        return reactor
    
    #8
    elif parameter == 'xa_k_s2':
        reactor.microbiome.subpopD['xa'].feedingTerms[1].monodKs[reactor.microbiome.subpopD['xa'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #9
    elif parameter == 'xb_k_s3':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glucose')] = new_value
        return reactor
    
    #10
    elif parameter == 'xb_k_s4':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glutamate')] = new_value
        return reactor
    
    #11
    elif parameter == 'xb_k_s2':
        reactor.microbiome.subpopD['xb'].feedingTerms[1].monodKs[reactor.microbiome.subpopD['xb'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #12
    elif parameter == 'xa_g_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].yields[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('trehalose')] = new_value
        return reactor
    
    #13
    elif parameter == 'xa_g_s2':
        reactor.microbiome.subpopD['xa'].feedingTerms[1].yields[reactor.microbiome.subpopD['xa'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #14
    elif parameter == 'xa_g_s6_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].yields[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('acetate')] = new_value
        return reactor
    
    #15
    elif parameter == 'xa_g_s5_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].yields[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('lactate')] = new_value
        return reactor
    
    #16
    elif parameter == 'xa_g_s6_s2':
        reactor.microbiome.subpopD['xa'].feedingTerms[1].yields[reactor.microbiome.subpopD['xa'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    #17
    elif parameter == 'xb_g_s3':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].yields[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glucose')] = new_value
        return reactor
    
    #18
    elif parameter == 'xb_g_s4':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].yields[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glutamate')] = new_value
        return reactor
    
    #19
    elif parameter == 'xb_g_s2':
        reactor.microbiome.subpopD['xb'].feedingTerms[1].yields[reactor.microbiome.subpopD['xb'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #20
    elif parameter == 'xb_g_s6_s3_s4':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].yields[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('acetate')] = new_value
        return reactor
    
    #21
    elif parameter == 'xb_g_s6_s2':
        reactor.microbiome.subpopD['xb'].feedingTerms[1].yields[reactor.microbiome.subpopD['xb'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    
    #22
    elif parameter == 'z1_r':
        reactor.microbiome.bacteria['bh'].connections['xa'][0] = ('xb', getTransitionFunction(f"(4.572019894399926e-08**1/(4.572019894399926e-08 + metObj.metD['trehalose'].concentration**1.0))"), new_value)
        return reactor
    
    #23
    elif parameter == 'z2_r':
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = ('xa', getTransitionFunction(f"(metObj.metD['trehalose'].concentration**50.0)/(metObj.metD['trehalose'].concentration**50.0 + 0.21**50.0)"), new_value)
        return reactor
    
    #24
    elif parameter == 'z3_r':
        reactor.microbiome.bacteria['bh'].connections['xa'][1] = ('xc', getTransitionFunction('""'), new_value)
        return reactor
    
    #25
    elif parameter == 'z4_r':
        reactor.microbiome.bacteria['bh'].connections['xb'][1] = ('xc', getTransitionFunction(f"(0.0001293346**2.9055456987/(0.0001293346**2.9055456987 + ((metObj.metD['glucose'].concentration + metObj.metD['glutamate'].concentration - ((metObj.metD['glucose'].concentration - metObj.metD['glutamate'].concentration)**2)**0.5)/2)**2.905545698700455))"), new_value)
        return reactor
    
    #26
    elif parameter == 'z5_r':
        reactor.microbiome.bacteria['bh'].connections['xc'][0] = ('xd', getTransitionFunction('""'), new_value)
        return reactor
    
    #27
    elif parameter == 'z1_l_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xa'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xa'][0] = (current[0], getTransitionFunction(f"({new_value}**1/({new_value} + metObj.metD['trehalose'].concentration**1.0))"), current[2])
        return reactor
    
    #28
    elif parameter == 'z2_l_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = (current[0], getTransitionFunction(f"(metObj.metD['trehalose'].concentration**50.0)/(metObj.metD['trehalose'].concentration**50.0 + {new_value}**50.0)"), current[2])
        return reactor
    
    #29
    elif parameter == 'z4_l_s3_s4':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][1][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][1] = (current[0], getTransitionFunction(f"({new_value}**2.9055456987/({new_value}**2.9055456987 + ((metObj.metD['glucose'].concentration + metObj.metD['glutamate'].concentration - ((metObj.metD['glucose'].concentration - metObj.metD['glutamate'].concentration)**2)**0.5)/2)**2.905545698700455))"), current[2])
        return reactor
    
    #30
    elif parameter == 'z1_h_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xa'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xa'][0] = (current[0], getTransitionFunction(f"(4.572019894399926e-08**{new_value}/(4.572019894399926e-08 + metObj.metD['trehalose'].concentration**{new_value}))"), current[2])
        return reactor
    
    #31
    elif parameter == 'z2_h_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = (current[0], getTransitionFunction(f"(metObj.metD['trehalose'].concentration**{new_value})/(metObj.metD['trehalose'].concentration**{new_value} + 0.21**{new_value})"), current[2])
        return reactor
    
    #32
    elif parameter == 'z4_h_s3_s4':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = (current[0], getTransitionFunction(f"(metObj.metD['trehalose'].concentration**{new_value})/(metObj.metD['trehalose'].concentration**{new_value} + 0.21**{new_value})"), current[2])
        return reactor
   


def update_bt_Params(parameter, new_value, reactor):
    
    #33
    if parameter == 'xe_mumax':
        reactor.microbiome.subpopD['xe'].mumax = new_value
        return reactor
    
    #34
    if parameter == 'xf_mumax':
        reactor.microbiome.subpopD['xf'].mumax = new_value
        return reactor
    
    #35
    if parameter == 'xe_pHopt':
        reactor.microbiome.subpopD['xe'].pHopt = new_value
        return reactor
    
    #36
    if parameter == 'xf_pHopt':
        reactor.microbiome.subpopD['xf'].pHopt = new_value
        return reactor
    
    #37
    if parameter == 'xe_pHalpha':
        reactor.microbiome.subpopD['xe'].pHalpha = new_value
        return reactor
    
    #38
    if parameter == 'xf_pHalpha':
        reactor.microbiome.subpopD['xf'].pHalpha = new_value
        return reactor
    
    #39
    if parameter == 'xe_k_s3':
        reactor.microbiome.subpopD['xe'].feedingTerms[1].monodKs[reactor.microbiome.subpopD['xe'].feedingTerms[1].metIDs.index('glucose')] = new_value
        return reactor
    
    #40
    if parameter == 'xe_k_s2':
        reactor.microbiome.subpopD['xe'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xe'].feedingTerms[0].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #41
    if parameter == 'xf_k_s7':
        reactor.microbiome.subpopD['xf'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xf'].feedingTerms[0].metIDs.index('mannose')] = new_value
        return reactor
    
    #42
    if parameter == 'xe_g_s2':
        reactor.microbiome.subpopD['xe'].feedingTerms[0].yields[reactor.microbiome.subpopD['xe'].feedingTerms[0].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #43
    if parameter == 'xe_g_s3':
        reactor.microbiome.subpopD['xe'].feedingTerms[1].yields[reactor.microbiome.subpopD['xe'].feedingTerms[1].metIDs.index('glucose')] = new_value
        return reactor
    
    #44
    if parameter == 'xe_g_s5_s2':
        reactor.microbiome.subpopD['xe'].feedingTerms[0].yields[reactor.microbiome.subpopD['xe'].feedingTerms[0].metIDs.index('lactate')] = new_value
        return reactor
    
    #45
    if parameter == 'xe_g_s6_s2':
        reactor.microbiome.subpopD['xe'].feedingTerms[0].yields[reactor.microbiome.subpopD['xe'].feedingTerms[0].metIDs.index('acetate')] = new_value
        return reactor
    
    #46
    if parameter == 'xe_g_s6_s3':
        reactor.microbiome.subpopD['xe'].feedingTerms[1].yields[reactor.microbiome.subpopD['xe'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    #47
    if parameter == 'xe_g_s8_s3':
        reactor.microbiome.subpopD['xe'].feedingTerms[1].yields[reactor.microbiome.subpopD['xe'].feedingTerms[1].metIDs.index('succinate')] = new_value
        return reactor
    
    #48
    if parameter == 'xe_g_s9_s2':
        reactor.microbiome.subpopD['xe'].feedingTerms[0].yields[reactor.microbiome.subpopD['xe'].feedingTerms[0].metIDs.index('formate')] = new_value
        return reactor
    
    #49
    if parameter == 'xf_g_s7':
        reactor.microbiome.subpopD['xf'].feedingTerms[0].yields[reactor.microbiome.subpopD['xf'].feedingTerms[0].metIDs.index('mannose')] = new_value
        return reactor
    
    #50
    if parameter == 'xf_g_s6_s7':
        reactor.microbiome.subpopD['xf'].feedingTerms[0].yields[reactor.microbiome.subpopD['xf'].feedingTerms[0].metIDs.index('acetate')] = new_value
        return reactor
    
    #51
    if parameter == 'xf_g_s8_s7':
        reactor.microbiome.subpopD['xf'].feedingTerms[0].yields[reactor.microbiome.subpopD['xf'].feedingTerms[0].metIDs.index('succinate')] = new_value
        return reactor
    
    #52
    if parameter == 'z6_r':
        reactor.microbiome.bacteria['bt'].connections['xe'][0] = ('xf', getTransitionFunction(f"(0.0099187739**1.0006644305/(0.0099187739**1.0006644305 + (metObj.metD['glucose'].concentration**1.0006644305))) * (metObj.metD['mannose'].concentration**1.0022092044/(metObj.metD['mannose'].concentration**1.0022092044 + 0.50994613706**1.0022092044))"), new_value)
        return reactor
    
    #53
    if parameter == 'z7_r':
        reactor.microbiome.bacteria['bt'].connections['xe'][1] = ('xg', getTransitionFunction(f"(0.0065072254**1.4599162310/(0.0065072254**1.4599162310 + (metObj.metD['glucose'].concentration**1.4599162310)))*(5.5**10)/((metObj.pH**10)+(5.5**10))"), new_value)
        return reactor
    
    #54
    if parameter == 'z8_r':
        reactor.microbiome.bacteria['bt'].connections['xf'][0] = ('xg', getTransitionFunction(f"(0.0000007162**29.999999996778822/(0.0000007162**29.999999996778822 + (metObj.metD['mannose'].concentration**29.9999999967)))*(5.5**10)/((metObj.pH**10)+(5.5**10))"), new_value)
        return reactor
    
    #55
    if parameter == 'z9_r':
        reactor.microbiome.bacteria['bt'].connections['xg'][0] = ('xh', getTransitionFunction('""'), new_value)
        return reactor
    
    #56
    if parameter == 'z10_r':
        reactor.microbiome.bacteria['bt'].connections['xf'][1] = ('xe', getTransitionFunction(f"((metObj.metD['glucose'].concentration)**10/((metObj.metD['glucose'].concentration)**10 + 0.5**10))"), new_value)
        return reactor
    
    #57
    if parameter == 'z6_l_s3':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][0] = (current[0], getTransitionFunction(f"({new_value}**1.0006644305/({new_value}**1.0006644305 + (metObj.metD['glucose'].concentration**1.0006644305))) * (metObj.metD['mannose'].concentration**1.0022092044/(metObj.metD['mannose'].concentration**1.0022092044 + 0.50994613706**1.0022092044))"), current[2])
        return reactor
    
    #58
    if parameter == 'z6_l_s7':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][0] = (current[0], getTransitionFunction(f"({new_value}**1.0006644305/({new_value}**1.0006644305 + (metObj.metD['glucose'].concentration**1.0006644305))) * (metObj.metD['mannose'].concentration**1.0022092044/(metObj.metD['mannose'].concentration**1.0022092044 + {new_value}**1.0022092044))"), current[2])
        return reactor
    
    #59
    if parameter == 'z7_l_s3':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][1][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][1] = (current[0], getTransitionFunction(f"({new_value}**1.4599162310/({new_value}**1.4599162310 + (metObj.metD['glucose'].concentration**1.4599162310)))*(5.5**10)/((metObj.pH**10)+(5.5**10))"), current[2])
        return reactor
    
    
    #60
    if parameter == 'z7_l_pH':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][1][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][1] = (current[0], getTransitionFunction(f"(0.0065072254**1.4599162310/(0.0065072254**1.4599162310 + (metObj.metD['glucose'].concentration**1.4599162310)))*({new_value}**10)/((metObj.pH**10)+({new_value}**10))"), current[2])
        return reactor
    
    #61
    if parameter == 'z8_l_s7':
        current = reactor.microbiome.bacteria['bt'].connections['xf'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xf'][0] = (current[0], getTransitionFunction(f"({new_value}**29.999999996778822/({new_value}**29.999999996778822 + (metObj.metD['mannose'].concentration**29.9999999967)))*(5.5**10)/((metObj.pH**10)+(5.5**10))"), current[2])
        return reactor
    
    #62
    if parameter == 'z8_l_pH':
        current = reactor.microbiome.bacteria['bt'].connections['xf'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xf'][0] = (current[0], getTransitionFunction(f"(0.0000007162**29.999999996778822/(0.0000007162**29.999999996778822 + (metObj.metD['mannose'].concentration**29.9999999967)))*({new_value}**10)/((metObj.pH**10)+({new_value}**10))"), current[2])
        return reactor
    
    #63
    if parameter == 'z10_l_s3':
        current = reactor.microbiome.bacteria['bt'].connections['xf'][1][:]
        reactor.microbiome.bacteria['bt'].connections['xf'][1] = (current[0], getTransitionFunction(f"((metObj.metD['glucose'].concentration)**10/((metObj.metD['glucose'].concentration)**10 + {new_value}**10))"), current[2])
        return reactor
    
    
    #64
    if parameter == 'z6_h_s3':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][0] = (current[0], getTransitionFunction(f"(0.0099187739**{new_value}/(0.0099187739**{new_value} + (metObj.metD['glucose'].concentration**{new_value}))) * (metObj.metD['mannose'].concentration**1.0022092044/(metObj.metD['mannose'].concentration**1.0022092044 + 0.50994613706**1.0022092044))"), current[2])
        return reactor
    
    #65
    if parameter == 'z6_h_s7':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][0] = (current[0], getTransitionFunction(f"(0.0099187739**1.0006644305/(0.0099187739**1.0006644305 + (metObj.metD['glucose'].concentration**1.0006644305))) * (metObj.metD['mannose'].concentration**{new_value}/(metObj.metD['mannose'].concentration**{new_value} + 0.50994613706**{new_value}))"), current[2])
        return reactor
    
    #66
    if parameter == 'z7_h_s3':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][1][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][1] = (current[0], getTransitionFunction(f"(0.0065072254**{new_value}/(0.0065072254**{new_value} + (metObj.metD['glucose'].concentration**{new_value})))*(5.5**10)/((metObj.pH**10)+(5.5**10))"), current[2])
        return reactor
    
    #67
    if parameter == 'z7_h_pH':
        current = reactor.microbiome.bacteria['bt'].connections['xe'][1][:]
        reactor.microbiome.bacteria['bt'].connections['xe'][1] = (current[0], getTransitionFunction(f"(0.0065072254**1.4599162310/(0.0065072254**1.4599162310 + (metObj.metD['glucose'].concentration**1.4599162310)))*(5.5**{new_value})/((metObj.pH**{new_value})+(5.5**{new_value}))"), current[2])
        return reactor
    
    #68
    if parameter == 'z8_h_s7':
        current = reactor.microbiome.bacteria['bt'].connections['xf'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xf'][0] = (current[0], getTransitionFunction(f"(0.0000007162**{new_value}/(0.0000007162**{new_value} + (metObj.metD['mannose'].concentration**{new_value})))*(5.5**10)/((metObj.pH**10)+(5.5**10))"), current[2])
        return reactor
    
    #69
    if parameter == 'z8_h_pH':
        current = reactor.microbiome.bacteria['bt'].connections['xf'][0][:]
        reactor.microbiome.bacteria['bt'].connections['xf'][0] = (current[0], getTransitionFunction(f"(0.0000007162**29.999999996778822/(0.0000007162**29.999999996778822 + (metObj.metD['mannose'].concentration**29.9999999967)))*(5.5**{new_value})/((metObj.pH**{new_value})+(5.5**{new_value}))"), current[2])
        return reactor
    
    #70
    if parameter == 'z10_h_s3':
        current = reactor.microbiome.bacteria['bt'].connections['xf'][1][:]
        reactor.microbiome.bacteria['bt'].connections['xf'][1] = (current[0], getTransitionFunction(f"((metObj.metD['glucose'].concentration)**{new_value}/((metObj.metD['glucose'].concentration)**{new_value} + 0.5**{new_value}))"), current[2])
        return reactor


def update_ri_Params(parameter, new_value, reactor):
    
    #71
    if parameter == 'xi_mumax':
        reactor.microbiome.subpopD['xi'].mumax = new_value
        return reactor
    
    #72
    if parameter == 'xj_mumax':
        reactor.microbiome.subpopD['xj'].mumax = new_value
        return reactor
        
    #73
    if parameter == 'xi_pHopt':
        reactor.microbiome.subpopD['xi'].pHopt = new_value
        return reactor
    
    #74
    if parameter == 'xj_pHopt':
        reactor.microbiome.subpopD['xj'].pHopt = new_value
        return reactor
    
    #75
    if parameter == 'xi_pHalpha':
        reactor.microbiome.subpopD['xi'].pHalpha = new_value
        return reactor
    
    #76
    if parameter == 'xj_pHalpha':
        reactor.microbiome.subpopD['xj'].pHalpha = new_value
        return reactor
    
    #77
    if parameter == 'xi_k_s2':
        reactor.microbiome.subpopD['xi'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xi'].feedingTerms[0].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #78
    if parameter == 'xi_k_s3':
        reactor.microbiome.subpopD['xi'].feedingTerms[1].monodKs[reactor.microbiome.subpopD['xi'].feedingTerms[1].metIDs.index('glucose')] = new_value
        return reactor
    
    #79
    if parameter == 'xj_k_s5':
        reactor.microbiome.subpopD['xj'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xj'].feedingTerms[0].metIDs.index('lactate')] = new_value
        return reactor
    
    #80
    if parameter == 'xj_k_s6':
        reactor.microbiome.subpopD['xj'].feedingTerms[1].monodKs[reactor.microbiome.subpopD['xj'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    #81
    if parameter == 'xi_g_s2':
        reactor.microbiome.subpopD['xi'].feedingTerms[0].yields[reactor.microbiome.subpopD['xi'].feedingTerms[0].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #82
    if parameter == 'xi_g_s3':
        reactor.microbiome.subpopD['xi'].feedingTerms[1].yields[reactor.microbiome.subpopD['xi'].feedingTerms[1].metIDs.index('glucose')] = new_value
        return reactor
    
    
    #83
    if parameter == 'xj_g_s5':
        reactor.microbiome.subpopD['xj'].feedingTerms[0].yields[reactor.microbiome.subpopD['xj'].feedingTerms[0].metIDs.index('lactate')] = new_value
        return reactor
    
    
    #84
    if parameter == 'xj_g_s6':
        reactor.microbiome.subpopD['xj'].feedingTerms[1].yields[reactor.microbiome.subpopD['xj'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    
    #85
    if parameter == 'xi_g_s5_s3':
        reactor.microbiome.subpopD['xi'].feedingTerms[1].yields[reactor.microbiome.subpopD['xi'].feedingTerms[1].metIDs.index('lactate')] = new_value
        return reactor
    
    
    #86
    if parameter == 'xi_g_s6_s2':
        reactor.microbiome.subpopD['xi'].feedingTerms[0].yields[reactor.microbiome.subpopD['xi'].feedingTerms[0].metIDs.index('acetate')] = new_value
        return reactor
    
    #87
    if parameter == 'xi_g_s6_s3':
        reactor.microbiome.subpopD['xi'].feedingTerms[1].yields[reactor.microbiome.subpopD['xi'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    #88
    if parameter == 'xi_g_s10_s2':
        reactor.microbiome.subpopD['xi'].feedingTerms[0].yields[reactor.microbiome.subpopD['xi'].feedingTerms[0].metIDs.index('butyrate')] = new_value
        return reactor
    
    #89
    if parameter == 'xi_g_s10_s3':
        reactor.microbiome.subpopD['xi'].feedingTerms[1].yields[reactor.microbiome.subpopD['xi'].feedingTerms[1].metIDs.index('butyrate')] = new_value
        return reactor
    
    
    #90
    if parameter == 'xj_g_s10_s5':
        reactor.microbiome.subpopD['xj'].feedingTerms[0].yields[reactor.microbiome.subpopD['xj'].feedingTerms[0].metIDs.index('butyrate')] = new_value
        return reactor
    
    #91
    if parameter == 'xj_g_s10_s6':
        reactor.microbiome.subpopD['xj'].feedingTerms[1].yields[reactor.microbiome.subpopD['xj'].feedingTerms[1].metIDs.index('butyrate')] = new_value
        return reactor
    
    #92
    if parameter == 'z11_r':
        reactor.microbiome.bacteria['ri'].connections['xi'][0] = ('xj', getTransitionFunction(f"((metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**29.9999999958/((metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**29.9999999958 + 2.3311915035**29.9999999958))"), new_value)
        return reactor
    
    #93
    if parameter == 'z12_r':
        reactor.microbiome.bacteria['ri'].connections['xi'][1] = ('xl', getTransitionFunction(f"(0.0034909802**1.1123778083/(0.0034909802**1.1123778083 + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**1.1123778083))"), new_value)
        return reactor
    
    #94
    if parameter == 'z13_r':
        reactor.microbiome.bacteria['ri'].connections['xj'][0] = ('xk', getTransitionFunction('""'), new_value)
        return reactor
    
    #95
    if parameter == 'z14_r':
        reactor.microbiome.bacteria['ri'].connections['xk'][0] = ('xl', getTransitionFunction('""'), new_value)
        return reactor
    
    #96
    if parameter == 'z15_r':
        reactor.microbiome.bacteria['ri'].connections['xj'][1] = ('xi', getTransitionFunction(f"((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**50/((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**50 + 8.0**50))"), new_value)
        return reactor
    
    
    #97
    if parameter == 'z11_l_s5_s6':
        current = reactor.microbiome.bacteria['ri'].connections['xi'][0][:]
        reactor.microbiome.bacteria['ri'].connections['xi'][0] = (current[0], getTransitionFunction(f"((metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**29.9999999958/((metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**29.9999999958 + {new_value}**29.9999999958))"), current[2])
        return reactor
    
    #98
    if parameter == 'z12_l_s3_s2':
        current = reactor.microbiome.bacteria['ri'].connections['xi'][1][:]
        reactor.microbiome.bacteria['ri'].connections['xi'][1] = (current[0], getTransitionFunction(f"({new_value}**1.1123778083/({new_value}**1.1123778083 + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**1.1123778083))"), current[2])
        return reactor
    
    
    #99
    if parameter == 'z15_l_s3_s2':
        current = reactor.microbiome.bacteria['ri'].connections['xj'][1][:]
        reactor.microbiome.bacteria['ri'].connections['xj'][1] = (current[0], getTransitionFunction(f"((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**50/((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**50 + {new_value}**50))"), current[2])
        return reactor
    
    #100
    if parameter == 'z11_h_s5_s6':
        current = reactor.microbiome.bacteria['ri'].connections['xi'][0][:]
        reactor.microbiome.bacteria['ri'].connections['xi'][0] = (current[0], getTransitionFunction(f"((metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**{new_value}/((metObj.metD['lactate'].concentration+metObj.metD['acetate'].concentration)**{new_value} + 2.3311915035**{new_value}))"), current[2])
        return reactor
    
    #101
    if parameter == 'z12_h_s3_s2':
        current = reactor.microbiome.bacteria['ri'].connections['xi'][1][:]
        reactor.microbiome.bacteria['ri'].connections['xi'][1] = (current[0], getTransitionFunction(f"(0.0034909802**{new_value}/(0.0034909802**{new_value} + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**{new_value}))"), current[2])
        return reactor
    
    #102
    if parameter == 'z15_h_s3_s2':
        current = reactor.microbiome.bacteria['ri'].connections['xj'][1][:]
        reactor.microbiome.bacteria['ri'].connections['xj'][1] = ('xi', getTransitionFunction(f"((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**{new_value}/((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration)**{new_value} + 8.0**{new_value}))"), current[2])
        return reactor


