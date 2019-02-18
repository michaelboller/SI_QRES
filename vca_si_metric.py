import os, sys

import csv
import math
import numpy as np
import pandas as pd

import pfnet
import gridopt
import epri.gat as gat


def main():
    # Power Flow Case
    rawFile = r"D:\Box Sync\Box Sync\VCA\RTS_33.raw"
    #rawFile = r"D:\Box Sync\Box Sync\VCA\2383wp_psse.raw"
    #rawFile = r"D:\Box Sync\Box Sync\VCA\2383_gen_fix.raw"
    
    # Output CSV file
    outFile = r"D:\Box Sync\Box Sync\VCA\SI_QRES\si_rts_redo.csv"
    #outFile = r"D:\Box Sync\Box Sync\VCA\SI_QRES\si_2383fix.csv"
    
    # Initalizing gridopt solver
    method = gridopt.power_flow.new_method('ACPF')
    method.set_parameters({'solver': 'nr', 'limit_vars': True, 'quiet': True})
    
    # Parse RAW file to PFNET network
    parser_raw = pfnet.ParserRAW()
    net_o = parser_raw.parse(rawFile)
    
    # Load VCA Solutions
    vcaFile = r"D:\Box Sync\Box Sync\VCA\RTS_33_sklearn-spectral_3_8_500.csv"
    #vcaFile = r"D:\Box Sync\Box Sync\VCA\PolandWinterVCASample.csv"
    #solNum = 5 # Remember zero indexing!
    
    vcaSol = pd.read_csv(vcaFile)
    
    #for solNum in range(1):
    for solNum in range(len(vcaSol.index)):
        # Create a list of buses in each VCAaa
        numVcas = vcaSol['num vcas'][solNum]
        vcaBuses = []
        for vca_k in range(numVcas):
            vcas_k = []
            for bus_i in net_o.buses:
                iNum = bus_i.number
                if vcaSol[str(iNum)][solNum] == vca_k:
                    vcas_k.append(bus_i.index)
            vcaBuses.append(vcas_k)
        
        # Solve original network
        method.solve(net_o)
        method.update_network(net_o)
        
        # Calculate the generator reactive power injections (original)
        Q_o = []
        for gen_i in net_o.generators:
                Q_o.append(gen_i.Q)
        
        # Determining initial paramerters for binary search
        maxStress = 0;
        for gen_i in net_o.generators:
            maxStress = maxStress + gen_i.P_max
        
        lambda_l = 1.0
        mu_l = 1.0
        lambda_g = 1.0
        
        Rmin = []
        R = []
        SI = []
        critBus = []
        
        for k in range(numVcas):
            # Initialize binary search variables
            S_u = maxStress
            S_l = 0
            S_err = 0.35/net_o.base_power
            
            # Find largest contingency (largest generator in VCA)
            max = 0;
            for gen_i in net_o.generators:
                if (vcaSol[str(gen_i.bus.number)][solNum] == k):
                    if gen_i.P_max > max:
                        max = gen_i.P_max
                        contInd = gen_i.index
            
            critBus.append(net_o.get_generator(contInd).bus.number)
            
            # Running binary search to find maximum stress with contingency
            while abs(S_u - S_l) > S_err:
                # Set S to the floored integer (S_u+S_l)/2
                S = (S_u + S_l) / 2;
                
                # Create a new case to run stresses
                net_c = parser_raw.parse(rawFile)
                
                # Update load/generator values based on stress (S)
                for load_i in net_o.loads:
                    if (vcaSol[str(load_i.bus.number)][solNum] == k):
                        iNum = load_i.index
                        cLoad = net_c.get_load(iNum)
                        cLoad.P = load_i.P + lambda_l * S
                        cLoad.Q = load_i.Q + mu_l * S
                
                for gen_i in net_o.generators:
                    if (vcaSol[str(gen_i.bus.number)][solNum] == k):
                        iNum = gen_i.index
                        cGen = net_c.get_generator(iNum)
                        cGen.P = gen_i.P + lambda_g * S
                        if cGen.P > cGen.P_max:
                            cGen.P = cGen.P_max
                        elif cGen.P < cGen.P_min:
                            cGen.P = cGen.P_min
                            
                # Create contingency
                contGen = net_c.get_generator(contInd)
                contGen.outage = True
                
                # Run contingency
                try:
                    method.solve(net_c)
                except:
                    S_u = S - S_err
                    success = 0
                else:   
                    if method.results['solver status'] == 'solved':
                        # If S converges, set S_l to (S_u+S_l)/2+1
                        method.update_network(net_c)
                        S_l = S + S_err
                        success = 1
                    else:
                        # If S does not converge, set S_u to (S_u+S_l)/2-1
                        S_u = S - S_err
                        success = 0
            
            # Continuing binary search while not converged
            while success == 0 and S_u > 0:
                # Set S to the floored integer (S_u+S_l)/2
                S = (S_u + S_l) / 2;
                
                # Create a new case to run stresses
                net_c = parser_raw.parse(rawFile)
                
                # Update load/generator values based on stress (S)
                for load_i in net_o.loads:
                    if (vcaSol[str(load_i.bus.number)][solNum] == k):
                        iNum = load_i.index
                        cLoad = net_c.get_load(iNum)
                        cLoad.P = load_i.P + lambda_l * S
                        cLoad.Q = load_i.Q + mu_l * S
                
                for gen_i in net_o.generators:
                    if (vcaSol[str(gen_i.bus.number)][solNum] == k):
                        iNum = gen_i.index
                        cGen = net_c.get_generator(iNum)
                        cGen.P = gen_i.P + lambda_g * S
                        if cGen.P > cGen.P_max:
                            cGen.P = cGen.P_max
                        elif cGen.P < cGen.P_min:
                            cGen.P = cGen.P_min
                    
                # Create contingency
                contGen = net_c.get_generator(contInd)
                contGen.outage = True
                
                # Run contingency
                try:
                    method.solve(net_c)
                except:
                    S_u = S - S_err
                    success = 0
                else:   
                    if method.results['solver status'] == 'solved':
                        # If S converges, set S_l to (S_u+S_l)/2+1
                        method.update_network(net_c)
                        S_l = S + S_err
                        success = 1
                    else:
                        # If S does not converge, set S_u to (S_u+S_l)/2-1
                        S_u = S - S_err
                        success = 0
                        
            if success == 1:
                # Calculating pre-contingency reactive power production
                # Create a new case to run stresses
                net_c = parser_raw.parse(rawFile)
                
                # Creating pre-contingency load/generator values based on stress (S)
                for load_i in net_o.loads:
                    if (vcaSol[str(load_i.bus.number)][solNum] == k):
                        iNum = load_i.index
                        cLoad = net_c.get_load(iNum)
                        cLoad.P = load_i.P + lambda_l * S
                        cLoad.Q = load_i.Q + mu_l * S
                
                for gen_i in net_o.generators:
                    if (vcaSol[str(gen_i.bus.number)][solNum] == k):
                        iNum = gen_i.index
                        cGen = net_c.get_generator(iNum)
                        cGen.P = gen_i.P + lambda_g * S
                        if cGen.P > cGen.P_max:
                            cGen.P = cGen.P_max
                        elif cGen.P < cGen.P_min:
                            cGen.P = cGen.P_min
                        
                method.solve(net_c)
                method.update_network(net_c)
                
                # Calculate the generator reactive power injections
                Q_pre = []
                for gen_i in net_c.generators:
                    Q_pre.append(gen_i.Q)
                
                # Calculating post-contingency reactive power production
                contGen = net_c.get_generator(contInd)
                contGen.outage = True
                
                method.solve(net_c)
                method.update_network(net_c)
                
                # Calculate the generator reactive power injections
                Q_post = []
                for gen_i in net_c.generators:
                    Q_post.append(gen_i.Q)
                
                # Set the tolerance of change that the generators must have
                cTol = 5.0/net_o.base_power # Previous 20.0
                
                # Sum reactive reserves pre and post contingency to find Rmin
                Rmin_k = 0
                m = 0
                for gen_i in net_o.generators:
                    if (vcaSol[str(gen_i.bus.number)][solNum] == k):
                        if (Q_post[m] - Q_pre[m]) > cTol:
                            Rmin_k = Rmin_k + (Q_post[m] - Q_pre[m])
                    m = m + 1
                Rmin.append(Rmin_k)
                
                # Find reactive reserves (R) for original case
                R_k = 0
                m = 0
                for gen_i in net_o.generators:
                    if (vcaSol[str(gen_i.bus.number)][solNum] == k):
                        if (Q_post[m] - Q_pre[m]) > cTol:
                            R_k = R_k + (Q_post[m] - Q_o[m])
                    m = m + 1
                R.append(R_k)
                
                # Finding overall security index using R and Rmin
                try:
                    SI.append((R[k]-Rmin[k])/Rmin[k])
                except:
                    SI.append('NAN')
            else:
                Rmin.append(0)
                R.append(0)
                SI.append('U')
        
    
        print(SI)
        with open(outFile, 'ab') as f:
            outLines = []
            outLines.append(critBus)
            outLines.append(Rmin)
            outLines.append(R)
            outLines.append(SI)
            writer = csv.writer(f)
            writer.writerows(outLines)
            writer.writerow("")
            
            
            
                
            
                
    
if __name__=="__main__":
	main()
    