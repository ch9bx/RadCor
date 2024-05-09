# import matplotlib
# matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import math as m
import numpy as np

energyStr = "10380"
Ebeam = 10.38
theta_values = [11, 18, 30]
theta = theta_values[0]
thetaStr = str(theta)
deg2rad = 0.0174533
Mp = 0.938
degree_sign = u'\N{DEGREE SIGN}'
unpolLabel = "_unpol"
unpolLabel2 = " Unpolarized"
longLabel = "_long"
longLabel2 = " Longitudinal"
transLabel = "_trans"
transLabel2 = " Transverse"
AparLabel = "_Apar"
AparLabel2 = r'$\mathrm{A}_{\parallel}^{^{3}He}$'
AperpLabel = "_Aperp"
AperpLabel2 = r'$\mathrm{A}_{\perp}^{^{3}He}$'

outputFileFormat = "png"
outputFileExtension = "."+outputFileFormat


def q2_calc(Ep, theta):
    return 2*Ebeam*Ep*(1 - m.cos(theta*deg2rad))

def x_calc(Ep, theta):
    Q2 = q2_calc(Ep, theta)
    nu = Ebeam - Ep
    denom = 2*Mp*nu
    return Q2/denom

def gpLabel(theta):
    if (theta != theta_values[2]): return r'$\mathrm{g}_1 += 0.05\mathrm{F}_1$'
    return r'$\mathrm{g}_1 += 0.1\mathrm{F}_1$'

def gmLabel(theta):
    if (theta != theta_values[2]): return r'$\mathrm{g}_1 -= 0.05\mathrm{F}_1$'
    return r'$\mathrm{g}_1 -= 0.1\mathrm{F}_1$'

def main():

    # Theta value loop
    for theta_value in theta_values:

        theta = theta_value
        thetaStr = str(theta)

        print("\nCreating Plots for theta = " + thetaStr)
        for i in range(30): print("-", end='')
        print()

        # Get the input filename from the user
        inUnpolFileName = "Radiated Data/Model_IQE_Smearing/radiated_model_" + thetaStr + "deg" + unpolLabel + ".dat"
        inLongFileName = "Radiated Data/Model_IQE_Smearing/radiated_model_" + thetaStr + "deg" + longLabel + ".dat"
        inTransFileName = "Radiated Data/Model_IQE_Smearing/radiated_model_" + thetaStr + "deg" + transLabel + ".dat"
        inUnpol2FileName = "Radiated Data/Model_Inelastic_Smearing/radiated_model_" + thetaStr + "deg" + unpolLabel + ".dat"
        inLong2FileName = "Radiated Data/Model_Inelastic_Smearing/radiated_model_" + thetaStr + "deg" + longLabel + ".dat"
        inTrans2FileName = "Radiated Data/Model_Inelastic_Smearing/radiated_model_" + thetaStr + "deg" + transLabel + ".dat"
        outUnpolFileName = "Plots/" + thetaStr + "deg" + unpolLabel + outputFileExtension
        outLongFileName = "Plots/" + thetaStr + "deg" + longLabel + outputFileExtension
        outTransFileName = "Plots/" + thetaStr + "deg" + transLabel + outputFileExtension
        outAparFileName = "Plots/" + thetaStr + "deg" + AparLabel + outputFileExtension
        outAperpFileName = "Plots/" + thetaStr + "deg" + AperpLabel + outputFileExtension
        outApardiffFileName = "Plots/" + thetaStr + "deg" + AparLabel + "diff" + outputFileExtension
        outAperpdiffFileName = "Plots/" + thetaStr + "deg" + AperpLabel + "diff" + outputFileExtension

        # Read data from the input file
        try:
            with open(inUnpolFileName, 'r') as infile:
                next(infile)
                dataUnpol = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataUnpol.append(vals)
        except IOError:
            print("Unable to open file:", inUnpolFileName)
            return

        try:
            with open(inLongFileName, 'r') as infile:
                next(infile)
                dataLong = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataLong.append(vals)
        except IOError:
            print("Unable to open file:", inLongFileName)
            return

        try:
            with open(inTransFileName, 'r') as infile:
                next(infile)
                dataTrans = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataTrans.append(vals)
        except IOError:
            print("Unable to open file:", inTransFileName)
            return

        try:
            with open(inUnpol2FileName, 'r') as infile:
                next(infile)
                dataUnpol2 = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataUnpol2.append(vals)
        except IOError:
            print("Unable to open file:", inUnpol2FileName)
            return

        try:
            with open(inLong2FileName, 'r') as infile:
                next(infile)
                dataLong2 = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataLong2.append(vals)
        except IOError:
            print("Unable to open file:", inLong2FileName)
            return

        try:
            with open(inTrans2FileName, 'r') as infile:
                next(infile)
                dataTrans2 = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataTrans2.append(vals)
        except IOError:
            print("Unable to open file:", inTrans2FileName)
            return


        if (theta == 11):
            xvals = np.array([0.394,0.504,0.663,0.912])
            aparvals = np.array([-0.0044,-0.0000,0.0042,-0.001])
            aparerrs = np.array([0.0016,0.0019,0.0028,0.006])
            aperpvals = np.array([0.0007,-0.0001,-0.0023,-0.0073])
            aperperrs = np.array([0.0006,0.0007,0.0011,0.0023])
        elif (theta == 18):
            xvals = np.array([0.542,0.637,0.747,0.885])
            aparvals = np.array([-0.002,-0.002,0.003,-0.008])
            aparerrs = np.array([0.005,0.006,0.010,0.020])
            aperpvals = np.array([-0.0004,-0.0014,0.0050,-0.018])
            aperperrs = np.array([0.0016,0.0020,0.0032,0.006])
        else:
            xvals = np.array([0.423,0.440,0.465,0.492,0.517,0.544,0.571,0.599,0.629,0.647,0.676,0.708,0.737,0.767,0.799,0.834,0.867,0.901])
            aparvals = np.array([-0.00615,-0.00233,-0.00462,0.00448,0.00402,0.00348,0.00753,0.00873,0.00217,0.00302,0.00387,0.00429,0.00777,0.008,0.00812,0.00501,0.00717,0.01247])
            aparerrstat = np.array([0.00711,0.003,0.00326,0.00356,0.00242,0.00262,0.00285,0.00323,0.00189,0.002,0.0018,0.00206,0.00239,0.00282,0.00346,0.00584,0.00907,0.01179])
            aparerrsyst = np.array([0.00078,0.00125,0.00138,0.0015,0.00106,0.00097,0.00095,0.00114,0.0008,0.00078,0.00072,0.00031,0.0005,0.00037,0.00055,0.00072,0.00111,0.00141])
            aparerrs = np.sqrt(aparerrstat**2 + aparerrsyst**2)
            aperpvals = np.array([-0.02459,0.00831,0.00093,-0.01353,0.00866,0.01762,0.01071,0.0152,0.00809,0.005,0.00936,0.00627,0.00723,0.01034,0.01797,-0.0043,-0.00852,0.00376])
            aperperrstat = np.array([0.01754,0.00729,0.00787,0.00864,0.00626,0.00666,0.00746,0.00823,0.0048,0.00463,0.0039,0.00448,0.00518,0.00607,0.00759,0.01359,0.02221,0.02985])
            aperperrsyst = np.array([0.00128,0.00066,0.00085,0.00154,0.00279,0.00281,0.00332,0.00517,0.00299,0.00149,0.00134,0.00174,0.00254,0.00309,0.00378,0.00615,0.00948,0.01257])
            aperperrs = np.sqrt(aperperrstat**2 + aperperrsyst**2)

        xsradUnpol = [[],[]]
        xsbornUnpol = [[],[]]
        xsradLong = [[],[]]
        xsbornLong = [[],[]]
        xsradTrans = [[],[]]
        xsbornTrans = [[],[]]

        xsradUnpol2 = [[],[]]
        xsbornUnpol2 = [[],[]]
        xsradLong2 = [[],[]]
        xsbornLong2 = [[],[]]
        xsradTrans2 = [[],[]]
        xsbornTrans2 = [[],[]]


        Apar_rad = [[],[]]
        Apar_born = [[],[]]
        Aperp_rad = [[],[]]
        Aperp_born = [[],[]]
        Apar_diff = [[],[]]
        Aperp_diff = [[],[]]

        Apar2_rad = [[],[]]
        Apar2_born = [[],[]]
        Aperp2_rad = [[],[]]
        Aperp2_born = [[],[]]
        Apar2_diff = [[],[]]
        Aperp2_diff = [[],[]]


        for d in dataUnpol:
            xsradUnpol[0].append(Ebeam - d[1]/1000)
            xsradUnpol[1].append(d[2])
            xsbornUnpol[0].append(Ebeam - d[1]/1000)
            xsbornUnpol[1].append(d[3])
        for d in dataLong:
            xsradLong[0].append(Ebeam - d[1]/1000)
            xsradLong[1].append(d[2])
            xsbornLong[0].append(Ebeam - d[1]/1000)
            xsbornLong[1].append(d[3])
        for d in dataTrans:
            xsradTrans[0].append(Ebeam - d[1]/1000)
            xsradTrans[1].append(d[2])
            xsbornTrans[0].append(Ebeam - d[1]/1000)
            xsbornTrans[1].append(d[3])

        for d in dataUnpol2:
            xsradUnpol2[0].append(Ebeam - d[1]/1000)
            xsradUnpol2[1].append(d[2])
            xsbornUnpol2[0].append(Ebeam - d[1]/1000)
            xsbornUnpol2[1].append(d[3])
        for d in dataLong2:
            xsradLong2[0].append(Ebeam - d[1]/1000)
            xsradLong2[1].append(d[2])
            xsbornLong2[0].append(Ebeam - d[1]/1000)
            xsbornLong2[1].append(d[3])
        for d in dataTrans2:
            xsradTrans2[0].append(Ebeam - d[1]/1000)
            xsradTrans2[1].append(d[2])
            xsbornTrans2[0].append(Ebeam - d[1]/1000)
            xsbornTrans2[1].append(d[3])

        for i in range(len(dataLong)):
            if (dataLong[i][3] == 0): continue
            tempEp = Ebeam - (dataLong[i][1] + dataUnpol[i][1])/2000
            tempx = x_calc(tempEp, theta)
            if (tempx > 1): continue
            Apar_rad[0].append(tempx)
            Apar_born[0].append(tempx)
            Apar_rad[1].append(dataLong[i][2]/(2.0*dataUnpol[i][2]))
            Apar_born[1].append(dataLong[i][3]/(2.0*dataUnpol[i][3]))
            Apar_diff[0].append(tempx)
            Apar_diff[1].append(dataLong[i][3]/(2.0*dataUnpol[i][3]) - dataLong[i][2]/(2.0*dataUnpol[i][2]))
            
        for i in range(len(dataTrans)):
            if (dataTrans[i][3] == 0): continue
            tempEp = Ebeam - (dataTrans[i][1] + dataUnpol[i][1])/2000
            tempx = x_calc(tempEp, theta)
            if (tempx > 1): continue
            Aperp_rad[0].append(tempx)
            Aperp_born[0].append(tempx)
            Aperp_rad[1].append(dataTrans[i][2]/(2.0*dataUnpol[i][2]))
            Aperp_born[1].append(dataTrans[i][3]/(2.0*dataUnpol[i][3]))
            Aperp_diff[0].append(tempx)
            Aperp_diff[1].append(dataTrans[i][3]/(2.0*dataUnpol[i][3]) - dataTrans[i][2]/(2.0*dataUnpol[i][2]))


        for i in range(len(dataLong2)):
            if (dataLong2[i][3] == 0): continue
            tempEp = Ebeam - (dataLong2[i][1] + dataUnpol2[i][1])/2000
            tempx = x_calc(tempEp, theta)
            if (tempx > 1): continue
            Apar2_rad[0].append(tempx)
            Apar2_born[0].append(tempx)
            Apar2_rad[1].append(dataLong2[i][2]/(2.0*dataUnpol2[i][2]))
            Apar2_born[1].append(dataLong2[i][3]/(2.0*dataUnpol2[i][3]))
            Apar2_diff[0].append(tempx)
            Apar2_diff[1].append(dataLong2[i][3]/(2.0*dataUnpol2[i][3]) - dataLong2[i][2]/(2.0*dataUnpol2[i][2]))
            
        for i in range(len(dataTrans2)):
            if (dataTrans2[i][3] == 0): continue
            tempEp = Ebeam - (dataTrans2[i][1] + dataUnpol2[i][1])/2000
            tempx = x_calc(tempEp, theta)
            if (tempx > 1): continue
            Aperp2_rad[0].append(tempx)
            Aperp2_born[0].append(tempx)
            Aperp2_rad[1].append(dataTrans2[i][2]/(2.0*dataUnpol2[i][2]))
            Aperp2_born[1].append(dataTrans2[i][3]/(2.0*dataUnpol2[i][3]))
            Aperp2_diff[0].append(tempx)
            Aperp2_diff[1].append(dataTrans2[i][3]/(2.0*dataUnpol2[i][3]) - dataTrans2[i][2]/(2.0*dataUnpol2[i][2]))

                
        horizLineUnpol = [[xsradUnpol[0][0],xsradUnpol[0][len(xsradUnpol[0])-1]],[0,0]]
        horizLineLong = [[xsradLong[0][0],xsradLong[0][len(xsradLong[0])-1]],[0,0]]
        horizLineTrans = [[xsradTrans[0][0],xsradTrans[0][len(xsradTrans[0])-1]],[0,0]]
        horizLineApar = [[Apar_born[0][0],Apar_born[0][len(Apar_rad[0])-1]],[0,0]]
        horizLineAperp = [[Aperp_rad[0][0],Aperp_rad[0][len(Aperp_rad[0])-1]],[0,0]]

        plt.figure(figsize=(10, 5))  # Adjust width and height as needed

        bbox_inches = 'tight' #Bbox.from_extents(0.1,0.1,0.9,0.9)

        # Plot the data
        plt.plot(xsradUnpol[0], xsradUnpol[1], color="blue", label='XS rad (I+QE)')
        plt.plot(xsbornUnpol[0], xsbornUnpol[1], color="green", label='XS born (I+QE)')
        plt.plot(xsradUnpol2[0], xsradUnpol2[1], color="blue", linestyle='--', label='XS rad (I)')
        plt.plot(xsbornUnpol2[0], xsbornUnpol2[1], color="green", linestyle='--', label='XS born (I)')       
        plt.plot(horizLineUnpol[0], horizLineUnpol[1], color='black', linestyle='--')
        plt.xlabel('E\' [GeV]')
        plt.ylabel('XS')
        plt.title(unpolLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign + ' Inelastic + Quasi-Elastic Smearing and Inelastic Smearing')
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outUnpolFileName, format=outputFileFormat, bbox_inches=bbox_inches)
        print("Plot saved to " + outUnpolFileName)

        plt.clf()

        # Plot the data
        plt.plot(xsradLong[0], xsradLong[1], color="blue", label='XS rad (I+QE)')
        plt.plot(xsbornLong[0], xsbornLong[1], color="green", label='XS born (I+QE)')
        plt.plot(xsradLong2[0], xsradLong2[1], color="blue", linestyle='--', label='XS rad (I)')
        plt.plot(xsbornLong2[0], xsbornLong2[1], color="green", linestyle='--', label='XS born (I)')
        plt.plot(horizLineLong[0], horizLineLong[1], color='black', linestyle='--')
        plt.xlabel('E\' [GeV]')
        plt.ylabel('XS')
        plt.title(longLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign + ' Inelastic + Quasi-Elastic Smearing and Inelastic Smearing')
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outLongFileName, format=outputFileFormat, bbox_inches=bbox_inches)
        print("Plot saved to " + outLongFileName)

        plt.clf()

        # Plot the data
        plt.plot(xsradTrans[0], xsradTrans[1], color="blue", label='XS rad (I+QE)')
        plt.plot(xsbornTrans[0], xsbornTrans[1], color="green", label='XS born (I+QE)')
        plt.plot(xsradTrans2[0], xsradTrans2[1], color="blue", linestyle='--', label='XS rad (I)')
        plt.plot(xsbornTrans2[0], xsbornTrans2[1], color="green", linestyle='--', label='XS born (I)')
        plt.plot(horizLineTrans[0], horizLineTrans[1], color='black', linestyle='--')
        plt.xlabel('E\' [GeV]')
        plt.ylabel('XS')
        plt.title(transLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign + ' Inelastic + Quasi-Elastic Smearing and Inelastic Smearing')
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outTransFileName, format=outputFileFormat, bbox_inches=bbox_inches)
        print("Plot saved to " + outTransFileName)

        plt.clf()

        # Plot the data
        plt.plot(Apar_rad[0], Apar_rad[1], color="blue", label=AparLabel2+' rad (I+QE)')
        plt.plot(Apar_born[0], Apar_born[1], color="green", label=AparLabel2+' born (I+QE)')
        plt.plot(Apar2_rad[0], Apar2_rad[1], color="blue", linestyle='--', label=AparLabel2+' rad (I)')
        plt.plot(Apar2_born[0], Apar2_born[1], color="green", linestyle='--', label=AparLabel2+' born (I)')
        plt.plot(horizLineApar[0], horizLineApar[1], color='black', linestyle='--')
        plt.errorbar(xvals, aparvals, yerr=aparerrs, fmt='o', color='red', markersize=5, capsize=5, label=AparLabel2+' data')
        plt.xlabel('x')
        plt.ylabel(AparLabel2)
        plt.title(AparLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign + ' Inelastic + Quasi-Elastic Smearing and Inelastic Smearing')
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outAparFileName, format=outputFileFormat, bbox_inches=bbox_inches)
        print("Plot saved to " + outAparFileName)

        plt.clf()

        # Plot the data
        plt.plot(Aperp_rad[0], Aperp_rad[1], color="blue", label=AperpLabel2+' rad (I+QE)')
        plt.plot(Aperp_born[0], Aperp_born[1], color="green", label=AperpLabel2+' born (I+QE)')
        plt.plot(Aperp2_rad[0], Aperp2_rad[1], color="blue", linestyle='--', label=AperpLabel2+' rad (I)')
        plt.plot(Aperp2_born[0], Aperp2_born[1], color="green", linestyle='--', label=AperpLabel2+' born (I)')
        plt.plot(horizLineAperp[0], horizLineAperp[1], color='black', linestyle='--')
        plt.errorbar(xvals, aperpvals, yerr=aperperrs, fmt='o', color='red', markersize=5, capsize=5, label=AperpLabel2+' data')
        plt.xlabel('x')
        plt.ylabel(AperpLabel2)
        plt.title(AperpLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign + ' Inelastic + Quasi-Elastic Smearing and Inelastic Smearing')
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outAperpFileName, format=outputFileFormat, bbox_inches=bbox_inches)
        print("Plot saved to " + outAperpFileName)

        plt.clf()

        # Plot the data
        plt.plot(Apar_diff[0], Apar_diff[1], color="blue", label=AparLabel2+'(born - rad) (I+QE)')
        plt.plot(Apar2_diff[0], Apar2_diff[1], color="blue", linestyle='--', label=AparLabel2+'(born - rad) (I)')
        plt.plot(horizLineAperp[0], horizLineAperp[1], color='black', linestyle='--')
        plt.xlabel('x')
        plt.ylabel(AparLabel2 + '(born - rand)')
        plt.title(AparLabel2 + '(born - rad) E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign + ' Inelastic + Quasi-Elastic Smearing and Inelastic Smearing')
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outApardiffFileName, format=outputFileFormat, bbox_inches=bbox_inches)
        print("Plot saved to " + outApardiffFileName)

        plt.clf()

        # Plot the data
        plt.plot(Aperp_diff[0], Aperp_diff[1], color="blue", label=AperpLabel2+'(born - rad) (I+QE)')
        plt.plot(Aperp2_diff[0], Aperp2_diff[1], color="blue", linestyle='--', label=AperpLabel2+'(born - rad) (I)')
        plt.plot(horizLineAperp[0], horizLineAperp[1], color='black', linestyle='--')
        plt.xlabel('x')
        plt.ylabel(AperpLabel2 + '(born - rand)')
        plt.title(AperpLabel2 + '(born - rad) E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign + ' Inelastic + Quasi-Elastic Smearing and Inelastic Smearing')
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outAperpdiffFileName, format=outputFileFormat, bbox_inches=bbox_inches)
        print("Plot saved to " + outAperpdiffFileName)

        plt.clf()
    
    for i in range(30): print("-", end='')
    print()




if __name__ == "__main__":
    main()