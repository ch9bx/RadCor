import csv
import math as m

deg2rad = 0.0174533     # degrees to radians
Mp = 938                # MeV

def quad(x,y):
    return m.sqrt(x*x + y*y)

def diff_errcalc(dx, dy):
    return m.sqrt(dx*dx + dy*dy)

def asymm_errcalc(x, y, delta_x, delta_y):
    partial_z_x = 1 / (2 * y)
    partial_z_y = -x / (2 * y**2)
    delta_z = m.sqrt((delta_x / (2 * y))**2 + (x * delta_y / (2 * y**2))**2)
    return delta_z

def q2_calc(Ebeam, Ep, theta):
    return 2*Ebeam*Ep*(1 - m.cos(theta*deg2rad))

def x_calc(Ebeam, Ep, theta):
    Q2 = q2_calc(Ebeam, Ep, theta)
    nu = Ebeam - Ep
    denom = 2*Mp*nu
    return Q2/denom

class XS:
    def __init__(self, unpol, unpolerr, para, paraerr, perp, perperr):
        self.unpol = unpol
        self.unpolerr = unpolerr
        self.para = para
        self.paraerr = paraerr
        self.perp = perp
        self.perperr = perperr

class Asymm:
    def __init__(self, xs):
        if (xs.unpol == 0):
            self.para = ''
            self.paraerr = ''
            self.perp = ''
            self.perperr = ''
        else:
            self.para = xs.para/(2*xs.unpol)
            self.paraerr = asymm_errcalc(xs.para, xs.unpol, xs.paraerr, xs.unpolerr)
            self.perp = xs.perp/(2*xs.unpol)
            self.perperr = asymm_errcalc(xs.perp, xs.unpol, xs.perperr, xs.unpolerr)

class AsymmDiff:
    def __init__(self, asymm1, asymm2):
        if (asymm1.para == '' or asymm2.para == ''):
            self.para = ''
            self.paraerr = ''
        else:
            self.para = asymm1.para - asymm2.para
            self.paraerr = diff_errcalc(asymm1.paraerr,asymm2.paraerr)
        if (asymm1.perp == '' or asymm2.para == ''):
            self.perp = ''
            self.perperr = ''
        else:
            self.perp = asymm1.perp - asymm2.perp
            self.perperr = diff_errcalc(asymm1.perperr,asymm2.perperr)


class Data:
    def __init__(self, E, nu, xsrad, xsborn, errstat, errsyst):
        self.E = E
        self.nu = nu
        self.xsrad = xsrad
        self.xsborn = xsborn
        self.errstat = errstat
        self.errsyst = errsyst

class Output:
    def __init__(self, E, nu, x, q2, xsborn, asymmborn, xsrad, asymmrad, asymmdiff):
        self.E = E
        self.nu = nu
        self.x = x
        self.q2 = q2
        self.xsborn = xsborn
        self.asymmborn = asymmborn
        self.xsrad = xsrad
        self.asymmrad = asymmrad
        self.asymmdiff = asymmdiff

def read_data_from_file(filename):
    data = []
    try:
        with open(filename, 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) == 6:
                    E,nu,xsrad,xsborn,errstat,errsyst = parts
                    if float(E) != 10380: continue
                    data.append(Data(float(E), float(nu), float(xsrad), float(xsborn), float(errstat), float(errsyst)))
    except FileNotFoundError:
        print(f"Error: Unable to open file {filename}")
    return data

def combine_data(unpol, para, perp, theta):
    combined = []
    for i in range(len(unpol)):
        if (unpol[i].xsborn == 0 or para[i].xsborn == 0 or perp[i].xsborn == 0): continue
        xsborn = XS(unpol[i].xsborn,0,para[i].xsborn,0,perp[i].xsborn,0)
        xsrad = XS(unpol[i].xsrad, quad(unpol[i].errstat,unpol[i].errsyst), para[i].xsrad, 
        quad(para[i].errstat,para[i].errsyst), perp[i].xsrad, quad(perp[i].errstat,perp[i].errsyst))
        q2 = q2_calc(unpol[i].E, unpol[i].E - unpol[i].nu, theta)
        x = x_calc(unpol[i].E, unpol[i].E - unpol[i].nu, theta)
        asymmborn = Asymm(xsborn)
        asymmrad = Asymm(xsrad)
        asymmdiff = AsymmDiff(asymmborn,asymmrad)
        combined.append(Output(unpol[i].E, unpol[i].nu, x, q2, xsborn, asymmborn, xsrad, asymmrad, asymmdiff))
    return combined

def write_data_to_csv(data, filename):
    try:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Ebeam","nu","x","Q2","XSborn_unpol","XSborn_unpol_err","XSborn_para","XSborn_para_err",
            "XSborn_perp","XSborn_perp_err", "Asymmborn_para", "Asymmborn_para_err", "Asymmborn_perp", 
            "Asymmborn_perp_err", "XSrad_unpol","XSrad_unpol_err","XSrad_para","XSrad_para_err",
            "XSrad_perp","XSrad_perp_err", "Asymmrad_para", "Asymmrad_para_err", "Asymmrad_perp", 
            "Asymmrad_perp_err","Asymmdiff_para","Asymmdiff_para_err","Asymmdiff_perp","Asymmdiff_perp_err"])
            for d in data:
                writer.writerow([d.E, d.nu, d.x, d.q2, d.xsborn.unpol, d.xsborn.unpolerr, d.xsborn.para, d.xsborn.paraerr,
                d.xsborn.perp, d.xsborn.perperr, d.asymmborn.para, d.asymmborn.paraerr, d.asymmborn.perp, 
                d.asymmborn.perperr, d.xsrad.unpol, d.xsrad.unpolerr, d.xsrad.para, d.xsrad.paraerr, d.xsrad.perp,
                d.xsrad.perperr, d.asymmrad.para, d.asymmrad.paraerr, d.asymmrad.perp, d.asymmrad.perperr,
                d.asymmdiff.para, d.asymmdiff.paraerr, d.asymmdiff.perp, d.asymmdiff.perperr])
    except IOError:
        print(f"Error: Unable to create file {filename}")

def main():
    theta_vals = [11,18,30]

    for theta in theta_vals:
        thetastr = str(theta)

        # Change these filenames as per your file paths
        # unpolFile = "CAnalyzer-master/example/output/radiated_model_"+thetastr+"deg_unpol.dat"
        # paraFile = "CAnalyzer-master/example/output/radiated_model_"+thetastr+"deg_long.dat"
        # perpFile = "CAnalyzer-master/example/output/radiated_model_"+thetastr+"deg_trans.dat"
        unpolFile = "Radiated Data/Model_IQE_Smearing_g1-F1/radiated_model_"+thetastr+"deg_unpol.dat"
        paraFile = "Radiated Data/Model_IQE_Smearing_g1-F1/radiated_model_"+thetastr+"deg_long.dat"
        perpFile = "Radiated Data/Model_IQE_Smearing_g1-F1/radiated_model_"+thetastr+"deg_trans.dat"

        combined_data = []

        # Read data from each file
        unpol = read_data_from_file(unpolFile)
        para = read_data_from_file(paraFile)
        perp = read_data_from_file(perpFile)

        # Combine data from all files
        combined = combine_data(unpol, para, perp, theta)

        # Write combined data to CSV
        outputFile = "Radiated Data/radiated_data_"+thetastr+"deg.csv"
        write_data_to_csv(combined, outputFile)

        print("Combined data has been written to " + outputFile)

if __name__ == "__main__":
    main()
