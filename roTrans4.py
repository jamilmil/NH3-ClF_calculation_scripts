import psi4
import numpy as np
import rotation_matrices
import os
#from termcolor import colored

#defining variables for computation / user-friendly / Change these values to suite your study
#   define the number of threads that must be used
num_threads=4
#   define the basis set
basis_set="jun-cc-pVDZ"
#   define SAPT's level of theory
method="sapt0"
#   define the starting angle in degree
theta_0=0
#   define the ending angle in degree
theta_f=120
#   define the step for the angle change
rot_step=1
#   define the minimal distance between the 2 monomers
d_0=2.4
#define the maximal distance between the monomers
d_f=2.5
#   define the step for the distance change
trans_step=0.1
#   In this field, past the XYZ file of the acceptor molecule
acceptor="""
2

F	0.0000000	0.0000000	-1.629471
Cl	0.0000000	0.0000000	0.0000000
"""
#   In this field, past the XY file of the donor molecule
donor="""
4

N	0.0000000	0.0000000	 0.0000000
H	0.0000000	0.9380680	0.380024
H	0.8123910	-0.4690340	0.380024
H	-0.8123910	-0.4690340	0.380024
"""

#   From here the main code for the computation starts / not user-friendly / Please don't change unless you know what you're doing
#function for extracting informations from xyz file
def get(x_file):
    x=open(x_file,"r")
    x_header=int(x.readline())
    a=x.readlines()[1:]
    x_elements=[]; x_coord=[]
    for i in range(len(a)):
        set=[]
        x_elements.append(a[i].split()[0])
        set.append(float(a[i].split()[1]))
        set.append(float(a[i].split()[2]))
        set.append(float(a[i].split()[3]))
        x_coord.append(set[:])
        set.clear()
    return x_header,x_elements, x_coord

def get_data(acceptor):
    file=acceptor.split(sep="\n")
    x_header=int(file[1])
    x_elements=[]; x_coord=[]
    for i in range(len(file)-4):
        set=[]
        x_elements.append(file[i+3].split()[0])
        set.append(float(file[i+3].split()[1]))
        set.append(float(file[i+3].split()[2]))
        set.append(float(file[i+3].split()[3]))
        x_coord.append(set[:])
        set.clear()
    return x_header,x_elements, x_coord
#function for seperating the 2 fragments by a user defined intitial distance
def seperate_fragments(x_header,x_elements,x_coord,d_min):
    file=open("tmp.xyz","w")
    file.write(str(x_header)+"\n\n")
    for i in range(len(x_elements)):
        file.write(str(x_elements[i]+"\t"+str(x_coord[i][0]))+"\t"+str(x_coord[i][1])+"\t"+str(x_coord[i][2]+d_min)+"\n")
    file.close()
    return 0
#function for creating shifts according to a computed increment from a user defined number of step, returns a list of shifts
def incrementation(a,b,c):
    x=(b-a)/c
    y=[a]
    z=a+c
    while len(y)<=x:
        y.append(z)
        z+=c
    return y

psi4.set_num_threads(num_threads)

acceptor_elements=get_data(acceptor)[1]; acceptor_coord=get_data(acceptor)[2];acceptor_header=get_data(acceptor)[0]
donor_elements=get_data(donor)[1]; donor_coord=get_data(donor)[2];donor_header=get_data(donor)[0]

translation_incs=incrementation(d_0,d_f,trans_step)
rotation_incs=incrementation(theta_0,theta_f,rot_step)


os.mkdir("calc_test")
os.chdir("calc_test")
result=open("summary.dat","w")
result.write("""d(ang)\ttheta(deg)\tEint\tEelst\tExr\tEind\tEdisp
""")
for d in translation_incs:
    os.mkdir(str(d)+"_trans")
    os.chdir(str(d)+"_trans")
    seperate_fragments(donor_header,donor_elements,donor_coord,d)
    donor_coord_sep=get("tmp.xyz")[2]
    os.remove("tmp.xyz")
    #print(colored("\nCalculations for a separation distance of d="+str(d)+" angstrom\n",on_color="on_red"))
    for theta in rotation_incs:
        a=open("tmp2.xyz","w")
        R=rotation_matrices.x_axis_rot(theta)
        #R=rotation_matrices.y_axis_rot(theta)
        #R=rotation_matrices.z_axis_rot(theta)
        
        for i in range(len(donor_elements)):
            a_new=np.zeros(3)
            a_new=np.matmul(R,donor_coord_sep[i])
            a.write(str(donor_elements[i])+"\t"+str(a_new[0])+"\t"+str(a_new[1])+"\t"+str(a_new[2])+"\n")
        a.close()

        b=open("tmp3.xyz","w"); dimer=open("dimer_"+str(d)+"_ang_"+str(theta)+"_deg.xyz","w"); c=open("tmp2.xyz","r")
        dimer.write(str(acceptor_header+donor_header)+"\n\n")
        r=c.read()
        b.write("0 1\n")
        for j in range(len(acceptor_elements)):
            b.write(str(acceptor_elements[j])+"\t"+str(acceptor_coord[j][0])+"\t"+str(acceptor_coord[j][1])+"\t"+str(acceptor_coord[j][2])+"\n")
            dimer.write(str(acceptor_elements[j])+"\t"+str(acceptor_coord[j][0])+"\t"+str(acceptor_coord[j][1])+"\t"+str(acceptor_coord[j][2])+"\n")
        dimer.write(str(r)); b.write("--\n0 1\n"); b.write(str(r))
        b.close(); c.close(); dimer.close()
        os.remove("tmp2.xyz")

        e=open("tmp3.xyz","r"); f=e.read()
        #print(colored("\nCalculations for a rotation of theta="+str(theta)+" deg\n",on_color="on_green"))
        psi4.geometry(f,name="dimer")
        e.close()
        os.remove("tmp3.xyz")

        psi4.energy(method+"/"+basis_set)
        etot=psi4.variable('SAPT TOTAL ENERGY')
        eelst=psi4.variable('SAPT ELST ENERGY')
        exr=psi4.variable('SAPT EXCH ENERGY')
        eind=psi4.variable('SAPT IND ENERGY')
        edisp=psi4.variable('SAPT DISP ENERGY')
        result.write(str(d)+"\t"+str(theta)+"\t"+str(etot)+"\t"+str(eelst)+"\t"+str(exr)+"\t"+str(eind)+"\t"+str(edisp)+"\n")
    result.write("\n")
    os.chdir("../")
result.close()