#pylint:disable=C0114,C0103
import subprocess
import numpy as np
import pandas as pd
import json 


program_type = "hp_3d"
file_name = "/".join(["./fortran",program_type,"input.dat"])

with open(file_name,"w") as input_file:
    ib=10	    # Pannelli chordwise
    jb=10	    # Pannelli spanwise
    nst=60	    # Numero step temporali
    vt=1.	    # Velocità
    b=4.	    # apertura alare	    
    amed=0.	    # angolo di attacco medio in gradi 
    temp=0.5    # divisore temporale (default è 4)
    aa=5.       # ampiezza pitch in gradi
    oma=2.	    # omega pitch
    sza=0.1	    # ampiezza heave
    omh=2.	    # omega heave
    fase=0.	    # fase heave in gradi (fase di pitch == 0)
    thic=.04    # spessore NACA (rispetto alla corda unitaria)

    input_file.write(str(ib) + "\n")
    input_file.write(str(jb) + "\n")
    input_file.write(str(nst) + "\n")
    input_file.write(str(vt) + "\n")
    input_file.write(str(b) + "\n")
    input_file.write(str(amed) + "\n")
    input_file.write(str(temp) + "\n")
    input_file.write(str(aa) + "\n")
    input_file.write(str(oma) + "\n")
    input_file.write(str(sza) + "\n")
    input_file.write(str(omh) + "\n")
    input_file.write(str(fase) + "\n")
    input_file.write(str(thic) + "\n")

subprocess.run(["./p17"],cwd="./fortran/hp_3d",capture_output=True, text=True)

with open("./fortran/hp_3d/out/output.csv","r") as output_file:
    results = pd.read_csv(output_file,dtype=np.float64)

parameters = {}
parameters["chordwise_panels"] = ib
parameters["spanwise_panels"] = jb
parameters["total_steps"] = nst
parameters["speed"] = vt
parameters["wing_span"] = b
parameters["mean_alpha0"] = amed
parameters["temp_division"] = temp
parameters["pitch_amplitude"] = aa
parameters["pitch_omega"] = oma
parameters["heave_amplitude"] = sza
parameters["heave_omega"] = omh
parameters["heave_phase"] = fase
parameters["naca_thickness"] = thic

results_dict = {}
results_dict["lift_coefficient"] = results["CL"].tolist()
results_dict["drag_coefficient"] = results["CD"].tolist()
results_dict["moment_coefficient"] = results["CM"].tolist()

output_dict = {}
output_dict["parameters"] = parameters
output_dict["data"] = results_dict

with open("./output.json","w") as json_file:
    json.dump(output_dict,json_file)