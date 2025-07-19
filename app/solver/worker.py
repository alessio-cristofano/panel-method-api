#pylint:disable=C0114,C0103
import subprocess
import numpy as np
import pandas as pd
from flask import jsonify

def solve(solve_input:dict):
    program_type = "hp_3d"
    file_name = "/".join(["./app/solver",program_type,"input.dat"])

    with open(file_name,"w") as input_file:
        ib      = solve_input["ib"]	    # Pannelli chordwise
        jb      = solve_input["jb"]	    # Pannelli spanwise
        nst     = solve_input["nst"]    # Numero step temporali
        vt      = solve_input["vt"]     # Velocità
        b       = solve_input["b"]      # apertura alare	    
        amed    = solve_input["amed"]   # angolo di attacco medio in gradi 
        temp    = solve_input["temp"]   # divisore temporale (default è 4)
        aa      = solve_input["aa"]     # ampiezza pitch in gradi
        oma     = solve_input["oma"]    # omega pitch
        sza     = solve_input["sza"]    # ampiezza heave
        omh     = solve_input["omh"]    # omega heave
        fase    = solve_input["fase"]   # fase heave in gradi (fase di pitch == 0)
        thic    = solve_input["thic"]   # spessore NACA (rispetto alla corda unitaria)

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

    subprocess.run(["./p17"],cwd="./app/solver/hp_3d",capture_output=True, text=True)

    with open("./app/solver/hp_3d/out/output.csv","r") as output_file:
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

    # with open("./output.json","w") as json_file:
    #     json.dump(output_dict,json_file)
    return output_dict