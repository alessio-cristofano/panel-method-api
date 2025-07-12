#pylint:disable=C0114,C0103
import subprocess

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

result = subprocess.run(["./p17"],cwd="./fortran/hp_3d",capture_output=True, text=True)

print(result)