#!/bin/bash
#
# Compilazione. 
# Rimuovere i flag -freal e -fdefault se si vuole precisione singola
# ===============================================================================
gfortran -c moduhp.f90  -g  -fcheck=all		-freal-4-real-8 -fdefault-real-8
gfortran -c p17.f90   -g  -fcheck=all		-freal-4-real-8 -fdefault-real-8
gfortran -o p17 p17.o moduhp.o -g -fcheck=all	-freal-4-real-8 -fdefault-real-8
# ================================================================================
#
# INPUT
# ------------------------------
#
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
#
scriviinput () {
echo $ib    >  input.dat
echo $jb    >> input.dat
echo $nst   >> input.dat 
echo $b	    >> input.dat 
echo $vt    >> input.dat 
echo $amed  >> input.dat
echo $temp  >> input.dat
echo $aa    >> input.dat
echo $oma   >> input.dat
echo $sza   >> input.dat
echo $omh   >> input.dat
echo $fase  >> input.dat
echo $thic  >> input.dat
}
# ------------------------------
# Esecuzione comandi
# ------------------------------
mkdir out/
scriviinput 
./p17
