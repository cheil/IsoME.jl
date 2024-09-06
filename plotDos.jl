using IsoME
using Plots

# dos
dos_file = "/afs/itp.tugraz.at/user/spathd/Downloads/Nb/dos.dat"
energies, dos, ef, unit = IsoME.readIn_Dos(dos_file, 0, 1, 1)
ef = 0

dos_file2 = "/temp/kogler_e/materials/Nb/bands/Nb.dos"
energies2, dos2, ef2, unit2 = IsoME.readIn_Dos(dos_file2, 1)

default(show=true)
p= plot(energies, dos, label="Antonio")
p=plot(p, energies2, dos2, label="wir")
savefig("dos.pdf")


# a2f
a2f_file = "/afs/itp.tugraz.at/user/spathd/Downloads/Nb/elias_ph.in"
omega, a2f = IsoME.readIn_a2f(a2f_file) 

a2f_file2 = "/temp/kogler_e/materials/Nb/phonons/Nb.a2F"
omega2, a2f2 = IsoME.readIn_a2f(a2f_file2, 1) 

default(show=true)
p= plot(omega, a2f, label="Antonio")
p=plot(p, omega2, a2f2, label="wir")
savefig("a2f.pdf")
