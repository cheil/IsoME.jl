using IsoME
using Plots
using CSV
using DataFrames
using LinearAlgebra

# dos
dos_file = "/temp/spathd/MasterThesis/materials/Nb_Antonio/dos.dat"
en_antonio, dos_antonio, ef, unit = IsoME.readIn_Dos(dos_file, 0, 1, 1)


dos_file2 = "/temp/kogler_e/materials/Nb/bands/Nb.dos"
energies2, dos2, ef2, unit2 = IsoME.readIn_Dos(dos_file2)

default(show=true)
plot(en_antonio, dos_antonio)
plot!(energies2, dos2)


#=
# a2f
a2f_file = "/temp/spathd/MasterThesis/materials/Nb_Antonio/elias_ph.in"
omega, a2f = IsoME.readIn_a2f(a2f_file) 

a2f_file2 = "/temp/kogler_e/materials/Nb/phonons/Nb.a2F"
omega2, a2f2 = IsoME.readIn_a2f(a2f_file2, 1) 

default(show=true)
p= plot(omega, a2f, label="Antonio")
p=plot(p, omega2, a2f2, label="wir")
#savefig("a2f.pdf")

# Weep
Weep_file           = "/temp//spathd/MasterThesis/materials/Nb_Antonio/KC.OUT"
Wen_file            = "/temp/spathd/MasterThesis/materials/Nb_Antonio/Wen.dat"


Weep, unit  = IsoME.readIn_Weep(Weep_file)
Wen         = IsoME.readIn_Wen(Wen_file)
idx_ef = findmin(abs.(Wen))
idx_ef = idx_ef[2]


Wen_file            = "/temp/kogler_e/materials/Nb/10.DOS.dat"
Weep_file           = "/temp/kogler_e/materials/Nb/10.Weep.dat"

Weep2, unit2  = IsoME.readIn_Weep(Weep_file)
Wen2         = IsoME.readIn_Wen(Wen_file)
Wen2 = Wen2 .- ef2
idx_ef2 = findmin(abs.(Wen2))
idx_ef2 = idx_ef2[2]

#=
### RPA & KO from Antonio
Ry2meV = 13605.662285137
RPA = CSV.read("/afs/itp.tugraz.at/user/spathd/Downloads/Weep_RPA.csv", DataFrame)
en_RPA = RPA[:,2].*Ry2meV
Weep_RPA = RPA[:,1].*Ry2meV
KO = CSV.read("/afs/itp.tugraz.at/user/spathd/Downloads/Weep_KO.csv", DataFrame)
en_KO = KO[:,2].*Ry2meV
Weep_KO = KO[:,1].*Ry2meV
=#



### weep/bad_dos*dos
en_weep, dos_weep, ef_weep, unit_weep = IsoME.readIn_Dos(Wen_file, 1,0,1)

=#

#=
# interpolate weep on dos
en_interval = [en_weep[findfirst(en_weep .> energies2[1])]; en_weep[findlast(en_weep .< energies2[end])]]
# number of points overlapping
idxOverlap = [findfirst(energies2 .> en_weep[1]), findlast(energies2 .< en_weep[end])]
Nitp = idxOverlap[2] - idxOverlap[1] + 1
energies2_itp = energies2[idxOverlap[1]:idxOverlap[2]]
dos2_itp = dos2[idxOverlap[1]:idxOverlap[2]]

#dos_en, dos_antonio = interpolateDos(dos_en, dos_antonio, en_interval, Nitp)
Weep2 = IsoME.interpolateWeep(en_weep, Weep2, en_interval, Nitp)
en_weep, dos_weep = IsoME.interpolateDos(en_weep, dos_weep, en_interval, Nitp)

default(show=true)
p= plot(Wen, diag(Weep), label="Antonio")
p= plot(p,energies2_itp, diag(Weep2.*(dos_weep*dos_weep')./(dos2_itp*dos2_itp')), label="Wir pp")
p= plot(p,Wen2, diag(Weep2), label="Wir")
savefig("Weep_dos2.pdf")
=#

# dos
#=
default(show=true)
p= plot(en_antonio, dos_antonio, label="Antonio")
p=plot(p, energies2, dos2, label="wir")
p = plot(p,en_weep, dos_weep, label="weep")
savefig("dos.pdf")
=#

#=
### Al
# ef
dos_file2 = "/temp/kogler_e/materials/Al/qe_2016/bands/Al.dos"
energies2, dos2, ef2, unit2 = IsoME.readIn_Dos(dos_file2, 1)

Weep_file           = "/temp//spathd/MasterThesis/Al_Antonio/KC_Al.OUT"
Wen_file            = "/temp/spathd/MasterThesis/Al_Antonio/Wen.dat"


Weep, unit  = IsoME.readIn_Weep(Weep_file)
Wen         = IsoME.readIn_Wen(Wen_file)
idx_ef = findmin(abs.(Wen))
idx_ef = idx_ef[2]

Wen_file            = "/temp/kogler_e/materials/Al/qe_2016/10.DOS.dat"
Weep_file           = "/temp/kogler_e/materials/Al/qe_2016/10.Weep.dat"

Weep2, unit2  = IsoME.readIn_Weep(Weep_file)
Wen2         = IsoME.readIn_Wen(Wen_file)
Wen2 = Wen2 .- ef2
idx_ef2 = findmin(abs.(Wen2))
idx_ef2 = idx_ef2[2]

default(show=true)
p= plot(Wen, diag(Weep), label="Antonio")
p= plot(p,Wen2, diag(Weep2), label="Wir")
savefig("Weep_Al.pdf")
=#