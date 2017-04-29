#/bin/bash
gfortran numerovL.f90

#Triplet
energytr=( "-2.224" "-11.5" "-19.5" )
energysi=( "0" "-7.4" "-15.9" )
lamda=( 2 4 6 8 )
lec=( -150 -400 -900 -1500 )
mass=( 140 510 805)

#for j in {0..2}
#do
#	for i in {0..3}
#	do
#		echo ./a.out ${energytr[$j]} ${mass[$j]} ${lamda[$i]} ${lec[$i]}
#		./a.out ${energy[$j]} ${mass[$j]} ${lamda[$i]} ${lec[$i]}
#	done
#done

echo "Triplet"
./a.out -2.224 140 2 -150
./a.out -2.224 140 4 -500
./a.out -2.224 140 6 -1000
./a.out -2.224 140 8 -1900
./a.out -11.5 510 2 -150
./a.out -11.5 510 4 -400
./a.out -11.5 510 6 -900
./a.out -11.5 510 8 -1500
./a.out -19.5 805 2 -150
./a.out -19.5 805 4 -400
./a.out -19.5 805 6 -800
./a.out -19.5 805 8 -1300

#pion mass = 140 MeV uses a(nn)=22.9/fm
echo "Singlet"
./a.out -7.4 510 2 -150
./a.out -7.4 510 4 -400
./a.out -7.4 510 6 -900
./a.out -7.4 510 8 -1500
./a.out -15.9 805 2 -150
./a.out -15.9 805 4 -400
./a.out -15.9 805 6 -750
./a.out -15.9 805 8 -1300