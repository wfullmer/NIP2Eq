#!/usr/local/bin/bash
#
# 	Back up
#
echo '   Backing up files ... '
#
cd ./HOLD/
mv * ../HOLD_OLD/
cd ../
cp *.f90 ./HOLD/
cp *.in *.sh ./HOLD
#
#	clear
#
echo '   Cleaning up old files ... '
#
rm output*.dat
rm *.mod
rm nip2eq.x
#
echo '   Building mod_subs.f90 ... '
cat header.f90 > mod_subs.f90
cat ghost.f90 >> mod_subs.f90
cat grid1d.f90 >> mod_subs.f90
cat init.f90 >> mod_subs.f90
cat update.f90 >> mod_subs.f90
cat bndry.f90 >> mod_subs.f90
cat exact.f90 >> mod_subs.f90
cat source.f90 >> mod_subs.f90
cat fou.f90 >> mod_subs.f90
cat mma.f90 >> mod_subs.f90
cat mmu.f90 >> mod_subs.f90
cat edfou2.f90 >> mod_subs.f90
cat getrhs.f90 >> mod_subs.f90
cat rk3s1.f90 >> mod_subs.f90
cat rk3s2.f90 >> mod_subs.f90
cat rk3s3.f90 >> mod_subs.f90
#
#cat predictor.f90 >> mod_subs.f90
#cat corrector.f90 >> mod_subs.f90
echo 'END MODULE subs ' >> mod_subs.f90
#
echo '   Compiling mod_prec ... '
ifort -c mod_prec.f90
#
echo '   Compiling mod_global ... '
ifort -c mod_global.f90
#
echo '   Compiling mod_inout ... '
ifort -c mod_inout.f90
#
echo '   Compiling mod_subs ... '
ifort -c mod_subs.f90
#
#
#
echo '   Building nip2eq.x ... '
#
ifort -g -CA -CB -CU -o nip2eq.x main.f90 mod_prec.o mod_global.o mod_inout.o mod_subs.o 
#
#
#	Clean up Objects 
#
rm *.o
#
