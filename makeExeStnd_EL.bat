echo ***
echo *** Compiling the EcoLight code
del *.obj
 copy ..\common\incfiles_default.for ..\common\incfiles_user.for 
 lf95 -dbl -nco -nlst -nap -ndal -nchk -ntrace -inln -npca -nsav -stchk -o1 -nw -nwo -c ..\common\incfiles_user.for > link_stnd.log
 lf95 -dbl -chk -nco -nlst -pca -nsav -stchk -ntrace -ml winapi -win -nvsw -nw -c ..\common\w*.f90 >> link_stnd.log
 lf95 -dbl -chk -nco -nlst -pca -nsav -stchk -ntrace -ml winapi -win -nvsw -nw -c ..\common\*.f90 >> link_stnd.log
 lf95 -dbl -nco -nlst -nap -ndal -nchk -trace -inln -npca -nsav -stchk -o1 -nw -nwo -c -ml msvb *.f >> link_stnd.log
 lf95 -dbl -co -nlst -nchk -pca -nsav -stchk -nw -nwo -c -ml msvb *.for >> link_stnd.log
echo ***
echo *** Linking the EcoLight code
 lf95 *.obj -ml msvb -lib %WinDir%\system32\HE5info.lib -nomap -winconsole -out mainEL_stnd.exe >> link_stnd.log
echo ***
echo *** Moving the EcoLight executable
 move mainEL_stnd.exe ..
del *.obj
del *.mod