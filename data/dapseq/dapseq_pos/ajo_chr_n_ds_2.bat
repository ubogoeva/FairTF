del *.sta
del *.overlap
call ajo_chr_n_ds.bat  %1 1
call ajo_chr_n_ds.bat  %1 2
call ajo_chr_n_ds.bat  %1 3
call ajo_chr_n_ds.bat  %1 4
call ajo_chr_n_ds.bat  %1 5
call cp3.bat %1
call cp4.bat %1