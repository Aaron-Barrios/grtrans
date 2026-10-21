       program call_grtrans
       implicit none

       character(len=300) :: ifile, ofile
       character(len=40) :: readfile

       namelist /files/ ifile, ofile

       readfile='files.in'
       open(unit=8,file=readfile)
       read(8,nml=files)
       close(unit=8)
       call grtrans_run(ifile,ofile)
       end program
