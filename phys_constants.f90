      module phys_constants
      
      ! CGS, CODATA 2018 / IAU 2015; the same values as the pipeline's
      ! constants.py (change both together).  msun was 1.998d33 upstream, a
      ! typo for 1.989d33 that made r_g, and so every path length, 0.45% too
      ! large (history.tex, Stage 21).
      real(kind=8), parameter :: h=6.62607015d-27,k=1.380649d-16, &
       c=2.99792458d10,e=4.80320471d-10,G=6.67430d-8, &
       m=9.1093837015d-28,mp=1.67262192369d-24,pi=acos(-1d0),c2=c*c,sigb=5.670374d-5, &
       msun=1.98841d33, sigt=6.6524587d-25

      contains
        
        subroutine phys_constants_dummy()
        end subroutine phys_constants_dummy

      end module
