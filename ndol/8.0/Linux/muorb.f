      subroutine muorb (NA,IA,M)
      include 'ndoldim.inc'
      COMMON /N11/ NAT(NATMAX)
      dimension IA(*), M(*)
c MU es el contador de orbitales atomicos generado in situ
c M(MU) es el identificador del tipo de OA:
c                        1 = s, 2 = px, 3 = py, 4 = pz
c IA(MU) es el numero del atomo al que pertenece el OA MU
      MU = 1
      do I=1,NA
        M(MU) = 1
        IA(MU) = I
        MU = MU + 1
        if (NAT(I).gt.2) then
          M(MU) = 2
          M(MU+1) = 3
          M(MU+2) = 4
          IA(MU) = I
          IA(MU+1) = I
          IA(MU+2) = I
          MU = MU + 3
        endif
      enddo
      return
      end
      
