      subroutine readint (NATOMS, x, y, z, title, QENDMMH)
      include 'ndoldim.inc'
c Los siguientes commons son para casos especiales en la subrutina GMETRY
c Si el numero atomico del ultimo atomo es 107, se procesa como un solido y
c es preciso evaluar NUMCAL.Ver la subrutina GMETRY.
      COMMON /KEYWRD/ KEYWRD
      CHARACTER KEYWRD*241, TITLE*162
      common /geos/ geo(3,NATMAX)
      DIMENSION x(*), y(*), z(*),
     &          NA(NATMAX),NB(NATMAX),NC(NATMAX)
      call gettxt (title,QENDMMH)
      call getgeo (geo,na,nb,nc,NATOMS,qaint)
      call geome (NATOMS, x, y, z, na, nb, nc)
      return
      end
