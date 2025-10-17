      j = 0
      igm = ( mmmax - 1 ) * 2 + 1

c region name
      if(nmnr_flag) then
        idum = 0
        do im = 1, nm
          do ir = 1, nr
            idum = idum + 1

            j = j + 1
            ntrn = kr(j)
            j = j + 1
            mtrn = kr(j)
            call echrg2(mtrn,kr(j+1),dum1_n(idum),lng1_n(idum),
     &                  icmb,igm)
            j = j + mtrn
            j = j + 1
            ntrn = kr(j)
            j = j + 1
            mtrn = kr(j)
            call echrg2(mtrn,kr(j+1),dum2_n(idum),lng2_n(idum),
     &                  icmb,igm)
            j = j + mtrn
          end do
        end do
      end if
