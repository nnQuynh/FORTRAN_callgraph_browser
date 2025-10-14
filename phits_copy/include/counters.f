!-------------------------------counters.f------------------------------
! Version: 051227-1200
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
      integer*8 iannih,iaphi,ibhabha,ibrems,icollis,
     z icompt,iedgbin,ieii,ielectr,ihardx,
     z ihatch,ikauger,ikshell,ikxray,ilauger,
     z ilshell,ilxray,imoller,imscat,ipair,
     z iphoto,iphoton, iraylei,ishower,iuphi,
     z itmxs,noscat,iblock
      common/COUNTERS/                       ! Subroutine-entry counters
     * iannih,
     * iaphi,
     * ibhabha,
     * ibrems,
     * icollis,
     * icompt,
     * iedgbin,
     * ieii,
     * ielectr,
     * ihardx,
     * ihatch,
     * ikauger,
     * ikshell,
     * ikxray,
     * ilauger,
     * ilshell,
     * ilxray,
     * imoller,
     * imscat,
     * ipair,
     * iphoto,
     * iphoton,
     * iraylei,
     * ishower,
     * iuphi,
     * itmxs,
     * noscat,
     * iblock
!$OMP THREADPRIVATE(/COUNTERS/)
      common/COUNTERS_sum/
     z iannih_sum, iaphi_sum, ibhabha_sum, ibrems_sum, icollis_sum,
     z icompt_sum, iedgbin_sum, ieii_sum, ielectr_sum, ihardx_sum,
     z ihatch_sum, ikauger_sum, ikshell_sum, ikxray_sum, ilauger_sum,
     z ilshell_sum, ilxray_sum, imoller_sum, imscat_sum, ipair_sum,
     z iphoto_sum, iphoton_sum, iraylei_sum, ishower_sum, iuphi_sum,
     z itmxs_sum, noscat_sum, iblock_sum
c
      integer*8
     z iannih_sum,iaphi_sum,ibhabha_sum,ibrems_sum,icollis_sum,
     z icompt_sum,iedgbin_sum,ieii_sum,ielectr_sum,ihardx_sum,
     z ihatch_sum,ikauger_sum,ikshell_sum,ikxray_sum,ilauger_sum,
     z ilshell_sum,ilxray_sum,imoller_sum,imscat_sum,ipair_sum,
     z iphoto_sum,iphoton_sum,iraylei_sum,ishower_sum,iuphi_sum,
     z itmxs_sum,noscat_sum,iblock_sum
c
!-----------------------last line of counters.f-------------------------
