diff --git a/code/inputparam.f90 b/code/inputparam.f90
index e342931..38c54ec 100644
--- a/code/inputparam.f90
+++ b/code/inputparam.f90
@@ -14,7 +14,7 @@ module inputparam
   integer,parameter:: imagn_default=0,ianiso_default=0,ipop3_default=0,ibasnet_default=0,iopac_default=3,&
     ikappa_default=5,istati_default=0,igamma_default=0,nndr_default=1,iledou_default=0,idifcon_default=0,&
     iunder_default=0,nbchx_default=200,nrband_default=1,icncst_default=0,iprn_default=99,&
-    iout_default=0,itmin_default=5,idebug_default=0,itests_default=0
+    iout_default=0,itmin_default=5,idebug_default=0,itests_default=0 , idcirch_default=1
   real(kindreal),parameter:: fenerg_default=1.0d0,richac_default=1.0d0,zsol_default=1.40d-2,frein_default=0.0d0,&
     K_Kawaler_default=0.d0,Omega_saturation_default=14.d0,vwant_default=0.0d0,xfom_default=1.0d0, &
     dunder_default=0.0d0,dgro_default=0.010d0,dgr20_default=0.010d0,binm2_default=0.d0,periodini_default=0.d0,&
@@ -52,14 +52,15 @@ module inputparam
 !-----------------------------------------------------------------------
 
 ! **** Rotation-linked parameters
-  integer,save:: idiff,iadvec,istati=istati_default,icoeff,igamma=igamma_default,idialo,idialu
+  integer,save:: idiff,iadvec,istati=istati_default,icoeff,igamma=igamma_default,idialo,idialu, &
+  idcirch=idcirch_default
   real(kindreal),save:: fenerg=fenerg_default,richac=richac_default,frein=frein_default,K_Kawaler=K_Kawaler_default, &
                           Omega_saturation=Omega_saturation_default,rapcrilim,vwant=vwant_default,&
                           xfom=xfom_default,omega,xdial,B_initial=B_initial_default,add_diff=add_diff_default
   logical,save:: Add_Flux=Add_Flux_default,diff_only=diff_only_default
 !-----------------------------------------------------------------------
   namelist /RotationParams/idiff,iadvec,istati,icoeff,fenerg,richac,igamma,frein,K_Kawaler,Omega_saturation,rapcrilim, &
-                           vwant,xfom,omega,xdial,idialo,idialu,Add_Flux,diff_only,B_initial,add_diff
+                           vwant,xfom,omega,xdial,idialo,idialu,Add_Flux,diff_only,B_initial,add_diff,idcirch
 !-----------------------------------------------------------------------
 
 ! **** Surface parameters
@@ -111,7 +112,7 @@ module inputparam
     icncst_default,iprn_default,iout_default,itmin_default,fenerg_default,richac_default,zsol_default, &
     frein_default,K_Kawaler_default,Omega_saturation_default,vwant_default,xfom_default,dunder_default,dgr20_default, &
     xyfiles_default,idebug_default,bintide_default,binm2_default,periodini_default,const_per_default, &
-    var_rates_default,verbose_default,stop_deg_default
+    var_rates_default,verbose_default,stop_deg_default,idcirch_default
 
 contains
 !=======================================================================
@@ -225,6 +226,7 @@ subroutine Write_namelist(Unit,nwseqnew,modanfnew,nzmodnew,xcnwant)
   call Write_param(Unit,"diff_only=",diff_only,diff_only_default)
   call Write_param(Unit,"B_initial=",B_initial,B_initial_default)
   call Write_param(Unit,"add_diff=",add_diff,add_diff_default)
+  call Write_param(Unit,"idcirch=" , idcirch, idcirch_default)
   write(Unit,'("&END"/)')
 
   write(Unit,'(a)') "&SurfaceParams"
