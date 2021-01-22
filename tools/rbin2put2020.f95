program rbin2010

implicit none

integer, parameter:: ldi=4999,nbelx=8,npondcouche=200,isugi=1,inetburn=3,headlength=22
real, parameter:: um=2.3025850929940457d0

integer:: i,ii,l,nwseq,modanf,ialflu,m,nsugi,ierror=0

real(8):: gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,&
        ab,dm_lost,drl,drte,dk,drp,drt,drr,rlp,rlt,rlc,rrp,rrt,&
        rrc,rtp,rtt,rtc,tdiff,suminenv,xltotbeg,dlelexprev,zams_radius,vsuminenv

real(8), dimension(npondcouche):: CorrOmega

real(8), dimension(ldi)::q,p,t,r,s,x,y,xc12,vp,vt,vr,vs,xo16,vx, &
vy,vxc12,vxo16,y3,xc13,xn14,xn15,xo17,xo18,vy3,vxc13,vxn14,&
vxn15,vxo17,vxo18,xne20,xne22,xmg24,xmg25,xmg26,vxne20,vxne22,&
vxmg24,vxmg25,vxmg26,omegi,vomegi,xf19,xne21,xna23,xal26,xal27,&
xsi28,vxf19,vxne21,vxna23,vxal26g,vxal27,vxsi28,xneut,xprot,&
xc14,xf18,xbid,xbid1,vxneut,vxprot,vxc14,vxf18,vxbid,vxbid1

real(8), dimension(nbelx,ldi):: abelx,vabelx


character(120):: fnameB,fnameout
character(256):: LineString,formstring,format_out
!======================================================================
format_out = '(" $INPUT"/" GMS=",d13.6,", ALTER=",1pd22.15,", GLS=",&
    1pd12.5,",  TEFF=",0pf8.0,","/47x,",GLSV=",1pd12.5,", TEFFV=",&
    0pf8.0,","/" DZEITJ=",1pd12.5,", DZEIT=",1pd12.5,", DZEITV=",&
    1pd12.5,","/" SUMMAS=",0pf8.4,", AB=",1pd9.3,", M=",i4,",")'
write(fnameout,'(a14)') 'inimodel.start'
write(*,*) 'Enter b file: '
read(5,*) fnameB
write(*,*)'input .b file: ',fnameB
write(*,*)'output file: ',fnameout

open(12,file=fnameB,form='unformatted')
do while (ierror == 0)
 read(12,iostat=ierror) gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,&
   ab,dm_lost,m,(q(i),p(i),t(i),r(i),s(i),x(i),y(i),xc12(i), &
   vp(i),vt(i),vr(i),vs(i),xo16(i),vx(i),vy(i),vxc12(i),vxo16(i),&
   i=1,m),drl,drte,dk,drp,drt,drr,rlp,rlt,rlc,rrp,rrt,rrc,rtp, &
   rtt,rtc,tdiff,vsuminenv,(CorrOmega(i),i=1,npondcouche),&
   xLtotbeg,dlelexprev,zams_radius
 if (ierror == 0) then
   read(12)(y3(i),xc13(i),xn14(i),xn15(i),xo17(i),xo18(i),vy3(i),&
     vxc13(i),vxn14(i),vxn15(i),vxo17(i),vxo18(i),xne20(i), &
     xne22(i),xmg24(i),xmg25(i),xmg26(i),vxne20(i),vxne22(i),vxmg24(i),&
     vxmg25(i),vxmg26(i),omegi(i),vomegi(i),i=1,m)
   read(12)(xf19(i),xne21(i),xna23(i),xal26(i),xal27(i),xsi28(i),vxf19(i),&
     vxne21(i),vxna23(i),vxal26g(i),vxal27(i),vxsi28(i), xneut(i), &
     xprot(i),xc14(i),xf18(i),xbid(i),xbid1(i),vxneut(i),vxprot(i),&
     vxc14(i),vxf18(i),vxbid(i),vxbid1(i),i=1,m)

   do ii=1,nbelx
    read(12) (abelx(ii,i),vabelx(ii,i),i=1,m)
   enddo

   if (isugi >= 1) then
     read(12) nsugi
   endif
 endif
enddo

write(*,*)'End of reading'

write(formstring,'(i4)')m
formstring='('//trim(formstring)//'(0pf15.10,","))'
write(*,*) trim(formstring)

do i=1,m
 q(i)=q(i)/um
 p(i)=p(i)/um
 t(i)=t(i)/um
 r(i)=r(i)/um
 s(i)=s(i)/um
 vp(i)=vp(i)/um
 vt(i)=vt(i)/um
 vr(i)=vr(i)/um
 vs(i)=vs(i)/um
enddo

! Ecriture du modele initial
open(13,file=fnameout,form='formatted')
write(13,format_out) gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv, &
  xmini,ab,m
write(13,*) ' Q= '
write(13,trim(formstring)) (q(i),i=1,m)
write(13,*) ' P= '
write(13,trim(formstring)) (p(i),i=1,m)
write(13,*) ' T= '
write(13,trim(formstring)) (t(i),i=1,m)
write(13,*) ' R= '
write(13,trim(formstring)) (r(i),i=1,m)
write(13,*) ' S= '
write(13,trim(formstring)) (s(i),i=1,m)

write(13,*) ' VP= '
write(13,trim(formstring)) (vp(i),i=1,m)
write(13,*) ' VT= '
write(13,trim(formstring)) (vt(i),i=1,m)
write(13,*) ' VR= '
write(13,trim(formstring)) (vr(i),i=1,m)
write(13,*) ' VS= '
write(13,trim(formstring)) (vs(i),i=1,m)

write(13,*) ' x= '
write(13,trim(formstring)) (x(i),i=1,m)
write(13,*) ' y3= '
write(13,trim(formstring)) (y3(i),i=1,m)
write(13,*) ' y= '
write(13,trim(formstring)) (y(i),i=1,m)
write(13,*) ' xc12= '
write(13,trim(formstring)) (xc12(i),i=1,m)
write(13,*) ' xc13= '
write(13,trim(formstring)) (xc13(i),i=1,m)
write(13,*) ' xn14= '
write(13,trim(formstring)) (xn14(i),i=1,m)
write(13,*) ' xn15= '
write(13,trim(formstring)) (xn15(i),i=1,m)
write(13,*) ' xo16= '
write(13,trim(formstring)) (xo16(i),i=1,m)
write(13,*) ' xo17= '
write(13,trim(formstring)) (xo17(i),i=1,m)
write(13,*) ' xo18= '
write(13,trim(formstring)) (xo18(i),i=1,m)
write(13,*) ' xne20= '
write(13,trim(formstring)) (xne20(i),i=1,m)
write(13,*) ' xne22= '
write(13,trim(formstring)) (xne22(i),i=1,m)
write(13,*) ' xmg24= '
write(13,trim(formstring)) (xmg24(i),i=1,m)
write(13,*) ' xmg25= '
write(13,trim(formstring)) (xmg25(i),i=1,m)
write(13,*) ' xmg26= '
write(13,trim(formstring)) (xmg26(i),i=1,m)
write(13,*) ' xc14= '
write(13,trim(formstring)) (xc14(i),i=1,m)
write(13,*) ' xf18= '
write(13,trim(formstring)) (xf18(i),i=1,m)
write(13,*) ' xf19= '
write(13,trim(formstring)) (xf19(i),i=1,m)
write(13,*) ' xne21= '
write(13,trim(formstring)) (xne21(i),i=1,m)
write(13,*) ' xna23= '
write(13,trim(formstring)) (xna23(i),i=1,m)
write(13,*) ' xal26= '
write(13,trim(formstring)) (xal26(i),i=1,m)
write(13,*) ' xal27= '
write(13,trim(formstring)) (xal27(i),i=1,m)
write(13,*) ' xsi28= '
write(13,trim(formstring)) (xsi28(i),i=1,m)
write(13,*) ' xneut= '
write(13,trim(formstring)) (xneut(i),i=1,m)
write(13,*) ' xprot= '
write(13,trim(formstring)) (xprot(i),i=1,m)
write(13,*) ' xbid= '
write(13,trim(formstring)) (xbid(i),i=1,m)
write(13,*) ' xbid1= '
write(13,trim(formstring)) (xbid1(i),i=1,m)
write(13,*) ' xabel1= '
write(13,trim(formstring)) (abelx(1,i),i=1,m)
write(13,*) ' xabel2= '
write(13,trim(formstring)) (abelx(2,i),i=1,m)
write(13,*) ' xabel3= '
write(13,trim(formstring)) (abelx(3,i),i=1,m)
write(13,*) ' xabel4= '
write(13,trim(formstring)) (abelx(4,i),i=1,m)
write(13,*) ' xabel5= '
write(13,trim(formstring)) (abelx(5,i),i=1,m)
write(13,*) ' xabel6= '
write(13,trim(formstring)) (abelx(6,i),i=1,m)
write(13,*) ' xabel7= '
write(13,trim(formstring)) (abelx(7,i),i=1,m)
write(13,*) ' xabel8= '
write(13,trim(formstring)) (abelx(8,i),i=1,m)
write(13,*) ' omegi= '
write(13,trim(formstring)) (omegi(i),i=1,m)
write(13,*) ' $END '

close(12)
close(13)
end program rbin2010
