	module array	
	  complex*16,allocatable::H(:,:),work(:)
	  double precision, allocatable::rwork(:),evl(:)
	  integer,allocatable::px(:),py(:),pz(:)
	  double precision, allocatable::ion_pot(:),n_total(:)
	  double precision, allocatable::m_s(:),th_s(:),ph_s(:)
	  double precision, allocatable::op_fl(:,:),dos(:),ttemp(:),mu(:)
	end module array
	
	module global
	  integer::d,dc,Tgrid_max,MCSW,intrval
	  double precision::t1,t2,U1,filling,gama,strnth,ds
	  double precision::pi,T,mu_sys
	end module global
	
	program  dos_and_energy
	  use array
	  use global
	  implicit none
	  integer::l,x,y,z,i
	  integer::temp,count1,config
	  integer::mc_count,p,q,iw
	  double precision::get_mu,mc_mu,w,dw
	  double precision::inp
	  double precision::fs,sys_energy,mc_energy
	  character(len=8 )::grnsfn
	  
	  open (21,file='dos_input.dat',status='unknown')
	  do i=1,11
	    read (21,*)inp
	    if (i.eq.1)d=int(inp)
	    if (i.eq.2)dc=int(inp)
	    if (i.eq.3)t1=dble(inp)
	    if (i.eq.4)t2=dble(inp)
	    if (i.eq.5)Tgrid_max=int(inp)
	    if (i.eq.6)MCSW=int(inp)
	    if (i.eq.7)intrval=int(inp)
	    if (i.eq.8)U1=dble(inp)
	    if (i.eq.9)filling=dble(inp)
	    if (i.eq.10)gama=dble(inp)
	    if (i.eq.11)strnth=dble(inp)
	  enddo

	  
	  Print*,"System size,  d=", d
	  Print*,"cluster size, dc=", dc
	  Print*,"NN hopping,   t1=", t1
	  Print*,"NNN hopping,  t2=", t2
	  Print*,"Total temperature grid,   Tgrid_max=", Tgrid_max
	  Print*,"Total system sweep,  MCSW=", MCSW
	  Print*,"Left system sweep, intrval=", intrval
	  print*,"the interaction value, U1=", U1
	  print*,"filling in the system=", filling
	  print*,"the broadening of the lorentzian for cal. of DOS", gama
	  print*,"strength of the ionic potential=", strnth  
	  
	  ds=filling*d**3 
	  config=(MCSW/(2*intrval))+1	  
	  allocate(px(d**3),py(d**3),pz(d**3),H(2*d**3,2*d**3))
	  allocate(evl(2*d**3),work(2*(2*d**3)-1),rwork(3*(2*d**3)-2))
	  allocate(m_s(d**3),th_s(d**3),ph_s(d**3))
	  allocate(ion_pot(d**3),n_total(d**3))
	  allocate(op_fl(Tgrid_max*config*d**3,5),dos(2000),ttemp(Tgrid_max),mu(Tgrid_max))
	  12 format('fort.',I3, '')
	  
      pi=acos(-1.0d0)

	 do l=1,d**3
	   do z=1,d
	     do y=1,d
	       do x=1,d
	         if((((z-1)*d**2)+d*(y-1)+x).eq.l)then
	           px(l)=x
	           py(l)=y
	           pz(l)=z
	         endif
	       enddo
	     enddo
	   enddo
	 enddo          

	  l=1
	  do z=1,d
        do y=1,d
          do x=1,d
	        ion_pot(l)=strnth*(-1.0d0)**(x+y+z)
	         l=l+1
	      enddo
	    enddo
	  enddo

	 
       do temp=1,Tgrid_max
         write(grnsfn,12) 300+temp
         open(unit=7,file=grnsfn,status='unknown')
         do q=1,config*d**3
           read(7,*)op_fl(q+((temp-1)*config*d**3),1:5)
         enddo
         close(7)
       enddo
       
       open(unit=27,file="fort.18",status="unknown")
         do temp=1,Tgrid_max
           read(27,*)ttemp(temp),mu(temp)
         enddo
       close(27)
       
!       do temp=1,Tgrid_max
!         print*,ttemp(temp),mu(temp)
 !      enddo


	T=1.10d0

	do temp=1,Tgrid_max
      if(temp.le.7)T=T-0.10d0
      if((temp.gt.7).and.(temp.le.19))T=T-0.025d0
      if((temp.gt.19).and.(temp.le.28))T=T-0.01d0
      if(temp.gt.28)T=0.005d0
      print*,temp,T
      
      mc_mu=0.0d0
      mc_energy=0.0d0
      dos=0.0d0   
      mc_count=0
	  do count1=1,config
	    mc_count=mc_count+1
1919    format (1x,i4,4f16.8)
	    do l=1,d**3
          p=(temp-1)*config*d**3+(count1-1)*d**3+l
          m_s(l)=op_fl(p,2)
          th_s(l)=op_fl(p,3)
          ph_s(l)=op_fl(p,4)
          n_total(l)=op_fl(p,5)
        enddo   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! site 
        call system_matrix_gen
 !       mu_sys=get_mu(dble(ds))
         mu_sys=mu(temp)
        call energy(mu_sys,sys_energy)
        call cal_dos
        
        mc_energy=mc_energy+sys_energy
  !      mc_mu=mc_mu+mu_sys
        
      enddo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! mc 
      
 !     fs=0.0d0
!	  do i=1,2*d**3
!	    fs=fs+(1.0d0/(exp((evl(i)-(mc_mu/dble(mc_count)))/T)+1.0d0))
!	  end do

      fs=0.0d0
	  do i=1,2*d**3
	    fs=fs+(1.0d0/(exp((evl(i)-mu_sys)/T)+1.0d0))
	  end do
	
      write(23,*)T,mc_energy/dble(mc_count),fs,mu_sys
      flush(23)
      
!      write(24,*)T,mc_mu/dble(mc_count),fs
!      flush(24)
      
      w=-10.0d0
      dw=0.02
      do iw=1,1000
        w=w+dw
        write(800+temp,*)w,dos(iw)/dble(mc_count)
        flush(800+temp)     
      enddo
        
    enddo  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! temp        
        
        
        
    end
!....................................................................
!subroutine for system matrix generation
!....................................................................        
    subroutine system_matrix_gen
    use array
    use global
	implicit none
	integer::l,k,x,y,z,xi,yi,zi,xd,yd,zd,a,b,i,j,info

	H=cmplx(0.0d0,0.0d0)
		
	do l=1,d**3
	  xi=1
	  xd=-1
	  yi=1
	  yd=-1
	  zi=1
	  zd=-1
	  
	  if(px(l).eq.d)xi=1-d
	  if(py(l).eq.d)yi=1-d
	  if(pz(l).eq.d)zi=1-d
	  
	  if(px(l).eq.1)xd=d-1
	  if(py(l).eq.1)yd=d-1
	  if(pz(l).eq.1)zd=d-1
	  do k=1,d**3
	    if(l.eq.k)then
	      a=2*l-1
	      b=2*k-1
	      H(a,b)=(-abs(m_s(k))*cos(th_s(k))*(U1/2.0d0))+((n_total(k))*U1/2.0d0)+ ion_pot(k)
	      H(a+1,b+1)=(abs(m_s(k))*cos(th_s(k))*(U1/2.0d0))+((n_total(k))*U1/2.0d0)+ ion_pot(k)
	      H(a,b+1)=-abs(m_s(k))*sin(th_s(k))*cmplx(cos(ph_s(k)),-sin(ph_s(k)))*(U1/2.0d0)
	      H(a+1,b)=conjg(H(a,b+1))    
	    endif
	    if(((px(k).eq.px(l)+xi).and.(py(k).eq.py(l)).and.(pz(k).eq.pz(l))).or.&
	      &((px(k).eq.px(l)+xd).and.(py(k).eq.py(l)).and.(pz(k).eq.pz(l))))then
	      a=2*l-1
	      b=2*k-1
	      H(a,b)=t1
	      H(a+1,b+1)=t1
	    endif
	    if(((px(k).eq.px(l)).and.(py(k).eq.py(l)+yi).and.(pz(k).eq.pz(l))).or.&
	      &((px(k).eq.px(l)).and.(py(k).eq.py(l)+yd).and.(pz(k).eq.pz(l))))then
	      a=2*l-1
	      b=2*k-1
	      H(a,b)=t1
	      H(a+1,b+1)=t1
	    endif
	    if(((px(k).eq.px(l)).and.(py(k).eq.py(l)).and.(pz(k).eq.pz(l)+zi)).or.&
	      &((px(k).eq.px(l)).and.(py(k).eq.py(l)).and.(pz(k).eq.pz(l)+zd)))then
	      a=2*l-1
	      b=2*k-1
	      H(a,b)=t1
	      H(a+1,b+1)=t1
	    endif
	    
	    if(((px(k).eq.px(l)+xi).and.(py(k).eq.py(l)+yi).and.pz(k).eq.pz(l)).or.&
	      &((px(k).eq.px(l)+xd).and.(py(k).eq.py(l)+yd).and.pz(k).eq.pz(l)))then
	      a=2*l-1
	      b=2*k-1
	      H(a,b)=t2
	      H(a+1,b+1)=t2
	    endif
	    if(((px(k).eq.px(l)+xi).and.(py(k).eq.py(l)).and.pz(k).eq.pz(l)+zi).or.&
	      &((px(k).eq.px(l)+xd).and.(py(k).eq.py(l)).and.pz(k).eq.pz(l)+zd))then
	      a=2*l-1
	      b=2*k-1
	      H(a,b)=t2
	      H(a+1,b+1)=t2
	    endif	
	    if(((px(k).eq.px(l)).and.(py(k).eq.py(l)+yi).and.pz(k).eq.pz(l)+zi).or.&
	      &((px(k).eq.px(l)).and.(py(k).eq.py(l)+yd).and.pz(k).eq.pz(l)+zd))then
	      a=2*l-1
	      b=2*k-1
	      H(a,b)=t2
	      H(a+1,b+1)=t2
	    endif    
	  enddo	    
	enddo
        
    do i=1,2*d**3
      do j=1,2*d**3
        write(25,*)i,j,real(H(i,j))!,j=1,2*d**3)
      enddo       
    enddo
    
    call zheev ('N','U',2*d**3,H,2*d**3,evl,work,2*(2*d**3)-1,rwork,info)
    if(info.ne.0) print*,info
  
    do i=1,2*d**3
      write(26,*)i,evl(i)
    enddo
        
	end
	
!.......................................................................
!calculation of chemical potential for the system
!.......................................................................
    double precision function get_mu(fill)
    use array
    use global
	implicit none
	double precision f0, f, fL2, fR, mR, mL, rtmp,m_d
	integer i
	double precision fill
	
	mR = maxval(evl)       !right-side chemical potential
	fr=0.0d0
	do i=1,2*d**3
	  fr=fr+(1.0d0/(exp((evl(i)-mR)/T)+1.0d0))
	end do
	
    mL = minval(evl)       !left-side chemical potential
	fL2=0.0d0
 	do i=1,2*d**3
 	  fL2=fL2+(1.0d0/(exp((evl(i)-mL)/T)+1.0d0))
	end do
	
	m_d = 0.5d0*(mL+mR)    !middle chemical potential
	f=0.0d0
	do i=1,2*d**3
	  f=f+(1.0d0/(exp((evl(i)-m_d)/T)+1.0d0))
	end do
	
	!print*,f,fill
	do while(abs(f-fill).ge.1e-8)
	  m_d = 0.5d0*(mL+mR)
	  f=0.0d0
	  do i=1,2*d**3
	    f=f+(1.0d0/(exp((evl(i)-m_d)/T)+1.0d0))
	  end do
	  if(f.gt.fill)then
	    !if middle filling is above target, make it the new right bound.
	    mR = m_d
	    fR = f
	 elseif(f.lt.fill)then
	   !if middle filling is below target, make it the new left bound.
	   mL = m_d
	   fR = f
	 endif
	enddo
	!Return the middle value
	get_mu = m_d
	return
	end function get_mu
	

!.....................................................................
!subroutine for calculation of energy
!.....................................................................
	subroutine energy(EF,Ein)
	use array
	use global
	implicit none
	integer:: i
	double precision ::sum,sum1,Ein,EF   	
        
   	sum=0.0d0
    do i=1,2*d**3
   	  sum=sum+evl(i)*0.5d0*(1.0d0+tanh((EF-evl(i))/(2.0d0*T)))
   	  !/((exp(evl(i)-EF)/T)+1.0d0)
   	enddo 	
   	sum1=0.0d0
  	do i=1,d**3
   	  sum1=sum1+(m_s(i)**2)*(U1/4.0d0)-((n_total(i))**2)*(U1/4.0d0)
  	enddo
  	Ein=sum+sum1
	return
	end
	
	
!......................................................................
!calculation of DOS
!......................................................................
	subroutine cal_dos
	use array
	use global
	implicit none
	integer::iw,i
	double precision::w,dw,dos_sum
    w=-10.0d0
    dw=0.02
    do iw=1,1000
      w=w+dw
      dos_sum=0.0d0
      do i=1,2*d**3       
        dos_sum=dos_sum+((gama/pi)/((w-(evl(i)-mu_sys))**2+gama**2))
      enddo
      dos(iw)=dos(iw)+(dos_sum/dble(2*d**3))
    enddo    
    return
    end
