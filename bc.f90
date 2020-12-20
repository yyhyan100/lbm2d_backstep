subroutine bc()
!call non_eq_extra()
call inlet_vel()
call inner_wall()
call outlet_pressure()
call upper_wall()
call bc_macro()
end subroutine
!----------------------------------------------------------------------------
subroutine inlet_vel() ! Zou-He boundary condition for velocity inlet
	use vars
	i=1
	do j=1,jed
	if (ph(i,j)==2) then
		rho(i,j)=(f(0,i,j)+f(2,i,j)+f(4,i,j)+2*(f(3,i,j)+f(6,i,j)+f(7,i,j)))/(1-u(i,j))
		f(1,i,j)=f(3,i,j)+2*rho(i,j)*u(i,j)/3.0
		f(5,i,j)=f(7,i,j)-0.5*(f(2,i,j)-f(4,i,j))+rho(i,j)*u(i,j)/6.0
		f(8,i,j)=f(6,i,j)+0.5*(f(2,i,j)-f(4,i,j))+rho(i,j)*u(i,j)/6.0
	endif
	enddo
end subroutine
!----------------------------------------------------------------------------
subroutine upper_open() ! open boundary condition
	use vars
	j=jed
	do i=1,ied
		f(4,i,j)=2*f(4,i,j-1)-f(4,i,j-2)
		f(7,i,j)=2*f(7,i,j-1)-f(7,i,j-2)
		f(8,i,j)=2*f(8,i,j-1)-f(8,i,j-2)
!		v(i,j)=0.0
!		u(i,j)=2*u(i,j-1)-u(i,j-2)
	enddo
end subroutine
!----------------------------------------------------------------------------
subroutine upper_wall() ! bounce back boundary condition
	use vars
	j=jed
	do i=1,ied
		f(4,i,j)=f(2,i,j)
		f(7,i,j)=f(5,i,j)
		f(8,i,j)=f(6,i,j)
!		v(i,j)=0.0
!		u(i,j)=0.0
	enddo 
end subroutine
!----------------------------------------------------------------------------
subroutine outlet_pressure() ! pressure boundary condition
	use vars
	i=ied
	do j=2,jed-1
	if (ph(i,j)==2) then
		u(i,j)=(f(0,i,j)+f(2,i,j)+f(4,i,j)+2.0*(f(1,i,j)+f(5,i,j)+f(8,i,j)))/rho_out-1.0
		f(3,i,j)=f(1,i,j)-2*rho_out*u(i,j)/3
		f(6,i,j)=f(8,i,j)-0.5*(f(2,i,j)-f(4,i,j))-rho_out*u(i,j)/6.0
		f(7,i,j)=f(5,i,j)+0.5*(f(2,i,j)-f(4,i,j))-rho_out*u(i,j)/6.0
		v(i,j)=0.0
	endif
	enddo
end subroutine
!----------------------------------------------------------------------------
subroutine inner_wall()  ! bounce back boundary condition for inner bounds
	use vars
	integer i,j,k
	do i=2,ied-1
	do j=2,jed-1
		if(ph(i,j)==1) then
			if(ph(i+1,j)==2) f(1,i,j)=f(3,i,j)
			if(ph(i,j+1)==2) f(2,i,j)=f(4,i,j)
			if(ph(i-1,j)==2) f(3,i,j)=f(1,i,j)
			if(ph(i,j-1)==2) f(4,i,j)=f(2,i,j)
			if(ph(i+1,j+1)==2) f(5,i,j)=f(7,i,j)
			if(ph(i-1,j+1)==2) f(6,i,j)=f(8,i,j)
			if(ph(i-1,j-1)==2) f(7,i,j)=f(5,i,j)
			if(ph(i+1,j-1)==2) f(8,i,j)=f(6,i,j)
		endif
	enddo
	enddo
	j=1
	do i=1,ied
		if(ph(i,j)==1) then
			if(ph(i,j+1)==2) f(2,i,j)=f(4,i,j)
			if(ph(i+1,j+1)==2) f(5,i,j)=f(7,i,j)
			if(ph(i-1,j+1)==2) f(6,i,j)=f(8,i,j)
		endif
	enddo
!	j=jed
!	do i=1,ied
!		if(ph(i,j)==1) then
!			if(ph(i,j-1)==2) f(4,i,j)=f(2,i,j)
!			if(ph(i-1,j-1)==2) f(7,i,j)=f(5,i,j)
!			if(ph(i+1,j-1)==2) f(8,i,j)=f(6,i,j)
!		endif
!	enddo
	i=1
	do j=1,jed
		if(ph(i,j)==1) then
			if(ph(i+1,j)==2) f(1,i,j)=f(3,i,j)
			if(ph(i,j+1)==2) f(2,i,j)=f(4,i,j)
			if(ph(i,j-1)==2) f(4,i,j)=f(2,i,j)
			if(ph(i+1,j+1)==2) f(5,i,j)=f(7,i,j)
			if(ph(i+1,j-1)==2) f(8,i,j)=f(6,i,j)
		endif
	enddo
	i=ied
	do j=1,jed
		if(ph(i,j)==1) then
			if(ph(i,j+1)==2) f(2,i,j)=f(4,i,j)
			if(ph(i-1,j)==2) f(3,i,j)=f(1,i,j)
			if(ph(i,j-1)==2) f(4,i,j)=f(2,i,j)
			if(ph(i-1,j+1)==2) f(6,i,j)=f(8,i,j)
			if(ph(i-1,j-1)==2) f(7,i,j)=f(5,i,j)
		endif
	enddo
end subroutine
