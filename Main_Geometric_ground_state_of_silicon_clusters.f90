program silicon
	!Author: @Philipp Bronner, 12.2020
	!This programm tries to find the energedicall lowest constellation of a silicon cluster.
	!The silicon atoms interact with a Bazant force field and for movement a velocity verlet algorithm is used.
	!The optimiztion is dome with the methode of of simulated annealing and furthermore every 100 steps a
	!steepest descent optimiztion with energy feedback is applied.

	! Initial Variables
	Implicit none
	Integer :: iat, jj,ll, nat, kk = 1
	real*8 :: t1, t2, t3
	real*8, dimension(3,16) :: rxyz= 0, gxyz, fxyz, vxyz, particles
	real*8, dimension(3) :: alat
	real*8 :: pener, coord, ener_var, coord_var, count
	real*8 :: ref_ekin, act_ekin, etot, etot_sd
	real*8, dimension(3,16,10) :: minpos_rxyz 
	real*8, dimension(10) :: minpos_E
	
	Character(10) :: filename
	Character(1) :: fn
	real*8 :: dt = 5.d-2 ! Time step of the verlocity verlet algorythm
	real*8, dimension(10) :: x 
	
	! Dummy list
	minpos_E = (/10,9,8,7,6,5,4,3,2,1 /)
	! Import and create Files for storing the data
	open(unit = 11, file = 'initial_config.ascii')
	open(unit=22,file='minilal_energies.dat',status='unknown')
	open(unit=44,file='actual_energies.dat',status='unknown')
	open(unit=33,file='tenminimalenergies.dat',status='unknown')
	
	! dummy readout
	read(11,*) nat
	read(11,*) t1, t2, t3 
	read(11,*) t1, t2, t3
	
	! create initial line of position file
	write(33,*) nat, 'Positions'
	write(33,*) 0.d0, 0.d0, 0.d0
	write(33,*) 25.d0, 25.d0, 25.d0
	
	! set box size
	alat = (/25.d0, 25.d0, 25.d0/)
	! read initial position
	do iat = 1,nat
		read(11,*) (rxyz(jj,iat), jj = 1,3)
	end do
	
	!First forces
	!Routine with bazant force field is called
	call bazant(nat,alat,rxyz,gxyz,pener, coord, ener_var, coord_var,count)
	
	!set velocities according to the maxwell distribution
	call gausdist(nat,vxyz)
	ref_ekin = 5.d0 ! Initial temperature
	call velnorm(nat,ref_ekin,vxyz)

	
	! Begin calculation. The algorythm is applied 
	!initialization of a loop
1000 continue 
	ref_ekin = ref_ekin *0.99995d0 ! Decreasing factor of the Temperature
	if (ref_ekin .le. 1.d-2) goto 2000 !end condition
	
	!Apply a verlocity verlet step
	rxyz = rxyz + dt*vxyz + dt**2/2*gxyz  ! Update positions
	call bazant(nat,alat,rxyz,fxyz,pener,coord,ener_var,coord_var,count) ! Update force
	vxyz = vxyz + (fxyz+gxyz)/2*dt ! Update velocity
	
	act_ekin = 0.5d0*sum(vxyz**2)/(3*nat-6) ! Acuall kinetic energy
	
	! Store total energy
	etot = act_ekin + pener 
	write(44,*) etot
	
	! Compare with reference energy
	!Velocities are increased or decreased dependent
	!on an incresig or decreasing kinetic energy
	if (act_ekin .gt. ref_ekin) then
		vxyz = vxyz*0.99d0
	else
		vxyz = vxyz*1.01d0
	end if
	
	
	!Every 100 steps, application of a steepest descent optimiztion
	if (mod(kk,100) == 0) then
		particles = rxyz
		!Call of steepest descent routine with energy feedback.
		call energyfeedback(particles,1.d-2,alat,etot_sd)
		! store minimal energy after stepest descent
		write(22,*) etot_sd
		call optimalconstelation(etot_sd, particles,minpos_rxyz,minpos_E)
	end if
	if (mod(kk,1000) == 0) then
		print*, kk, etot_sd, ref_ekin
	end if
	kk = kk+ 1
	
	gxyz = fxyz 
	goto 1000
	
2000 continue

	! Store ten minimal constellations
	do kk = 0,9
		write(33,*) minpos_E(kk+1) !Write down the ten lowest energies

		! Write down the contributed positions
		write(fn,'(i1)') kk
		filename = 'pos'//fn//'.ascii'
		open(unit = 44, file = filename, status = 'unknown')
		write(44,*) 16, minpos_E(kk+1)
		write(44,*) 0.d0, 0.d0, 0.d0
		write(44,*) 25.d0, 25.d0, 25.d0
		do ll = 1,nat
			write(44,*) minpos_rxyz(1,ll,kk), minpos_rxyz(2,ll,kk), minpos_rxyz(3,ll,kk)
		end do
	end do
	
end program


! This subroutine findes the 10 lowest from a list 
!and overwrites the highes entry if the new entry is lower.
subroutine optimalconstelation(epot, rxyz, minpos_rxyz, minpos_E)
	implicit none
	integer :: ii
	real*8 :: epot
	real*8, dimension(1) :: x
	real*8, dimension(3,16) :: rxyz
	real*8, dimension(3,16,10) :: minpos_rxyz 
	real*8, dimension(10) :: minpos_E
	epot = float(nint(epot))
	
	if (epot < maxval(minpos_E)) then
		if (Any (minpos_E == epot)) then
		
		else
			x = int(maxloc(minpos_E))
			ii = int(x(1))
			minpos_E(ii) = epot
			minpos_rxyz(:,:,ii) = rxyz
		end if
		
	end if	
end

! This subroutine does an enrgy feedback steepest descent minimization starting from an given constellation
subroutine energyfeedback(particles,a,alat, etot)
	implicit none
	real*8 :: etotpre, norm, etot, abar,a
	real*8 :: coord_var, ener_var, coord, count
	integer :: n
    real(8), dimension(3,16) :: particles, fxyz, fxyz2
	real*8, dimension(3) :: alat
	n = 0
	etotpre = 0
	etot = 0
	abar = a
	call bazant(16,alat,particles,fxyz,etot,coord,ener_var,coord_var,count)
	
	norm = sqrt(sum(fxyz**2))
    do while (norm > 1.d-3)
		if (n == 0) then
			abar = abar
		else if ((n > 0) .and. (etot > etotpre)) then
			abar = abar/2
		else if ((n > 0) .and. (etot <= etotpre)) then
			abar = abar*1.05
		end if
		etotpre = etot
        particles = particles + abar*fxyz
		call bazant(16,alat,particles,fxyz,etot,coord,ener_var,coord_var,count)
		norm = sqrt(sum(fxyz**2))
        n = n+1
    end do
end