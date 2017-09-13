!********************************************************************************
!   Cassandra - An open source atomistic Monte Carlo software package
!   developed at the University of Notre Dame.
!   http://cassandra.nd.edu
!   Prof. Edward Maginn <ed@nd.edu>
!   Copyright (2013) University of Notre Dame du Lac
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*******************************************************************************

subroutine virialMC_Driver

  USE Run_Variables
  USE Random_Generators
  USE File_Names
  USE Energy_Routines
  USE Cluster_Routines
  USE Read_Write_Checkpoint
  USE Rotation_Routines
  USE io_utilities

  IMPLICIT NONE

  INTEGER :: idist, iconf, jconf, is, js, im, jm, ifrag
  INTEGER :: ndist, irot, jrot, n_pos
  REAL(DP) :: en_min, dist, b2, e_min, e_max
  LOGICAL :: overlap
  CHARACTER(240) :: virialfile

  ! local fragment declarations
  INTEGER :: itotal_frags, jtotal_frags  ! total number of conformations for this fragment in the library
  INTEGER :: this_fragment, frag_start   ! random fragment to start growing from
  INTEGER :: ifrag_type, jfrag_type, this_atom
  REAL(DP):: total_n

  ndist = INT( (mcvirial%max_dist - mcvirial%min_dist) / mcvirial%dist_step ) + 1

  ALLOCATE( mcvirial%effective(ndist) )
  ALLOCATE( mcvirial%coefficient(ndist) )
  mcvirial%effective = 0.0
  mcvirial%coefficient = 0.0
  !second_Bn = 0

  CALL Name_Files(run_name,'.virial',virialfile)
  OPEN(991,file=virialfile)
  write(991,"(A)") '# distance   effective_potential (kJ/mol) virial_coefficient '

  energy(1)%total = 0.0
  energy(1)%intra = 0.0_DP
  energy(1)%bond  = 0.0_DP
  energy(1)%angle = 0.0_DP
  energy(1)%dihedral = 0.0_DP
  energy(1)%improper = 0.0_DP
  energy(1)%intra_vdw = 0.0_DP
  energy(1)%intra_q = 0.0_DP
  energy(1)%erf_self = 0.0_DP

  en_min = 10.0

  total_n = DBLE(mcvirial%nrotations(1) * mcvirial%nrotations(2) * mcvirial%nconfs(1) * mcvirial%nconfs(2))

  im = 1
  jm = 1
  is = mcvirial%species(1)
  js = mcvirial%species(2)
  frag_start = 1
  ifrag_type = frag_list(frag_start,is)%type
  jfrag_type = frag_list(frag_start,js)%type
  itotal_frags = frag_list(frag_start,is)%nconfig
  jtotal_frags = frag_list(frag_start,js)%nconfig

  ! calculate at different separation distances
  dist = mcvirial%max_dist
  dist_loop: DO idist = 1, ndist
    
    e_min =  100000000.0
    e_max = -100000000.0
    n_pos = 0
    ! sample different configurations of molecule 1
    iconf_loop: DO iconf = 1, mcvirial%nconfs(1)
        
        ! Pull from the fragment reservoir with uniform probability
        this_fragment = INT(rranf() * itotal_frags) + 1
    
        ! Read the coordinates for every atom
        DO ifrag = 1, frag_list(frag_start,is)%natoms 
            this_atom = frag_list(frag_start,is)%atoms(ifrag)
            
            atom_list(this_atom,im,is)%rxp = frag_coords(ifrag,this_fragment,ifrag_type)%rxp
            atom_list(this_atom,im,is)%ryp = frag_coords(ifrag,this_fragment,ifrag_type)%ryp
            atom_list(this_atom,im,is)%rzp = frag_coords(ifrag,this_fragment,ifrag_type)%rzp
    
            atom_list(this_atom,im,is)%exist = .TRUE.
        END DO
    
        ! Turn on the molecule and its individual atoms
        molecule_list(im,is)%which_box = 1
        molecule_list(im,is)%molecule_type = int_normal
        molecule_list(im,is)%live = .TRUE.
        ! Recompute the COM in case the molecule was wrapped
        CALL Get_COM(im,is)
    
        ! Compute the distance of the atom farthest from COM
        CALL Compute_Max_COM_Distance(im,is)
        atom_list(:,im,is)%rxp = atom_list(:,im,is)%rxp - molecule_list(im,is)%xcom
        atom_list(:,im,is)%ryp = atom_list(:,im,is)%ryp - molecule_list(im,is)%ycom
        atom_list(:,im,is)%rzp = atom_list(:,im,is)%rzp - molecule_list(im,is)%zcom

        ! rotate species 1
        DO irot = 1, mcvirial%nrotations(1)
            ! sample different configurations of molecule 2
            DO jconf = 1, mcvirial%nconfs(2)
        
                ! Pull from the fragment reservoir with uniform probability
                this_fragment = INT(rranf() * itotal_frags) + 1
            
                ! Read the coordinates for every atom
                DO ifrag = 1, frag_list(frag_start,js)%natoms 
                    this_atom = frag_list(frag_start,js)%atoms(ifrag)
                    
                    atom_list(this_atom,jm,js)%rxp = frag_coords(ifrag,this_fragment,jfrag_type)%rxp
                    atom_list(this_atom,jm,js)%ryp = frag_coords(ifrag,this_fragment,jfrag_type)%ryp
                    atom_list(this_atom,jm,js)%rzp = frag_coords(ifrag,this_fragment,jfrag_type)%rzp
            
                    atom_list(this_atom,jm,js)%exist = .TRUE.
                END DO
            
                CALL Get_COM(jm,js)
                atom_list(:,jm,js)%rxp = atom_list(:,jm,js)%rxp - molecule_list(jm,js)%xcom
                atom_list(:,jm,js)%ryp = atom_list(:,jm,js)%ryp - molecule_list(jm,js)%ycom
                atom_list(:,jm,js)%rzp = atom_list(:,jm,js)%rzp - molecule_list(jm,js)%zcom

                ! Turn on the molecule and its individual atoms
                molecule_list(jm,js)%which_box = 1
                molecule_list(jm,js)%molecule_type = int_normal
                molecule_list(jm,js)%live = .TRUE.
            
                molecule_list(jm,js)%xcom = dist
                molecule_list(jm,js)%ycom = 0.0
                molecule_list(jm,js)%zcom = 0.0
         
                atom_list(:,jm,js)%rxp = atom_list(:,jm,js)%rxp + molecule_list(jm,js)%xcom
            
                DO jrot = 1, mcvirial%nrotations(2)
                    ! Recompute the COM in case the molecule was wrapped
                    !CALL Get_COM(im,is)
                
                    ! Compute the distance of the atom farthest from COM
                    !CALL Compute_Max_COM_Distance(im,is)
             
                    ! do not compute intramolecular energy for B2, only use it for the configuration generation
                    ! See Harismiadis & Szleifer, Mol. Phys., 1993
                    CALL Compute_Total_System_Energy(1,.false.,overlap)
                    
    !write(*,"(5F15.3)") dist, energy(1)%total, energy(1)%intra, energy(1)%intra_vdw, energy(1)%inter_vdw
                    mcvirial%coefficient(idist) = mcvirial%coefficient(idist) + dexp(-1.0_DP * energy(1)%total * beta(1))
                    IF (overlap .eqv. .FALSE.) THEN
                        mcvirial%effective(idist) = mcvirial%effective(idist) + energy(1)%total
                        if (e_max < energy(1)%total) then
                            e_max = energy(1)%total
                            CALL Write_Coords(1)
                        endif
                        if (e_min > energy(1)%total) e_min = energy(1)%total
                        n_pos = n_pos+1
                    !ELSE
                    !    WRITE(logunit,*) 'Overlaping atoms at', dist, 'stopping computation'
                    !    EXIT dist_loop
                    END IF
    
                    CALL Rotate_Molecule_Eulerian(1,js)
                
                END DO
        
                IF (energy(1)%intra_vdw > 0 ) THEN
                    CALL Write_Coords(1)
                    Exit iconf_loop
                END IF
            END DO
        END DO
    END DO iconf_loop
    
    write(991,"(F11.4,F23.5,F15.8)") dist, mcvirial%effective(idist) * atomic_to_kJmol / total_n, &
                                           mcvirial%coefficient(idist) / total_n

    dist = dist - mcvirial%dist_step

  END DO dist_loop

  CLOSE(991)

  mcvirial%coefficient(:) = mcvirial%coefficient(:) / total_n
  mcvirial%coefficient(:) = mcvirial%coefficient(:) - 1.0_DP

  ! integrate trapezoidal
  b2 = 0.0
  dist = mcvirial%min_dist
  DO idist = 1, ndist-1
    b2 = b2 + ( (mcvirial%coefficient(idist) + mcvirial%coefficient(idist+1)) * dist**2.0)
    dist = dist + mcvirial%dist_step
    !mcvirial%coefficient(idist) =  dexp(-1.0_DP * mcvirial%effective(idist) / total_n * beta(1))
  END DO
  !mcvirial%coefficient(ndist) =  dexp(-1.0_DP * mcvirial%effective(ndist) / total_n * beta(1))

  b2 = -twoPI * (mcvirial%dist_step / 2.0_DP) * b2 
  WRITE(logunit,'(A,2x,F20.7)') 'VIRIAL B2 <exp[-U*beta]-1>: ', b2
  WRITE(logunit, *)

!  mcvirial%coefficient(:) = mcvirial%coefficient(:) - 1.0_DP
!  b2 = 0.0
!  dist = mcvirial%min_dist
!  DO idist = 1, ndist-1
!    b2 = b2 + ( (mcvirial%coefficient(idist) + mcvirial%coefficient(idist+1)) * dist**2.0)
!    dist = dist + mcvirial%dist_step
!  END DO
!
!  b2 = -twoPI * (mcvirial%dist_step / 2.0_DP) * b2 
!  WRITE(logunit,'(A,2x,F20.7)') 'VIRIAL B2 exp[<-U*beta>]-1: '
!  WRITE(logunit,'(2x,F20.7)') b2

end subroutine virialMC_Driver
