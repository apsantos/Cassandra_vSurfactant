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

SUBROUTINE TEST_Driver
  !****************************************************************************
  !
  ! The subroutine performs GEMC moves. 
  ! 
  ! Called by
  !
  !   main.f90
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
!*******************************************************************************

  USE Run_Variables
  USE Random_Generators
  USE File_Names
  USE Energy_Routines
  USE Read_Write_Checkpoint
  USE Cluster_Routines

  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER :: i, this_box, ibox, is, which_step, nfailure, ntests
  INTEGER :: ioldN,  delta_n  ! old molecule number and change in molecule number

  REAL(DP) :: rand_no
  REAL(DP) :: time_start, now_time

  LOGICAL :: overlap, complete

  LOGICAL, DIMENSION(:), ALLOCATABLE :: next_write, next_rdf_write

  ! The total number of trial move array may not have been set if this
  ! is a fresh run i.e. start_type == make_config. Otherwise this array
  ! is set in read_checkpoint subroutine in the module Read_Write_Checkpoint
  nspecies = 1
  nbr_boxes = 1
  IF(.NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF(.NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))

  ntests = 4
  nfailure = 0
  this_box = 1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test 1                                              !
  !  - cluster joining move rejection                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nbr_atomtypes = 1
  vdw_param1_table(1,1) = 1.0
  vdw_param2_table(1,1) = 1.0
  vdw_param9_table(1,1) = 0.5
  vdw_param10_table(1,1) = 0.5

  !ALLOCATE( nmolecules(nspecies), natoms(nspecies), Stat = AllocateStatus )
  !nmolecules(1) = 3
  !natoms(1) = 3
  !ALLOCATE( atom_list(MAXVAL(natoms), MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
  !ALLOCATE( nonbond_list(MAXVAL(natoms), nspecies), Stat = AllocateStatus )
  !ALLOCATE( molecule_list(MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
  molecule_list(:,:)%live = .false.
  molecule_list(1:4,1)%live = .true.

  ! make positions                   !
  nonbond_list(1,1)%element = 'C'
  atom_list(1,1,1)%rxp = 0
  atom_list(1,1,1)%rxp = 0
  atom_list(1,2,1)%rxp = 0
  atom_list(1,4,1)%rxp = -8
  atom_list(1,3,1)%rxp = 3
  atom_list(1,1,1)%ryp = 0
  atom_list(1,2,1)%ryp = 1
  atom_list(1,3,1)%ryp = 0
  atom_list(1,4,1)%ryp = -8
  atom_list(1,1:4,1)%rzp = 0

  molecule_list(1,1)%xcom = 0
  molecule_list(2,1)%xcom = 0
  molecule_list(3,1)%xcom = 3
  molecule_list(4,1)%xcom = -8
  molecule_list(1,1)%ycom = 0
  molecule_list(2,1)%ycom = 1
  molecule_list(3,1)%ycom = 0
  molecule_list(4,1)%ycom = -8
  molecule_list(1:4,1)%zcom = 0

  CALL Find_Clusters(this_box,1)
  IF (cluster%M(1) /= 2 .or. &
      cluster%M(2) /= 1) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster Finding failed'
  ELSE
    write(*,*) 'Cluster Finding passed'
  ENDIF

  CALL Translate_Cluster(this_box)
  CALL Compute_Total_System_Energy(this_box,.TRUE.,overlap)
  IF (cluster%M(1) /= 4 .or. &
      cluster%M(2) /= 2 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster energy acceptance failed - wrong # clusters found'
  ELSE IF( abs(energy(this_box)%total - 0.371681304684064) > 0.001) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster energy acceptance failed - wrong energy calculation'
  ELSEIF (reject_type /= 0) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster energy acceptance failed - the move was not accepted'
  ELSE
    write(*,*) 'Cluster energy acceptance passed'
  ENDIF
    
  CALL Translate_Cluster(this_box)
  IF (reject_type /= -1) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster joining move rejection test failed'

  ELSE
    write(*,*) 'Cluster joining move rejection test passed'
  ENDIF

  atom_list(1,1,1)%rxp = -0.7
  atom_list(1,2,1)%rxp = -0.7
  atom_list(1,4,1)%rxp = 3
  atom_list(1,4,1)%ryp = 2
  molecule_list(1,1)%xcom = -0.7
  molecule_list(2,1)%xcom = -0.7
  molecule_list(4,1)%xcom = 3
  molecule_list(4,1)%ycom = 2
  vdw_param9_table(1,1) = 10
  CALL Translate_Cluster(this_box)
  IF (reject_type /= -3) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster move energy rejection test failed'
  ELSE
    write(*,*) 'Cluster move energy rejection test passed'
  ENDIF
    
  vdw_param9_table(1,1) = 0.5
  ! let us check if at the end of the simulation, the energies are properly updated
  write(*,*) 'fail rate, nfailures/ntests', nfailure / FLOAT(ntests)

  CALL Compute_Total_System_Energy(this_box,.TRUE.,overlap)
  write(*,'(A)') '---Energies'
  write(*,'(A)') '          Total          intra-molecular           intra-vdw         intra-elec'
  write(*,'(4F20.3)') energy(this_box)%total, energy(this_box)%intra, energy(this_box)%intra_vdw, energy(this_box)%intra_q
  write(*,'(A)') ''
  write(*,'(A)') '          bond                angle                  dihedral            improper'
  write(*,'(4F20.3)') energy(this_box)%bond, energy(this_box)%angle, energy(this_box)%dihedral, energy(this_box)%improper
  write(*,'(A)') ''
  write(*,'(A)') '          Inter-vdw          longrange-corr           inter-elec           reci-elec         self-ewald'
  write(*,'(5F20.3)') energy(this_box)%inter_vdw, energy(this_box)%lrc, energy(this_box)%inter_q, energy(this_box)%ewald_reciprocal, energy(this_box)%ewald_self

  END SUBROUTINE TEST_Driver
