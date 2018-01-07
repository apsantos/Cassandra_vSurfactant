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

  INTEGER :: this_box, nfailure, ntests

  LOGICAL :: overlap, get_energy

  ! The total number of trial move array may not have been set if this
  ! is a fresh run i.e. start_type == make_config. Otherwise this array
  ! is set in read_checkpoint subroutine in the module Read_Write_Checkpoint
  IF(.NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF(.NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))

  get_energy = .false.
  ntests = 17
  nfailure = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Force Routine Tests
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL Test_Pair_Force_QQ(nfailure)

  CALL Test_Pair_Force_QQcorr(nfailure)

  CALL Test_Pair_Force_Hydr(nfailure)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Cluster routine Tests
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cluster joining move rejection 
  CALL Test_Translate_Cluster(nfailure)

  ! cluster finding              
  CALL Test_Find_Cluster(nfailure)

  ! Test neighbor finding routine
  CALL Test_Neighbor_Cluster(nfailure)

  ! Test update labels
  CALL Test_Update_Labels(nfailure)

  ! Test cluster energy calc
  CALL Test_Compute_Cluster_Nonbond_Inter_Energy(nfailure)

  ! let us check if at the end of the simulation, the energies are properly updated
  write(*,'(A,I6,A,I6)') 'nfailures', nfailure , ' out of ntests', ntests

  IF (get_energy ) THEN
     CALL Compute_Total_System_Energy(this_box,.TRUE.,overlap)
     write(*,'(A)') '---Energies'
     write(*,'(A)') '          Total          intra-molecular           intra-vdw         intra-elec'
     write(*,'(4F20.3)') energy(this_box)%total, energy(this_box)%intra, energy(this_box)%intra_vdw, energy(this_box)%intra_q
     write(*,'(A)') ''
     write(*,'(A)') '          bond                angle                  dihedral            improper'
     write(*,'(4F20.3)') energy(this_box)%bond, energy(this_box)%angle, energy(this_box)%dihedral, energy(this_box)%improper
     write(*,'(A)') ''
     write(*,'(A)') '          Inter-vdw          longrange-corr      '
     write(*,'(2F20.3)') energy(this_box)%inter_vdw, energy(this_box)%lrc
     write(*,'(A)') ''
     write(*,'(A)') '          inter-elec           reci-elec         self-ewald'
     write(*,'(3F20.3)') energy(this_box)%inter_q, energy(this_box)%ewald_reciprocal, energy(this_box)%ewald_self
  END IF

END SUBROUTINE TEST_Driver

SUBROUTINE Test_Pair_Force_Hydr(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  REAL(DP) :: negrForce_vdw, negrForce_qq
  INTEGER, INTENT(INOUT) :: nfailure

  int_vdw_style_mix(:,:,:) = .false.
  int_vdw_style_mix(1,1,vdw_hydra) = .true.
  nonbond_list(1,1)%charge = 0.0_DP
  vdw_param3_table(1,1) = 10.0
  vdw_param4_table(1,1) = 4.0
  vdw_param5_table(1,1) = 0.5

  molecule_list(:,:)%live = .false.
  molecule_list(1:2,1)%live = .true.
  atom_list(:,:,:)%exist = .false.
  atom_list(1,1:2,1)%exist = .true.
  static_perm(1) = 78.0_DP
  charge_factor(1) = charge_factor_vacuum / static_perm(1)
  rcut_coul(1) = 22.9
  !                ewald_tol
  !                  ^^^^^
  ewald_p(1) = -DLOG( 0.001_DP )
  ewald_p_sqrt(1) = DSQRT(ewald_p(1))
  alpha_ewald(1) = ewald_p_sqrt(1) / rcut_coul(1)
  h_ewald_cut(1) = 2.0_DP * ewald_p(1) / rcut_coul(1)

  CALL Pair_Force(9.0_DP,1,1,1,1,2,1,.false.,.true.,negrForce_vdw,negrForce_qq)

  IF( abs(negrForce_vdw + 12.958) > 0.001) THEN
    nfailure = nfailure + 1
    write(*,*) 'Hydr force acceptance FAILED'
  ELSE
    write(*,*) 'Hydr force acceptance PASSED'
  ENDIF
END SUBROUTINE Test_Pair_Force_Hydr

SUBROUTINE Test_Pair_Force_QQ(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  REAL(DP) :: negrForce_vdw, negrForce_qq
  INTEGER, INTENT(INOUT) :: nfailure

  int_vdw_style_mix(:,:,:) = .false.
  nonbond_list(1,1)%charge = -1.0_DP

  molecule_list(:,:)%live = .false.
  molecule_list(1:2,1)%live = .true.
  atom_list(:,:,:)%exist = .false.
  atom_list(1,1:2,1)%exist = .true.
  static_perm(1) = 78.0_DP
  charge_factor(1) = charge_factor_vacuum / static_perm(1)
  rcut_coul(1) = 22.9
  !                ewald_tol
  !                  ^^^^^
  ewald_p(1) = -DLOG( 0.001_DP )
  ewald_p_sqrt(1) = DSQRT(ewald_p(1))
  alpha_ewald(1) = ewald_p_sqrt(1) / rcut_coul(1)
  h_ewald_cut(1) = 2.0_DP * ewald_p(1) / rcut_coul(1)

  CALL Pair_Force(9.0_DP,1,1,1,1,2,1,.false.,.true.,negrForce_vdw,negrForce_qq)

  IF( abs(negrForce_qq - 0.323796) > 0.001) THEN
    nfailure = nfailure + 1
    write(*,*) 'QQ coulomb force acceptance FAILED'
  ELSE
    write(*,*) 'QQ coulomb force acceptance PASSED'
  ENDIF
END SUBROUTINE Test_Pair_Force_QQ

SUBROUTINE Test_Pair_Force_QQcorr(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  REAL(DP) :: negrForce_vdw, negrForce_qq
  INTEGER, INTENT(INOUT) :: nfailure

  int_vdw_style_mix(:,:,:) = .false.
  int_vdw_style_mix(1,1,vdw_corr) = .true.
  nonbond_list(1,1)%charge = -1.0_DP
  int_charge_sum_style(1) = charge_ewald
  vdw_param6_table(1,1) = 2.0
  vdw_param7_table(1,1) = 0.5

  molecule_list(:,:)%live = .false.
  molecule_list(1:2,1)%live = .true.
  atom_list(:,:,:)%exist = .false.
  atom_list(1,1:2,1)%exist = .true.
  static_perm(1) = 78.0_DP
  charge_factor(1) = charge_factor_vacuum / static_perm(1)
  rcut_coul(1) = 22.9
  !                ewald_tol
  !                  ^^^^^
  ewald_p(1) = -DLOG( 0.001_DP )
  ewald_p_sqrt(1) = DSQRT(ewald_p(1))
  alpha_ewald(1) = ewald_p_sqrt(1) / rcut_coul(1)
  h_ewald_cut(1) = 2.0_DP * ewald_p(1) / rcut_coul(1)

  CALL Pair_Force(9.0_DP,1,1,1,1,2,1,.false.,.true.,negrForce_vdw,negrForce_qq)

  IF( abs(negrForce_qq - 0.37205) > 0.001) THEN
    nfailure = nfailure + 1
    write(*,*) 'QQ correction force acceptance FAILED'
  ELSE
    write(*,*) 'QQ correction force acceptance PASSED'
  ENDIF
END SUBROUTINE Test_Pair_Force_QQcorr

SUBROUTINE Test_Translate_Cluster(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  INTEGER :: i, this_box
  LOGICAL :: overlap
  INTEGER, INTENT(INOUT) :: nfailure
  cluster%M = 0
  cluster%N = 0
  cluster%clabel = 0

  this_box = 1

  int_vdw_style_mix(:,:,:) = .false.
  int_vdw_style_mix(:,:,vdw_yukawa) = .true.
  int_vdw_style_mix(:,:,vdw_lj) = .true.
  int_vdw_sum_style_mix(:,:) = vdw_cut_shift
  vdw_param1_table(1,1) = 1.0
  vdw_param2_table(1,1) = 1.0
  vdw_param8_table(1,1) = 0.5
  vdw_param9_table(1,1) = 0.5

  cluster%min_distance_sq(:,:,:, 0, 0 ) = 0
  cluster%min_distance_sq(:,1,1, 0, 0 ) = 1.8096**2.0
  molecule_list(:,:)%live = .false.
  molecule_list(1:4,1)%live = .true.
  atom_list(:,:,:)%exist = .false.
  atom_list(1,1:4,1)%exist = .true.

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

  DO i = 1, 8
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  CALL Translate_Cluster(this_box)
  CALL Compute_Total_System_Energy(this_box,.TRUE.,overlap)
  IF (cluster%M(1) /= 2 .or. &
      cluster%M(2) /= 1 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster energy acceptance FAILED - wrong # clusters found'
  ELSEIF (reject_type /= 0) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster energy acceptance FAILED - the move was not accepted'
  ELSE
    write(*,*) 'Cluster energy acceptance PASSED'
  ENDIF
    
  CALL Translate_Cluster(this_box)
  IF (reject_type /= -1) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster joining move rejection test FAILED'

  ELSE
    write(*,*) 'Cluster joining move rejection test PASSED'
  ENDIF

  atom_list(1,1,1)%rxp = -0.7
  atom_list(1,2,1)%rxp = -0.7
  atom_list(1,4,1)%rxp = 3
  atom_list(1,4,1)%ryp = 2
  DO i = 1, 8
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  vdw_param8_table(1,1) = 10
  CALL Translate_Cluster(this_box)
  IF (reject_type /= -3) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster move energy rejection test FAILED'
  ELSE
    write(*,*) 'Cluster move energy rejection test PASSED'
  ENDIF
    
  vdw_param8_table(1,1) = 0.5
END SUBROUTINE Test_Translate_Cluster

SUBROUTINE Test_Find_Cluster(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  INTEGER :: i, this_box
  INTEGER, INTENT(INOUT) :: nfailure

  this_box = 1

  cluster%M = 0
  cluster%N = 0
  cluster%clabel = 0

  cluster%min_distance_sq(:,:,:, 0, 0 ) = 0
  cluster%min_distance_sq(:,1,1, 0, 0 ) = 1.8096**2.0

  vdw_param1_table(1,1) = 1.0
  vdw_param2_table(1,1) = 1.0
  vdw_param8_table(1,1) = 0.5
  vdw_param9_table(1,1) = 0.5

  box_list(this_box)%length(1,1) = 17.0
  box_list(this_box)%length(2,2) = 17.0
  box_list(this_box)%length(3,3) = 17.0

  box_list(this_box)%hlength(1,1) = 0.5_DP * box_list(this_box)%length(1,1)
  box_list(this_box)%hlength(2,2) = 0.5_DP * box_list(this_box)%length(2,2)
  box_list(this_box)%hlength(3,3) = 0.5_DP * box_list(this_box)%length(3,3)

  nmolecules(1) = 8
  nmols(1,1) = 8
  DO i = 1, 8
    locate(i,1) = i
  ENDDO
  atom_list(:,:,:)%exist = .false.
  atom_list(1,1:8,1)%exist = .true.
  molecule_list(:,:)%live = .false.
  molecule_list(1:8,1)%live = .true.
  molecule_list(:,1)%which_box = 1

  ! make positions                   !
  ! AT THE VERTICIES
  !   __
  !  /_/|
  !  |_|/
  nonbond_list(1,1)%element = 'C'
  atom_list(1,1:4,1)%rxp = -8
  atom_list(1,5:8,1)%rxp = 8
  atom_list(1,1:2,1)%ryp = 8
  atom_list(1,5:6,1)%ryp = 8
  atom_list(1,3:4,1)%ryp = -8
  atom_list(1,7:8,1)%ryp = -8
  atom_list(1,1,1)%rzp = -8
  atom_list(1,2,1)%rzp = 8
  atom_list(1,3,1)%rzp = -8
  atom_list(1,4,1)%rzp = 8
  atom_list(1,5,1)%rzp = -8
  atom_list(1,6,1)%rzp = 8
  atom_list(1,7,1)%rzp = -8
  atom_list(1,8,1)%rzp = 8

  DO i = 1, 8
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  CALL Find_Clusters(this_box,1)
  IF ( cluster%M(8) /= 1 ) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster Finding FAILED - 8 corners'
  ELSE
    write(*,*) 'Cluster Finding PASSED - 8 corners'
  ENDIF

  CALL Write_Coords(1)

  ! make positions                   !
  cluster%M = 0
  cluster%N = 0
  cluster%clabel = 0

  nonbond_list(1,1)%element = 'C'
  atom_list(1,1,1)%rxp = -5
  atom_list(1,2,1)%rxp = -6
  atom_list(1,3,1)%rxp = -7
  atom_list(1,4,1)%rxp = -8
  atom_list(1,5:8,1)%rxp = 8
  atom_list(1,1:4,1)%ryp = 8
  atom_list(1,5:8,1)%ryp = -8
  atom_list(1,5,1)%rzp = 8
  atom_list(1,6,1)%rzp = 7
  atom_list(1,7,1)%rzp = 6
  atom_list(1,8,1)%rzp = 5
  atom_list(1,1:4,1)%rzp = -8

  DO i = 1, 8
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  CALL Find_Clusters(this_box,1)
  IF ( cluster%M(8) /= 1 ) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster Finding FAILED - 2 snakes on 2 corners'
  ELSE
    write(*,*) 'Cluster Finding PASSED - 2 snakes on 2 corners'
  ENDIF

  CALL Write_Coords(1)
END SUBROUTINE Test_Find_Cluster

SUBROUTINE Test_Neighbor_Cluster(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  INTEGER :: i, this_box
  INTEGER, INTENT(INOUT) :: nfailure

  this_box = 1
  cluster%criteria(:,:) = .false.
  cluster%criteria(:, int_com) = .true.

  nonbond_list(1,1)%element = 'C'
  atom_list(1,1,1)%rxp = 0
  atom_list(1,2,1)%rxp = 1
  atom_list(1,1,1)%ryp = 0
  atom_list(1,2,1)%ryp = 1
  atom_list(1,1,1)%rzp = 0
  atom_list(1,2,1)%rzp = 1
  DO i = 1, 2
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  IF ( .not. Neighbor(1, 1, 2, 1, 1) ) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster Neighbor com FAILED'
  ELSE
    write(*,*) 'Cluster Neighbor com PASSED'
  ENDIF

  atom_list(1,1,1)%rxp = 0
  atom_list(1,2,1)%rxp = 1
  atom_list(1,1,1)%ryp = 0
  atom_list(1,2,1)%ryp = 8
  atom_list(1,1,1)%rzp = 0
  atom_list(1,2,1)%rzp = 1
  DO i = 1, 2
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  IF ( Neighbor(1, 1, 2, 1, 1) ) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster Neighbor com FAILED'
  ELSE
    write(*,*) 'Cluster Neighbor com PASSED'
  ENDIF

  ! type criteria
  cluster%criteria(:,:) = .false.
  cluster%criteria(:, int_type) = .true.
  cluster%min_distance_sq(:,1,2, 1, 2 ) = 1.8096**2.0
  molecule_list(1,2)%live = .true.
  atom_list(:,1,2)%exist = .true.

  nonbond_list(1,1)%element = 'C'
  nmolecules(2) = 1
  nmols(2,1) = 1
  atom_list(1,1,1)%rxp = 0
  atom_list(1,1,1)%ryp = 0
  atom_list(1,1,1)%rzp = 0
  molecule_list(1,1)%xcom = 0
  molecule_list(1,1)%ycom = 0
  molecule_list(1,1)%zcom = 0

  atom_list(1,1,2)%rxp = 0.5
  atom_list(1,1,2)%ryp = 1
  atom_list(1,1,2)%rzp = 1
  atom_list(2,1,2)%rxp = -0.5
  atom_list(2,1,2)%ryp = 1
  atom_list(2,1,2)%rzp = 1
  molecule_list(1,2)%xcom = 0
  molecule_list(1,2)%ycom = 1
  molecule_list(1,2)%zcom = 1

  IF ( .not. Neighbor(1, 1, 1, 1, 2) ) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster Neighbor type FAILED'
  ELSE
    write(*,*) 'Cluster Neighbor type PASSED'
  ENDIF

  ! test the type criteria
  atom_list(1,1,1)%rxp = -1
  atom_list(1,1,1)%ryp = -1
  atom_list(1,1,1)%rzp = -1
  molecule_list(1,1)%xcom = -1
  molecule_list(1,1)%ycom = -1
  molecule_list(1,1)%zcom = -1
  IF ( Neighbor(1, 1, 1, 1, 2) ) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster Neighbor type FAILED'
  ELSE
    write(*,*) 'Cluster Neighbor type PASSED'
  ENDIF

END SUBROUTINE Test_Neighbor_Cluster

SUBROUTINE Test_Update_Labels(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  LOGICAL, DIMENSION(nmolecules(1)) :: neigh_list
  INTEGER, INTENT(INOUT) :: nfailure

  cluster%M = 0
  cluster%N = 0
  cluster%clabel = 0
  cluster%clusmax = 0

  ! graph
  ! __| | | | | |
  ! __ 1 5 3
  ! __     6
  ! __   2 4   7 <- unconnected

  ! stage one
  ! 1, 5, 6 have been set
  ! can it join cluster 1 (1,5) and cluster 2 (6) 
  cluster%clusmax = 2
  cluster%N(1) = 2
  cluster%N(2) = 1
  cluster%clabel(1,1) = 1
  cluster%clabel(5,1) = 1
  cluster%clabel(6,1) = 2
  neigh_list = .false.
  neigh_list(5) = .true.
  neigh_list(6) = .true.
  CALL Update_Labels(3, 1, 1, neigh_list )
  IF ( cluster%N(1) /= 4 .or. &
      cluster%N(2) /= -1 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster update labels FAILED - N '
  ELSEIF ( cluster%clabel(3,1) /= 1 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster update labels FAILED - clabel '
  ELSE
    write(*,*) 'Cluster update labels PASSED'
  ENDIF

  ! can it now join 4?
  neigh_list = .false.
  neigh_list(2) = .true.
  neigh_list(6) = .true.
  CALL Update_Labels(4, 1, 1, neigh_list )
  IF ( cluster%N(1) /= 6 .or. &
      cluster%N(2) /= -1 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster update labels FAILED - N '
  ELSEIF ( cluster%clabel(4,1) /= 1 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster update labels FAILED - clabel '
  ELSE
    write(*,*) 'Cluster update labels PASSED'
  ENDIF

  neigh_list = .false.
  CALL Update_Labels(7, 1, 1, neigh_list )
  IF ( cluster%N(4) /= 1 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster update labels FAILED - N '
  ELSEIF ( cluster%clabel(7,1) /= 4 ) THEN 
    nfailure = nfailure + 1
    write(*,*) 'Cluster update labels FAILED - clabel '
  ELSE
    write(*,*) 'Cluster update labels PASSED'
  ENDIF

END SUBROUTINE Test_Update_Labels

SUBROUTINE Test_Compute_Cluster_Nonbond_Inter_Energy(nfailure)
  USE Run_Variables
  USE Energy_Routines
  USE Cluster_Routines

  IMPLICIT NONE

  INTEGER :: i, this_box, nclus_mol
  REAL(DP) :: E_inter_vdw, E_inter_qq
  INTEGER, DIMENSION(2) :: clus_mol, clus_is
  INTEGER, INTENT(INOUT) :: nfailure
  LOGICAL :: overlap

  cluster%M = 0
  cluster%N = 0
  cluster%clabel = 0

  this_box = 1

  vdw_param1_table(1,1) = 1.0
  vdw_param2_table(1,1) = 1.0
  vdw_param8_table(1,1) = 0.5
  vdw_param9_table(1,1) = 0.5

  cluster%min_distance_sq(:,:,:, 0, 0 ) = 0
  cluster%min_distance_sq(:,1,1, 0, 0 ) = 1.8096**2.0
  molecule_list(:,:)%live = .false.
  molecule_list(1:4,1)%live = .true.
  atom_list(:,:,:)%exist = .false.
  atom_list(1,1:4,1)%exist = .true.

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

  DO i = 1, 8
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  nclus_mol = 2
  clus_is(:) = 1
  clus_mol(1) = 1
  clus_mol(2) = 2
  
  CALL Compute_Cluster_Nonbond_Inter_Energy(nclus_mol,clus_mol,clus_is,E_inter_vdw,E_inter_qq,overlap)

  IF( abs(E_inter_vdw - 0.057984436) > 0.001) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster energy acceptance FAILED - wrong energy calculation'
  ELSE
    write(*,*) 'Cluster energy acceptance PASSED'
  ENDIF

  atom_list(1,1,1)%rxp = 0
  atom_list(1,1,1)%rxp = 0
  atom_list(1,2,1)%rxp = 0
  atom_list(1,4,1)%rxp = -8
  atom_list(1,3,1)%rxp = 2
  atom_list(1,1,1)%ryp = 0
  atom_list(1,2,1)%ryp = 1
  atom_list(1,3,1)%ryp = 0
  atom_list(1,4,1)%ryp = -8
  atom_list(1,1:4,1)%rzp = 0

  DO i = 1, 8
    molecule_list(i,1)%xcom = atom_list(1,i,1)%rxp
    molecule_list(i,1)%ycom = atom_list(1,i,1)%ryp
    molecule_list(i,1)%zcom = atom_list(1,i,1)%rzp
  ENDDO

  nclus_mol = 2
  clus_is(:) = 1
  clus_mol(1) = 1
  clus_mol(2) = 2
  
  CALL Compute_Cluster_Nonbond_Inter_Energy(nclus_mol,clus_mol,clus_is,E_inter_vdw,E_inter_qq,overlap)

  IF( abs(E_inter_vdw - 0.06954544) > 0.001) THEN
    nfailure = nfailure + 1
    write(*,*) 'Cluster energy acceptance FAILED - wrong energy calculation'
  ELSE
    write(*,*) 'Cluster energy acceptance PASSED'
  ENDIF

END SUBROUTINE Test_Compute_Cluster_Nonbond_Inter_Energy

