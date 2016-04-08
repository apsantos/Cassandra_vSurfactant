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

SUBROUTINE PP_Driver
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
  USE File_Names
  USE Energy_Routines
  USE Read_Write_Checkpoint
  USE Cluster_Routines
  USE Excluded_Volume
  USE Degree_Association

  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER :: i, this_box, ibox, is, which_step
  INTEGER :: ioldN,  delta_n  ! old molecule number and change in molecule number

  REAL(DP) :: time_start, now_time

  LOGICAL :: overlap, complete, lopen

  LOGICAL, DIMENSION(:), ALLOCATABLE :: next_write, next_rdf_write

  ALLOCATE(next_write(nbr_boxes))
  ALLOCATE(next_rdf_write(nbr_boxes))
  next_write(:) = .false.
  next_rdf_write(:) = .false.
  complete = .FALSE.
 
  this_box = 1

  IF(.NOT. openmp_flag) THEN
     CALL cpu_time(time_start)
  ELSE
!$  time_start = omp_get_wtime()
  END IF

  i = 0

  DO WHILE (.NOT. complete)

     i = i + 1

     ! We will select a move from Golden Sampling scheme
  
     which_step = i
     ioldN = nmols(1,1)  ! store beginning molecule number 

     IF (i == n_mcsteps) complete = .TRUE.

     CALL Read_XYZ(i)
 
     next_write(this_box) = .true.
     next_rdf_write(this_box) = .true.


     INQUIRE(file=xyz_config_file(this_box),opened=lopen)
     IF (.not. lopen) THEN
        WRITE(*,*) 'Only, '//TRIM(Int_To_String(i-1))//', steps in xyz file'
        IF (i < n_equilsteps) THEN
            WRITE(*,*) 'Too many equil steps given'
        ELSE
            WRITE(*,*) 'Too many steps given'
        END IF
        complete = .TRUE.
        nthermo_freq = 0
     ENDIF

     IF (i < n_equilsteps) CYCLE

     CALL Accumulate(this_box)
     
     IF ( .NOT. block_average ) THEN

        ! instantaneous values are to be printed   

        IF ( ncoord_freq /= 0 ) THEN
           IF ( MOD(i,ncoord_freq) == 0 ) THEN
              
              DO ibox = 1, nbr_boxes
                 
                 CALL Write_Coords(ibox)
                 
              END DO
              
           END IF
        END IF

        IF ( ncluster_freq /= 0 ) THEN
           IF ( MOD(i,ncluster_freq) == 0 ) THEN
           
              DO ibox = 1, nbr_boxes
              
                 CALL Find_Clusters(ibox)
                 CALL Write_Cluster(ibox)
              
              END DO
           
           END IF
        END IF
        
        IF ( nexvol_freq /= 0 ) THEN
           IF ( MOD(i,nexvol_freq) == 0 ) THEN
           
              DO ibox = 1, nbr_boxes
                 IF ( MOD(i,ncluster_freq) /= 0 ) THEN
                    CALL Find_Clusters(ibox)
                 END IF
              
                 CALL Calculate_Excluded_Volume(ibox)
              
              END DO
           
           END IF
        END IF

        IF ( nalpha_freq /= 0 ) THEN
           IF ( MOD(i,nalpha_freq) == 0 ) THEN
           
              DO ibox = 1, nbr_boxes
                 IF ( MOD(i,ncluster_freq) /= 0 ) THEN
                    CALL Find_Clusters(ibox)
                 END IF
              
                 CALL Calculate_Degree_Association(ibox)
              
              END DO
           
           END IF
        END IF
        
        IF ( nthermo_freq /= 0) THEN
           IF ( MOD(i,nthermo_freq) == 0) THEN
   
              DO ibox = 1, nbr_boxes
              
                 CALL Write_Properties(i,ibox)
                 CALL Reset(ibox)
                 
              END DO
              
           END IF
        END IF
        
     ELSE
        
        DO ibox = 1, nbr_boxes
           
           IF (tot_trials(ibox) /= 0) THEN
              IF(MOD(tot_trials(ibox),nthermo_freq) == 0) THEN
                 IF (next_write(ibox)) THEN
                    CALL Write_Properties(tot_trials(ibox),ibox)
                    CALL Reset(ibox)
                    next_write(ibox) = .false.
                 END IF
              END IF
              
              IF (MOD(tot_trials(ibox), ncoord_freq) == 0) THEN
                 IF (next_rdf_write(ibox)) THEN
                    CALL Write_Coords(ibox)
                    next_rdf_write(ibox) = .false.
                 END IF
              END IF
              
           END IF
           
        END DO
     
     END IF

     IF ( ncoord_freq /= 0 ) THEN
        IF(MOD(i,ncoord_freq) == 0) THEN
           CALL Write_Checkpoint(i)
        END IF
     END IF

     DO is = 1,nspecies

        IF(species_list(is)%int_insert == int_igas) THEN
           IF(mod(i,n_igas_moves(is)) == 0) CALL Update_Reservoir(is)
        END IF

     END DO
     
  END DO

  CLOSE(50)

  ! let us check if at the end of the simulation, the energies are properly updated

  write(logunit,*) '*********** Ending simulation *****************'
  write(logunit,*)
  write(logunit,*)

    ! Display the components of the energy.
  WRITE(logunit,*) '*****************************************'
  WRITE(logunit,'(A36,2X,I2)') ' Starting energy components for box', this_box
  WRITE(logunit,*) ' Atomic units-Extensive'
  WRITE(logunit,*) '*****************************************'
  WRITE(logunit,*)

  write(logunit,'(A,T30,F20.3)') 'Total system energy is' , energy(this_box)%total
  write(logunit,'(A,T30,F20.3)') 'Intra molecular energy is', energy(this_box)%intra
  write(logunit,'(A,T30,F20.3)') 'Bond energy is', energy(this_box)%bond
  write(logunit,'(A,T30,F20.3)') 'Angle energy is', energy(this_box)%angle
  write(logunit,'(A,T30,F20.3)') 'Dihedral enregy is', energy(this_box)%dihedral
  WRITE(logunit,'(A,T30,F20.3)') 'Improper angle energy is', energy(this_box)%improper
  write(logunit,'(A,T30,F20.3)') 'Intra nonbond vdw is', energy(this_box)%intra_vdw
  write(logunit,'(A,T30,F20.3)') 'Intra nonbond elec is', energy(this_box)%intra_q
  write(logunit,'(A,T30,F20.3)') 'Inter molecule vdw is', energy(this_box)%inter_vdw
  write(logunit,'(A,T30,F20.3)') 'Long range correction is', energy(this_box)%lrc
  write(logunit,'(A,T30,F20.3)') 'Inter molecule q is', energy(this_box)%inter_q
  write(logunit,'(A,T30,F20.3)') 'Reciprocal ewald is', energy(this_box)%ewald_reciprocal
  write(logunit,'(A,T30,F20.3)') 'Self ewald is', energy(this_box)%ewald_self
  
  write(logunit,*) '**************************************************'


  CALL Compute_Total_System_Energy(this_box,.TRUE.,overlap)

    ! Display the components of the energy.
  write(logunit,*)
  write(logunit,*)
  WRITE(logunit,*) '*****************************************'
  write(logunit,'(A52,2X,I2)') 'Components of energy from total energy call for box', this_box
  WRITE(logunit,*) 'Atomic units-Extensive'
  WRITE(logunit,*) '*****************************************'
  WRITE(logunit,*)

  write(logunit,'(A,T30,F20.3)') 'Total system energy is' , energy(this_box)%total
  write(logunit,'(A,T30,F20.3)') 'Intra molecular energy is', energy(this_box)%intra
  write(logunit,'(A,T30,F20.3)') 'Bond energy is', energy(this_box)%bond
  write(logunit,'(A,T30,F20.3)') 'Angle energy is', energy(this_box)%angle
  write(logunit,'(A,T30,F20.3)') 'Dihedral enregy is', energy(this_box)%dihedral
  WRITE(logunit,'(A,T30,F20.3)') 'Improper angle energy is', energy(this_box)%improper
  write(logunit,'(A,T30,F20.3)') 'Intra nonbond vdw is', energy(this_box)%intra_vdw
  write(logunit,'(A,T30,F20.3)') 'Intra nonbond elec is', energy(this_box)%intra_q
  write(logunit,'(A,T30,F20.3)') 'Inter molecule vdw is', energy(this_box)%inter_vdw
  write(logunit,'(A,T30,F20.3)') 'Long range correction is', energy(this_box)%lrc
  write(logunit,'(A,T30,F20.3)') 'Inter molecule q is', energy(this_box)%inter_q
  write(logunit,'(A,T30,F20.3)') 'Reciprocal ewald is', energy(this_box)%ewald_reciprocal
  write(logunit,'(A,T30,F20.3)') 'Self ewald is', energy(this_box)%ewald_self

  IF(int_run_style == run_test) THEN

    OPEN(75,FILE='compare.dat',POSITION='APPEND')
    WRITE(75,"(T20,A,A)") testname, 'in the gcmc ensemble'
    WRITE(75,"(A,F24.12)") 'The total system energy is:', energy(1)%total
    WRITE(75,"(A,I10)") 'The total number of molecules is:', nmols(1,1)
    WRITE(75,*)
    CLOSE(75)

  END IF

  END SUBROUTINE PP_Driver
