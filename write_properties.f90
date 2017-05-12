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
!********************************************************************************
SUBROUTINE Write_Properties(this_mc_step,this_box)
  ! The subroutine will write desired properties to the property files. It is
  ! called by respective drivers such as.
  !
  ! CALLED BY
  !
  !        gcmc_driver
  !        gemc_driver
  !        nptmc_driver
  !        nvtmc_driver
  ! CALLS
  !
  !   None
  !
  ! 08/12/13 : Created beta version
!*********************************************************************************

  USE Run_Variables
  USE File_Names
  USE Energy_Routines

  IMPLICIT NONE

  CHARACTER(24) :: write_str
  INTEGER :: i, this_box, this_unit, this_mc_step


  DO i = 1, nbr_prop_files(this_box)
     
     !-- check to see if the file is open or not
     this_unit = (this_box-1)*MAXVAL(nbr_prop_files) + i + propunit
     
     IF (first_open(i,this_box)) THEN
        
        OPEN(unit=this_unit,file=prop_files(i,this_box))
        
        ! write the header information that indicates the properties contained
        ! in the file
        CALL Write_Header(i)

        first_open(i,this_box) = .FALSE.

     END IF
        
     CALL Write_Properties_Buffer(i)
     
  END DO


CONTAINS

  SUBROUTINE Write_Header(file_number)

    IMPLICIT NONE

    INTEGER :: file_number, ii
    CHARACTER(120) :: prop_to_write
    CHARACTER(120), ALLOCATABLE :: prop_unit(:)

    IF (block_average) THEN
       WRITE(this_unit,'(A)') '# Block averages'
    ELSE
       WRITE(this_unit,'(A)') '# Instantaneous properties'
    END IF

    write_str = ""
    write_str = "# MC_STEP"
    
    WRITE(this_unit,'(A12,2X)',ADVANCE='NO') ADJUSTL(write_str)
    
    ! Now write strings for the rest of fields.
    
    DO ii = 1, prop_per_file(file_number,this_box)-1
       
       WRITE(this_unit,'(A16,2X)',ADVANCE='NO') (TRIM(prop_output(ii,file_number,this_box)))
       
    END DO

    WRITE(this_unit,'(A16,2X)') (TRIM(prop_output(ii,file_number,this_box)))


    ALLOCATE(prop_unit(prop_per_file(file_number,this_box)+1))

    prop_unit(:) = ""
    prop_unit(1) ='# '

    WRITE(this_unit, '(A12,2X)',ADVANCE='NO') ADJUSTL(prop_unit(1))

    DO ii = 1, prop_per_file(file_number,this_box) - 1

       prop_to_write = prop_output(ii,file_number,this_box)

       IF (prop_to_write(1:6) == 'Energy') THEN

          prop_unit(ii) = '(kJ/mol)-Ext '

       ELSE IF (prop_to_write == 'Pressure') THEN
          
          prop_unit(ii) = '(bar)'

       ELSE IF (prop_to_write == 'Volume') THEN

          prop_unit(ii) = '(A^3)'

       ELSE IF (prop_to_write == 'Density') THEN

          prop_unit(ii) = '(molec/A^3)'

       END IF

       WRITE(this_unit,'(A16,2X)',ADVANCE='NO') (TRIM(prop_unit(ii)))

    END DO

    prop_to_write = prop_output(ii,file_number,this_box)
    
    IF (prop_to_write(1:6) == 'Energy') THEN
       
       prop_unit(ii) = '(kJ/mol)-Ext '
       
    ELSE IF (prop_to_write == 'Pressure') THEN
       
       prop_unit(ii) = '(bar)'
       
    ELSE IF (prop_to_write == 'Volume') THEN
       
       prop_unit(ii) = '(A^3)'
       
    ELSE IF (prop_to_write == 'Density') THEN
       
       prop_unit(ii) = '(molec/A^3)'

    END IF
    
    WRITE(this_unit,'(A16,2X)') (TRIM(prop_unit(ii)))

    DEALLOCATE(prop_unit)
    
  END SUBROUTINE Write_Header

 SUBROUTINE Write_Properties_Buffer(file_number)
   !************************************************************************
   ! The subroutine fills in a line buffer based on which properties are to
   ! be written and then write the buffer to a file.
   !
   ! Writte by Jindal Shah
   !
   ! Revision History
   !
   !*************************************************************************

   USE Simulation_Properties
   
   INTEGER :: file_number, ii, is, is_dens, is_cp, is_lambda
   REAL(DP),DIMENSION(:), ALLOCATABLE :: write_buff
   CHARACTER(FILENAME_LEN) :: prop_written

   ALLOCATE(write_buff(prop_per_file(file_number,this_box)+1))

   write_buff(1) = this_mc_step
  
   !***********************************************************************
   ! Fill the elements of write_buff with each of the properties
   ! Note that these are average properties over the frequency interval
   !***********************************************************************

!   DO ii = 1, prop_per_file(file_number,this_box)

   ii = 1
   is = 1
   is_dens = 1
   is_cp = 1
   is_lambda = 1

   DO WHILE ( ii <= prop_per_file(file_number,this_box))

      prop_written = prop_output(ii,file_number,this_box)

      IF (prop_written == 'Energy_Total') THEN

         IF ( block_average) THEN
            write_buff(ii+1) = ac_energy(this_box)%total / REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%total
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_LJ') THEN

         IF ( block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%inter_vdw + ac_energy(this_box)%intra_vdw) / &
                 REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_vdw + energy(this_box)%intra_vdw
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Elec') THEN

         IF ( block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%inter_q + ac_energy(this_box)%intra_q) / &
                 REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_q + energy(this_box)%intra_q
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Intra') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%intra
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Pressure') THEN

         CALL Compute_Forces(this_box)

         Pressure_tensor(:,:,this_box) = W_tensor_total(:,:,this_box) / box_list(this_box)%volume
         P_inst(this_box) = ((Pressure_tensor(1,1,this_box) + Pressure_tensor(2,2,this_box) + &
                             Pressure_tensor(3,3,this_box)) / 3.0_DP) * atomic_to_bar
    
         IF(int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
            P_inst(this_box) = P_inst(this_box) + ((virial(this_box)%lrc / box_list(this_box)%volume) * atomic_to_bar) 
         END IF

         P_ideal(this_box) = SUM(nmols(:,this_box)) / box_list(this_box)%volume * temperature(this_box) * p_const

         write_buff(ii+1) = P_ideal(this_box) + P_inst(this_box)

      ELSE IF (prop_written == 'Energy_Intra_VDW') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%intra_vdw
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Inter_VDW') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_vdw
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Intra_Q') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%intra_q
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Inter_Q') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_q
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Recip') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%ewald_reciprocal
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Self') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%ewald_self
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Volume') THEN
         
         IF (block_average) THEN
            write_buff(ii+1) = (ac_volume(this_box))/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = box_list(this_box)%volume
         END IF

      ELSE IF (prop_written == 'Enthalpy') THEN
         IF (block_average) THEN
            write_buff(ii+1) = (ac_enthalpy(this_box))/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%total + pressure(this_box) * box_list(this_box)%volume
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Nmols') THEN

         IF (block_average) THEN
            write_buff(ii+1) = ac_nmols(is,this_box) / REAL(nthermo_freq,DP)
         ELSE
               write_buff(ii+1) = nmols(is,this_box)
         END IF

         ! increment the species index by 1 so that if there is
         ! another species and if nmols is to be output for that
         ! species, we will have correct index
         is = is + 1
         
         
      ELSE IF (prop_written == 'Density') THEN

         IF (block_average) THEN
            write_buff(ii+1) = ac_density(is_dens,this_box)/ REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = REAL(nmols(is_dens,this_box),DP) / box_list(this_box)%volume
         END IF
         ! increment the species index by 1 for the same reason as species
         ! in 'Nmols' was incremented
         is_dens = is_dens + 1

      ELSE IF (prop_written == 'Chemical_Potential') THEN
         write_buff(ii+1) = chpot(is_cp,this_box) / REAL(ntrials(is_cp,this_box)%cpcalc)
         is_cp = is_cp + 1
         
      ELSE IF (prop_written == 'Noligomers') THEN

         IF (block_average) THEN
            write_buff(ii+1) = cluster%n_oligomers / REAL(ncluster_freq,DP)
         ELSE
            write_buff(ii+1) = cluster%n_oligomers
         END IF

      END IF
      
      ! At the end increment property counter by 1

      ii = ii + 1

   END DO

   ! write the line buffer to the property file

   WRITE(this_unit,'(I12,2X)',ADVANCE='NO') this_mc_step
   DO ii = 1, prop_per_file(file_number,this_box)-1

      WRITE(this_unit,'(E16.8,2X)',ADVANCE='NO') write_buff(ii+1)

   END DO
   WRITE(this_unit,'(E16.8,2X)') write_buff(prop_per_file(file_number,this_box)+1)
   
   DEALLOCATE(write_buff)

 END SUBROUTINE Write_Properties_Buffer
 
END SUBROUTINE Write_Properties
!**************************************************************************************
SUBROUTINE Write_Coords(this_box)
  !************************************************************************************
  ! The subroutine writes coordinates of simulation box for later analyis of
  ! RDFs. It gets called by driver routines.
  !
  ! CALLED BY
  !
  !        gcmc_driver
  !        gemc_driver
  !        nptmc_driver
  !        nvtmc_driver
  !
  ! 08/12/13 (JS) : Created beta version
  !*************************************************************************************
  
  USE Run_Variables
  USE Simulation_Properties
  USE File_Names, ONLY : movie_header_unit,movie_xyz_unit 

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box

 !************************************************************************************
  
  INTEGER :: ii, jj, is, nmolecules_is, im, this_im, ia 
  INTEGER :: M_XYZ_unit,MH_unit,Num_Atoms

  ! Write the information about volume

  MH_unit = movie_header_unit + this_box
  M_XYZ_unit = movie_xyz_unit + this_box

  Num_Atoms = 0

  WRITE(MH_unit,*) box_list(this_box)%volume
  ! The cell matrix
  DO ii = 1, 3
     WRITE(MH_unit,*)(box_list(this_box)%length(ii,jj), jj=1,3)
  END DO

  WRITE(MH_unit,*)
  WRITE(MH_unit,*) nspecies
  
  !-- Number of molecules of each of the species
  DO is = 1, nspecies
     CALL Get_Nmolecules_Species(this_box,is,nmolecules_is)
     WRITE(MH_unit,*) is,nmolecules_is
     Num_Atoms = Num_Atoms + nmolecules_is*natoms(is)
  END DO

  !--- Write the coordinates of molecules in this box
  WRITE(M_XYZ_unit,*) Num_Atoms
  WRITE(M_XYZ_unit,*)
  DO is = 1,nspecies
     DO im = 1,nmolecules(is)
        this_im = locate(im,is)
        IF(molecule_list(this_im,is)%live  .AND. &
           molecule_list(this_im,is)%which_box == this_box ) THEN
           DO ia = 1, natoms(is)
!FSL Write is and im values           
              WRITE(M_XYZ_unit,'(A3,F20.13,F20.13,F20.13,I6,I6)') &
              nonbond_list(ia,is)%atom_name, & !element, &
                   atom_list(ia,this_im,is)%rxp, &
                   atom_list(ia,this_im,is)%ryp, &
                   atom_list(ia,this_im,is)%rzp, &
                   is, im
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE Write_Coords

SUBROUTINE Write_Cluster(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        *_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE Cluster_Routines
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: iM, box_unit

  cluster_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.clu'
  box_unit = cluster_file_unit + this_box
  OPEN(unit=box_unit, file=cluster_file)

  WRITE(box_unit,*) '#   M    pop'
  
  DO iM = 1, SIZE(cluster%M)
     IF (cluster%M(iM) > 0) THEN
        WRITE(box_unit,'(I6, I10)') iM, cluster%M(iM)
     END IF
  END DO
  CLOSE(unit=box_unit)
 
END SUBROUTINE Write_Cluster

SUBROUTINE Write_Histogram(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        *_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE Cluster_Routines
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER             :: this_unit, is, imol, istart, iend, i, namph, hist_index
  REAL(DP)            :: e_total

  e_total = energy(this_box)%total

  DO is = 1, nspecies

      namph = nmols(is,this_box)
      if (energy_hist(is,-2,namph) == energy_hist(is,0,namph)) then    ! first energy observed for this namph - set
                                                              ! at mid-point
          energy_hist(is,0,namph) = REAL( energy_hist_width * int(e_total/energy_hist_width - n_energy_hist/2), SP)
          energy_hist(is,-2,namph) = REAL( energy_hist(is,0,namph) - &
                                     energy_hist_width * (int((energy_hist(is,0,namph)-e_total)/energy_hist_width)+1), SP)
          energy_hist(is,-1,namph) = REAL( energy_hist(is,0,namph) + &
                                     energy_hist_width * (int((-energy_hist(is,0,namph)+e_total)/energy_hist_width)+1), SP)
      endif

      if (e_total < energy_hist(is,-2,namph)) then    ! found lower value than min
          energy_hist(is,-2,namph) = REAL( energy_hist(is,0,namph) - &
                                     energy_hist_width * (int((energy_hist(is,0,namph)-e_total)/energy_hist_width)+1), SP)
          if (energy_hist(is,-2,namph) <= energy_hist(is,0,namph) + energy_hist_width ) then
              write (logunit,*) 'Maximum width of energy histograms exceeded (low)'
              write (logunit,*) 'namph, energy(-2),energy(-1),energy(0),e_total', namph,energy_hist(is,-2:0,namph),e_total
              stop
          endif

      else if (e_total > energy_hist(is,-1,namph)) then    ! found greater value than max
          energy_hist(is,-1,namph) = REAL( energy_hist(is,0,namph) + &
                                     energy_hist_width * (int((-energy_hist(is,0,namph)+e_total)/energy_hist_width)+1), SP)
          if (energy_hist(is,-1,namph) > energy_hist(is,0,namph)+energy_hist_width*(n_energy_hist-2)) then
              write (logunit,*) 'Maximum width of energy histograms exceeded (high)'
              write (logunit,*) 'namph, energy(-2),energy(-1),energy(0),e_total', namph,energy_hist(is,-2:0,namph),e_total
              stop
          endif

      endif
      ! value ok
      hist_index = int((e_total-energy_hist(is,0,namph))/energy_hist_width + 0.5) 
      energy_hist(is,hist_index,namph) = energy_hist(is,hist_index,namph) + 1

      this_unit = histogram_file_unit + is 
      histogram_file = 'his_'//TRIM(run_name)// '.spec' // TRIM(Int_To_String(is)) //'.dat'
      OPEN(unit=this_unit, file=histogram_file)

      WRITE (this_unit,'(A)') '    kBT       mu          width     x- y- z-dim  (Atomistic units)'
      WRITE (this_unit,'(3f11.3,3f8.3)') 1./beta(this_box), species_list(is)%chem_potential, energy_hist_width, &
                                        box_list(this_box)%length(1,1), &
                                        box_list(this_box)%length(2,2), &
                                        box_list(this_box)%length(3,3)

      DO imol = 0, nmolecules(is)
          istart = INT( (energy_hist(is, -2,imol) - energy_hist(is, 0, imol)) / energy_hist_width )
          iend   = INT( (energy_hist(is, -1,imol) - energy_hist(is, 0, imol)) / energy_hist_width )

          IF (istart /= iend) THEN 
              WRITE (this_unit,'(2i7,f13.5)') imol, iend-istart+1, energy_hist(is, -2,imol)
              WRITE (this_unit,'(8f9.0)') (energy_hist(is, i, imol), i = istart, iend)
          ENDIF

      END DO

      CLOSE (this_unit) 
  END DO

END SUBROUTINE Write_Histogram
