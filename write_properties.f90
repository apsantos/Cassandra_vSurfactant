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
    CHARACTER(240) :: prop_to_write
    CHARACTER(240), ALLOCATABLE :: prop_unit(:)

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

       prop_to_write = TRIM( prop_output(ii,file_number,this_box) )

       IF (prop_to_write(1:6) == 'Energy') THEN

          prop_unit(ii) = '(kJ/mol)-Ext '

       ELSE IF (prop_to_write == 'Pressure') THEN
          
          prop_unit(ii) = '(bar)'

       ELSE IF (prop_to_write == 'Volume') THEN

          prop_unit(ii) = '(A^3)'

       ELSE IF (prop_to_write == 'Density') THEN

          prop_unit(ii) = '(molec/A^3)'

       ELSE IF (prop_to_write(1:8) == 'Excluded') THEN

          prop_unit(ii) = '(% excluded)'

       ELSE IF (prop_to_write(1:6) == 'Degree') THEN

          prop_unit(ii) = '(% associated)'

       END IF

       WRITE(this_unit,'(A16,2X)',ADVANCE='NO') (TRIM(prop_unit(ii)))

    END DO

    prop_to_write = TRIM( prop_output(ii,file_number,this_box) )
    
    IF (prop_to_write(1:6) == 'Energy') THEN
       
       prop_unit(ii) = '(kJ/mol)-Ext '
       
    ELSE IF (prop_to_write == 'Pressure') THEN
       
       prop_unit(ii) = '(bar)'
       
    ELSE IF (prop_to_write == 'Volume') THEN
       
       prop_unit(ii) = '(A^3)'
       
    ELSE IF (prop_to_write == 'Density') THEN
       
       prop_unit(ii) = '(molec/A^3)'

    ELSE IF (prop_to_write(1:8) == 'Excluded') THEN

       prop_unit(ii) = '(% excluded)'

    ELSE IF (prop_to_write(1:6) == 'Degree') THEN

       prop_unit(ii) = '(% associated)'

    ELSE IF (prop_to_write(1:6) == 'Virial') THEN

       prop_unit(ii) = '(kJ/mol)-Ext '

    ELSE IF (prop_to_write(1:6) == 'Effect') THEN

       prop_unit(ii) = '(kJ/mol)-Ext '

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
   
   INTEGER :: file_number, ii, is, is_dens, is_cp
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

      ELSE IF (prop_written == 'Energy_Intra_Angle') THEN

         write_buff(ii+1) = energy(this_box)%angle
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Intra_Dihedral') THEN

         write_buff(ii+1) = energy(this_box)%dihedral
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

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

      ELSE IF (prop_written == 'NclustersOlig') THEN

         IF (block_average) THEN
            write_buff(ii+1) = cluster%n_olig_clus / REAL(ncluster_freq,DP)
         ELSE
            write_buff(ii+1) = cluster%n_olig_clus
         END IF

      ELSE IF (prop_written == 'NclustersMicelle') THEN

         IF (block_average) THEN
            write_buff(ii+1) = cluster%n_mic_clus / REAL(ncluster_freq,DP)
         ELSE
            write_buff(ii+1) = cluster%n_mic_clus
         END IF

      ELSE IF (prop_written == 'MicelleSize') THEN

         IF (block_average) THEN
            write_buff(ii+1) = cluster%Mave / REAL(ncluster_freq,DP)
         ELSE
            write_buff(ii+1) = cluster%Mave
         END IF

      ELSE IF (prop_written == 'OligNNdist') THEN

         IF (block_average) THEN
            write_buff(ii+1) = cluster%olig_nn_dist / REAL(noligdist_freq,DP)
         ELSE
            write_buff(ii+1) = cluster%olig_nn_dist
         END IF

      ELSE IF (prop_written == 'Excluded_Volume') THEN

         IF (block_average) THEN
            write_buff(ii+1) = exvol%excluded / REAL(nexvol_freq * exvol%n_iter,DP)
         ELSE
            write_buff(ii+1) = exvol%excluded / REAL(exvol%n_iter)
            IF (lattice_sim) write_buff(ii+1) = exvol%excluded / REAL(exvol%n_iter * box_list(1)%volume)

         END IF

      ELSE IF (prop_written == 'Degree_Association') THEN

         IF (cluster%n_clusters == 0) THEN
            write_buff(ii+1) = 0.0
         ELSE IF (block_average) THEN
            write_buff(ii+1) = alpha%n_assoc / REAL(nalpha_freq * cluster%n_clusters,DP)
         ELSE
            write_buff(ii+1) = alpha%n_assoc / REAL(cluster%n_clusters, DP)
         END IF

!      ELSE IF (prop_written == 'Virial_Coefficient') THEN
!
!         IF (block_average) THEN
!            write_buff(ii+1) = virial%coefficient/REAL(nthermo_freq,DP)
!         ELSE
!            write_buff(ii+1) = virial%coefficient
!         END IF
!         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol
!
!      ELSE IF (prop_written == 'Effective_Potential') THEN
!
!         IF (block_average) THEN
!            write_buff(ii+1) = virial%effective/REAL(nthermo_freq,DP)
!         ELSE
!            write_buff(ii+1) = virial%effective
!         END IF
!         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

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
              WRITE(M_XYZ_unit,'(A3,F20.13,F20.13,F20.13,I5,I5)') &
              nonbond_list(ia,is)%element, &
                   atom_list(ia,this_im,is)%rxp, &
                   atom_list(ia,this_im,is)%ryp, &
                   atom_list(ia,this_im,is)%rzp, &
                   is, im
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE Write_Coords

SUBROUTINE Write_MSD(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        pp_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE Transport_Properties
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: it, box_unit, is
  DOUBLE PRECISION :: dstep, norm

  DO is = 1, nspecies
     IF (.not. trans%msd_species(is)) CYCLE

     msd_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.species' // TRIM(Int_To_String(is)) // '.msd'
     box_unit = msd_file_unit + this_box
     OPEN(unit=msd_file_unit+this_box, file=msd_file)
   
     WRITE(box_unit,'(A)') '#  time(ps)    msd(A^2) -     x           y           z'
           
     dstep = trans%sim_step * trans%sim_freq

     DO it = 1, trans%ntel
        norm = DBLE(nmolecules(is))*DBLE(trans%ntime(it))
        WRITE(box_unit,'(5E12.5)') DBLE(it-1) * dstep, &
                                   (trans%x_msd(is, it) + trans%y_msd(is, it) + trans%z_msd(is, it)) / norm, &
                                   trans%x_msd(is, it) / norm, &
                                   trans%y_msd(is, it) / norm, &
                                   trans%z_msd(is, it) / norm
                                   !trans%msd(is, it) / norm, &
   
     END DO

     CLOSE(unit=msd_file_unit+this_box)
  END DO

END SUBROUTINE Write_MSD

SUBROUTINE Write_VACF(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        pp_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE Transport_Properties
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: it, box_unit, is
  DOUBLE PRECISION :: dstep, norm

  DO is = 1, nspecies
     IF (.not. trans%vacf_species(is)) CYCLE

     vacf_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.species' // TRIM(Int_To_String(is)) // '.vacf'
     box_unit = vacf_file_unit + this_box
     OPEN(unit=vacf_file_unit+this_box, file=vacf_file)
   
     WRITE(box_unit,'(A)') '#  time(ps)    vacf(A/ps) -     x           y           z'
           
     dstep = trans%sim_step * trans%sim_freq

     DO it = 1, trans%ntel
        norm = DBLE(nmolecules(is))*DBLE(trans%ntime(it))
        WRITE(box_unit,'(5E14.5)') DBLE(it-1) * dstep, &
                                   trans%vacf(is, it) / norm, &
                                   trans%x_vacf(is, it) / norm, &
                                   trans%y_vacf(is, it) / norm, &
                                   trans%z_vacf(is, it) / norm
                                   !trans%vacf(is, it) / norm, &
   
     END DO

     CLOSE(unit=vacf_file_unit+this_box)
  END DO

END SUBROUTINE Write_VACF

SUBROUTINE Write_Cluster(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        pp_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE Cluster_Routines
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: iM, box_unit, is

  cluster_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.clu'
  box_unit = cluster_file_unit + this_box
  OPEN(unit=cluster_file_unit+this_box, file=cluster_file)

  WRITE(box_unit,'(A)', ADVANCE='NO') '#    M       pop'

  IF (nalphaclus_freq > 0) THEN
     WRITE(box_unit,'(A)', ADVANCE='NO') '   bound-counterions'
  END IF

  IF (ncluslife_freq > 0) THEN
     WRITE(box_unit,'(5A)', ADVANCE='NO') '  lifetime       n_birth   n_death  ',&
                                          'fusion M->1   fusion M->2  fusion M->3  fusion M->4  ',&
                                          'fusion M->5  fusion M->6  fusion M->7  fusion M->8  ', &
                                          'fission M->1   fission M->2  fission M->3  fission M->4  ',&
                                          'fission M->5  fission M->6  fission M->7  fission M->8'
  END IF

  IF (nendclus_freq > 0) THEN
     WRITE(box_unit,'(A)', ADVANCE='NO') '   End-to-end distance (A) for species:'
     DO is = 1, nspecies
        IF (.not. measure_mol%end2end_spec(is) ) CYCLE
        WRITE(box_unit,'(I10)', ADVANCE='NO') is
     END DO
  END IF
        
  DO iM = 1, SIZE(cluster%M)
     IF (cluster%M(iM) > 0) THEN
        WRITE(box_unit,'(/)', ADVANCE='NO')
        WRITE(box_unit,'(I6, I10)', ADVANCE='NO') iM, cluster%M(iM)

        IF (nalphaclus_freq > 0) THEN
            WRITE(box_unit,'(I10)', ADVANCE='NO') alpha%n_assoc_clus(iM)
        END IF
        
        IF (ncluslife_freq > 0) THEN
            WRITE(box_unit,'(19I10)', ADVANCE='NO') cluster%lifetime(iM), cluster%n_clus_birth(iM), cluster%n_clus_death(iM), &
                                      cluster%fusion(iM,1), cluster%fusion(iM,2), cluster%fusion(iM,3), cluster%fusion(iM,4), &
                                      cluster%fusion(iM,5), cluster%fusion(iM,6), cluster%fusion(iM,7), cluster%fusion(iM,8), &
                                      cluster%fission(iM,1), cluster%fission(iM,2), cluster%fission(iM,3), cluster%fission(iM,4), &
                                      cluster%fission(iM,5), cluster%fission(iM,6), cluster%fission(iM,7), cluster%fission(iM,8) 
        END IF
        
        IF (nendclus_freq > 0) THEN
            DO is = 1, nspecies
                IF (.not. measure_mol%end2end_spec(is) ) CYCLE
                WRITE(box_unit,'(E12.5)', ADVANCE='NO') measure_mol%end2end(iM, is) / cluster%M(iM)
            END DO
        END IF
        
     END IF
  END DO
  CLOSE(unit=cluster_file_unit+this_box)

END SUBROUTINE Write_Cluster

SUBROUTINE Write_Bond(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        pp_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: box_unit, is, i, im, ib, ib_bin
  CHARACTER(5) :: Tatom
  REAL(DP)  :: length, bin_width

  DO is = 1, nspecies
     IF( .NOT. ANY(measure_mol%bond_spec(:,is) )) CYCLE
     bond_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.species' // TRIM(Int_To_String(is)) // '.bond'
     box_unit = bond_file_unit + this_box
     OPEN(unit=bond_file_unit+this_box, file=bond_file)

     WRITE(box_unit,'(A)', ADVANCE='NO') '# species molecule bond:'
     DO ib = 1, nbonds(is)
        IF (.not. measure_mol%bond_spec(ib, is) ) CYCLE
        write (Tatom, '(I5)') bond_list(ib,is)%atom2
        WRITE(box_unit,'(I5,A,A)', ADVANCE='NO') bond_list(ib,is)%atom1, '-', adjustl(Tatom)
     END DO
        
     WRITE(box_unit,'(/)', ADVANCE='NO')
  
     DO i = 1, nmolecules(is)
        im = locate(i, is)
        IF( .NOT. molecule_list(im, is)%live ) CYCLE
        IF( .NOT. ANY(measure_mol%bond_spec(:,is) )) CYCLE

        WRITE(box_unit,'(I3, I6)', ADVANCE='NO') is, im

        DO ib = 1, nbonds(is)
           IF (.not. measure_mol%bond_spec(ib, is) ) CYCLE
           WRITE(box_unit,'(A,E11.3)', ADVANCE='NO') ' ', measure_mol%bond(ib, im, is) / FLOAT(measure_mol%nbondcall)
        END DO
        IF (i < nmolecules(is)) WRITE(box_unit,'(/)', ADVANCE='NO')
     END DO
     CLOSE(unit=bond_file_unit+this_box)

     bond_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.species' // TRIM(Int_To_String(is)) // '.bond_his'
     box_unit = bond_file_unit + this_box
     OPEN(unit=bond_file_unit+this_box, file=bond_file)

     WRITE(box_unit,'(A)', ADVANCE='NO') '# bond length  bond:'
     DO ib = 1, nbonds(is)
        IF (.not. measure_mol%bond_spec(ib, is) ) CYCLE
        write (Tatom, '(I5)') bond_list(ib,is)%atom2
        WRITE(box_unit,'(I5,A,A)', ADVANCE='NO') bond_list(ib,is)%atom1, '-', adjustl(Tatom)
     END DO
        
     bin_width = (8.0/3.0) * measure_mol%l0ave(is) / FLOAT(measure_mol%nb_bins)
     DO ib_bin = 1, measure_mol%nb_bins

        WRITE(box_unit,'(/)', ADVANCE='NO')
        length = (measure_mol%l0ave(is)/3.0) + (ib_bin * bin_width)
        WRITE(box_unit,'(E11.3)', ADVANCE='NO') length
        DO ib = 1, nbonds(is)
           IF (.not. measure_mol%bond_spec(ib, is) ) CYCLE
           WRITE(box_unit,'(I10)', ADVANCE='NO') measure_mol%bond_his(ib_bin, ib, is) 
        END DO
     END DO
     CLOSE(unit=bond_file_unit+this_box)
  END DO


END SUBROUTINE Write_Bond

SUBROUTINE Write_Angle(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        pp_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: box_unit, is, i, im, ia, ia_bin
  CHARACTER(5) :: Tatom
  REAL(DP)  :: bin_width

  DO is = 1, nspecies
     IF( .NOT. ANY(measure_mol%angle_spec(:,is) )) CYCLE
     angle_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.species' // TRIM(Int_To_String(is)) // '.angle'
     box_unit = angle_file_unit + this_box
     OPEN(unit=angle_file_unit+this_box, file=angle_file)

     WRITE(box_unit,'(A)', ADVANCE='NO') '# species molecule angle:'
     DO ia = 1, nangles(is)
        IF (.not. measure_mol%angle_spec(ia, is) ) CYCLE
        write (Tatom, '(I5)') angle_list(ia,is)%atom3
        WRITE(box_unit,'(I5,A,I5,A,A)', ADVANCE='NO') angle_list(ia,is)%atom1, '-', angle_list(ia,is)%atom2, '-', adjustl(Tatom)
     END DO
        
     WRITE(box_unit,'(/)', ADVANCE='NO')
  
     DO i = 1, nmolecules(is)
        im = locate(i, is)
        IF( .NOT. molecule_list(im, is)%live ) CYCLE
        IF( .NOT. ANY(measure_mol%angle_spec(:,is) )) CYCLE

        WRITE(box_unit,'(I3, I6)', ADVANCE='NO') is, im

        DO ia = 1, nangles(is)
           IF (.not. measure_mol%angle_spec(ia, is) ) CYCLE
           WRITE(box_unit,'(A,E11.4)', ADVANCE='NO') ' ', measure_mol%angle(ia, im, is) / FLOAT(measure_mol%nanglecall)
        END DO
        IF (i < nmolecules(is)) WRITE(box_unit,'(/)', ADVANCE='NO')
     END DO
     CLOSE(unit=angle_file_unit+this_box)

     angle_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.species' // TRIM(Int_To_String(is)) // '.angle_his'
     box_unit = angle_file_unit + this_box
     OPEN(unit=angle_file_unit+this_box, file=angle_file)

     WRITE(box_unit,'(A)', ADVANCE='NO') '# angle length  angle:'
     DO ia = 1, nangles(is)
        IF (.not. measure_mol%angle_spec(ia, is) ) CYCLE
        write (Tatom, '(I5)') angle_list(ia,is)%atom3
        WRITE(box_unit,'(I5,A,I5,A,A)', ADVANCE='NO') angle_list(ia,is)%atom1, '-', angle_list(ia,is)%atom2, '-', adjustl(Tatom)
     END DO
        
     bin_width = twoPI / FLOAT(measure_mol%na_bins)
     DO ia_bin = 1, measure_mol%na_bins

        WRITE(box_unit,'(/)', ADVANCE='NO')
        WRITE(box_unit,'(E11.3)', ADVANCE='NO') (ia_bin * bin_width)
        DO ia = 1, nangles(is)
           IF (.not. measure_mol%angle_spec(ia, is) ) CYCLE
           WRITE(box_unit,'(I10)', ADVANCE='NO') measure_mol%angle_his(ia_bin, ia, is) 
        END DO
     END DO
     CLOSE(unit=angle_file_unit+this_box)
  END DO


END SUBROUTINE Write_Angle

SUBROUTINE Write_Dihedral(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        pp_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: box_unit, is, i, im, id, id_bin
  CHARACTER(5) :: Tatom
  REAL(DP)  :: bin_width

  DO is = 1, nspecies
     IF( .NOT. ANY(measure_mol%dihedral_spec(:,is) )) CYCLE
     dihedral_file = TRIM(run_name)//'.box'//TRIM(Int_To_String(this_box))//'.species'//TRIM(Int_To_String(is))//'.dihedral'
     box_unit = dihedral_file_unit + this_box
     OPEN(unit=dihedral_file_unit+this_box, file=dihedral_file)

     WRITE(box_unit,'(A)', ADVANCE='NO') '# species molecule dihedral:'
     DO id = 1, ndihedrals(is)
        IF (.not. measure_mol%dihedral_spec(id, is) ) CYCLE
        write (Tatom, '(I5)') dihedral_list(id,is)%atom4
        WRITE(box_unit,'(I5,A,I5,A,I5,A,A)', ADVANCE='NO') dihedral_list(id,is)%atom1, '-', &
                                                           dihedral_list(id,is)%atom2, '-', &
                                                           dihedral_list(id,is)%atom3, '-', &
                                                           adjustl(Tatom)
     END DO
        
     WRITE(box_unit,'(/)', ADVANCE='NO')
  
     DO i = 1, nmolecules(is)
        im = locate(i, is)
        IF( .NOT. molecule_list(im, is)%live ) CYCLE
        IF( .NOT. ANY(measure_mol%dihedral_spec(:,is) )) CYCLE

        WRITE(box_unit,'(I3, I6)', ADVANCE='NO') is, im

        DO id = 1, ndihedrals(is)
           IF (.not. measure_mol%dihedral_spec(id, is) ) CYCLE
           WRITE(box_unit,'(A,E11.3)', ADVANCE='NO') ' ', measure_mol%dihedral(id, im, is) / FLOAT(measure_mol%ndihedralcall)
        END DO
        IF (i < nmolecules(is)) WRITE(box_unit,'(/)', ADVANCE='NO')
     END DO
     CLOSE(unit=dihedral_file_unit+this_box)

     dihedral_file = TRIM(run_name)//'.box'//TRIM(Int_To_String(this_box))//'.species'//TRIM(Int_To_String(is))//'.dihedral_his'
     box_unit = dihedral_file_unit + this_box
     OPEN(unit=dihedral_file_unit+this_box, file=dihedral_file)

     WRITE(box_unit,'(A)', ADVANCE='NO') '# dihedral length  dihedral:'
     DO id = 1, ndihedrals(is)
        IF (.not. measure_mol%dihedral_spec(id, is) ) CYCLE
        write (Tatom, '(I5)') dihedral_list(id,is)%atom4
        WRITE(box_unit,'(I5,A,I5,A,I5,A,A)', ADVANCE='NO') dihedral_list(id,is)%atom1, '-', &
                                                           dihedral_list(id,is)%atom2, '-', &
                                                           dihedral_list(id,is)%atom3, '-', &
                                                           adjustl(Tatom)
     END DO
        
     bin_width = twoPI / FLOAT(measure_mol%nd_bins)
     DO id_bin = 1, measure_mol%nd_bins

        WRITE(box_unit,'(/)', ADVANCE='NO')
        WRITE(box_unit,'(E11.3)', ADVANCE='NO') (id_bin * bin_width) - PI
        DO id = 1, ndihedrals(is)
           IF (.not. measure_mol%dihedral_spec(id, is) ) CYCLE
           WRITE(box_unit,'(I10)', ADVANCE='NO') measure_mol%dihedral_his(id_bin, id, is) 
        END DO
     END DO
     CLOSE(unit=dihedral_file_unit+this_box)
  END DO


END SUBROUTINE Write_Dihedral

SUBROUTINE Write_Atom_Distribution(this_box)
  !************************************************************************************
  ! The subroutine writes the cluster vists in the simulation box
  !
  ! CALLED BY
  !
  !        pp_driver
  !
  !************************************************************************************

  USE Run_Variables
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  INTEGER :: box_unit
  INTEGER :: iap, is, ia, ja, i, im, iap_bin
  REAL(DP)  :: bin_width

  a_dist_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.atomdist'
  box_unit = a_dist_file_unit + this_box
  OPEN(unit=a_dist_file_unit+this_box, file=a_dist_file)

  WRITE(box_unit,'(A)') '# atom-pair-type ia ja:'
  DO iap = 1, measure_mol%natom_dists
      ia = measure_mol%a_dist_pairs(iap, 2)
      ja = measure_mol%a_dist_pairs(iap, 4)
      WRITE(box_unit,'(A,3I3)') '# ', iap, ia, ja
  END DO
  WRITE(box_unit,'(A)', ADVANCE='NO') '# is im js jm distance'
  
  DO is = 1, nspecies
      IF ( .NOT. ANY(measure_mol%a_dist_pairs(:, 1) == is) ) CYCLE
      DO i = 1, nmolecules(is)
          im = locate(i, is)
        IF (.TRUE.) THEN
            print*, im, cluster%N(cluster%clabel(im, is))
        END IF 
          IF( .NOT. molecule_list(im, is)%live ) CYCLE
          WRITE(box_unit,'(/)', ADVANCE='NO')
          WRITE(box_unit,'(I3, I6)', ADVANCE='NO') is, im 
          DO iap = 1, measure_mol%natom_dists
              WRITE(box_unit,'(E13.5)', ADVANCE='NO') ( measure_mol%a_dist_sq(iap, im, is) )**(0.5) / FLOAT(measure_mol%nadistcall)
              !WRITE(box_unit,'(E11.3)', ADVANCE='NO') ( measure_mol%a_dist_sq(iap, im, is) )**(0.5) / FLOAT(measure_mol%nadistcall * nmolecules(is))
          END DO
     END DO
  END DO
  CLOSE(unit=a_dist_file_unit+this_box)

  a_dist_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(this_box)) // '.atomdist_his'
  box_unit = a_dist_file_unit + this_box
  OPEN(unit=a_dist_file_unit+this_box, file=a_dist_file)

  WRITE(box_unit,'(A)', ADVANCE='NO') '# distance atom pair:'
  DO iap = 1, measure_mol%natom_dists
     WRITE(box_unit,'(I5)', ADVANCE='NO') iap
  END DO
     
  
  bin_width =  measure_mol%a_dist_max_sq / FLOAT(measure_mol%nad_bins)
  DO iap_bin = 1, measure_mol%nad_bins

     WRITE(box_unit,'(/)', ADVANCE='NO')
     WRITE(box_unit,'(E11.4)', ADVANCE='NO') (iap_bin * bin_width)
     DO iap = 1, measure_mol%natom_dists
        WRITE(box_unit,'(I10)', ADVANCE='NO') measure_mol%a_dist_his(iap_bin, iap)
     END DO
  END DO
  CLOSE(unit=a_dist_file_unit+this_box)

END SUBROUTINE Write_Atom_Distribution

