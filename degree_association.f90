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

MODULE Degree_Association

  !******************************************************************************
  ! 
  ! - Andrew P. Santos
  !
  !*******************************************************************************

  USE Run_Variables
  USE IO_Utilities
  USE Cluster_Routines

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Calculate_Degree_Association(this_box)

    !*********************************************************************************
    !
    ! Calculate the Degree of ion Association
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: i, j, am, cm
    INTEGER :: ierr, line_nbr, nbr_entries
    REAL(DP) :: rxij, ryij, rzij, rijsq, rxijp, ryijp, rzijp

    alpha%n_assoc = 0
    DO i = 1, nmolecules(alpha%assoc_species)
        am = locate(i, alpha%assoc_species)
        IF( .NOT. molecule_list(am, alpha%assoc_species)%live ) CYCLE

        DO j = 1, nmolecules(alpha%clus_species)
            cm = locate(j, alpha%clus_species)
            IF( .NOT. molecule_list(cm, alpha%clus_species)%live ) CYCLE

            ! IF the molecule is not in a cluster of the specified size
            IF (cluster%N( cluster%clabel(cm, alpha%clus_species) ) <= cluster%M_olig(alpha%clus_species)) CYCLE

            ! Get the positions of the COM of the two molecule species
            rxijp = atom_list(alpha%atype(alpha%assoc_species), am, alpha%assoc_species)%rxp - &
                    atom_list(alpha%atype(alpha%clus_species), cm, alpha%clus_species)%rxp
            ryijp = atom_list(alpha%atype(alpha%assoc_species), am, alpha%assoc_species)%ryp - &
                    atom_list(alpha%atype(alpha%clus_species), cm, alpha%clus_species)%ryp
            rzijp = atom_list(alpha%atype(alpha%assoc_species), am, alpha%assoc_species)%rzp - &
                    atom_list(alpha%atype(alpha%clus_species), cm, alpha%clus_species)%rzp
            
            ! Now get the minimum image separation 
            CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)
    
            rijsq = rxij*rxij + ryij*ryij + rzij*rzij
            
            IF (rijsq < alpha%cutoff_sq) THEN
                alpha%n_assoc = alpha%n_assoc + 1
                EXIT
            END IF

        END DO
    END DO
                    
  END SUBROUTINE Calculate_Degree_Association

END MODULE Degree_Association
