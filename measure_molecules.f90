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

MODULE Measure_Molecules

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

  SUBROUTINE Calculate_End_To_End_Distance(this_box)

    !*********************************************************************************
    !
    ! Calculate the End_To_End Distance of species as a function of cluster size
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: i, j, is, im
    INTEGER :: clus_size
    REAL(DP) :: rxij, ryij, rzij, rijsq, rxijp, ryijp, rzijp

    DO is = 1, nspecies
        IF (.not. measure_mol%end2end_spec(is) ) CYCLE

        DO i = 1, nmolecules(is)
            im = locate(i, is)
            IF( .NOT. molecule_list(im, is)%live ) CYCLE

            ! Get the positions of the COM of the two molecule species
            rxijp = atom_list(1, im, is)%rxp - atom_list(natoms(is), im, is)%rxp
            ryijp = atom_list(1, im, is)%ryp - atom_list(natoms(is), im, is)%ryp
            rzijp = atom_list(1, im, is)%rzp - atom_list(natoms(is), im, is)%rzp
            
            ! Now get the minimum image separation 
            CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)
    
            rijsq = (rxij*rxij) + (ryij*ryij) + (rzij*rzij)
            
            clus_size = cluster%N( cluster%clabel(im, is) )
            measure_mol%end2end(clus_size, is) = measure_mol%end2end(clus_size, is) + ( SQRT(rijsq) / DBLE(clus_size) )

        END DO
    END DO
                    
  END SUBROUTINE Calculate_End_To_End_Distance

  SUBROUTINE Calculate_Bond_His(this_box)

    !*********************************************************************************
    !
    ! Calculate the End_To_End Distance of species as a function of cluster size
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: is, i, im, ib, ib_bin
    REAL(DP) :: length

    measure_mol%nbondcall = measure_mol%nbondcall + 1
    DO is = 1, nspecies
        DO ib = 1, nbonds(is)
            IF( .NOT. measure_mol%bond_spec(ib,is) ) CYCLE
            
            DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE
                CALL Get_Bond_Length(ib,im,is,length)

                measure_mol%bond(ib, im, is) = measure_mol%bond(ib, im, is) + length
                ! calculate histogram bin number
                ib_bin = FLOOR( (measure_mol%bondpre * (length / measure_mol%l0ave(is))) - measure_mol%bondadd) + 1
                ib_bin = MIN(ib_bin, measure_mol%nb_bins)
                measure_mol%bond_his(ib_bin, ib, is) = measure_mol%bond_his(ib_bin, ib, is) + 1

            END DO
        END DO
    END DO
                    
  END SUBROUTINE Calculate_Bond_His

  SUBROUTINE Calculate_Angle_His(this_box)

    !*********************************************************************************
    !
    ! Calculate the End_To_End Distance of species as a function of cluster size
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: is, i, im, ia, ia_bin
    REAL(DP) :: theta

    measure_mol%nanglecall = measure_mol%nanglecall + 1
    DO is = 1, nspecies

        DO ia = 1, nangles(is)
            IF( .NOT. measure_mol%angle_spec(ia,is) ) CYCLE
            
            DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE
                CALL Get_Bond_Angle(ia, im, is, theta)

                measure_mol%angle(ia, im, is) = measure_mol%angle(ia, im, is) + theta
                ! calculate histogram bin number
                ia_bin = FLOOR((theta / twoPI) * measure_mol%na_bins) + 1
                ia_bin = MIN(ia_bin, measure_mol%na_bins)
                measure_mol%angle_his(ia_bin, ia, is) = measure_mol%angle_his(ia_bin, ia, is) + 1

            END DO
        END DO
    END DO
                    
  END SUBROUTINE Calculate_Angle_His

  SUBROUTINE Calculate_Dihedral_His(this_box)

    !*********************************************************************************
    !
    ! Calculate the End_To_End Distance of species as a function of cluster size
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: is, i, im, id, id_bin
    REAL(DP) :: phi

    measure_mol%ndihedralcall = measure_mol%ndihedralcall + 1
    DO is = 1, nspecies

        DO id = 1, ndihedrals(is)
            IF( .NOT. measure_mol%dihedral_spec(id,is) ) CYCLE
            
            DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE
                CALL Get_Dihedral_Angle(id, im, is, phi)

                measure_mol%dihedral(id, im, is) = measure_mol%dihedral(id, im, is) + phi
                ! calculate histogram bin number
                id_bin = FLOOR(((PI + phi) / twoPI) * measure_mol%na_bins) + 1
                id_bin = MIN(id_bin, measure_mol%nd_bins)
                measure_mol%dihedral_his(id_bin, id, is) = measure_mol%dihedral_his(id_bin, id, is) + 1

            END DO
        END DO
    END DO
                    
  END SUBROUTINE Calculate_Dihedral_His

  SUBROUTINE Calculate_Atom_Distribution(this_box)

    !*********************************************************************************
    !
    ! Calculate the End_To_End Distance of species as a function of cluster size
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: iap, is, ia, js, ja, i, im, j, jm, iap_bin
    REAL(DP) :: rxij, ryij, rzij, rijsq, rxijp, ryijp, rzijp, min_dist_sq

    measure_mol%nadistcall = measure_mol%nadistcall + 1
    DO iap = 1, measure_mol%natom_dists
        is = measure_mol%a_dist_pairs(iap, 1)
        ia = measure_mol%a_dist_pairs(iap, 2)
        js = measure_mol%a_dist_pairs(iap, 3)
        ja = measure_mol%a_dist_pairs(iap, 4)

        DO i = 1, nmolecules(is)
            im = locate(i, is)
            IF( .NOT. molecule_list(im, is)%live ) CYCLE
            min_dist_sq = 1000000.0
            DO j = 1, nmolecules(js)
                jm = locate(j, js)
                IF( .NOT. molecule_list(jm, js)%live ) CYCLE
                IF( is == js .AND. im == jm ) CYCLE

                ! Get the positions of the COM of the two molecule species
                rxijp = atom_list(ia, im, is)%rxp - atom_list(ja, jm, js)%rxp
                ryijp = atom_list(ia, im, is)%ryp - atom_list(ja, jm, js)%ryp
                rzijp = atom_list(ia, im, is)%rzp - atom_list(ja, jm, js)%rzp
                
                ! Now get the minimum image separation 
                CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)
                rijsq = (rxij*rxij) + (ryij*ryij) + (rzij*rzij)
    
                IF (rijsq < min_dist_sq) min_dist_sq = rijsq
                ! calculate histogram bin number
                iap_bin = FLOOR( ((rijsq)**(0.5) / measure_mol%a_dist_max_sq) * measure_mol%nad_bins) + 1
                iap_bin = MIN(iap_bin, measure_mol%nad_bins)
                measure_mol%a_dist_his(iap_bin, iap) = measure_mol%a_dist_his(iap_bin, iap) + 1

            END DO
            measure_mol%a_dist_sq(iap, im, is) = measure_mol%a_dist_sq(iap, im, is) + min_dist_sq 
        END DO
    END DO
                    
  END SUBROUTINE Calculate_Atom_Distribution

END MODULE Measure_Molecules
