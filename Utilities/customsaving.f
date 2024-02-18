            ! Author: Konstantinos
            ! Date: 8 Feb 2024
            ! Add to run_star_extras_evol_Mors0.f

         ! ## Konstantinos' Custom Saving ##############
         ! Declare variables for custom saving
         ! Needs to be on top of the function for some godforsaken reason
          logical :: flag
          real, parameter :: tolerance = 0.1  ! Set your desired tolerance
          real :: agenow
          real, dimension(100) :: wanted_ages  ! Adjust the array size as needed
          integer :: i
          
          ! ## Konstantinos' Custom Saving ##############
          ! Initialize flag to false
          flag = .false.
          
          ! linspace in log10(yrs)
          agenow = s% star_age
          ! write (*,*) agenow
            wanted_ages = [1.00000000e+07, 1.10909091e+08, 2.11818182e+08, 3.12727273e+08, &
       4.13636364e+08, 5.14545455e+08, 6.15454545e+08, 7.16363636e+08, &
       8.17272727e+08, 9.18181818e+08, 1.01909091e+09, 1.12000000e+09, &
       1.22090909e+09, 1.32181818e+09, 1.42272727e+09, 1.52363636e+09, &
       1.62454545e+09, 1.72545455e+09, 1.82636364e+09, 1.92727273e+09, &
       2.02818182e+09, 2.12909091e+09, 2.23000000e+09, 2.33090909e+09, &
       2.43181818e+09, 2.53272727e+09, 2.63363636e+09, 2.73454545e+09, &
       2.83545455e+09, 2.93636364e+09, 3.03727273e+09, 3.13818182e+09, &
       3.23909091e+09, 3.34000000e+09, 3.44090909e+09, 3.54181818e+09, &
       3.64272727e+09, 3.74363636e+09, 3.84454545e+09, 3.94545455e+09, &
       4.04636364e+09, 4.14727273e+09, 4.24818182e+09, 4.34909091e+09, &
       4.45000000e+09, 4.55090909e+09, 4.65181818e+09, 4.75272727e+09, &
       4.85363636e+09, 4.95454545e+09, 5.05545455e+09, 5.15636364e+09, &
       5.25727273e+09, 5.35818182e+09, 5.45909091e+09, 5.56000000e+09, &
       5.66090909e+09, 5.76181818e+09, 5.86272727e+09, 5.96363636e+09, &
       6.06454545e+09, 6.16545455e+09, 6.26636364e+09, 6.36727273e+09, &
       6.46818182e+09, 6.56909091e+09, 6.67000000e+09, 6.77090909e+09, &
       6.87181818e+09, 6.97272727e+09, 7.07363636e+09, 7.17454545e+09, &
       7.27545455e+09, 7.37636364e+09, 7.47727273e+09, 7.57818182e+09, &
       7.67909091e+09, 7.78000000e+09, 7.88090909e+09, 7.98181818e+09, &
       8.08272727e+09, 8.18363636e+09, 8.28454545e+09, 8.38545455e+09, &
       8.48636364e+09, 8.58727273e+09, 8.68818182e+09, 8.78909091e+09, &
       8.89000000e+09, 8.99090909e+09, 9.09181818e+09, 9.19272727e+09, &
       9.29363636e+09, 9.39454545e+09, 9.49545455e+09, 9.59636364e+09, &
       9.69727273e+09, 9.79818182e+09, 9.89909091e+09, 1.00000000e+10]
          
          ! Check if A is within tolerance of any element in B
          do i = 1, size(wanted_ages)
             if (  abs( agenow - wanted_ages(i) ) / agenow <= tolerance) then
                flag = .true.
                exit  ! Exit the loop once a match is found
             end if
          end do

          ! if flag is true, then save a profile
          if (flag) then
            !write(*,*) 'flag is true'
            s% need_to_save_profiles_now = .true.
          end if 
