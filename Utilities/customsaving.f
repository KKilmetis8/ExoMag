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
            wanted_ages = [1.00000000e+07, 1.09854114e+07, 1.20679264e+07, 1.32571137e+07, &
       1.45634848e+07, 1.59985872e+07, 1.75751062e+07, 1.93069773e+07, &
       2.12095089e+07, 2.32995181e+07, 2.55954792e+07, 2.81176870e+07, &
       3.08884360e+07, 3.39322177e+07, 3.72759372e+07, 4.09491506e+07, &
       4.49843267e+07, 4.94171336e+07, 5.42867544e+07, 5.96362332e+07, &
       6.55128557e+07, 7.19685673e+07, 7.90604321e+07, 8.68511374e+07, &
       9.54095476e+07, 1.04811313e+08, 1.15139540e+08, 1.26485522e+08, &
       1.38949549e+08, 1.52641797e+08, 1.67683294e+08, 1.84206997e+08, &
       2.02358965e+08, 2.22299648e+08, 2.44205309e+08, 2.68269580e+08, &
       2.94705170e+08, 3.23745754e+08, 3.55648031e+08, 3.90693994e+08, &
       4.29193426e+08, 4.71486636e+08, 5.17947468e+08, 5.68986603e+08, &
       6.25055193e+08, 6.86648845e+08, 7.54312006e+08, 8.28642773e+08, &
       9.10298178e+08, 1.00000000e+09, 1.00000000e+09, 1.18367347e+09, &
       1.36734694e+09, 1.55102041e+09, 1.73469388e+09, 1.91836735e+09, &
       2.10204082e+09, 2.28571429e+09, 2.46938776e+09, 2.65306122e+09, &
       2.83673469e+09, 3.02040816e+09, 3.20408163e+09, 3.38775510e+09, &
       3.57142857e+09, 3.75510204e+09, 3.93877551e+09, 4.12244898e+09, &
       4.30612245e+09, 4.48979592e+09, 4.67346939e+09, 4.85714286e+09, &
       5.04081633e+09, 5.22448980e+09, 5.40816327e+09, 5.59183673e+09, &
       5.77551020e+09, 5.95918367e+09, 6.14285714e+09, 6.32653061e+09, &
       6.51020408e+09, 6.69387755e+09, 6.87755102e+09, 7.06122449e+09, &
       7.24489796e+09, 7.42857143e+09, 7.61224490e+09, 7.79591837e+09, &
       7.97959184e+09, 8.16326531e+09, 8.34693878e+09, 8.53061224e+09, &
       8.71428571e+09, 8.89795918e+09, 9.08163265e+09, 9.26530612e+09, &
       9.44897959e+09, 9.63265306e+09, 9.81632653e+09, 1.00000000e+10]
          
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
