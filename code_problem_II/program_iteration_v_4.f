      program iteration
c     The sign of kappa*u1 has been changed 24 Oct 2023
c     This program computes only kappa and u1. created 25 Oct 2023
c     Print of u.z_o has been added on 28-04-24

      implicit none
c     vectors
      doubleprecision beta(25),ev(25),ef1(25),const_norm(25)
      doubleprecision u1(150),g(150),akappa(150),cap_a(150)
      doubleprecision sum_ev_zj1(150),sum_ev_zj1xi(150)
       doubleprecision u_dot_zj(150),u_dot_zjxi(150),u_dot_z_zero(150)
        doubleprecision conc(150)

c     scalars
       doubleprecision alpha,akappa_in,u1_in,cap_lambda,a,delta_tau
       doubleprecision tau,tau_max,xi,cap_lambda_e
      doubleprecision sum_ev_zj1_zero,sum_ev_zj1xi_zero
       doubleprecision u1_it,akappa_it,error_u1,error_kappa
       

c     Integers
       integer i,j,k,j_ev_max,iter_max,i_tau_max,iter

c     USE i FOR TIME INCREMENTS j FOR EIGEN VALUES
       
c     Opening data file to read alpha and roots of characteristic equation
       open (unit=12,file='input_test_ev.dat')
       read (12,*) alpha,(beta(j), j=1, 20) !beta are the eigenvalues of Carslaw!assuming that max ev are 20 

c     opening file to read dat required for execuation of program
       open (unit=11,file='input_iteration.dat') 

c     Try with eigen values less than 20
       read (11,*)  j_ev_max,akappa_in,u1_in,
     ccap_lambda,a,delta_tau,tau_max,iter_max

c     opening a file suitable for plotting data
       open (unit=16,file='output_plot_tau_u1.dat')

c     opening file for writing concentrations in a row into a common file
c       open (unit=17,file='output_all_xi_tau.dat',status='old',
c     cposition='append')

c     opening file for recording output
       open (unit=15,file='output_iteration_tau_u1.dat')
       write (15,*) 'alpha  =  ', alpha,'   j_ev_max  =  ',j_ev_max
c     calculating other values
       cap_lambda_e=a*cap_lambda*alpha
       do j=1,j_ev_max  
          ev(j)=beta(j)*beta(j)
          const_norm(j)=dsqrt(2.*(ev(j)+alpha*alpha)/
     c         (ev(j)+alpha+alpha*alpha))
          ef1(j)=const_norm(j)*dcos(beta(j))
       enddo
c     calculating sum_ev_zj1. This is sum of Zji**2*exp-ev(j)*tau
c     calculating sum_ev_zj1xi. This is sum of Zji zj(xi)*exp-ev(j)*tau
c       It is needed only upto max tau of the program
       i_tau_max=tau_max/delta_tau
       do i=1, i_tau_max
         sum_ev_zj1(i)=0.
c         sum_ev_zj1xi(i)=0.
       do j=1,j_ev_max
          sum_ev_zj1(i)=sum_ev_zj1(i)+ef1(j)*ef1(j)*
     c             dexp(-ev(j)*dble(i)*delta_tau)
c          sum_ev_zj1xi(i)=sum_ev_zj1xi(i)+ef1(j)*const_norm(j)*
c     c      dcos(xi*beta(j))*dexp(-ev(j)*dble(i)*delta_tau)
          enddo
          enddo

c we need these at tau = 0

          sum_ev_zj1_zero= 0.
       do j=1,j_ev_max
          sum_ev_zj1_zero=sum_ev_zj1_zero+ef1(j)*ef1(j)
c   sum_ev_zj1_zero=alpha*alpha/(alpha+1)!, this is the theoretical value
          enddo
c          write (*,*) sum_ev_zj1xi(1),sum_ev_zj1xi_zero

c     Now we need to integrate with time

       i=1
c     This is the trial value. akappa_it is the result of trial
       akappa(1)=akappa_in
c     u1_it is the trial value. u1_it is the result of trial
       u1(1)=u1_in
       iter=0
 20    continue
       iter=iter+1
c     calculation of g 
        g(1)=1.d0*dexp(-delta_tau*0.5d00*
     c         (akappa_in+akappa(1)))

c     calculation of integral of ku+Lambda_e/g
           cap_a(1)=delta_tau*0.5d00*(
     c (cap_lambda_e-akappa_in*u1_in)+
     c (cap_lambda_e- akappa(1)*u1(1))/g(1))

           u_dot_z_zero(1) = g(1)*dsqrt((alpha+1)/alpha)- g(1)*cap_a(1)/
     c dsqrt(alpha*(alpha+1))
c        write (*,*) u_dot_z_zero(1)
c     calculation of term involving exp(lambda_j*tau)zj1-square
           u_dot_zj(1)=-(g(1)/alpha)*0.5d00*delta_tau*(
     c (-akappa_in*u1_in+cap_lambda_e)*sum_ev_zj1(1)+
     c (-akappa(1)*u1(1)+cap_lambda_e)*sum_ev_zj1_zero/g(1))

c     calculate the results of the assumption
           u1_it=dsqrt(alpha/(alpha+1))*u_dot_z_zero(1) +u_dot_zj(1)
           akappa_it= cap_lambda/(dsqrt(alpha/(alpha+1))*u_dot_z_zero(1) 
     c   -u_dot_zj(1)/alpha)

c     testing the accurqacy
           error_u1=dabs(u1_it-u1(1))/dabs(u1_it)
           error_kappa=dabs(akappa_it-akappa(1))/dabs(akappa_it)
c           write (*,*) 'error_u1',error_u1,'error_kappa', error_kappa
           if (error_u1.gt.1d-05) go to 30
           if (error_kappa.lt.1d-05) go to 41
c     we come here if accuracy is low. redefine variables. go back to
c     iteration if iter_max is not exceeded 
 30       u1(1)=u1_it
          akappa(1)=akappa_it
          if (iter.gt.iter_max)  then
             write (*,*) 'Iteration did not converge'
           write (*,*) 'i = ', i, ',  u1(i)  =  ',u1(i), 
     c    ', kappa(i)  =', akappa(i)
c     we come here it iter_max exceeded. quit
            go to 40
             else
c     we come here if we have to continue iteration
          go to 20
          endif
c     we come here if iteration converged. print and proceed with i=2

 41       continue
          tau=real(i)*delta_tau
        write (15,*) 'tau  =',tau,'i = ', i, ',  u1(i)  =  ',u1(i), 
     c    ', kappa(i)  =', akappa(i),'u_dot_z_zero(i)',u_dot_z_zero(i)
        write (16,*) tau,i,u1(i), akappa(i),u_dot_z_zero(i)
c     Now go on to further time intervals
 51     i=i+1
c     setting intial guesses to be the previous converged value
          akappa(i)=akappa(i-1)
          u1(i)=u1(i-1)
          iter=0
 50    continue
       iter=iter+1
c     calculation of g 
        g(i)=g(i-1)*dexp(-delta_tau*0.5d00*
     c         (akappa(i-1)+akappa(i)))

c     calculation of integral of ku+Lambda_e/g
           cap_a(i)=cap_a(i-1)+delta_tau*0.5d00*(
     c       (cap_lambda_e-akappa(i-1)*u1(i-1))/g(i-1)+
     c (cap_lambda_e-akappa(i)*u1(i))/g(i))          
c     c (cap_lambda_e+akappa(i)*u1(i))/g(i))

           u_dot_z_zero(i) = g(i)*dsqrt((alpha+1)/alpha)- g(i)*cap_a(i)/
     c dsqrt(alpha*(alpha+1))

           u_dot_zj(i)=0.5d00*delta_tau*(
     c (-akappa_in*u1_in+cap_lambda_e)*sum_ev_zj1(i)+
     c (-akappa(i)*u1(i)+cap_lambda_e)*sum_ev_zj1_zero/g(i))
        do k=1,i-1    
           u_dot_zj(i)=u_dot_zj(i)+delta_tau*
     c (-akappa(k)*u1(k)+cap_lambda_e)*sum_ev_zj1(i-k)/g(k)
           enddo
           u_dot_zj(i)=-g(i)*u_dot_zj(i)/alpha

c     calculate the results of the assumption
           u1_it=dsqrt(alpha/(alpha+1))*u_dot_z_zero(i) +u_dot_zj(i)
           akappa_it= cap_lambda/(dsqrt(alpha/(alpha+1))*u_dot_z_zero(i) 
     c   -u_dot_zj(i)/alpha)   

c     testing the accurqacy
           error_u1=dabs(u1_it-u1(i))/dabs(u1_it)
           error_kappa=dabs(akappa_it-akappa(i))/dabs(akappa_it)
c           write (*,*) 'error_u1',error_u1,'error_kappa', error_kappa
           if (error_u1.gt.1d-05) go to 31
           if (error_kappa.lt.1d-05) go to 42

c     we come here if accuracy is low. redefine variables. go back to
c     iteration if iter_max is not exceeded 
 31       u1(i)=u1_it
          akappa(i)=akappa_it 
          if (iter.gt.iter_max)  then
             write (*,*) iter, 'Iteration did not converge'
          write (*,*) 'i = ', i, ',  u1(i)  =  ',u1(i), ', kappa(i)  =',
     c  akappa(i)
c     we come here it iter_max exceeded. quit
             go to 40
             else
c     we come here if we have to continue iteration
          go to 50
          endif

c     we come here if iteration converged. print and proceed with next i
 42        continue
          tau=real(i)*delta_tau
        write (15,*) 'tau   =', tau,'i = ', i,'u1(i)  =  ',u1(i),
     c  'kappa(i)  =', akappa(i), 'u_dot_z_zero(i)', u_dot_z_zero(i)
        write (16,*) tau,i,u1(i),akappa(i), u_dot_z_zero(i)
c     check if we have reached maximum tau.
c        k=i-tau_max/delta_tau
       if ((i-tau_max/delta_tau).lt.0) go to 51

       write (15,*) 'a =', a, 'alpha = ',alpha, 'max-ev ='
     c  ,j_ev_max, 'cap Lambda = ', cap_lambda 
      write (16,*) 'a =', a, 'alpha = ',alpha, 'max-ev ='
     c  ,j_ev_max, 'cap Lambda = ', cap_lambda 
c         We come here only when we want to quit
 40   continue


           stop
          end
