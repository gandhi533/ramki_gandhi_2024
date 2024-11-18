      program iteration
c     The sign of kappa*u1 has been changed 24 Oct 2023
c     This takes input of tau, kappa, u1 and xi to calculate conc
c     This does calculations for different xi results 2 Nov 2023
c     input file has been made as input_conc_plot.dat to allow plotting for different times
      
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
       open (unit=11,file='input_conc_plot.dat') !This is same as input_iteration.dat except that it allows to compute upto different tau_max      
c       open  (unit=11,file='input_iteration.dat')
       
c     Try with eigen values less than 20
       read (11,*)  j_ev_max,akappa_in,u1_in,
     ccap_lambda,a,delta_tau,tau_max,iter_max  

c     opening a file for plotting data
c       open (unit=16,file='plot_conc_time.dat')

c     opening file for writing concentrations in a row into a common file
       open (unit=17,file='output_all_xi_tau.dat')

c     opening file for reading output of v_4.f
       open (unit=15,file='output_plot_tau_u1.dat')
        i_tau_max=tau_max/delta_tau
              do i=1,i_tau_max
       read (15,*) tau,k,u1(i),akappa(i),u_dot_z_zero(i)
       write (*,*) tau,k,u1(i), akappa(i),u_dot_z_zero(i)
              enddo

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
         xi=0.0
 100  xi=xi+5.d-02 
       do i=1, i_tau_max
         sum_ev_zj1(i)=0.
         sum_ev_zj1xi(i)=0.
       do j=1,j_ev_max
          sum_ev_zj1(i)=sum_ev_zj1(i)+ef1(j)*ef1(j)*
     c             dexp(-ev(j)*dble(i)*delta_tau)
          sum_ev_zj1xi(i)=sum_ev_zj1xi(i)+ef1(j)*const_norm(j)*
     c      dcos(xi*beta(j))*dexp(-ev(j)*dble(i)*delta_tau)
          enddo
          enddo

c we need these at tau = 0

          sum_ev_zj1_zero= 0.
          sum_ev_zj1xi_zero=0.
       do j=1,j_ev_max
          sum_ev_zj1_zero=sum_ev_zj1_zero+ef1(j)*ef1(j)
c   sum_ev_zj1_zero=alpha*alpha/(alpha+1)!, this is the theoretical value
          sum_ev_zj1xi_zero=sum_ev_zj1xi_zero+ef1(j)*const_norm(j)*
     c      dcos(xi*beta(j))
          enddo

c calculate concentrations          
c     calculating concentration. i=1 is special
c     calculation of g 
        g(1)=1.d0*dexp(-delta_tau*0.5d00*
     c         (akappa_in+akappa(1)))

c     calculation of integral of ku+Lambda_e/g
           cap_a(1)=delta_tau*0.5d00*(
     c (cap_lambda_e-akappa_in*u1_in)+
     c (cap_lambda_e- akappa(1)*u1(1))/g(1))

           u_dot_z_zero(1) = g(1)*dsqrt((alpha+1)/alpha)- g(1)*cap_a(1)/
     c dsqrt(alpha*(alpha+1))

          u_dot_zjxi(1)=-(g(1)/alpha)*0.5d00*delta_tau*(
     c         (-akappa_in*u1_in+cap_lambda_e)*sum_ev_zj1xi(1)+
     c         (-akappa(1)*u1(1)+cap_lambda_e)*sum_ev_zj1xi_zero/g(1))

          
         conc(1)=dsqrt(alpha/(alpha+1))*u_dot_z_zero(1)+u_dot_zjxi(1)   
         
c     calculating concentration for other i
       do i=2,i_tau_max
c     calculation of g 
        g(i)=g(i-1)*dexp(-delta_tau*0.5d00*
     c         (akappa(i-1)+akappa(i)))

c     calculation of integral of ku+Lambda_e/g
           cap_a(i)=cap_a(i-1)+delta_tau*0.5d00*(
     c       (cap_lambda_e-akappa(i-1)*u1(i-1))/g(i-1)+
     c (cap_lambda_e-akappa(i)*u1(i))/g(i))          

           u_dot_z_zero(i) = g(i)*dsqrt((alpha+1)/alpha)- g(i)*cap_a(i)/
     c dsqrt(alpha*(alpha+1))

          u_dot_zjxi(i)=0.5d00*delta_tau*(
     c         (-akappa_in*u1_in+cap_lambda_e)*sum_ev_zj1xi(i)+
     c         (-akappa(i)*u1(i)+cap_lambda_e)*sum_ev_zj1xi_zero/g(i))
        do k=1,i-1 
                u_dot_zjxi(i)= u_dot_zjxi(i)+delta_tau*
     c      (-akappa(k)*u1(k)+cap_lambda_e)*sum_ev_zj1xi(i-k)/g(k)
                enddo
              u_dot_zjxi(i)= -g(i)*u_dot_zjxi(i)/alpha  

        conc(i)=dsqrt(alpha/(alpha+1))*u_dot_z_zero(i) +u_dot_zjxi(i)
        enddo
         write (17,*) xi,(real(i)*delta_tau,i,conc(i), i=1,i_tau_max)
         if (xi.lt.0.94) go to 100
         write (17,*) '1.0',(real(i)*delta_tau,i,u1(i), i=1,i_tau_max)
         write (17,*) 'a =', a, 'alpha = ',alpha, 'cap Lambda',
     c  cap_lambda, 'max-ev =',j_ev_max
            stop
          end

