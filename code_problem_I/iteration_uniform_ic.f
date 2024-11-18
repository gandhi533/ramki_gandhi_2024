        program iteration
c     This is for iterating to find K_1 and K_2
c     This is for uniform concentration ic
      implicit none

        doubleprecision  delta1,delta2,delta3,epsilon1,epsilon2
        doubleprecision alpha,beta,gamma,dummy
       doubleprecision alambda(20),b1(20),a2(20),b2(20),a3(20),b3(20)
       doubleprecision z1_int(20),z2_int(20),z3_int(20)
       doubleprecision z1z1(20,20), z3z3(20,20)
       doubleprecision coscos,sinsin,sincos !These are needed to evaluate integrals
      integer i,k,m,n,j,num_ev,num_iter,num_iter_max,i_max

c     more data
       doubleprecision a,cap_gamma,cnot,delta_tau,tau,tau_max
       doubleprecision error_1,error_2,error_3
       doubleprecision error1(20),error2(20),error3(20)
c     these are variables related to actual iteration
       doubleprecision kay_1(150),kay_2(150) !These are K_1 and K_2
       doubleprecision d01_in,d02_in,d03_in !These are d_{01}(0),d_{02}(0),d_{03}(0)
       doubleprecision d01(150),d02(150),d03(150) ! These are d_{01}(tau),d_{02}(tau),d_{03}(tau) in notes
       doubleprecision d1_in(20),d2_in(20),d3_in(20) ! These are d_{j1}(0),d_{j2}(0),d_{j3}(0)  in notes
       doubleprecision d1(150,20),d2(150,20),d3(150,20) !first is time, These are d_{j1}(tau),d_{j2}(tau),d_{j3}(tau) in notes
       doubleprecision integral
       dimension integral(150,20) ! first is time, second is eigen value, this is integral in notes


c     These are for intial guesses
       doubleprecision g_d01(150),g_d02(150),g_d03(150)
       doubleprecision g_d1(150,20),g_d2(150,20),g_d3(150,20)

c     opening file to read inputs

      open(11,file='input_plotcheq.dat')
      read (11,*) delta1,delta2,delta3,epsilon1,epsilon2

c     reading other data
      open(15,file='operating.dat')
      read (15,*) num_ev,a,cap_gamma,cnot,delta_tau,num_iter_max,tau_max
      i_max=tau_max/delta_tau
      
c     reading data for constants
       open(12,file='output_normalized_vector.dat') ! this contains all data 

c     reading data of inner product of uniform ic with kappa free eigenfunctions
       open(16,file='inner_prod_uniform_ic_diff_5bc.dat') ! 
 
c     opening file for recording ouput
       open(unit=14,file='output_iteration_general.dat')
c     record input data

       write (14,*) 'num_ev = ',num_ev, 'delta_tau = ', delta_tau
       write(14,*) ' cap_gamma = ',cap_gamma,'cnot = ',cnot

c     read eigenvlaues and constants in eigenfunctions
       do k=1,num_ev
       read(12,*) alambda(k),b1(k),a2(k),b2(k),a3(k),b3(k)
c       print *, alambda(k),b1(k),a2(k),b2(k),a3(k),b3(k)
       enddo

c     evaluate integrals of z1 and ze
              do k=1,num_ev
          alpha=dsqrt(alambda(k)/delta1)
          z1_int(k) =b1(k)*dsin(alpha*epsilon1)/alpha
          gamma=dsqrt(alambda(k)/delta2)
          z2_int(k)=a2(k)*(dcos(gamma*epsilon1) -dcos(gamma*epsilon2))
     c     /gamma
     c  +b2(k)*(dsin(gamma*epsilon2) - dsin(gamma*epsilon1))/gamma     
          beta=dsqrt(alambda(k)/delta3)  
                 z3_int(k)=a3(k)*(dcos(beta*epsilon2) -dcos(beta))/beta
     c  +b3(k)*(dsin(beta) - dsin(beta*epsilon2))/beta
c            print*,  k, z1_int(k),z2_int(k),z3_int(k)
                 enddo

c     evaluate the z1z1 integrals
       do k=1,num_ev
          alpha=dsqrt(alambda(k)/delta1)
          do m=1,num_ev
         beta=dsqrt(alambda(m)/delta1)           
         z1z1(k,m)=b1(k)*b1(m)*coscos(alpha,beta,0.d0,epsilon1)
c         print *, k,m,z1z1(k,m)
         enddo
         enddo

c     evaluate the z3z3 integrals
       do k=1,num_ev
          alpha=dsqrt(alambda(k)/delta3)
          do m=1,num_ev
         beta=dsqrt(alambda(m)/delta3)          
         z3z3(k,m)=a3(k)*a3(m)*sinsin(alpha,beta,epsilon2,1.d0)
     c     +a3(k)*b3(m)*sincos(alpha,beta,epsilon2,1.d0)     
     c      +a3(m)*b3(k)*sincos(beta,alpha,epsilon2,1.d0)     
     c       +b3(k)*b3(m)*coscos(alpha,beta,epsilon2,1.d0)     
c        print *, k,m,z3z3(k,m)
            enddo
            enddo
c     Setting initial values
            read (16,*) dummy, d01_in,d02_in,d03_in
            do j=1,num_ev
            read (16,*) dummy,d1_in(j),d2_in(j),d3_in(j)
         enddo           

c     begin iteration variables

              i=0
c     starting with i=1
 2            i=i+1
          if (i.eq.1) then    
c     make guesses for zero eigen value at first time step
           d01(i)=d01_in
           d02(i)=d02_in
           d03(i)=d03_in
c     make guesses for other eigen values at first time step
          do j=1,num_ev
              d1(i,J)=d1_in(j)
              d2(i,j)=d2_in(j)
              d3(i,j)=d3_in(j)
              enddo
c ------------------------ guesses for t=0 ends
          else
c     make guesses for zero eigen value at other times
           d01(i)=d01(i-1)
           d02(i)=d02(i-1)
           d03(i)=d03(i-1)
c     make guesses for other eigen values  at other  times      
          do j=1,num_ev
              d1(i,J)=d1(i-1,j)
              d2(i,j)=d2(i-1,j)
              d3(i,j)=d3(i-1,j)
              enddo          
          endif    
              num_iter=1
c     this is where iteration starts
 1            continue

c     computing the integral
              if (i.eq.1) then
          do j=1,num_ev
             integral(i,j)=(0.5d0*delta_tau*cap_gamma/cnot)*
     c  (d1_in(j)/d01_in+a*d3_in(j)/d03_in
     c   +(d1(i,j)/d01(i)+a*d3(i,j)/d03(i))
     c   *dexp(alambda(j)*real(i)*delta_tau))
             enddo
           else
          do j=1,num_ev
              integral(i,j)=integral(i-1,j)+
     c (0.5d0*delta_tau*cap_gamma/cnot)*(
     c (d1(i-1,j)/d01(i-1)+a*d3(i-1,j)/d03(i-1))*
     c   dexp(alambda(j)*real(i-1)*delta_tau) +
     c (d1(i,j)/d01(i)+a*d3(i,j)/d03(i))*
     c   dexp(alambda(j)*real(i)*delta_tau))
           enddo
        endif

c     update guesses for zero eigen vlaue
          g_d01(i)=(d01_in+d02_in+d03_in-cap_gamma*(1.d0+a)*
     c  real(i)*delta_tau/cnot)*epsilon1 
          g_d03(i)=(d01_in+d02_in+d03_in-cap_gamma*(1.d0+a)
     c  *real(i)*delta_tau/cnot)*(1.d0-epsilon2) 
           do j=1,num_ev
           g_d01(i)=g_d01(i)
     c    +z1_int(j)*dexp(-alambda(j)*real(i)*delta_tau)*(
     c      d1_in(j)+d2_in(j)+d3_in(j)-integral(i,j))
           g_d03(i)=g_d03(i)
     c     +z3_int(j)*dexp(-alambda(j)*real(i)*delta_tau)*(
     c      d1_in(j)+d2_in(j)+d3_in(j)-integral(i,j))
           enddo
c     update guesses for other eigen values
          do k=1,num_ev
              g_d1(i,k)=(d01_in+d02_in+d03_in-cap_gamma*(1.d0+a)*
     c  real(i)*delta_tau/cnot)*z1_int(k) 
           do j=1,num_ev
              g_d1(i,k)=g_d1(i,k)
     c +z1z1(k,j)*dexp(-alambda(j)*real(i)*delta_tau)*(
     c      d1_in(j)+d2_in(j)+d3_in(j)-integral(i,j))
              enddo
              enddo
          do k=1,num_ev
              g_d3(i,k)=(d01_in+d02_in+d03_in-cap_gamma*(1.d0+a)*
     c  real(i)*delta_tau/cnot)*z3_int(k) 
           do j=1,num_ev
              g_d3(i,k)=g_d3(i,k)
     c +z3z3(k,j)*dexp(-alambda(j)*real(i)*delta_tau)*(
     c      d1_in(j)+d2_in(j)+d3_in(j)-integral(i,j))
              enddo
              enddo

c --------------------------
c     test error
              error_1=dabs(g_d01(i)-d01(i))
              error_3=dabs(g_d03(i)-d03(i))
             print*, 'errors', num_iter,error_1,error_3
              do k=1,num_ev
                 error1(k)=dabs(g_d1(i,k)-d1(i,k))
               error3(k)=dabs(g_d3(i,k)-d3(i,k))
             print*, error1(k),error3(k)
              enddo
c --------------------------------
c              reset values
             d01(i)= g_d01(i)
             d03(i)= g_d03(i)
          do k=1,num_ev
              d1(i,k)=g_d1(i,k)
              d3(i,k)=g_d3(i,k)
              enddo
              if (error_1.gt.1.d-07) go to 1
              if (error_3.gt.1.d-07) go to 1
              do k=1,num_ev
              if (error1(k).gt.1.d-07) go to 1
              if (error3(k).gt.1.d-07) go to 1
              enddo
              num_iter=num_iter+1
              if (num_iter.gt.num_iter_max) go to 100
             print*, 'i',i,'iteration converged'
             kay_1(i)=cap_gamma/cnot/d01(i)
             kay_2(i)=cap_gamma/cnot/d03(i)
c     print results
             tau=real(i)*delta_tau
c             if (i.eq.1) write (14,*) i, tau,(integral(i,j),j=1,num_ev),
c     c  kay_1(i),kay_2(i)
        write (14,*) i, tau,(integral(i,j),j=1,num_ev),
     c kay_1(i),kay_2(i)
             
c       if (mod(i,10)==0) write (14,*) i, tau,(integral(i,j),j=1,num_ev),
c     c kay_1(i),kay_2(i)

       
c     test if maximum time steps have been reached
             if (i.lt.i_max) go to 2
             go to 110
 100         print*, 'iteration did not converge'
 110         continue
       stop
       end

      doubleprecision  function coscos(y1,y2,x1,x2)
      doubleprecision y1,y2,x1,x2 !x1 and x2 are the bottom and top limits of integration
      if (y1.eq.y2) then
         coscos=0.5d0*(x2-x1)+(dsin(2.d0*y1*x2)-dsin(2.d0*y1*x1))
     c        /(4.d0*y1)
         else
      coscos=0.5d0*(dsin((y1+y2)*x2)-dsin((y1+y2)*x1))
     c /(y1+y2) +
     c 0.5d0*(dsin((y1-y2)*x2)-dsin((y1-y2)*x1))
     c           /(y1-y2)
      endif
c      write (*,*) coscos
      return
      end

      doubleprecision  function sinsin(y1,y2,x1,x2)
      doubleprecision y1,y2,x1,x2
c      write (*,*) y1,y2
c     write (*,*) x1,x2
      if (y1.eq.y2) then
         sinsin=0.5d0*(x2-x1)-(dsin(2.d0*y1*x2)-dsin(2.d0*y1*x1))
     c        /(4.d0*y1)
         else
      sinsin=0.5d0*(dsin((y1-y2)*x2)-dsin((y1-y2)*x1))
     c /(y1-y2) -
     c 0.5d0*(dsin((y1+y2)*x2)-dsin((y1+y2)*x1))
     c           /(y1+y2)
      endif
c      write (*,*) sinsin
      return
      end

      doubleprecision  function sincos(y1,y2,x1,x2)
      doubleprecision y1,y2,x1,x2 !NOTE THE ORDER y1 IS THE ARGUMENT OF SIIN
c      write (*,*) y1,y2
c     write (*,*) x1,x2
      if (y1.eq.y2) then
         sincos=-(dcos(2.d0*y1*x2)-dcos(2.d0*y1*x1))/
     c        (4.d0*y1)
         else
       sincos=-0.5d0*(dcos((y1+y2)*x2)-dcos((y1+y2)*x1))
     c /(y1+y2) -
     c 0.5d0*(dcos((y1-y2)*x2)-dcos((y1-y2)*x1))
     c           /(y1-y2)
       endif
c       write (*,*) sincos
       return
       end

