        program iteration
c     This is for iterating to find K_1 and K_2
c     This is for uniform concentration ic
c     This is for nonlinear kinetics

c     NOTE: UNLIKE EARLIER PROGRAMS, i=1 IS THE INITIAL CONDITION!!!!
      implicit none

        doubleprecision  delta1,delta2,delta3,epsilon1,epsilon2
        doubleprecision alpha,beta,gamma,dummy,xi,deltaxi
        doubleprecision cap_gamma,cnot,delta_tau,tau,tau_max
        doubleprecision cs1max,cs3max
        
c     data on ef and ev
        doubleprecision alambda(20),b1(20),a2(20),b2(20),a3(20),b3(20)
        doubleprecision z1(21,20),z2(21,20),z3(21,20) ! These are eigenfunctions (j,k)
c     data on initial inner products from older version
        doubleprecision d01_in,d03_in ! these are defined in notes on LINEAR version
       doubleprecision d1_in(20),d2_in(20),d3_in(20) ! These are d_{j1}(0),d_{j2}(0),d_{j3}(0)  in notes

c     These are definitions of error

       doubleprecision errorudotz(20)
       doubleprecision error1,error2  ! errors in kappa1 and kappa2
       
c     these are variables related to actual iteration
       doubleprecision udotz0(150)
       doubleprecision kay_1(150),kay_2(150) !These are Kappa_1 and Kappa_2
       doubleprecision d01n(150),d03n(150) ! These defined in notes. NOTE d01n_in=d01n(1) etc
       doubleprecision d1n(150,20),d2n(150,20),d3n(150,20) !first is time, These are defined in notes. NOTE d1n_in(k)=d01n(1,k) etc
       doubleprecision u1(150,21),u2(150,21),u3(150,21),u4(150,21) !These are u1(i,j) as defined in the notes
       doubleprecision u5(150,21),udotz(150,20) ! these are <u,z>(i,k)
       doubleprecision integral(150,20) ! first is time, second is eigen value, this is integral in notes


c     These are for  guesses
       doubleprecision gkay_1(150),gkay_2(150)
       doubleprecision gudotz(150,20)

       integer i,k,m,n,j,num_ev,num_iter,num_iter_max,i_max,nst


c     opening file to read inputs

      open(11,file='input_plotcheq.dat')
      read (11,*) delta1,delta2,delta3,epsilon1,epsilon2

c     reading other data
      open(15,file='operating_non_linear.dat')
      read (15,*) num_ev,nst,cap_gamma,cnot,delta_tau,num_iter_max,
     c     tau_max,cs1max,cs3max
      print*,  num_ev,nst,cap_gamma,cnot,delta_tau,num_iter_max,
     c     tau_max,cs1max,cs3max
      i_max=tau_max/delta_tau
      deltaxi=1.d0/3.d0/real(nst)
      
      
c    read eigenvlaues, constants in eigenfunctions and compute eigenfunctions
       open(12,file='output_normalized_vector.dat') ! this contains all data 
       do k=1,num_ev
       read(12,*) alambda(k),b1(k),a2(k),b2(k),a3(k),b3(k)
       do j=1,nst+1
c     first interval
          xi=real(j-1)*deltaxi
          z1(j,k)=b1(k)*dcos(dsqrt(alambda(k)/delta1)*xi)
c     second interval
          xi=epsilon1+real(j-1)*deltaxi
          z2(j,k)=a2(k)*dsin(dsqrt(alambda(k)/delta2)*xi)
     c          +b2(k)*dcos(dsqrt(alambda(k)/delta2)*xi)
c     third interval
          xi=epsilon2+real(j-1)*deltaxi           
          z3(j,k)=a3(k)*dsin(dsqrt(alambda(k)/delta3)*xi)
     c     +b3(k)*dcos(dsqrt(alambda(k)/delta3)*xi)          
       enddo
       enddo
       
c---------------------------       
c     Setting initial values
c     reading data of inner product of uniform ic with kappa free eigenfunctions
       open(16,file='inner_prod_uniform_ic_diff_5bc.dat') ! 
            read (16,*) dummy, d01_in,dummy,d03_in
            d01n(1)=0.5d0*d01_in
            d03n(1)=d03_in*dsqrt(3.d0)/4.d0
           udotz0(1)=1.0d0
           kay_1(1)=cap_gamma/cnot/d01n(1)
           kay_2(1)=cap_gamma/cnot/d03n(1)
           do k=1,num_ev
               read (16,*) dummy,d1_in(k),dummy,d3_in(k)
               udotz(1,k)=0.0
               integral(1,k)=0.0
               d1n(1,k)=0.5d0*d1_in(k)
               d3n(1,k)=d3_in(k)*dsqrt(3.d0)/4.d0
         enddo           
         do j=1,nst+1
            u1(1,j)=1.0d0
            u2(1,j)=1.0d0
            u3(1,j)=1.0d0
            u4(1,j)=0.5d0
            u5(1,j)=0.75d0
            enddo
c---------------------------------            
c     opening file for recording ouput
         open(unit=14,file='output_iteration_general.dat')
         open(unit=17,file='output_li_content.dat')   
         open(unit=18,file='output_conc_profile.dat')
         open(unit=19,file='output_integral.dat')
c     record input data
       write (14,*) 'num_ev = ',num_ev, 'delta_tau = ', delta_tau
       write (14,*) 'delta xi = ',deltaxi
       write(14,*) ' cap_gamma = ',cap_gamma,'cnot = ',cnot
       write (14,*) 'kappa1(0)= ', kay_1(1),'   kappa2(0)  = ', kay_2(1)   
c----------------------------------
       i=1
c     begin iteration variables
 2            i=i+1

c     make first guesses after success in the previous step for inner products
                 udotz0(i)=udotz0(i-1)
                 kay_1(i)=kay_1(i-1)
                 kay_2(i)=kay_2(i-1)
             do k=1,num_ev
                udotz(i,k)=udotz(i-1,k)
                integral(i,k)=integral(i-1,k)
             enddo
c     initial guesses for concentration variables              
          do j=1,nst+1
             u1(i,j)=u1(i-1,j)
             u2(i,j)=u2(i-1,j)
             u3(i,j)=u3(i-1,j)
            u4(i,j)=u4(i-1,j)
            u5(i,j)=u5(i-1,j)
          enddo    
c---------------------
          num_iter=0

c     this is where iteration starts
 1        continue
          num_iter=num_iter+1
c---------------------------
c     evaluating d0n, d1n etc
           d01n(i)=0.5*deltaxi*(dsqrt(u4(i,1)*(1.d0-u4(i,1))*u1(i,1))+
     c        +dsqrt(u4(i,nst+1)*(1.d0-u4(i,nst+1))*u1(i,nst+1)))
           d03n(i)=0.5*deltaxi*(dsqrt(u5(i,1)*(1.d0-u5(i,1))*u3(i,1))+
     c        +dsqrt(u5(i,nst+1)*(1.d0-u5(i,nst+1))*u3(i,nst+1)))
           do j=2,nst
          d01n(i)=d01n(i)+deltaxi*dsqrt(u4(i,j)*(1.d0-u4(i,j))*u1(i,j))
          d03n(i)=d03n(i)+deltaxi*dsqrt(u5(i,j)*(1.d0-u5(i,j))*u3(i,j))              
       enddo
       do k=1,num_ev
        d1n(i,k)=0.5*deltaxi*(dsqrt(u4(i,1)*(1-u4(i,1))*u1(i,1))*z1(1,k)
     c    +dsqrt(u4(i,nst+1)*(1-u4(i,nst+1))*u1(i,nst+1))*z1(nst+1,k))
        d3n(i,k)=0.5*deltaxi*(dsqrt(u5(i,1)*(1-u5(i,1))*u3(i,1))*z3(1,k)
     c    +dsqrt(u5(i,nst+1)*(1-u5(i,nst+1))*u3(i,nst+1))*z3(nst+1,k))
        do j=2,nst
           d1n(i,k)=d1n(i,k)+
     c          deltaxi*sqrt(u4(i,j)*(1-u4(i,j))*u1(i,j))*z1(j,k)   
           d3n(i,k)=d3n(i,k)+
     c          deltaxi*sqrt(u5(i,j)*(1-u5(i,j))*u3(i,j))*z3(j,k)          
        enddo
        enddo
       
c     computing the integral
          do k=1,num_ev
              integral(i,k)=integral(i-1,k)+
     c 0.5d0*delta_tau*(
     c (d1n(i-1,k)/d01n(i-1)-d3n(i-1,k)/d03n(i-1))*
     c   dexp(alambda(k)*real(i-1)*delta_tau) +
     c (d1n(i,k)/d01n(i)-d3n(i,k)/d03n(i))*
     c   dexp(alambda(k)*real(i)*delta_tau))
           enddo
           
c     update guesses
c     inner products
           do k=1,num_ev
           gudotz(i,k)=(udotz(1,k)-(cap_gamma/cnot)*integral(i,k))
     c            *dexp(-alambda(k)*real(i)*delta_tau)
           enddo
c     update valuesguesses for kappa
           gkay_1(i)=cap_gamma/cnot/d01n(i)
           gkay_2(i)=cap_gamma/cnot/d03n(i)
           
c     Evaluate variables with the updated values
c     concentrations
           do j=1,nst+1
c     first interval
              xi=real(j-1)*deltaxi
              u1(i,j)=udotz0(1)
              do k=1,num_ev
                 u1(i,j)=u1(i,j)+gudotz(i,k)*z1(j,k)
              enddo   
c     second interval
              xi=epsilon1+real(j-1)*deltaxi
              u2(i,j)=udotz0(1)
              do k=1,num_ev
                 u2(i,j)=u2(i,j)+gudotz(i,k)*z2(j,k)
              enddo   
c     third interval
              xi=epsilon2+real(j-1)*deltaxi
              u3(i,j)=udotz0(1)
              do k=1,num_ev
                 u3(i,j)=u3(i,j)+gudotz(i,k)*z3(j,k)
              enddo   
           enddo
c     update guesses for lithium content variables
           do j=1,nst+1
              u4(i,j)=u4(i-1,j)+(cnot/cs1max)*0.5d0*delta_tau*(
     c     kay_1(i-1)*dsqrt(u4(i-1,j)*(1.d0-u4(i-1,j))*u1(i-1,j))
     c             +  gkay_1(i)*dsqrt(u4(i,j)*(1.d0-u4(i,j))*u1(i,j)))
               u5(i,j)=u5(i-1,j)-(cnot/cs3max)*0.5d0*delta_tau*(
     c     kay_2(i-1)*dsqrt(u5(i-1,j)*(1.d0-u5(i-1,j))*u3(i-1,j))
     c             +  gkay_2(i)*dsqrt(u5(i,j)*(1.d0-u5(i,j))*u3(i,j)))
            enddo
c            --------------------------
c     test error
              error1=dabs(gkay_1(i)-kay_1(i))
              error2=dabs(gkay_2(i)-kay_2(i))
c             print*, 'errors', num_iter,error1,error2
              do k=1,num_ev
                 errorudotz(k)=dabs(gudotz(i,k)-udotz(i,k))
c             print*, errorudotz(k)
              enddo
c --------------------------------
c              reset values
             kay_1(i)= gkay_1(i)
             kay_2(i)= gkay_2(i)
          do k=1,num_ev
              udotz(i,k)=gudotz(i,k)
           enddo
c----------------------------------           
c     exit if iterations exceeded the maximum
              if (num_iter.gt.num_iter_max) go to 100
c-----------------------------
c     if errors are large, iterate
              if (error1.gt.1.d-07) go to 1
              if (error2.gt.1.d-07) go to 1
              do k=1,num_ev
              if (errorudotz(k).gt.1.d-07) go to 1
              enddo
             print*, 'i',i,'iteration converged'
c---------------------
c     print results
             tau=real(i-1)*delta_tau
c             if (i.eq.1) write (14,*) i, tau,(integral(i,j),j=1,num_ev),
c     c  kay_1(i),kay_2(i)

          write (14,*) i, tau,kay_1(i),kay_2(i),(udotz(i,k),k=1,num_ev)
          write (19,*)  i,(integral(i,k),k=1,num_ev)
           write (17,*) i,(u4(i,j), u5(i,j),j=1,nst+1)
           write (18,*) i,(j,u1(i,j),j=1,nst+1),(j,u2(i,j),j=1,nst+1),
     c      (j,u3(i,j),j=1,nst+1) 
             
c       if (mod(i,10)==0) write (14,*) i, tau,(integral(i,j),j=1,num_ev),
c     c kay_1(i),kay_2(i)

c       print*,i_max
c     test if maximum time steps have been reached
             if (i.lt.i_max) go to 2
             go to 110
 100         print*, 'iteration did not converge'
 110         continue
       stop
       end


