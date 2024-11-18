          program concentration constructs concentration profile
        implicit none

        doubleprecision  delta1,delta2,delta3,epsilon1,epsilon2
        doubleprecision alpha,beta,gamma,dummy,tau_max
       doubleprecision alambda,b1,a2,b2,a3,b3
       dimension alambda(20),b1(20),a2(20),b2(20),a3(20),b3(20)
       doubleprecision a,cap_gamma,cnot,delta_tau,tau,eta
       doubleprecision d01_in,d02_in,d03_in !These are d_{01}(0),d_{02}(0),d_{03}(0)
       doubleprecision d1_in(20),d2_in(20),d3_in(20) ! These are d_{j1}(0),d_{j2}(0),d_{j3}(0)  in notes       
       doubleprecision integral(150,20) ! first is time, second is eigen value, this is integral in notes
       doubleprecision u1(150,30),u2(150,30),u3(150,30)
      integer i,k,m,n,j,l,num_ev,num_iter,num_iter_max,i_max


c     opening file for recording output
       open(16,file='output_conc_profile.dat')

c     opening file to read inputs

      open(11,file='input_plotcheq.dat')
      read (11,*) delta1,delta2,delta3,epsilon1,epsilon2

c     opening data file to read operating parameters
      open(15,file='operating.dat')
      read (15,*) num_ev,a,cap_gamma,cnot,delta_tau,num_iter_max,tau_max
      i_max=tau_max/delta_tau
c     reading data for constants
       open(12,file='output_normalized_vector.dat') ! this contains all data 
c     read eigenvlaues and constants in eigenfunctions
       do k=1,num_ev
       read(12,*) alambda(k),b1(k),a2(k),b2(k),a3(k),b3(k)
      enddo

c     reading data of inner product of uniform ic with kappa free eigenfunctions
       open(13,file='inner_prod_uniform_ic_diff_5bc.dat') ! 
    
c     Setting initial values
            read (13,*) dummy, d01_in,d02_in,d03_in
            do j=1,num_ev
            read (13,*) dummy,d1_in(j),d2_in(j),d3_in(j)
         enddo
         
       open(unit=14,file='output_iteration_general.dat')
       read (14,*)  !skipping the first line of the file
       read (14,*)  !skipping the second line of the file
       do l=1,i_max
             read (14,*)  i, dummy,(integral(i,j),j=1,num_ev)
             print*, i, (integral(i,j),j=1,num_ev)
             tau=real(i)*delta_tau
               do  m=1,11  !10 spatial points
c     m=1
                  eta=epsilon1*0.1*real(m-1)
             u1(i,m)=(d01_in+d02_in+d03_in-cap_gamma*(1.d0+a)*
     c  real(i)*delta_tau/cnot)
               do k=1,num_ev
            u1(i,m)=u1(i,m)+(d1_in(k)+d2_in(k)+d3_in(k)-integral(i,k))
     c     *dexp(-alambda(k)*real(i)*delta_tau)
     c    *b1(k)*dcos(dsqrt(alambda(k)/delta1)*eta)
         enddo
            write (16,*) i, tau,eta, u1(i,m)
            enddo
               do  m=1,11  !10 spatial points
c     m=1
               eta=epsilon1+(epsilon2-epsilon1)*0.1*real(m-1)
               u2(i,m)=(d01_in+d02_in+d03_in-cap_gamma*(1.d0+a)*
     c   real(i)*delta_tau/cnot)
               do k=1,num_ev
            u2(i,m)=u2(i,m)+(d1_in(k)+d2_in(k)+d3_in(k)-integral(i,k))
     c     *dexp(-alambda(k)*real(i)*delta_tau)*(
     c      a2(k)*dsin(dsqrt(alambda(k)/delta2)*eta)
     c    +b2(k)*dcos(dsqrt(alambda(k)/delta2)
     c                 *eta))
            enddo
            write (16,*) i,tau, eta,u2(i,m)
            enddo
               do  m=1,11  !10 spatial points
c     m=1
              eta=epsilon2+(1.d0-epsilon2)*0.1*real(m-1)
             u3(i,m)=(d01_in+d02_in+d03_in-cap_gamma*(1.d0+a)*
     c  real(i)*delta_tau/cnot)
               do k=1,num_ev
            u3(i,m)=u3(i,m)+(d1_in(k)+d2_in(k)+d3_in(k)-integral(i,k))
     c     *dexp(-alambda(k)*real(i)*delta_tau)*(
     c      a3(k)*dsin(dsqrt(alambda(k)/delta3)*eta)
     c    +b3(k)*dcos(dsqrt(alambda(k)/delta3)*eta))
         enddo
             write (16,*) i,tau,eta,u3(i,m)
            enddo
            enddo
             stop
             end
